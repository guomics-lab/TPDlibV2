# Reference: DIALib-QC_RPlot.pl


# 1.Reset environment -----------------------------------------------------
pacman::p_unload(pacman::p_loaded(), character.only = T)
rm(list = ls())
libs <- c('RColorBrewer', 'tidyverse', 'magrittr', 'scales', 'ggpubr', 'pheatmap')
table(sapply(libs, require, character.only=TRUE))
rm(libs)


# 2.Figures ---------------------------------------------------------------
## 2.1 Data input ------
# read libraries
libFolder <- '//172.16.13.114/share/members/LiLu/L01_Project/L_TPDlibV2/library_FragPipe_log'
libPaths <- list.files(libFolder, 'library\\.tsv$', full.names = T) %>% str_subset('Isoform|Phospho', negate = T)
names(libPaths) <- str_extract(libPaths, '480|7600|QE|tims')
names(libPaths)[which(names(libPaths) == 'QE')] <- 'Q Exactive HF'
names(libPaths)[which(names(libPaths) == '480')] <- 'Orbitrap Exploris 480'
names(libPaths)[which(names(libPaths) == '7600')] <- 'ZenoTOF 7600'
names(libPaths)[which(names(libPaths) == 'tims')] <- 'timsTOF Pro'

df <- lapply(libPaths, rio::import) %>% plyr::ldply(.id = 'Instrument')

# read fasta
fasPath <- '//172.16.13.136/share/project-thy/TPDlibV2/FASTA/2023-06-17-decoys-contam-20230616_20422_Human_fasta_iRT_nonIsoform.fasta.fas'

# set color style
heat_colors <- viridis::viridis_pal(option = 'H')(4)
instrument_colors <- c("#4DB6AC", "#D4A24C", "#1E4E79", "#A64B2A") %>% setNames(c('ZenoTOF 7600', 'Q Exactive HF', 'Orbitrap Exploris 480', 'timsTOF Pro'))
# allColors <- c('#639d98', '#1e5d64', '#203d26', '#b49259', '#7d2a16')
# allFills <- c('#75bcb6', '#2b828c', '#2f5b38', '#e0b66e', '#a6381d')
# image(x = 1:5, y = 1, z = as.matrix(1:5), col = colorRampPalette(allFills)(5))
df$Instrument %<>% factor(levels = names(instrument_colors))

## 2.2 Draw plots -----
df_mz <- df %>%
  distinct(Instrument, ModifiedPeptideSequence, PrecursorCharge, PrecursorMz)
summary(df_mz$PrecursorMz)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 280.2   516.6   625.8   659.3   766.3  1697.3 
length(df_mz$PrecursorMz[df_mz$PrecursorMz > 400 & df_mz$PrecursorMz < 1200]) / nrow(df_mz) * 100
# 91.22416%

# Plot D: precursor m/z distribution
y_anno4 <- 550
p4 <- ggplot(df_mz, aes(x = PrecursorMz))+
  geom_density(aes(y = (..count..), color = Instrument, fill = Instrument), linewidth = 1, alpha = 0.1)+
  # geom_vline(xintercept = c(400, 1200, median(df_mz$PrecursorMz))) +
  # annotate('text', label = round(median(df_mz$PrecursorMz), 2), x = median(df_mz$PrecursorMz), y = y_anno4) +
  # annotate('text', label = 400, x = 400, y = y_anno4) +
  # annotate('text', label = 1200, x = 1200, y = y_anno4) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  scale_color_manual(values = instrument_colors) +
  scale_fill_manual(values = instrument_colors) +
  labs(x = "Precursor m/z (Th)", y = "Frequency")+
  theme_bw()+
  theme(text = element_text(size = 10), legend.position = 'top')
# p4 <- ggplot(df_mz)+
#   aes(x = Instrument, y = PrecursorMz, fill = Instrument) +
#   geom_violin(color = NA) +
#   geom_boxplot(color = '#666666', fill = '#FFFFFF', width = 0.1, linewidth = 1, outlier.shape = NA) +
#   scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
#   scale_fill_manual(values = instrument_colors) +
#   labs(x = '', y = 'Precursor m/z (Th)')+
#   theme_bw()+
#   theme(text = element_text(size = 15),
#         axis.text.x = element_blank(),
#         legend.position = 'top')

# Plot E: precursor charge distribution
# label_p5 <- str_c(round(table(df_mz$PrecursorCharge) / nrow(df_mz) * 100, 2), '% (', table(df_mz$PrecursorCharge), ')') %>% setNames(sort(unique(df_mz$PrecursorCharge)))
p5A <- ggplot(df_mz) +
  aes(x = PrecursorCharge, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_x_discrete(limits = factor(1:6)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Precursor m/z (Th)', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')

tbl5 <- df_mz %>%
  count(Instrument, PrecursorCharge) %>% 
  pivot_wider(id_cols = 'Instrument', names_from = 'PrecursorCharge', values_from = 'n', values_fill = 0) %>% 
  column_to_rownames('Instrument') %>%
  set_colnames(., str_c('Charge=', colnames(.)))
p5B <- pheatmap(tbl5, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p5 <- ggarrange(p5A, p5B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


df_rt <- df %>%
  distinct(Instrument, ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime) %>% 
  filter(PrecursorCharge %in% 2:3)
# Plot K: +2/+3 RT linear regression
p11_ls <- plyr::dlply(df_rt, 'Instrument', function(dfsub){
  dfcor <- dfsub %>%
    count(ModifiedPeptideSequence) %>%
    filter(n == 2) %>% 
    semi_join(dfsub, .) %>% 
    pivot_wider(id_cols = ModifiedPeptideSequence,
                names_from = 'PrecursorCharge',
                values_from = 'NormalizedRetentionTime')
  lbl <- paste0('n=', nrow(dfcor), ', R2=', sprintf("%.2f", cor(dfcor$`2`, dfcor$`3`, method = "pearson")^2))
  ggplot(dfcor, aes(x = `2`, y = `3`))+
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 500) +
    scale_fill_continuous(low = "white", high = instrument_colors[dfsub$Instrument[1]])+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    #geom_point(shape = 20, size = 1, color = "dodgerblue4", alpha = 0.1)+
    labs(x = "+2 iRT", y = "+3 iRT", subtitle = dfsub$Instrument[1])+
    annotate("text", x = min(dfcor$`2`), y = max(dfcor$`3`) * 1.05, label = lbl, hjust = 0, vjust = 0, size = 4)+
    theme_bw()+
    theme(text = element_text(size = 10), legend.position = 'none')
})
p11 <- ggarrange(plotlist = p11_ls)


df_len <- df %>% distinct(PeptideSequence) %>% 
  mutate(PeptideLength = str_count(PeptideSequence))
length(df_len$PeptideLength[df_len$PeptideLength >= 8 & df_len <= 20]) / length(df_len$PeptideLength)
# 0.8037977 (164158)
# Plot F: Peptide length distribution
tbl6 <- df_len %>%
  count(Instrument, PeptideLength)
p6 <- ggplot(tbl6, aes(x = PeptideLength))+
  geom_ribbon(aes(ymin = 0, ymax = n, fill = Instrument), alpha = 0.1)+
  geom_line(aes(y = n, color = Instrument), linewidth = 1) +
  geom_point(aes(y = n, color = Instrument)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  scale_color_manual(values = instrument_colors) +
  scale_fill_manual(values = instrument_colors) +
  labs(x = "Peptide length", y = "Frequency")+
  theme_bw()+
  theme(text = element_text(size = 10), legend.position = 'top')


df_mod <- df %>% distinct(Instrument, ModifiedPeptideSequence)
df_mod['[+42]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:1)', e, fixed = T)) }))
df_mod['[+57]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:4)', e, fixed = T)) }))
df_mod['[+16]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:35)', e, fixed = T)) }))
# consider_phospho <- F
# if(consider_phospho){
#   df_mod['[+80]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:21)', e, fixed = T)) }))
# }
df_mod <- df_mod %>%
  pivot_longer(cols = -c('Instrument', 'ModifiedPeptideSequence'),
               names_to = 'Modification', values_to = 'Flag') %>% 
  filter(Flag)

# Plot G: Modification counts
p7A <- ggplot(df_mod) +
  aes(x = Modification, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Mass shift of modification', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')

tbl7 <- df_mod %>%
  count(Instrument, Modification) %>%
  pivot_wider(id_cols = 'Instrument', names_from = 'Modification', values_from = 'n', values_fill = 0) %>%
  column_to_rownames('Instrument')
p7B <- pheatmap(tbl7, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p7 <- ggarrange(p7A, p7B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


df_peppro <- df %>%
  distinct(Instrument, ProteinId, PeptideSequence) %>%
  count(Instrument, ProteinId, name = 'PeptidePerProtein') %>% 
  mutate(PeptidePerProtein = ifelse(PeptidePerProtein < 8, PeptidePerProtein, '8 (>=8)'))
# Plot I: Peptides per protein
p9A <- ggplot(df_peppro) +
  aes(x = PeptidePerProtein, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Peptide per protein', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')
tbl9 <- df_peppro %>%
  count(Instrument, PeptidePerProtein) %>%
  pivot_wider(id_cols = 'Instrument', names_from = 'PeptidePerProtein', values_from = 'n', values_fill = 0) %>%
  column_to_rownames('Instrument')
p9B <- pheatmap(tbl9, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p9 <- ggarrange(p9A, p9B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


df$Precursor <- paste0(df$ModifiedPeptideSequence, '_', df$'PrecursorCharge')
df_frpr <- df %>%
  distinct(Instrument, Precursor, Annotation) %>%
  count(Instrument, Precursor, name = 'FragmentPerPrecursor') %>% 
  mutate(FragmentPerPrecursor = ifelse(FragmentPerPrecursor < 6, FragmentPerPrecursor, '6 (>=6)'))

# Plot C: Fragments per precursor
p3A <- ggplot(df_frpr) +
  aes(x = FragmentPerPrecursor, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Fragments per precursor', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')
tbl3 <- df_frpr %>%
  count(Instrument, FragmentPerPrecursor) %>%
  pivot_wider(id_cols = 'Instrument', names_from = 'FragmentPerPrecursor', values_from = 'n', values_fill = 0) %>%
  column_to_rownames('Instrument')
p3B <- pheatmap(tbl3, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p3 <- ggarrange(p3A, p3B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


df_frtype <- df %>%
  distinct(Instrument, Precursor, Annotation, FragmentType)
# Plot A: Fragment ion type
p1A <- ggplot(df_frtype) +
  aes(x = FragmentType, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Fragment ion type', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')
tbl1 <- df_frtype %>%
  count(Instrument, FragmentType) %>%
  pivot_wider(id_cols = 'Instrument', names_from = 'FragmentType', values_from = 'n', values_fill = 0) %>%
  column_to_rownames('Instrument')
p1B <- pheatmap(tbl1, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p1 <- ggarrange(p1A, p1B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


df_frz <- df %>% distinct(Instrument, Precursor, Annotation, FragmentCharge)
# Plot B: Fragment ion charge
p2A <- ggplot(df_frz) +
  aes(x = FragmentCharge, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Fragment ion charge', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')
tbl2 <- df_frz %>%
  count(Instrument, FragmentCharge) %>%
  pivot_wider(id_cols = 'Instrument', names_from = 'FragmentCharge', values_from = 'n', values_fill = 0) %>%
  column_to_rownames('Instrument')
p2B <- pheatmap(tbl2, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p2 <- ggarrange(p2A, p2B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


print('Calculating protein sequence coverage which will cost several minutes')
fastafile <- seqinr::read.fasta(file = fasPath, seqtype = 'AA', as.string = TRUE)
df_fas <- data.frame(protein = names(fastafile) %>% stringr::str_split(pattern = '\\|') %>% sapply(., function(e){e[2]}), sequence = unlist(fastafile))
df_fas <- df_fas[stringr::str_detect(rownames(df_fas), '^rev_', negate = T), ] # remove "rev_" decoys
protein_coverage_cal <- function(df, df_fas){
  # 20,000 proteins cost 13 minutes
  prots <- sort(unique(df$ProteinId))
  coverage <- rep(NaN, length(prots))
  for(i in seq_along(prots)){
    cat(i, '--', prots[i], '...\r')
    prot <- prots[i]
    peptides <- unique(df$PeptideSequence[df$ProteinId == prot])
    full_seq <- df_fas$sequence[df_fas$protein == prot]
    covered <- rep(0, nchar(full_seq))# recording matched peptide sequence with number 1
    loc <- na.omit(stringr::str_locate(full_seq, peptides))
    if(nrow(loc) == 0) next;
    loc <- lapply(seq_len(nrow(loc)), function(j) {loc[j, 1]:loc[j, 2]} )
    loc <- Reduce(base::union, loc)
    covered[loc] <- 1
    coverage[i] <- sum(covered) / length(covered)
  }
  return(coverage)
}
coverage_list <- list()
for(i in 1:4){
  cat('Instrument', i, '...\n')
  coverage_list[[i]] <- protein_coverage_cal(split(df, df$Instrument)[[i]], df_fas)
}
df_coverage <- coverage_list %>% setNames(unique(df$Instrument)) %>% 
  plyr::ldply(function(vec){
    data.frame(Coverage = vec)
  }, .id = 'Instrument')
# Plot J: protein sequence coverage
p10 <- ggplot(df_coverage, aes(x = Coverage))+
  geom_density(aes(y = (..count..), color = Instrument, fill = Instrument), linewidth = 1, alpha = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  scale_color_manual(values = instrument_colors) +
  scale_fill_manual(values = instrument_colors) +
  labs(x = "Protein coverage (%)", y = "Frequency")+
  theme_bw()+
  theme(text = element_text(size = 10), legend.position = 'top')


print('Calculating missed cleavage')
df_peps <- df %>% distinct(PeptideSequence)
df_peps$MissedCleavage <- 0
df_peps$MissedCleavage <- apply(df_peps, 1, function(row){
  Seq <- row['PeptideSequence']
  # exclude the last cleavage site
  if (stringr::str_sub(Seq, -1, -1) == 'P' & stringr::str_sub(Seq, -2, -2) %in% c('K', 'R')){
    Seq <- stringr::str_sub(Seq, 1, -3)
  }else{
    Seq <- stringr::str_sub(Seq, 1, -2)
  }
  #
  Seqs <- stringr::str_sub(Seq, 1:nchar(Seq), 1:nchar(Seq))
  is_KR <- Seqs %in% c('K', 'R')
  isnt_P <- Seqs != 'P'
  flag <- c(isnt_P[2:length(isnt_P)], T)
  mis_cleav <- is_KR & flag
  return(sum(mis_cleav))
})
df_missCleavage <- df %>%
  distinct(Instrument, PeptideSequence) %>%
  left_join(df_peps) %>% 
  mutate(MissedCleavage = factor(MissedCleavage))
# Plot H: missed cleavage
p8A <- ggplot(df_missCleavage) +
  aes(x = MissedCleavage, fill = Instrument) +
  geom_bar(stat = 'count', color = '#000000', width = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01), add = c(0, 0)))+
  scale_fill_manual(values = instrument_colors) +
  labs(x = 'Missed cleavage', y = 'Frequency')+
  theme_classic()+
  theme(text = element_text(size = 10), legend.position = 'none')

tbl8 <- df_missCleavage %>%
  count(Instrument, MissedCleavage) %>% 
  pivot_wider(id_cols = 'Instrument', names_from = 'MissedCleavage', values_from = 'n', values_fill = 0) %>% 
  column_to_rownames('Instrument') %>%
  set_colnames(., str_c(colnames(.)))
p8B <- pheatmap(tbl8, scale = 'none', cluster_cols = F, cluster_rows = F,
                show_rownames = F, angle_col = 0,
                na_col = '#999999', border_color = '#000000',
                color = '#FFFFFF', legend = F, annotation_legend = F,
                display_numbers = T, number_format = '%d', number_color = '#000000',
                fontsize = 10, fontsize_number = 10,
                annotation_row = data.frame(Instrument = names(instrument_colors), row.names = names(instrument_colors)),
                annotation_colors = list(Instrument = instrument_colors)) %>% 
  ggplotify::as.ggplot()
p8 <- ggarrange(p8A, p8B, nrow = 2, ncol = 1, heights = c(0.75, 0.25)) +
  theme(plot.margin = margin(20, 20, 20, 20))


# IM-RT plot
df_imrt <- df %>%
  filter(Instrument == 'timsTOF Pro') %>% 
  select(ModifiedPeptideSequence, PrecursorCharge, AverageExperimentalRetentionTime, PrecursorIonMobility) %>%
  distinct() %>% 
  rename(RT = AverageExperimentalRetentionTime, IM = PrecursorIonMobility) %>% 
  mutate(RT = RT / 60) %>% 
  rename(`RT (min)` = RT)
# Plot L: IM-RT for PASEF data
p12 <- ggplot(df_imrt, aes(x = `RT (min)`, y = IM))+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 300) +
  scale_fill_continuous(low = "white", high = instrument_colors[4])+ # dodgerblue
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  geom_point(shape = 20, size = 1, color = instrument_colors[4], alpha = 0.1, stroke = NA)+
  labs(x = "RT (min)", y = "Ion mobility (V·s·cm-2)", subtitle = 'Library modified precursors RT-IM plot')+
  theme_bw()+
  theme(text = element_text(size = 10), legend.position = 'none')



# save.image('DIALib-QC_RPlot_20240913combine.RData')
source_data <- list(
  Figure2A = tbl1 %>% rownames_to_column('Instrument'),
  Figure2B = tbl2 %>% rownames_to_column('Instrument'),
  Figure2C = tbl3 %>% rownames_to_column('Instrument'),
  Figure2D = df_mz,
  Figure2E = tbl5 %>% rownames_to_column('Instrument'),
  Figure2F = tbl6,
  Figure2G = tbl7 %>% rownames_to_column('Instrument'),
  Figure2H = tbl8 %>% rownames_to_column('Instrument'),
  Figure2I = tbl9 %>% rownames_to_column('Instrument'),
  Figure2J = df_coverage,
  Figure2K = df_rt,
  Figure2L = df_imrt
)
rio::export(source_data, 'DIALib-QC_RPlot_20240913combine_souceData.xlsx', row.names = T)
p <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, common.legend = T, labels = LETTERS[1:12], font.label = 'bold', nrow = 4, ncol = 3)
ggsave('DIALib-QC_RPlot_20240913combine.pdf', p, width = 5*3, height = 5*4)

pdf('DIALib-QC_RPlot_20240913combine_F2L_density_curve.pdf', width = 6, height = 3)
plot(density(df_imrt$`RT (min)`, bw = 0.1))
plot(density(df_imrt$IM, bw = 0.01))
graphics.off()



