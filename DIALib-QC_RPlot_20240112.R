# Reference: DIALib-QC_RPlot.pl

# --------- reset environment --------------------
rm(list = ls())
libs <- c('RColorBrewer', 'ggplot2', 'scales', 'ggpubr')
table(sapply(libs, require, character.only=TRUE))
rm(libs)

# ------- function: plotting ------------------------
plotting <- function(df, consider_phospho = FALSE, fasPath = '//172.16.13.136/share/project-thy/TPDlibV2/FASTA/2023-06-17-decoys-contam-20230616_20422_Human_fasta_iRT_nonIsoform.fasta.fas', pdfName, rlt_dir = './', allColor = c('#A2565B', '#6E8E84', '#1A476F', '#E37E00', '#90353A'),
                     allFill = c('#C79A9D', '#B7C7C2', '#8DA3B7', '#F1BE80', '#C89A9D')){
  # DIA-LibQC (For FragPipe-Easypqp library format)
  df_mz <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, PrecursorMz)
  
  # Plot A: precursor m/z distribution
  p1 <- ggplot(df_mz, aes(x = PrecursorMz))+
    geom_histogram(aes(y = (..count..)/sum(..count..)),
                   color = allColor[1],
                   fill = allFill[1])+
    geom_vline(xintercept = c(400, 1000, median(df_mz$PrecursorMz))) +
    annotate('text', label = median(df_mz$PrecursorMz), x = median(df_mz$PrecursorMz), y = 0.1) +
    annotate('text', label = 400, x = 400, y = 0.1) +
    annotate('text', label = 1000, x = 1000, y = 0.1) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = '', x = "Precursor m/z", y = "Frequency")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  # Plot B: precursor charge distribution
  label_p2 <- stringr::str_c(round(table(df_mz$PrecursorCharge) / nrow(df_mz) * 100, 2), '% (', table(df_mz$PrecursorCharge), ')') %>% setNames(sort(unique(df_mz$PrecursorCharge)))
  p2 <- ggplot(df_mz, aes(x = PrecursorCharge))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 0.6,
             color = allColor[2],
             fill = allFill[2],
    )+
    geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = '', x = "Precursor charge", y = "Frequency")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_rt <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime)
  df_rt <- df_rt[df_rt$PrecursorCharge %in% c(2, 3), ]
  tb_tmp <- table(df_rt$ModifiedPeptideSequence)
  pairs <- names(tb_tmp[tb_tmp == 2])
  df_rt <- df_rt[df_rt$ModifiedPeptideSequence %in% pairs, ]
  df_rt <- reshape2::dcast(df_rt, ModifiedPeptideSequence~PrecursorCharge, value.var = 'NormalizedRetentionTime')
  colnames(df_rt) <- c('modSeq', 'RT2value', 'RT3value')
  # Plot C: +2/+3 RT linear regression
  correlation <- sprintf("%1f", cor(df_rt$RT2value, df_rt$RT3value, method = "pearson"))
  N <- nrow(df_rt)
  lbl = paste0('n=', N, ', R2=', correlation)
  # scatter plot
  p3 <- ggplot(df_rt, aes(x = RT2value, y = RT3value))+
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 300) +
    scale_fill_continuous(low = "white", high = allFill[3])+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = 'Library +2/+3 pair iRT correlation', x = "+2 iRT", y = "+3 iRT")+
    annotate("text", x = min(df_rt$RT2value), y = max(df_rt$RT3value) * 1.05, label = lbl, hjust = 0, vjust = 0, size = 5)+
    theme_bw()+
    theme(
      legend.position = 'none',
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  df_len <- df %>% dplyr::distinct(PeptideSequence)
  df_len$PeptideLength <- nchar(df_len$PeptideSequence)
  length(df_len$PeptideLength[df_len$PeptideLength >= 8 & df_len <= 20]) / length(df_len$PeptideLength)
  # 0.8056372 (152416)
  # Plot D: Peptide length distribution
  p4 <- ggplot(df_len, aes(x = PeptideLength))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 1,
             color = allColor[4],
             fill = allFill[4])+
    geom_vline(xintercept = c(8, 20, median(df_len$PeptideLength))) +
    annotate('text', label = median(df_len$PeptideLength), x = median(df_len$PeptideLength), y = 0.1) +
    annotate('text', label = 8, x = 8, y = 0.1) +
    annotate('text', label = 20, x = 20, y = 0.1) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = '', x = 'Peptide length', y = 'Frequency')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_mod <- df %>% dplyr::distinct(ModifiedPeptideSequence)
  df_mod['[+42]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:1)', e, fixed = T)) }))
  df_mod['[+57]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:4)', e, fixed = T)) }))
  df_mod['[+16]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:35)', e, fixed = T)) }))
  if(consider_phospho){
    df_mod['[+80]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:21)', e, fixed = T)) }))
  }
  vec_mod <- unlist(lapply(df_mod[, -1], sum))
  df_mod <- data.frame(mod_type = names(vec_mod), freq = vec_mod)
  # Plot E: Modification counts
  p5 <- ggplot(df_mod, aes(mod_type, freq))+
    geom_col(color = allColor[5],
             fill = allFill[5],
             width = 0.5)+
    geom_text(aes(label = freq),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = '', x = 'Modification type', y = '# Modifications')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_peppro <- df %>% dplyr::distinct(ProteinId, PeptideSequence)
  df_peppro <- data.frame(table(table(df_peppro$ProteinId))); colnames(df_peppro) <- c('pep_num', 'count')
  if (sum(as.numeric(as.character(df_peppro$pep_num)) >= 8)){
    df_peppro$pep_num <- as.numeric(as.character(df_peppro$pep_num))
    over8 <- sum(df_peppro[df_peppro$pep_num >= 8, 'count'])
    df_peppro <- df_peppro[df_peppro$pep_num < 8, ]
    df_peppro <- rbind(df_peppro, c('>=8', over8))
    df_peppro$pep_num <- factor(df_peppro$pep_num, levels = c('1', '2', '3', '4', '5', '6', '7', '>=8'), ordered = T)
    df_peppro$count <- as.numeric(df_peppro$count)
  }
  # Plot F: Peptides per protein
  p6 <- ggplot(df_peppro, aes(pep_num, count))+
    geom_col(color = allColor[6], fill = allFill[6])+
    geom_text(aes(label = count), size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = '', x = 'Peptides per protein', y = '# Protiens')+
    scale_x_discrete("Peptides per protein",
                     #limits = 
    )+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  df_pie <- df_peppro %>%
    dplyr::arrange(desc(pep_num)) %>%
    dplyr::mutate(prop = count / sum(.$count) * 100) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5 * prop) %>%
    dplyr::mutate(pie_text = stringr::str_glue("{pep_num}\n{sprintf('%.2f', prop)}%\n({count})"))

  # Basic piechart
  p_pie <- ggplot(df_pie, aes(x = "", y = prop, fill = pep_num)) +
    geom_bar(stat = "identity", width = 1, linewidth = 1, color = allColor[6], fill= "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "right") +
    geom_text(aes(x = 1, y = ypos, label = pie_text), color = "black", size = 5)


  df$Precursor <- paste0(df$ModifiedPeptideSequence, '_', df$'PrecursorCharge')
  df_frpr <- df %>% dplyr::distinct(Precursor, Annotation)
  df_frpr <- data.frame(table(table(df_frpr$Precursor))); colnames(df_frpr) <- c('fr_num', 'count')
  if(sum(as.numeric(as.character(df_frpr$fr_num)) >= 6)){
    df_frpr$fr_num <- as.numeric(as.character(df_frpr$fr_num))
    over6 <- sum(df_frpr[df_frpr$fr_num >= 6, 'count'])
    df_frpr <- df_frpr[df_frpr$fr_num < 6, ]
    df_frpr <-  rbind(df_frpr, c('>=6', over6))
    df_frpr$fr_num <- factor(df_frpr$fr_num, levels = c('1', '2', '3', '4', '5', '>=6'), ordered = T)
    df_frpr$count <- as.numeric(df_frpr$count)
  }
  df_frpr$freq <- df_frpr$count / sum(df_frpr$count)
  # Plot G: Fragments per precursor
  p7 <- ggplot(df_frpr, aes(fr_num, freq))+
    geom_col(width = 0.5, color = allColor[7], fill = allFill[7])+
    geom_text(aes(label = paste0(sprintf('%.3f', freq * 100), '%')), size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = '', x = 'Fragments per precursor', y = 'Frequency')+
    scale_x_discrete("Fragments per precursor",
                     #limits = 
    )+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frtype <- df %>% dplyr::distinct(Precursor, Annotation, FragmentType)
  df_frtype <- data.frame(table(df_frtype$FragmentType)); colnames(df_frtype) <- c('fr_type', 'count')
  df_frtype$freq <- df_frtype$count / sum(df_frtype$count)
  # Plot H: Fragment ion type distribution
  p8 <- ggplot(df_frtype, aes(fr_type, freq))+ theme_bw()+
    geom_col(color = allColor[8], fill = allFill[8], width = 0.3)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = '', x = 'Fragment ion type', y = 'Frequency')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frz <- df %>% dplyr::distinct(Precursor, Annotation, FragmentCharge)
  df_frz <- data.frame(table(df_frz$FragmentCharge)); colnames(df_frz) <- c('fr_z', 'count')
  df_frz$freq <- df_frz$count / sum(df_frz$count)
  # Plot I: Fragment ion charge distribution
  p9 <- ggplot(df_frz, aes(fr_z, freq))+ theme_bw()+
    geom_col(color = allColor[9], fill = allFill[9], width = 0.5)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    labs(title = '', x = 'Fragment ion charge', y = 'Frequency')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  print('Calculating protein sequence coverage which will cost several minutes')
  fastafile <- seqinr::read.fasta(file = fasPath, seqtype = 'AA', as.string = TRUE)
  df_fas <- data.frame(protein = names(fastafile) %>% stringr::str_split(pattern = '\\|') %>% sapply(., function(e){e[2]}), sequence = unlist(fastafile))
  df_fas <- df_fas[stringr::str_detect(rownames(df_fas), '^rev_', negate = T), ] # remove "rev_" decoys
  prots <- sort(unique(df$ProteinId))
  
  # 20,000 proteins cost 13 minutes
  coverage <- rep(NaN, length(prots))
  for(i in seq_along(prots)){
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
  
  # Plot J: protein sequence coverage
  p10 <- ggplot(data.frame(coverage = coverage * 100), aes(x = coverage)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)),
                   color = allColor[10],
                   fill = allFill[10])+
    labs(title = '', x = "Protein coverage (%)", y = "Frequency")+
    geom_vline(xintercept = c(median(coverage * 100))) +
    annotate('text', label = median(coverage * 100), x = median(coverage * 100), y = 0.075) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  print('Calculating missed cleavage')
  df_peps <- data.frame(PeptideSequence = unique(df$PeptideSequence))
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
  # Plot K: missed cleavage
  df_mis <- as.data.frame(table(df_peps$MissedCleavage)) %>% apply(., c(1, 2), as.character)
  vec_mis <- c()
  for(i in 1:nrow(df_mis)){
    vec_mis <- append(vec_mis, rep(df_mis[i, 1], df_mis[i, 2]))
  }
  p11 <- ggplot(data.frame(V1 = vec_mis), aes(x=`V1`)) +
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 0.6,
             position="stack",
             color = allColor[11],
             fill = allFill[11],
    )+
    geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', size = 6, hjust = 0.5, vjust = -0.2, position = "stack", color = 'black')+
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
    labs(title = '', x = "Missed cleavage", y = "Frequency")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))

  p <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p_pie, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), font.label = list(size = 14, face = "bold"), ncol = 3, nrow = 4)
  pdf(file = paste0(rlt_dir, '/', pdfName), width = 8.3 *2 , height = 11.7 * 2)
  print(p)
  graphics.off()
}


# --------- main DIA-LibQC -------------
# non-isoform library
libPath <- '//172.16.13.114/share/members/LiLu/L01_Project/L_TPDlibV2/library_FragPipe_log/480_240DDA_FP_library.tsv'
df <- data.table::fread(libPath, sep = '\t', stringsAsFactors = F, check.names = F, data.table = F)
# read fasta
fasPath <- '//172.16.13.136/share/project-thy/TPDlibV2/FASTA/2023-06-17-decoys-contam-20230616_20422_Human_fasta_iRT_nonIsoform.fasta.fas'
# set color style
allColors <- c('#639d98', '#1e5d64', '#203d26', '#b49259', '#7d2a16')
allFills <- c('#75bcb6', '#2b828c', '#2f5b38', '#e0b66e', '#a6381d')
image(x = 1:5, y = 1, z = as.matrix(1:5), col = colorRampPalette(allFills)(5))

plotting(df, fasPath = fasPath, pdfName = '480_240DDA_FP_library_20231127.pdf',
         allColor = allColors[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)],
         allFill = allFills[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)])


# phospho
df_p <- data.table::fread('//172.16.13.136/share/members/lilu/000_Project/011_TC60/Library/TC60_Phospho_tims_FragPipe_library.tsv', sep = '\t', stringsAsFactors = F, check.names = F, data.table = F)
plotting(df_p, consider_phospho = T, fasPath = fasPath, pdfName = 'TC60_Phospho_tims_FragPipe_library_20231126.pdf',
         allColor = allColors[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)],
         allFill = allFills[c(1, 2, 3, 4, 5, 2, 3, 1, 4, 2, 4)])



# # --------- comparison with some database -------------
# use: 
# #transition groups
# #modified precursors
# #modified peptides
# #proteotypic proteins

df <- rio::import('//172.16.13.114/share/members/LiLu/L01_Project/L_TPDlibV2/library_FragPipe_log/QE_240DDA_FP_library.tsv')
dfprot <- rio::import('//172.16.13.136/share/project-thy/TPDlibV2/QE_Library_20230620/protein.tsv')
setequal(dfprot$`Protein ID`, df$ProteinId) # TRUE

df_old <- rio::import('//172.16.13.136/share/project-thy/TPD_library_2021_final/TPD_SPNlibrary_60min_46files_20210225.csv')

# unify the modification terms 
df_old$ModifiedPeptide %<>%
  str_remove_all('^_|_$') %>% 
  str_replace_all('\\[Oxidation \\(M\\)\\]', '[+16]') %>% 
  str_replace_all('\\[Acetyl \\(Protein N\\-term\\)\\]', '[+42]') %>% 
  str_replace_all('\\[Carbamidomethyl \\(C\\)\\]', '[+57]')
df$ModifiedPeptideSequence %<>% 
  str_replace_all('\\(UniMod:35\\)', '[+16]') %>% 
  str_replace_all('\\(UniMod:1\\)', '[+42]') %>% 
  str_replace_all('\\(UniMod:4\\)', '[+57]')

# add pr, tg
df_old %<>% mutate(
  pr = str_c(ModifiedPeptide, PrecursorCharge, sep = '_'),
  tg = str_c(pr, FragmentType, FragmentNumber, FragmentCharge, FragmentLossType, sep = '_')
)

df %<>% mutate(
  pr = str_c(ModifiedPeptideSequence, PrecursorCharge, sep = '_'),
  tg = str_c(pr, FragmentType, FragmentSeriesNumber, FragmentCharge, 'noloss', sep = '_')
)

tgv1 <- unique(df_old$tg)
tgv2 <- unique(df$tg)
prv1 <- unique(df_old$pr)
prv2 <- unique(df$pr)
pev1 <- unique(df_old$ModifiedPeptide)
pev2 <- unique(df$ModifiedPeptideSequence)
ppv1 <- df_old %>% filter(!str_detect(ProteinGroups, ';')) %>% pull(ProteinGroups) %>% unique()
ppv2 <- dfprot %>% filter(`Indistinguishable Proteins` == '') %>% pull(`Protein ID`) %>% intersect(df$ProteinId)

# TPDv1 - TPDv2
tg_list <- list(TPDv1 = tgv1, TPDv2 = tgv2)
pr_list <- list(TPDv1 = prv1, TPDv2 = prv2)
pep_list <- list(TPDv1 = pev1, TPDv2 = pev2)
prot_list <- list(TPDv1 = ppv1, TPDv2 = ppv2)


# VennDiagram
p1 <- VennDiagram::venn.diagram(x = tg_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pep_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white',
                                main=stringr::str_glue("Transition groups ({length(unique(unlist(tg_list)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
p2 <- VennDiagram::venn.diagram(x = pr_list,
                                resolution = 300,
                                alpha=rep(0.95, length(prot_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white',
                                main=stringr::str_glue("Precursors ({length(unique(unlist(pr_list)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
p3 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pep_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white',
                                main=stringr::str_glue("Peptides ({length(unique(unlist(pep_list)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
p4 <- VennDiagram::venn.diagram(x = prot_list,
                                resolution = 300,
                                alpha=rep(0.95, length(prot_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white',
                                main=stringr::str_glue("Proteotypic proteins ({length(unique(unlist(prot_list)))} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
graphics.off()
pdf('TPDv1_TPDv2_VennDiagram_20240112.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(p1)
grid::grid.newpage(); grid::grid.draw(p2)
grid::grid.newpage(); grid::grid.draw(p3)
grid::grid.newpage(); grid::grid.draw(p4)
graphics.off()

# IM-RT plot
df_imrt <- df %>%
  dplyr::select(ModifiedPeptideSequence, PrecursorCharge, AverageExperimentalRetentionTime, PrecursorIonMobility) %>%
  dplyr::distinct() %>% 
  rename(RT = AverageExperimentalRetentionTime, IM = PrecursorIonMobility) %>% 
  mutate(RT = RT / 60) %>% 
  rename(`RT (min)` = RT)

pdf('TPDlibv2_IM_RT.pdf', width = 10, height = 10)
plot(df_imrt$`RT (min)`, df_imrt$IM)
graphics.off()

pimrt1 <- ggplot(df_imrt, aes(x = `RT (min)`, y = IM))+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 300) +
  scale_fill_continuous(low = "white", high = 'dodgerblue')+ # dodgerblue
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  labs(title = 'Library modified precursors RT-IM plot', x = "RT (min)", y = "Ion mobility (V·s·cm-2)")+
  theme_bw()+
  theme(
    legend.position = 'none',
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))+
  theme(panel.grid =element_blank())+
  theme(axis.text = element_text(size = 15,color = "black"))+
  theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
  theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
  theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))

ggsave('TPDlibv2_modified_precursors_RT_IM_plot_v1.pdf', pimrt1, width = 6, height = 6)



pimrt2 <- ggplot(df_imrt, aes(x = `RT (min)`, y = IM))+
  stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, n = 300) +
  scale_fill_continuous(low = "white", high = 'dodgerblue')+ # dodgerblue
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1), add = c(0, 0)))+
  geom_point(shape = 20, size = 1, color = "dodgerblue4", alpha = 0.1)+
  labs(title = 'Library modified precursors RT-IM plot', x = "RT (min)", y = "Ion mobility (V·s·cm-2)")+
  theme_bw()+
  theme(
    legend.position = 'none',
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))+
  theme(panel.grid =element_blank())+
  theme(axis.text = element_text(size = 15,color = "black"))+
  theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
  theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
  theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))

ggsave('TPDlibv2_modified_precursors_RT_IM_plot_v2.pdf', pimrt2, width = 6, height = 6)

pdf('TPDlibv2_modified_precursors_RT_IM_density_curve.pdf', width = 6, height = 3)
plot(density(df_imrt$`RT (min)`))
plot(density(df_imrt$IM))
graphics.off()

