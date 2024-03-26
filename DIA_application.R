library(magrittr)
library(tidyverse)


# functions ---------------------------------------------------------------
removeRowsAllNa <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

# 1.Data prepare ----------------------------------------------------------
ordered_types <- c('N', 'MNG', 'FA', 'FTC', 'PTC', 'PDTC', 'ATC', 'MTC', 'LN')
type_colors <- c(N = '#969696',
                 MNG = '#8ECFC9',
                 FA = '#DB95B4',
                 FTC = '#82B0D2',
                 PTC = '#936E49',
                 PDTC = '#E8C76C',
                 ATC = '#E59E63',
                 MTC = '#7E71B2',
                 LN = '#6CA679')

info <- rio::import('TTTD_480_tpdlibV2.xlsx', sheet = 3)
info %<>% mutate(Histopathology_type = factor(Histopathology_type, levels = ordered_types)) %>% 
  arrange(Histopathology_type)

df_pro <- rio::import('report.pg_matrix.tsv') %>% filter(str_detect(Protein.Group, 'iRT', negate = T))
df_pr <- rio::import('//172.16.13.136/tttd/3MS_data/TPDlibV2_figure4_rawdata/TPDlibV2_45files_480lib_result/report.pr_matrix.tsv')

mat_pro <- df_pro %>% column_to_rownames('Protein.Group') %>% select(-(Protein.Ids:First.Protein.Description))
colnames(mat_pro) %<>% str_remove('^.+\\\\') %>% str_remove('\\.raw$')
mat_pro <- mat_pro[, info$MS_file_name]
mat_pro %<>% log2()
df <- mat_pro %>% t() %>% as.data.frame() %>% rownames_to_column('MS_file_name') %>% inner_join(info, .)

protinfo <- df_pro %>%
  select(Protein.Group:First.Protein.Description) %>%
  mutate(Label = str_c(Protein.Group, '_', Genes)) %>% set_rownames(.$Protein.Group)


# 2.protein/peptide number --------------------------------------------------------
## 2.1 protein -------
pronum <- plyr::ddply(df, 'Histopathology_type', function(dfsub){
  X <- dfsub %>% select(-MS_file_name) %>% removeColsAllNa()
  prot <- ncol(X)
  prot0.5 <- sum(apply(X, 2, function(y) sum(is.na(y)) / length(y)) < 0.5)
  data.frame(Protein = prot, Protein0.5 = prot0.5)
})
pronum$Histopathology_type %<>% factor(levels = ordered_types)
pronum %<>% arrange(Histopathology_type)

p <- ggplot(pronum, aes(x = Histopathology_type, y = Protein))+
  geom_col(color = '#000000', fill = '#CCCCCC')+
  geom_text(aes(label = Protein), vjust = -0.5)+
  labs(x = 'Histopathology', y = '# Proteins')+
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_protein_number_barplot.pdf', p, width = 6, height = 6)

pronum$y2 <- pronum$Protein0.5 # proteins; missing ratio < 50%
pronum$y1 <- pronum$Protein - pronum$Protein0.5
p <- pronum %>% pivot_longer(cols = c('y1', 'y2'), names_to = 'Type') %>% 
  ggplot()+
  geom_col(aes(x = Histopathology_type, y = value, fill = Type), color = '#000000')+
  geom_text(data = pronum, 
            aes(x = Histopathology_type, y = y2, label = y2),
            vjust = -0.5)+
  geom_text(data = pronum, 
            aes(x = Histopathology_type, y = Protein, label = Protein
                ),
            vjust = -0.5)+
  labs(x = 'Histopathology', y = '# Proteins')+
  scale_fill_manual(values = c(y1 = '#CCCCCC', y2 = '#AAAAAA'),
                    labels = c(y1 = '', y2 = 'NA ratio<0.5')) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_protein_number_barplot_v2.pdf', p, width = 7, height = 6)




set.seed(10)
p <- ggplot(pronum, aes(x = Histopathology_type, y = Protein, fill = Histopathology_type))+
  geom_col(color = '#000000', alpha = 0.5, width = 0.6)+
  geom_text(aes(label = Protein), vjust = -0.5)+
  geom_jitter(data = data.frame(Histopathology_type = df$Histopathology_type,
                                Protein = apply(df[, -(1:2)], 1, function(x) sum(!is.na(x)))),
              aes(x = Histopathology_type, y = Protein, color = Histopathology_type),
              size = 3, show.legend = F) +
  labs(x = 'Histopathology', y = '# Proteins')+
  scale_fill_manual(values = type_colors)+
  scale_color_manual(values = type_colors)+
  theme_bw()+
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_protein_number_barplot_v3.pdf', p, width = 10, height = 4.5)


## 2.2 peptide ------
#precursor matrix to peptide matrix
df_pep <- df_pr %>%
  filter(str_detect(Protein.Group, 'iRT', negate = T)) %>% 
  select(-Proteotypic, -Precursor.Charge) %>% 
  group_by(Stripped.Sequence) %>% summarise_if(is.numeric, mean, na.rm = T)

mat_pep <- df_pep %>% column_to_rownames('Stripped.Sequence')
colnames(mat_pep) %<>% str_remove('^.+\\\\') %>% str_remove('\\.raw$')

pepnum <- mat_pep %>% t() %>% as.data.frame() %>% rownames_to_column('MS_file_name') %>% inner_join(info, .) %>%
  plyr::ddply(., 'Histopathology_type', function(dfsub){
    X <- dfsub %>% select(-MS_file_name) %>% removeColsAllNa()
    pep <- ncol(X)
    pep0.5 <- sum(apply(X, 2, function(y) sum(is.na(y)) / length(y)) < 0.5)
    data.frame(Peptide = pep, Peptide0.5 = pep0.5)
  })
pepnum$Histopathology_type %<>% factor(levels = ordered_types)
pepnum %<>% arrange(Histopathology_type)

p <- ggplot(pepnum, aes(x = Histopathology_type, y = Peptide))+
  geom_col(color = '#000000', fill = '#CCCCCC')+
  geom_text(aes(label = Peptide), vjust = -0.5)+
  labs(x = 'Histopathology', y = '# Peptides')+
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_Peptide_number_barplot.pdf', p, width = 6, height = 6)

pepnum$y2 <- pepnum$Peptide0.5 # Peptides; missing ratio < 50%
pepnum$y1 <- pepnum$Peptide - pepnum$Peptide0.5
p <- pepnum %>% pivot_longer(cols = c('y1', 'y2'), names_to = 'Type') %>% 
  ggplot()+
  geom_col(aes(x = Histopathology_type, y = value, fill = Type), color = '#000000')+
  geom_text(data = pepnum, 
            aes(x = Histopathology_type, y = y2, label = y2),
            vjust = -0.5)+
  geom_text(data = pepnum, 
            aes(x = Histopathology_type, y = Peptide, label = Peptide
            ),
            vjust = -0.5)+
  labs(x = 'Histopathology', y = '# Peptides')+
  scale_fill_manual(values = c(y1 = '#CCCCCC', y2 = '#AAAAAA'),
                    labels = c(y1 = '', y2 = 'NA ratio<0.5')) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_Peptide_number_barplot_v2.pdf', p, width = 7, height = 6)



df_tmp <- mat_pep %>% t() %>% as.data.frame() %>% rownames_to_column('MS_file_name') %>% inner_join(info, .)
set.seed(10)
p <- ggplot(pepnum, aes(x = Histopathology_type, y = Peptide, fill = Histopathology_type))+
  geom_col(color = '#000000', alpha = 0.5, width = 0.6)+
  geom_text(aes(label = Peptide), vjust = -0.5)+
  geom_jitter(data = data.frame(Histopathology_type = df_tmp$Histopathology_type,
                                Peptide = apply(df_tmp[, -(1:2)], 1, function(x) sum(!is.na(x)))),
              aes(x = Histopathology_type, y = Peptide, color = Histopathology_type),
              size = 3, show.legend = F) +
  labs(x = 'Histopathology', y = '# Peptides')+
  scale_fill_manual(values = type_colors)+
  scale_color_manual(values = type_colors)+
  theme_bw()+
  theme(text = element_text(size = 15))
ggsave('TPDlibV2_peptide_number_barplot_v3.pdf', p, width = 10, height = 4.5)




# 3.Correlation -----------------------------------------------------------
library(corrplot)

identical(info$MS_file_name, colnames(mat_pro)) # TRUE
info %<>%
  with_groups(Histopathology_type, mutate,
              label = str_c(Histopathology_type, 1:length(Histopathology_type)))


cors <- cor(mat_pro %>% set_colnames(info$label), use = "pairwise.complete.obs", method = "pearson")
uptri <- cors[upper.tri(cors)]
median(uptri) # 0.7534793
min(uptri) # 0.5231689

# hist(uptri)
pdf("TPDlibV2_DIA_Pearson_correlation.pdf", width = 15, height = 15)
corrplot(cors, method = "square", type = "full", order ="original", cl.cex = 1.4, tl.pos = "lt", tl.col = "black", tl.cex = 1.5, tl.srt = 60,
         col = COL2('PuOr', 10),
         is.corr = F, col.lim = c(0.5, 1),
         )
graphics.off()

dfcor <- data.frame(Pearson = uptri)
p <- ggplot(dfcor) + 
  geom_density(aes(x = Pearson)) +
  geom_vline(aes(xintercept = median(uptri)), color = '#000000', linetype = "dashed", size = 1)+
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("TPDlibV2_DIA_Pearson_correlation_density.pdf", p, width = 4.5, height = 4.5)


# 4.ANOVA -----------------------------------------------------------------
## 4.1 ANOVA -----------
dfimp <- df

#missing
df_miss <- data.frame(
  `Missing ratio (%)` = seq(0.01, 0.99, 0.02) * 100,
  `# Proteins` = sapply(seq(0.01, 0.99, 0.02), function(i){
    tmp <- mat_pro[apply(mat_pro, 1, function(x) sum(is.na(x)) / length(x)) < i, ]
    nrow(tmp)
  }),
  check.names = F
)
p <- ggplot(df_miss)+
  geom_line(aes(x = `Missing ratio (%)`, y = `# Proteins`), size = 2)+
  theme_bw()
ggsave('TPDlibV2_DIA_missing.pdf', p, width = 5, height = 5)
rio::export(df_miss, 'TPDlibV2_DIA_missing.xlsx')


dfimp <- df[, apply(df, 2, function(y) sum(is.na(y)) / length(y)) < 0.5]
dfimp[is.na(dfimp)] <- min(mat_pro, na.rm = T) + log2(0.8)

aov_ls <- list()
for(i in 3:ncol(dfimp)){
  cat(i, '...\r')
  res_aov <- aov(dfimp[[i]] ~ dfimp$Histopathology_type)
  aov_ls[[i-2]] <- res_aov
}

p.anova <- sapply(aov_ls, function(res_aov){
  summary(res_aov)[[1]]$`Pr(>F)`[1]
})
df_anova <- data.frame(Protein = colnames(dfimp)[-(1:2)],
                       P.anova = p.anova)
df_anova$adj.P.anova <- p.adjust(df_anova$P.anova, 'BH')
rio::export(df_anova, 'TPDlibV2_DIA_ANOVA_50NA.xlsx')

## 4.2 heatmap -----------
library(pheatmap)

X <- dfimp %>% select(-Histopathology_type) %>% column_to_rownames('MS_file_name') %>%
  scale() %>% t()

sum(df_anova$adj.P.anova < 0.001) # 1246
pheatmap(X[df_anova$Protein[df_anova$adj.P.anova < 0.001], ],
         scale = 'none',
         cluster_rows = T, cluster_cols = T,
         # clustering_distance_rows = 'euclidean', clustering_method = 'ward.D2',
         annotation_col = dfimp %>% column_to_rownames('MS_file_name') %>% select(Histopathology_type),
         # annotation_colors = anno_color,
         show_rownames = F, show_colnames = F,
         fontsize = 9, main = '1246 proteins ANOVA adj.P < 0.001',
         filename = 'TPDlibV2_DIA_ANOVA_heatmap_50NA.pdf',  width = 10, height = 4,
)


sum(df_anova$adj.P.anova < 10e-8) # 59
XX <- X[df_anova$Protein[df_anova$adj.P.anova < 10e-8], ]
rownames(XX) <- df_pro %>% set_rownames(.$Protein.Group) %>%
  .[rownames(XX), ] %>% # here in the UNIPROT ID input
  unite(Label, Protein.Group, Genes) %>% pull(Label)
pheatmap(XX,
         scale = 'none',
         cluster_rows = T, cluster_cols = T,
         clustering_distance_rows = 'euclidean', clustering_method = 'ward.D2',
         annotation_col = dfimp %>% column_to_rownames('MS_file_name') %>% select(Histopathology_type),
         # annotation_colors = anno_color,
         show_rownames = T, show_colnames = F, border_color = NA,
         fontsize = 9, main = '59 proteins ANOVA adj.P < 10e-8',
         filename = 'TPDlibV2_DIA_ANOVA_heatmap_50NA_adjP.pdf',  width = 10, height = 8
)

## 4.3 specific proteins ------
library(ggstatsplot)
paletteer::palettes_d_names

pdf('TPDlibV2_DIA_ANOVA_targets.pdf', width = 5, height = 4)
# MTC high-expressed proteins
ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = A6NCE7,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['A6NCE7', 'Label']
)

ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = Q9UI12,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['Q9UI12', 'Label']
)

ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = O43488,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['O43488', 'Label']
)

ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = P31150,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['P31150', 'Label']
)

ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = P61421,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['P61421', 'Label']
)

ggbetweenstats(
  data = dfimp,
  x = Histopathology_type,
  y = Q8IXB1,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = T, # ANOVA or Welch ANOVA
  plot.type = "box", p.adjust.method = 'BH',
  pairwise.comparisons = F,
  pairwise.display = "significant",
  centrality.plotting = F,
  bf.message = F,
  palette = "royal",
  package = "basetheme",
  xlab = 'Histopathology',
  ylab = protinfo['Q8IXB1', 'Label']
)

graphics.off()


