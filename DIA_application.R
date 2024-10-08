library(magrittr)
library(tidyverse)


# functions ---------------------------------------------------------------
removeRowsAllNa <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
ulhen_class<-function(tissue_df, fct=3){
  # tissue_df should be structured as data.frame(SampleType_vector, ProteinMatrix_matrix)
  j=1
  clas<-apply(tissue_df[,-1], 2, function(num){
    t1<-as.data.frame(matrix(ncol=3,nrow = nrow(tissue_df)))
    colnames(t1)<-c("tissue_type","abundance","classification")
    t1$tissue_type<-tissue_df[, 1]
    t1$abundance<-num
    t1<-t1[order(t1$abundance,decreasing = T),]
    # class 1(not detected) and 2(tissue enriched)
    for (i in 1:nrow(tissue_df)){
      if (t1$abundance[i]<=min(tissue_df[,-1],na.rm = T)) {t1[i,3]<-"not detected"
      } else if(t1$abundance[i]>=(fct*max(t1$abundance[-i]))) {t1[i,3]<-"tissue enriched"
      } else {t1[i,3]<-NA}
    }
    c12<-t1[!is.na(t1$classification),]
    
    # class 3(group enriched)
    max7<-t1[is.na(t1$classification),]
    if (nrow(max7)>7){
      high=max7$abundance[1]
      for (n in 1:6){
        high=high+max7$abundance[n+1]
        nextp=max7$abundance[n+2]
        # fill=max7$classification[1:(n+1)]
        if (high/(n+1) >= fct*nextp)  {max7$classification[1:(n+1)]<-"group enriched"}
      } 
    } else if (nrow(max7)>=2 & nrow(max7)<=7) {max7$classification<-"group enriched"}
    
    #class 4 (expressed in all)
    if (sum(as.logical(which(t1$abundance<=min(tissue_df[,-c(1,ncol(tissue_df))],na.rm = T))))==0){
      max7$classification[is.na(max7$classification)]<-"Expressed in all tissues"}
    c34<-na.omit(max7)
    
    #class 5 (tissue enriched) and class 6(mixed)
    max7_na<-max7[is.na(max7$classification),]
    if (nrow(max7_na)>=1){
      for (i in 1:nrow(max7_na)){
        if(max7_na$abundance[i]>= (fct*mean(t1$abundance))) {max7_na$classification[i]<-"tissue enhanced"
        } else {max7_na$classification[i]<-"mixed"}
      }
    }
    com_max<-do.call(rbind,list(c12,c34,max7_na))
    j <<- j+1
    colnames(com_max)[3]<-colnames(tissue_df)[j]
    if (j %% 1000 ==0){
      cat(paste("Processed",j-1,"/",ncol(tissue_df)-1,sep = " "),"\n")}
    
    com_max[,c(1,3)]
    
    # colnames(class[[j]])[2]<-colnames(tissue_df)[j+1]
  })
  
  print("Finished classification and start merging classification matrix")
  
  clas.all<-clas%>% reduce_right(full_join, by = 'tissue_type')
  return(clas.all)
}

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


# 5.Specificity -----------------------------------------------------------
## 5.1 Classification using Ulhen's method -----
dfimp <- df[, apply(df, 2, function(y) sum(is.na(y)) / length(y)) < 0.5]
dfimp[is.na(dfimp)] <- min(mat_pro, na.rm = T) + log2(0.8)
dfimp_med <- dfimp %>% group_by(Histopathology_type) %>% summarise_if(is.numeric, function(y) median(2^y)) %>% 
  as.data.frame()

# cost ~1min using Intel Core i7-9700K
sample_classification<-list()
sample_classification[[1]]<-ulhen_class(dfimp_med[,1:50], fct = 5)
for (i in 2:floor(ncol(dfimp_med)/50)){
  sample_classification[[i]]<-ulhen_class(dfimp_med[,c(1,(50*(i-1)+1):(50*i))], fct = 5)
  cat(paste("Processed",i*50,"/",ncol(dfimp_med) ,sep = " "),"\n")
}
sample_classification[[ceiling(ncol(dfimp_med)/50)]]<-ulhen_class(dfimp_med[,c(1,(50*floor(ncol(dfimp_med)/50)+1):ncol(dfimp_med))])
class_all<-sample_classification%>% reduce_right(full_join, by = "tissue_type")
colnames(class_all)[1] <- 'Histopathology_type'
class_all %<>% arrange(Histopathology_type)

class_all_long <- class_all %>% pivot_longer(cols = -Histopathology_type, names_to = 'Protein', values_to = 'Classification')
tis_enrich <- class_all_long %>% filter(Classification == 'tissue enriched') %>% 
  arrange(Histopathology_type)
tis_enrich$`First/Second` <- dfimp_med %>% select(all_of(tis_enrich$Protein)) %>% 
  apply(2, function(y){
    y_sorted <- sort(y, decreasing = T)
    y_sorted[1] / y_sorted[2]
  })
tis_enrich %<>%
  inner_join(df_pro %>%
               select(Protein.Group, Genes) %>%
               rename(Protein = Protein.Group, Gene = Genes), .) %>% 
  arrange(Histopathology_type, desc(`First/Second`))
tis_enrich %<>% inner_join(df_anova)
# pro_top5 <- tis_enrich %>% group_by(Histopathology_type) %>% slice(1:5) %>% pull(Protein)
pro_top5 <- tis_enrich %>% 
  inner_join(df_anova) %>%
  filter(adj.P.anova < 0.05) %>% 
  group_by(Histopathology_type) %>%
  arrange(adj.P.anova) %>% 
  slice(1:5) %>% pull(Protein)
tis_enrich$IsTop5 <- tis_enrich$Protein %in% pro_top5

list(All = class_all_long, TissueEnriched = tis_enrich, ClassMatrix = class_all) %>%
  rio::export('TPDlibV2_DIA_classification.xlsx')

class_all_long %>% count(Classification)
# Classification               n
# Expressed in all tissues 24547
# group enriched           11354
# mixed                     3626
# not detected              5423
# tissue enriched             59



XX2 <- X[class_all_long %>% filter(Classification == 'tissue enriched') %>% arrange(Histopathology_type) %>% pull(Protein) %>%
           intersect(pro_top5), ]
rownames(XX2) <- df_pro %>% set_rownames(.$Protein.Group) %>%
  .[rownames(XX2), ] %>% # here in the UNIPROT ID input
  unite(Label, Protein.Group, Genes) %>% pull(Label)

anno_col <- dfimp %>% column_to_rownames('MS_file_name') %>% select(Histopathology_type)
anno_row <- class_all_long %>%
  filter(Classification == 'tissue enriched') %>%
  arrange(Histopathology_type) %>%
  left_join(protinfo, by = c(Protein = 'Protein.Group')) %>%
  filter(Protein %in% pro_top5) %>% 
  column_to_rownames('Label') %>% select(Histopathology_type)
mycolors <- RColorBrewer::brewer.pal(6,"PiYG")
bk <- unique(c(seq(-1.5, 1.5, length=50)))
pheatmap(XX2,
         scale = 'none',
         cluster_rows = F, cluster_cols = F,
         clustering_distance_rows = 'euclidean', clustering_method = 'ward.D2',
         annotation_col = anno_col,
         annotation_row = anno_row,
         annotation_colors = list(Histopathology_type = type_colors),
         show_rownames = T, show_colnames = F, border_color = NA,
         color = colorRampPalette(rev(mycolors))(50),
         breaks = bk,
         gaps_row = cumsum(table(anno_row$Histopathology_type)) %>% .[. != 0],
         gaps_col = cumsum(table(anno_col$Histopathology_type)),
         fontsize = 9, main = '59 histopathology-enriched proteins',
         filename = 'TPDlibV2_DIA_classification_heatmap_top5.pdf',  width = 10, height = 9
)
