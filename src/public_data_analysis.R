library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

# dataset Wu ------------------------------------
# Wu, T.D., Madireddi, S., de Almeida, P.E. et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. Nature 579, 274–278 (2020).

df_Wu <- read_csv('data/adult data from public/blood mCD8T from dataset_Wu.csv.gz')
panel_Wu <- colnames(df_Wu)[1:3734]

volcano_tmp <- function(df){
  tmp <- df %>% group_by(expansion) %>% summarise(across(all_of(panel_Wu), mean))
  od <- tmp$expansion
  tmp <- t(tmp[, panel_Wu])
  colnames(tmp) <- od
  resFC <- tmp %>% as.data.frame() %>% mutate(logFC = expanded - non_expanded)
  resFC$gene <- rownames(resFC)
  tmp <- df %>% tidyr::pivot_longer(values_to = 'level', names_to = 'gene', cols = all_of(panel_Wu))
  resP <- tmp %>% 
    group_by(gene) %>% 
    t_test(level ~ expansion) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  res <- left_join(resFC, resP, by = 'gene')
  top_p_genes <- (res %>% top_n(-50, p.adj))$gene
  interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
  lab = unique(c(top_p_genes, interested_genes))
  
  EnhancedVolcano(res, x='logFC', y='p.adj', lab = res$gene, 
                  pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE, arrowheads = F, 
                  selectLab = lab,
                  xlim=c(-1,1),
                  col= c('grey', 'black', 'red', 'orange'),
                  subtitle = paste0('expanded(', table(df$expansion)['expanded'],  
                                    ') vs non-expanded(', table(df$expansion)['non_expanded'], ') cells')) + 
    theme_minimal()
}

volcano_tmp(df_Wu) + ggtitle('memory CD8 T cells in Blood (all)')
ggsave('figures/DE/public dataset/volcano mCD8T in blood (data_Wu) expanded vs non-expanded (all).pdf', width = 10, height = 10)

df_Wu %>% filter(tumor_type == 'Renal') %>% volcano_tmp() + ggtitle('memory CD8 T cells in Blood (Renal)')
ggsave('figures/DE/public dataset/volcano mCD8T in blood (data_Wu) expanded vs non-expanded (Renal).pdf', width = 10, height = 10)

df_Wu %>% filter(tumor_type == 'Lung') %>% volcano_tmp() + ggtitle('memory CD8 T cells in Blood (Lung)')
ggsave('figures/DE/public dataset/volcano mCD8T in blood (data_Wu) expanded vs non-expanded (Lung).pdf', width = 10, height = 10)

# dataset Guo ------------------------------------
# Guo, X., Zhang, Y., Zheng, L. et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing. Nat Med 24, 978–985 (2018). 

obj_Guo <- LoadH5Seurat('data/adult data from public/dataset_Guo.h5Seurat')

obj_Guo_sub <- subset(obj_Guo, subset = Clone.Status != 'NA' &
                        sampleType %in% c('PTC', 'PTH', 'PTR', 'PTY') &
                        majorCluster %in% c('CD8_C2-CD28', 'CD8_C3-CX3CR1', 'CD8_C4-GZMK', 'CD8_C5-ZNF683', 'CD8_C6-LAYN', 'CD8_C7-SLC4A10')# peripheral blood memory CD8T cells
)
obj_Guo_sub[[]]
Idents(obj_Guo_sub) <- 'Clone.Status'
res <- FindMarkers(obj_Guo_sub, ident.1 = 'Clonal', ident.2 = 'NoClonal',
                   logfc.threshold = 0, min.pct = 0.5)
top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
lab = unique(c(top_p_genes, interested_genes))

EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = lab,
                title = 'memory CD8 T cells in Blood (Lung)',
                # xlim=c(-3,3),
                col= c('grey', 'black', 'red', 'orange'),
                arrowheads=F,
                subtitle = paste0('expanded(', table(obj_sub[[]]$Clone.Status)['Clonal'],  
                                  ') vs non-expanded(', table(obj_sub[[]]$Clone.Status)['NoClonal'], ') cells')) + 
  theme_minimal()
ggsave('figures/DE/public dataset/volcano mCD8T in blood (data_Guo) expanded vs non-expanded (Lung).pdf', width = 10, height = 10)

# dataset Zhang ---------------------------------
# Zhang, L., Yu, X., Zheng, L. et al. Lineage tracking reveals dynamic relationships of T cells in colorectal cancer. Nature 564, 268–272 (2018).

obj_Zhang <- LoadH5Seurat('data/adult data from public/dataset_Zhang.h5Seurat')

obj_Zhang_sub <- subset(obj_Zhang, subset = Clonal.status != 'NA' &
                          sampleType %in% c('PTC', 'PTH', 'PTR', 'PTY', 'PP7') & # peripheral blood cells
                          majorCluster %in% c('CD8_C02-GPR183', 'CD8_C03-CX3CR1', 'CD8_C04-GZMK', 'CD8_C05-CD6',
                                              'CD8_C06-CD160', 'CD8_C07-LAYN', 'CD8_C08-SLC4A10') # exclude CD8 T cells with LEF1 expression
)
obj_Zhang_sub[[]]
Idents(obj_Zhang_sub) <- 'Clonal.status'
res <- FindMarkers(obj_Zhang_sub, ident.1 = 'Clonal', ident.2 = 'NoClonal',
                   logfc.threshold = 0, min.pct = 0.5)
top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
lab = unique(c(top_p_genes, interested_genes))

EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = lab,
                title = 'memory CD8 T cells in Blood (Colorectal)',
                xlim=c(-3,3),
                col= c('grey', 'black', 'red', 'orange'),
                arrowheads=F,
                subtitle = paste0('expanded(', table(obj_sub[[]]$Clonal.status)['Clonal'],  
                                  ') vs non-expanded(', table(obj_sub[[]]$Clonal.status)['NoClonal'], ') cells')) + 
  theme_minimal()
ggsave('figures/DE/public dataset/volcano mCD8T in blood (data_Zhang) expanded vs non-expanded (Colorectal).pdf', width = 10, height = 10)

# compare the fraction of expanded clones ----------------------------------

df <- rbind(df_Wu %>% group_by(patient, expansion, tumor_type) %>% tally() %>% group_by(patient) %>% mutate(sum=sum(n), fraction = n/sum),
            obj_Guo_sub[[]] %>% group_by(Patient, Clone.Status) %>% tally() %>% group_by(Patient) %>% mutate(sum=sum(n), fraction = n/sum) %>% rename(patient = Patient, expansion = Clone.Status) %>% mutate(tumor_type = 'NSCLC'),
            obj_Zhang_sub[[]] %>% group_by(Patient_ID, Clonal.status) %>% tally() %>% group_by(Patient_ID) %>% mutate(sum=sum(n), fraction = n/sum) %>% rename(patient = Patient_ID, expansion = Clonal.status) %>% mutate(tumor_type = 'Colorectal')
) %>% filter(expansion %in% c('expanded', 'Clonal'))

clone_id_map <- read_csv('data/clone_id_map.csv.gz')
data_scTCR <- read_csv('data/scTCR_data_merge_paired_chain.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz')
virus_specific_list <- read_csv(file = 'data/virus_specific_clone_id.csv')
patient_info <- read_csv(file='data/ISAC_TCR_sample_info.csv')

obj <- LoadH5Seurat('data/seurat_results_update.h5Seurat')

expand_id <- (clone_exp %>% filter(clone_count >= 2,
                                   # grepl('(_v1$)|(_1$)', Sample_Name), # baseline only
                                   !clone_id %in% virus_specific_list$clone_id) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
expand_sc <- data_scTCR %>% filter(clone_id %in% expand_id) %>% 
  tibble::add_column(expansion = 'expanded') %>%
  select(unique_index, expansion)

obj_sub <- subset(obj, subset = 
                    cell_type == 'CD8T' & 
                    !(clone_id %in% virus_specific_list$clone_id) &
                    TCR_Paired_Chains)
tmp <- obj_sub[[]] %>% left_join(expand_sc, by = 'unique_index') %>%
  tidyr::replace_na(list(expansion='non-expanded'))
tmp2 <- tmp$expansion
names(tmp2) <- tmp$unique_index
obj_sub$expansion <- tmp2

df <- rbind(df %>% mutate(group = 'adults'),
      obj_sub[[]] %>% group_by(patient_id, expansion, tumor_type) %>% tally() %>% group_by(patient_id) %>% mutate(sum=sum(n), fraction = n/sum, group='isac') %>% rename(patient = patient_id) %>% filter(expansion == 'expanded')
)

write_csv(df, 'data/adult_isac_expansion_ratio.csv')

df <- df %>% filter(tumor_type != 'Other Neuroblastoma')

ggplot(df, aes(x=factor(tumor_type, levels=unique(df$tumor_type)), y=fraction, color=group)) + geom_jitter() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('tumor_type') + ylab('Fraction of expanded cells in blood memory CD8T cells')

ggsave('figures/DE/public dataset/Fraction of expanded cells.pdf', width = 5, height = 5)

library(rstatix)
df %>% ungroup() %>% t_test(fraction ~ group)

# update grouping of 02 and 46, used to be in the 'other neuroblastoma group but we have new mutation results of them'
df <- read_csv('data/adult_isac_expansion_ratio.csv')
df <- df %>% mutate(tumor_type = case_when(
  patient == 'ISAC02' ~ 'High mutation Neuroblastoma',
  patient == 'ISAC46' ~ 'Low mutation Neuroblastoma',
  TRUE ~ tumor_type
)) %>% filter(tumor_type != 'Other Neuroblastoma')

# ggplot(df, aes(x=factor(tumor_type, levels=unique(df$tumor_type)), y=fraction, color=group)) + 
#   geom_bar(stat = 'iden')
#   geom_jitter()
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   xlab('tumor_type') + ylab('Fraction of expanded cells in blood memory CD8T cells')
# ggsave('figures/DE/update 02 and 46/Fraction of expanded cells.pdf', width = 5, height = 5)

ggplot(df, aes(x=factor(tumor_type, levels=unique(tumor_type)), y=fraction, color=group)) + 
  geom_bar(stat = 'summary', fun=mean, aes(fill=group)) +
  geom_jitter(color='black', width=0.2) +
  theme_bw()
ggsave('figures/DE/update 02 and 46/Fraction of expanded cells_barplot.pdf', width = 5, height = 4)

library(rstatix)
df %>% ungroup() %>% t_test(fraction ~ group)

# compare chi-square
tmp <- df %>% mutate(expanded = n,
              non_expanded = sum-n) %>%
  select(expanded, non_expanded, group) %>%
  group_by(group) %>%
  mutate(expanded_sum = sum(expanded),
         non_expanded_sum = sum(non_expanded)) %>%
  select(expanded_sum, non_expanded_sum) %>%
  distinct()

xtab <- t(tmp[,2:3])
dimnames(xtab) <- list(
  Expanded = c("Expanded", "Non-expanded"),
  Group = c("Adults", "Pediatric")
)

chisq_test(xtab)
fisher_test(xtab)

