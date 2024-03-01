library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(EnhancedVolcano)
library(Seurat)
library(SeuratDisk)

`%nin%` = Negate(`%in%`)

# load data
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


volcano_from_obj <- function(oj, exp_sc, title){
  
  tmp <- oj[[]] %>% left_join(exp_sc, by = 'unique_index') %>%
    tidyr::replace_na(list(expansion='non-expanded'))
  tmp2 <- tmp$expansion
  names(tmp2) <- tmp$unique_index
  oj$expansion <- tmp2
  Idents(oj) <- 'expansion'
  res <- FindMarkers(oj, ident.1 = 'expanded', ident.2 = 'non-expanded',
                     logfc.threshold = 0, min.pct = 0)
  top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
  interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1')
  lab = unique(c(top_p_genes, interested_genes))
  
  EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                  pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                  selectLab = lab,
                  title = title,
                  xlim=c(-1.5,1.5),
                  col= c('grey', 'black', 'orange', 'red'),
                  arrowheads=F,
                  subtitle = paste0('expanded(', table(tmp2)['expanded'],  
                                    ') vs non-expanded(', table(tmp2)['non-expanded'], ') cells')) + 
    theme_minimal()
}

cur_subpop = 'CD8T'
for (cur_tumor_type in unique(obj[[]]$tumor_type)){
  # cur_tumor_type = 'Wilms'
  obj_sub <- subset(obj, subset = 
                      cell_type == cur_subpop & 
                      tumor_type == cur_tumor_type &
                      clone_id %nin% virus_specific_list$clone_id &
                      TCR_Paired_Chains)
  volcano_from_obj(obj_sub, expand_sc, title = paste0(cur_subpop, ' of ', cur_tumor_type))
  ggsave(paste0('figures/DE/update 02 and 46/', cur_subpop, '_', cur_tumor_type ,'_expansion.pdf'), width = 15, height = 15)
}

obj_sub <- subset(obj, subset = 
                    cell_type == cur_subpop & 
                    clone_id %nin% virus_specific_list$clone_id &
                    TCR_Paired_Chains)
volcano_from_obj(obj_sub, expand_sc, title = paste0(cur_subpop, ' of all samples'))
ggsave(paste0('figures/DE/CD8T_all_samples_expansion.pdf'), width = 15, height = 15)

# in the first 8 isac patients, DE of expanded vs non-expanded if virus-specific excluded.
obj_sub <- subset(obj, subset = 
                    patient_id %in% c('ISAC02', 'ISAC31', 'ISAC35', 'ISAC62', 'ISAC77', 'ISAC81', 'ISAC99', 'ISAC125') &
                    cell_type == cur_subpop & 
                    clone_id %nin% virus_specific_list$clone_id &
                    TCR_Paired_Chains)
volcano_from_obj(obj_sub, expand_sc, title = paste0(cur_subpop, ' of all samples'))
ggsave(paste0('figures/DE/CD8T_all_samples_expansion_first8_virus_specific_excluded.pdf'), width = 15, height = 15)

# in the first 8 isac patients, DE of expanded vs non-expanded if virus-specific not excluded.
obj_sub <- subset(obj, subset = 
                    patient_id %in% c('ISAC02', 'ISAC31', 'ISAC35', 'ISAC62', 'ISAC77', 'ISAC81', 'ISAC99', 'ISAC125') &
                    cell_type == cur_subpop & 
                    # clone_id %nin% virus_specific_list$clone_id &
                    TCR_Paired_Chains)
expand_id_wvs <- (clone_exp %>% filter(clone_count >= 2,
                                       # !clone_id %in% virus_specific_list$clone_id # do not exclude virus specific clones
) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
expand_sc_wvs <- data_scTCR %>% filter(clone_id %in% expand_id_wvs) %>% 
  tibble::add_column(expansion = 'expanded') %>%
  select(unique_index, expansion)

volcano_from_obj(obj_sub, expand_sc_wvs, title = paste0(cur_subpop, ' of all samples'))
ggsave(paste0('figures/DE/CD8T_all_samples_expansion_first8_virus_specific_not_excluded.pdf'), width = 15, height = 15)






