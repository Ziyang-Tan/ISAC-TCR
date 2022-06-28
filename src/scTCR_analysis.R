library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(EnhancedVolcano)
source("src/clone_expansion_plots.R")

# load data
clone_id_map <- read_csv('data/clone_id_map.csv.gz')
data_scTCR <- read_csv('data/scTCR_data_merge_paired_chain.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz')
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types_manual.csv.gz')
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type

# construct sample list

tmp <- tidyr::separate(data_scTCR %>% select(Sample_Name), Sample_Name, into=c('patient', NA), remove = F)
sample_list <-lapply(unique(tmp$patient), function(x){
  tmp %>% filter(patient == x) %>% select(Sample_Name) %>% distinct() %>% unlist(use.names = F) %>% gtools::mixedsort()
})
names(sample_list) <- unique(tmp$patient)

# donut chart 

for (patient in names(sample_list)){
  fig_dir_clone <- file.path('figures', 'clonal expansion', patient)
  dir.create(fig_dir_clone, showWarnings = FALSE, recursive = TRUE)
  for (sub_name in c('CD4T', 'CD8T', 'gdT')){
    
    clone_exp_sub <- data_scTCR %>%
      filter(cell_type == sub_name) %>%
      # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
      group_by(CDR3_concat, Sample_Name) %>%
      summarise(clone_count = n()) %>%
      ungroup() %>%
      inner_join(clone_id_map, by='CDR3_concat')
    g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
    ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
      ggexport(filename = file.path(fig_dir_clone, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
               width = 10, height = 20)
    #g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_gdT'))
    #ggsave(plot = g, filename = file.path(fig_dir_clone, paste0(patient, '_top_clone_changes_', sub_name, '.pdf')))
  }
}

# public clone id

clone_exp <- data_scTCR %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')

public_clone_id <- clone_exp %>% 
  tidyr::separate(Sample_Name, into = c('individual', NA), remove = F) %>% 
  select(individual, clone_id) %>% 
  distinct() %>% 
  group_by(clone_id) %>% 
  tally() %>%
  filter(n>1) %>%
  select(clone_id) %>% unlist(use.names = F)

public_clone <- clone_exp %>% filter(clone_id %in% public_clone_id)
write_csv(public_clone, file = 'data/public_clone.csv')

# DE analysis highly expanded vs others

for (patient in names(sample_list)){
  fig_dir_DE <- file.path('figures', 'DE highly expanded vs others', patient)
  dir.create(fig_dir_DE, showWarnings = FALSE, recursive = TRUE)
  for (cur_subpop in c('CD4T', 'CD8T', 'gdT')){
    # cur_subpop <- 'CD8T'
    # cur_sample <- 'ISAC99_6'
    g_list2 <- lapply(sample_list[[patient]], function(cur_sample){
      obj_sub <- subset(obj, subset = Sample_Name == cur_sample & cell_type == cur_subpop & TCR_Paired_Chains)
      high_expand_id <- (clone_exp %>% filter(clone_count >= 10, Sample_Name == cur_sample) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
      high_expand_sc <- data_scTCR %>% filter(Sample_Name == cur_sample, cell_type == cur_subpop, clone_id %in% high_expand_id) %>% 
        tibble::add_column(highly_expand = 'highly_expanded') %>%
        select(unique_index, highly_expand)
      
      if (dim(high_expand_sc)[1] < 3) return(NA) # no highly expanded clones in this sample x cell type
      
      tmp <- obj_sub[[]] %>% left_join(high_expand_sc, by = 'unique_index') %>%
        tidyr::replace_na(list(highly_expand='others'))
      tmp2 <- tmp$highly_expand
      names(tmp2) <- tmp$unique_index
      obj_sub$highly_expand <- tmp2
      Idents(obj_sub) <- 'highly_expand'
      res <- FindMarkers(obj_sub, ident.1 = 'highly_expanded', ident.2 = 'others',
                         logfc.threshold = 0)
      return(
        EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
                        pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
                        selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
                        title = paste0(cur_subpop, ' of ', cur_sample),
                        subtitle = 'highly_expanded cells vs others')
      )
    })
    g_list2 <- g_list2[!is.na(g_list2)]
    if (length(g_list2) == 0) next
    ggarrange(plotlist = g_list2, ncol = 3, nrow = 2) %>%
      ggexport(filename = file.path(fig_dir_DE, paste0(patient, '_', cur_subpop, '.pdf')), width = 20, height = 15)
    #ggsave(paste0('figures/gene_DE/', cur_sample, '_', cur_subpop, '.pdf'))
  }
}



