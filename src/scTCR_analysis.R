library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(EnhancedVolcano)
library(Seurat)
library(SeuratDisk)
source("src/clone_expansion_plots.R")

# load data
clone_id_map <- read_csv('data/clone_id_map.csv.gz')
data_scTCR <- read_csv('data/scTCR_data_merge_paired_chain.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz')
virus_specific_list <- read_csv(file = 'data/virus_specific_clone_id.csv')

obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types_manual.csv.gz') %>% as.data.frame()
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type
clone_id <- data_scTCR %>% select(clone_id) %>% mutate(clone_id=as.character(clone_id)) %>% as.data.frame()
rownames(clone_id) <- data_scTCR$unique_index
obj$clone_id <- clone_id
obj$patient_id <- sub('_.*$', '', obj$Sample_Name)

# construct sample list

tmp <- tidyr::separate(data_scTCR %>% select(Sample_Name), Sample_Name, into=c('patient', NA), remove = F)
sample_list <-lapply(unique(tmp$patient), function(x){
  tmp %>% filter(patient == x) %>% select(Sample_Name) %>% distinct() %>% unlist(use.names = F) %>% gtools::mixedsort()
})
names(sample_list) <- unique(tmp$patient)

# donut chart 
top_clones <- data.frame()
for (patient in names(sample_list)){
  fig_dir_clone <- file.path('figures', 'clonal expansion', patient)
  dir.create(fig_dir_clone, showWarnings = FALSE, recursive = TRUE)
  for (sub_name in c('CD4T', 'CD8T', 'gdT')){
    patient <- 'ISAC02'
    sub_name <- 'CD8T'
    clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
      inner_join(clone_id_map, by='CDR3_concat')
    g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
    ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
      ggexport(filename = file.path(fig_dir_clone, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
               width = 10, height = 20)
    g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_', sub_name))
    ggsave(plot = g, filename = file.path(fig_dir_clone, paste0(patient, '_top_clone_changes_', sub_name, '.pdf')))
    g <- clone_expansion_alluvium(patient, clone_exp_sub, clone_label = TRUE) + labs(title = paste0(patient, '_', sub_name))
    ggsave(plot = g, filename = file.path(fig_dir_clone, paste0(patient, '_top_clone_changes_', sub_name, '_with_label.pdf')))
    tmp <- data.frame(clone_id = unique(g$data$clone_id), sub_name = sub_name, patient = patient) %>%
      left_join(ggplot_build(g)$data[[1]] %>% select(fill, stratum) %>% distinct() %>% rename(clone_id = stratum), by='clone_id')
    top_clones <- rbind(top_clones, tmp)
  }
}
write_csv(top_clones, file = 'data/top_clones_all_timepoints.csv')

# for ISAC02
# trend determination

patient <- 'ISAC02'
sub_name <- 'CD8T'
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')
g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_', sub_name)) + guides(fill='none', color='none')

g2 <- trend_determination_plot('ISAC02', clone_exp_sub, top_mod='union')
increase_clone <- c(31042, 24360, 1709, 29062, 33659, 21963)
decrease_clone <- c(4834, 19762, 32933, 20882, 2490, 26712, 11642, 32213, 26805, 
                    30159, 12573, 11556, 30108, 20878, 27993, 20975, 34986)

# flow by trend
df <- clone_exp_alluvium_preparation('ISAC02', clone_exp_sub, top_mod='union')
df <- df %>% mutate(
  trend = case_when(
    clone_id %in% virus_specific_list$clone_id ~ 'virus_specific',
    clone_id %in% increase_clone ~ 'increase',
    clone_id %in% decrease_clone ~ 'decrease',
    TRUE ~ 'other'
  ))

g3 <- ggplot(df, aes(x = time_point, stratum = clone_id, alluvium = clone_id, 
                     #y= clone_ratio, 
                     y = relative_clone_ratio,
                     fill = trend, color=clone_id)) +
  geom_alluvium() +
  geom_stratum(size = 0.1) +
  guides(color = "none") + # hide the legend of 'color'
  #scale_fill_manual(values = wes_palette("Rushmore1", length(unique(df$trend)), type = "continuous"))+
  #theme(legend.position = "none")+
  labs(title = 'ISAC02')

obj_sub <- subset(obj, subset = patient_id == 'ISAC02' & 
                    clone_id %in% union(increase_clone, decrease_clone) &
                    cell_type == 'CD8T')
tmp <- obj_sub[[]] %>% left_join(df %>% select(clone_id, trend) %>% distinct(), by='clone_id') %>% select(trend)
rownames(tmp) <- obj_sub[[]]$unique_index
obj_sub$trend <- tmp
Idents(obj_sub) <- 'trend'
res <- FindMarkers(obj_sub, ident.1 = 'increase', ident.2 = 'decrease',
                   logfc.threshold = 0)

g4 <- EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
                      pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
                      selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
                      title = paste0('CD8T', ' of ', 'ISAC02'),
                      subtitle = 'increase vs decrease')

ggarrange(g, g3, g4, ncol = 3, common.legend=T) %>% ggexport(filename = 'figures/ISAC02_trend_RNA.pdf', width = 20, height = 8)

# add some GO?



 # PCA of top clones (by patient)

sub_name <- 'CD8T'
top_clones <- read_csv(file = 'data/top_clones_all_timepoints.csv')

#
patient <- 'ISAC02'
chosen_clones <- c(21963, 20878, 27993, 33659, 26805, 1709)
obj_sub <- subset(obj, subset = clone_id %in% chosen_clones & cell_type == sub_name & patient_id == patient)
obj_sub <- RunPCA(obj_sub, features = VariableFeatures(object = obj_sub))
#PCAPlot(obj=obj_sub, group.by = 'clone_id')
#cols = (top_clones %>% filter(clone_id %in% c(26805, 20878, 21963, 33659)))$fill)
g1<- dittoSeq::dittoDimPlot(obj_sub, reduction.use = 'pca', var = 'clone_id', do.ellipse = T)
#color.panel = (top_clones %>% filter(clone_id %in% chosen_clones))$fill)
g2<- VizDimLoadings(obj_sub, dims = 1:2, reduction = "pca")
ggarrange(g1, g2) %>% ggexport(filename = 'figures/ISAC02_CD8T_clone_RNA.pdf', width=12, height = 6)
dittoSeq::dittoDimPlot(obj_sub, reduction.use = 'pca', var = 'Sample_Name', do.ellipse = T)
ggsave(filename = 'figures/ISAC02_CD8T_clone_RNA_by_time.pdf')


# public clone id

# clone_exp <- data_scTCR %>%
#   group_by(CDR3_concat, Sample_Name) %>%
#   summarise(clone_count = n()) %>%
#   ungroup() %>%
#   inner_join(clone_id_map, by='CDR3_concat')
# 
# public_clone_id <- clone_exp %>% 
#   tidyr::separate(Sample_Name, into = c('individual', NA), remove = F) %>% 
#   select(individual, clone_id) %>% 
#   distinct() %>% 
#   group_by(clone_id) %>% 
#   tally() %>%
#   filter(n>1) %>%
#   select(clone_id) %>% unlist(use.names = F)
# 
# public_clone <- clone_exp %>% filter(clone_id %in% public_clone_id)
# write_csv(public_clone, file = 'data/public_clone.csv')

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



