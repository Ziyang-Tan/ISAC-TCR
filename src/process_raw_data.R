library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(ggpubr)
source('src/load_BD_scTCR.R')

# parameters

bulk_data_dir <- '/Users/tan/Library/CloudStorage/OneDrive-KI.SE/TCR_processed_data/bulk'
bulk_proj_list <- c('P23556', 'P24864', 'P25060', 'P25755')
sc_data_dir <- '/Users/tan/Library/CloudStorage/OneDrive-KI.SE/TCR_processed_data/single cell'
sc_proj_list <- c('P23359_1001', 'P23359_1002', 'P23359_1003',
                  'P24851_1001', 'P24851_1002',
                  'P25158_1001',
                  'P25651_1001', 'P25651_1002')

# bulk TCR--------------------------------------------------------------------------------------------------

# read data

sample_info_bulk <- lapply(bulk_proj_list, function(proj_id){
  return(
    read_delim(
      Sys.glob(file.path(bulk_data_dir, proj_id, '*sample_info.txt')),
      delim = '\t',
      show_col_types = FALSE
    ) %>%
      tibble::add_column(proj_id = proj_id)
  )
}) %>% do.call(what = rbind) %>%
  tidyr::separate(`User ID`, into = c('Sample_Name', 'timepoint'), remove = F)

data_bulk <- lapply(sample_info_bulk$`NGI ID`, function(ngi_id){
  return(
    read_delim(
      Sys.glob(file.path(bulk_data_dir, '*', 'data', paste0(ngi_id, '*TRB*'))),
      delim = '\t',
      show_col_types = FALSE
    ) %>%
      tibble::add_column(ngi_id = ngi_id)
  )
}) %>% do.call(what = rbind)

# write to file
write_csv(data_bulk, file = file.path('.', 'data', 'bulk_data.csv.gz'))
write_csv(sample_info_bulk, file = file.path('.', 'data', 'bulk_info.csv.gz'))
# BD Rhap--------------------------------------------------------------------------------------------------

# read data

raw_tcr <- lapply(sc_proj_list, BD_load_VDJ, dir_path = sc_data_dir) %>% do.call(what = rbind)
sample_tag <- lapply(sc_proj_list, BD_load_sample_tag, dir_path = sc_data_dir) %>% do.call(what = rbind)
raw_gene <- lapply(sc_proj_list, BD_load_gene_exp, dir_path = sc_data_dir, norm_method = 'DBEC') %>% do.call(what = rbind)

# scRNA--------------------------------------------------------------------------------------------------
# process scRNA, determine major subsets
raw_gene_merge <- raw_gene %>% 
  left_join(sample_tag, by = 'unique_index') %>% 
  left_join(raw_tcr %>% select(unique_index, TCR_Paired_Chains, at_least_one_chain, is_gdT, proj_id), by = 'unique_index') %>%
  mutate(Sample_Name = case_when(
    is.na(Sample_Name) ~ proj_id,
    TRUE ~ Sample_Name
  ))
data_gene <- raw_gene_merge %>% 
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
counts_seurat <- data_gene %>% 
  select(-c('Sample_Tag', 'at_least_one_chain', 'proj_id',
            'Sample_Name', 'TCR_Paired_Chains', 'unique_index', 'is_gdT')) %>%
  t()
colnames(counts_seurat) <- data_gene$unique_index

meta_seurat <- data_gene %>%
  #tidyr::separate(proj_id, into = c('batch', NA), remove=F) %>%
  mutate(batch = proj_id) %>%
  select('Sample_Name', 'TCR_Paired_Chains', 'at_least_one_chain', 'is_gdT', 'unique_index', 'batch') %>%
  as.data.frame()
rownames(meta_seurat) <- data_gene$unique_index

obj <- CreateSeuratObject(counts = counts_seurat,
                          meta.data = meta_seurat)
obj <- SCTransform(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE, k.param = 25)
obj <- FindClusters(obj, algorithm = 3) # leiden
g1 <- DimPlot(obj, label = TRUE) + NoLegend()
g2 <- DimPlot(obj, group.by='TCR_Paired_Chains')
g3 <- DimPlot(obj, group.by='at_least_one_chain')
g4 <- DimPlot(obj, group.by='is_gdT')
g5 <- DimPlot(obj, group.by = 'orig.ident')
ggarrange(g1,g2,g3,g4,g5) %>% ggexport(width = 15, height = 10, filename = 'figures/overview_all_cells.pdf')
FeaturePlot(obj, features = c('CD3E', 'CD4', 'CD8A')) %>%
  ggsave(filename = 'figures/CD4_CD8_expression.pdf')
DoHeatmap(AverageExpression(obj, return.seurat = T), features = c('CD3E', 'CD4', 'CD8A')) %>%
  ggsave(filename = 'figures/CD4_CD8_heatmap.pdf')

# manual CD4T/CD8T discrimination
df <- obj[[]] %>%
  as_tibble() %>%
  mutate(cell_type = case_when(
    is_gdT ~ 'gdT', # should come first, or some gdT cells will be assigned into CD4T/CD8T, 
    # since they are in those seurat clusters
    seurat_clusters %in% c('1', '2', '5', '7', '9', '10', '16', '18', '19', '21', '22', '26') ~ 'CD4T',
    seurat_clusters %in% c('0', '3', '4', '6', '12', '14') ~ 'CD8T',
    TRUE ~ 'others'
  )) %>%
  cbind(obj@reductions$umap[[]])
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.5)
ggsave(filename = 'figures/CD4_CD8_types_manual.pdf')

cell_type <- df %>% 
  select(unique_index, cell_type, UMAP_1, UMAP_2, seurat_clusters)
write_csv(cell_type, file = 'data/cell_types_manual.csv.gz')

SaveH5Seurat(obj,filename = 'data/seurat_results.h5Seurat', overwrite = T)

# scTCR--------------------------------------------------------------------------------------------------
# clonal info by scTCR
cell_type <- read_csv(file = 'data/cell_types_manual.csv.gz')
# process
raw_tcr_merge <- left_join(raw_tcr, sample_tag, by = 'unique_index') %>%
  left_join(cell_type, by = 'unique_index') %>%
  mutate(Sample_Name = case_when(
    is.na(Sample_Name) ~ proj_id,
    TRUE ~ Sample_Name
  ))

# summarize 
df <- raw_tcr_merge %>% filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
table_summary <- df %>% group_by(proj_id) %>% summarise(demultiplexed = n()) %>% left_join(
  df %>% filter(at_least_one_chain) %>% group_by(proj_id) %>% summarise(at_least_one_CDR3 = n())
) %>% left_join(
  df %>% filter(TCR_Paired_Chains) %>% group_by(proj_id) %>% summarise(paired = n())
)
write_csv(table_summary, file = 'data/exp_quality_summary.csv')

# at least Vb chain (for GLIPH)
data_scTCR_vb <- raw_tcr_merge %>%
  filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
write_csv(data_scTCR_vb, file = 'data/scTCR_data_merge_vb.csv.gz')
# clonal expansion
data_scTCR <- raw_tcr_merge %>%
  filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
  filter(!is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant)) %>%
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined')) %>%
  mutate(CDR3_concat = paste0(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, '_', 
                              TCR_Beta_Delta_CDR3_Nucleotide_Dominant),
         CDR3aa_concat = paste0(TCR_Alpha_Gamma_CDR3_Translation_Dominant, '_', 
                                TCR_Beta_Delta_CDR3_Translation_Dominant)) %>%
  mutate(clone_id = as.character(as.numeric(as.factor(CDR3_concat))))
write_csv(data_scTCR, file = 'data/scTCR_data_merge_paired_chain.csv.gz')
clone_id_map <- data_scTCR %>% select(CDR3_concat, clone_id) %>% unique()
clone_exp <- data_scTCR %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
#inner_join(data %>% select(clone_id, CDR3aa_concat) %>% unique(), by='clone_id')
write_csv(clone_id_map, file = 'data/clone_id_map.csv.gz')
write_csv(clone_exp, file = 'data/clone_expansion.csv.gz')


