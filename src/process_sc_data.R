library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(ggpubr)
source("src/load_BD_scTCR.R")

# parameters
sc_data_dir <- "/Users/tan/Library/CloudStorage/OneDrive-KI.SE/TCR_processed_data/single cell"
sc_proj_list <- c(
    "P23359_1001", "P23359_1002", "P23359_1003", "P24851_1001", "P24851_1002",
    "P25158_1001", "P25651_1001", "P25651_1002", "P29012_1001", "P29012_1002",
    "P29455_1001", "P29552_1001"
)
prior_label <- read_csv(file = "data/cell_types_manual.csv.gz")

# read data
raw_tcr <- lapply(sc_proj_list, BD_load_VDJ, dir_path = sc_data_dir) %>% do.call(what = rbind)
sample_tag <- lapply(sc_proj_list, BD_load_sample_tag, dir_path = sc_data_dir) %>% do.call(what = rbind)
raw_gene_list <- lapply(sc_proj_list, BD_load_gene_exp, dir_path = sc_data_dir, norm_method = "DBEC")
common_panel <- lapply(raw_gene_list, function(x) colnames(x)) %>% Reduce(f = intersect)
raw_gene <- lapply(raw_gene_list, function(x) x[common_panel]) %>% do.call(rbind, .)

# process scRNA, determine major subsets
raw_gene_merge <- raw_gene %>%
    left_join(sample_tag, by = "unique_index") %>%
    left_join(
        raw_tcr %>% select(unique_index, TCR_Paired_Chains, at_least_one_chain, is_gdT, proj_id),
        by = "unique_index"
    ) %>%
    mutate(Sample_Name = case_when(
        is.na(Sample_Name) ~ proj_id,
        TRUE ~ Sample_Name
    ))
data_gene <- raw_gene_merge %>%
    filter(!Sample_Name %in% c("Multiplet", "Undetermined"))
counts_seurat <- data_gene %>%
    select(-c(
        "Sample_Tag", "at_least_one_chain", "proj_id",
        "Sample_Name", "TCR_Paired_Chains", "unique_index", "is_gdT"
    )) %>%
    t()
colnames(counts_seurat) <- data_gene$unique_index
meta_seurat <- data_gene %>%
    # tidyr::separate(proj_id, into = c('batch', NA), remove=F) %>%
    mutate(batch = proj_id) %>%
    select("Sample_Name", "TCR_Paired_Chains", "at_least_one_chain", "is_gdT", "unique_index", "batch") %>%
    as.data.frame()
rownames(meta_seurat) <- data_gene$unique_index
obj <- CreateSeuratObject(
    counts = counts_seurat,
    meta.data = meta_seurat
)
obj <- SCTransform(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE, k.param = 25)
obj <- FindClusters(obj, algorithm = 3) # leiden
g1 <- DimPlot(obj, label = TRUE) + NoLegend()
g2 <- DimPlot(obj, group.by = "TCR_Paired_Chains")
g3 <- DimPlot(obj, group.by = "at_least_one_chain")
g4 <- DimPlot(obj, group.by = "is_gdT")
g5 <- DimPlot(obj, group.by = "orig.ident")
ggarrange(g1, g2, g3, g4, g5) %>% ggexport(width = 15, height = 10, filename = "figures/overview_all_cells_3.pdf")
FeaturePlot(obj, features = c("CD3E", "CD4", "CD8A")) %>%
    ggsave(filename = "figures/CD4_CD8_expression_3.pdf")
DoHeatmap(AverageExpression(obj, return.seurat = T), features = c("CD3E", "CD4", "CD8A")) %>%
    ggsave(filename = "figures/CD4_CD8_heatmap_3.pdf")

df <- obj[[]] %>%
    as_tibble() %>%
    left_join(prior_label[c("unique_index", "cell_type")]) %>%
    cbind(obj@reductions$umap[[]])

# ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
#   geom_point(size = 0.5)

dic <- df %>%
    filter(!is.na(cell_type), cell_type != "gdT") %>%
    group_by(cell_type, seurat_clusters) %>%
    tally() %>%
    group_by(seurat_clusters) %>%
    filter(n == max(n))

# manual CD4T/CD8T discrimination
# df <- obj[[]] %>%
#   as_tibble() %>%
#   mutate(cell_type = case_when(
#     is_gdT ~ 'gdT', # should come first, or some gdT cells will be assigned into CD4T/CD8T,
#     # since they are in those seurat clusters
#     seurat_clusters %in% c('1', '2', '5', '7', '9', '10', '16', '18', '19', '21', '22', '26') ~ 'CD4T',
#     seurat_clusters %in% c('0', '3', '4', '6', '12', '14') ~ 'CD8T',
#     TRUE ~ 'others'
#   )) %>%
#   cbind(obj@reductions$umap[[]])

# auto labeling by prior labels
df <- df %>%
    mutate(cell_type = case_when(
        !is.na(cell_type) ~ cell_type, # do not change the prior labels
        is_gdT ~ "gdT", # should come first, or some gdT cells will be assigned into CD4T/CD8T,
        # since they are in those seurat clusters
        seurat_clusters %in% (dic %>% filter(cell_type == "CD4T"))$seurat_clusters ~ "CD4T",
        seurat_clusters %in% (dic %>% filter(cell_type == "CD8T"))$seurat_clusters ~ "CD8T",
        seurat_clusters %in% (dic %>% filter(cell_type == "others"))$seurat_clusters ~ "others",
        # TRUE ~ 'others'
    ))

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.5)
ggsave(filename = "figures/CD4_CD8_types_auto.pdf")

cell_type <- df %>%
    select(unique_index, cell_type, UMAP_1, UMAP_2, seurat_clusters)
write_csv(cell_type, file = "data/cell_types_auto.csv.gz")

# SaveH5Seurat(obj,filename = 'data/seurat_results.h5Seurat', overwrite = T)

# scTCR--------------------------------------------------------------------------------------------------
# clonal info by scTCR
cell_type <- read_csv(file = "data/cell_types_auto.csv.gz")
# process
raw_tcr_merge <- left_join(raw_tcr, sample_tag, by = "unique_index") %>%
    left_join(cell_type, by = "unique_index") %>%
    mutate(Sample_Name = case_when(
        is.na(Sample_Name) ~ proj_id,
        TRUE ~ Sample_Name
    ))

# summarize
df <- raw_tcr_merge %>% filter(!Sample_Name %in% c("Multiplet", "Undetermined"))
table_summary <- df %>%
    group_by(proj_id) %>%
    summarise(demultiplexed = n()) %>%
    left_join(
        df %>% filter(at_least_one_chain) %>% group_by(proj_id) %>% summarise(at_least_one_CDR3 = n())
    ) %>%
    left_join(
        df %>% filter(TCR_Paired_Chains) %>% group_by(proj_id) %>% summarise(paired = n())
    )
write_csv(table_summary, file = "data/exp_quality_summary.csv")

# at least Vb chain (for GLIPH)
data_scTCR_vb <- raw_tcr_merge %>%
    filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
    filter(!Sample_Name %in% c("Multiplet", "Undetermined"))
write_csv(data_scTCR_vb, file = "data/scTCR_data_merge_vb.csv.gz")
# clonal expansion
data_scTCR <- raw_tcr_merge %>%
    filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
    filter(!is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant)) %>%
    filter(!Sample_Name %in% c("Multiplet", "Undetermined")) %>%
    mutate(
        CDR3_concat = paste0(
            TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, "_",
            TCR_Beta_Delta_CDR3_Nucleotide_Dominant
        ),
        CDR3aa_concat = paste0(
            TCR_Alpha_Gamma_CDR3_Translation_Dominant, "_",
            TCR_Beta_Delta_CDR3_Translation_Dominant
        )
    ) %>%
    mutate(clone_id = as.character(as.numeric(as.factor(CDR3_concat))))
write_csv(data_scTCR, file = "data/scTCR_data_merge_paired_chain.csv.gz")
clone_id_map <- data_scTCR %>%
    select(CDR3_concat, clone_id) %>%
    unique()
clone_exp <- data_scTCR %>%
    group_by(CDR3_concat, Sample_Name) %>%
    summarise(clone_count = n()) %>%
    ungroup() %>%
    inner_join(clone_id_map, by = "CDR3_concat")
# inner_join(data %>% select(clone_id, CDR3aa_concat) %>% unique(), by='clone_id')
write_csv(clone_id_map, file = "data/clone_id_map.csv.gz")
write_csv(clone_exp, file = "data/clone_expansion.csv.gz")

# write Seurat obj with all info
cell_type <- read.csv("data/cell_types_auto.csv.gz") %>% as.data.frame()
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type
clone_id <- data_scTCR %>%
    select(clone_id) %>%
    mutate(clone_id = as.character(clone_id)) %>%
    as.data.frame()
rownames(clone_id) <- data_scTCR$unique_index
obj$clone_id <- clone_id
# fix an ISAC id
tmp <- obj$Sample_Name
tmp[tmp == "ISAC_4"] <- "ISAC134_4"
obj$Sample_Name <- tmp
obj$patient_id <- sub("_.*$", "", obj$Sample_Name)
# assign tumor type
tmp <- obj[[]] %>%
    select(patient_id) %>%
    mutate(
        tumor_type = case_when(
            patient_id %in% c("ISAC35", "ISAC99") ~ "Wilms",
            patient_id %in% c("ISAC77") ~ "Osteosarcoma",
            patient_id %in% c("ISAC100", "ISAC141", "ISAC31", "ISAC02") ~ "High mutation Neuroblastoma",
            patient_id %in% c("ISAC112", "ISAC134", "ISAC62", "ISAC125", "ISAC46") ~ "Low mutation Neuroblastoma",
            patient_id %in% c("ISAC81") ~ "Other Neuroblastoma"
        )
    )
obj$tumor_type <- tmp$tumor_type


SaveH5Seurat(obj, filename = "data/seurat_results_update.h5Seurat", overwrite = T)

# integrate more information into the Seurat obj
obj <- LoadH5Seurat("data/seurat_results_update.h5Seurat")

clone_id_map <- read_csv("data/clone_id_map.csv.gz")
data_sctcr <- read_csv("data/scTCR_data_merge_paired_chain.csv.gz") %>%
    mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = "data/clone_expansion.csv.gz")
virus_specific_list <- read_csv(file = "data/virus_specific_clone_id.csv")
patient_info <- read_csv(file = "data/ISAC_TCR_sample_info.csv")

# is expanded?
expand_id <- (clone_exp %>% filter(
    clone_count >= 2,
    # !clone_id %in% virus_specific_list$clone_id
) %>% select(clone_id) %>% unique())[[1]] # can include expansion from different sub type
expand_sc <- data_sctcr %>%
    filter(clone_id %in% expand_id) %>%
    tibble::add_column(expansion = "expanded") %>%
    select(unique_index, expansion)
tmp <- obj[[]] %>%
    filter(TCR_Paired_Chains == 1) %>%
    left_join(expand_sc, by = "unique_index") %>%
    tidyr::replace_na(list(expansion = "non-expanded"))
tmp2 <- tmp$expansion
names(tmp2) <- tmp$unique_index
obj$expansion <- tmp2

# is virus specific?
virus_specific_sc <- data_sctcr %>%
    filter(clone_id %in% virus_specific_list$clone_id) %>%
    tibble::add_column(virus_specific = "virus_specific") %>%
    select(unique_index, virus_specific)
tmp <- obj[[]] %>%
    filter(TCR_Paired_Chains == 1) %>%
    left_join(virus_specific_sc, by = "unique_index") %>%
    tidyr::replace_na(list(virus_specific = "non_virus_specific"))
tmp2 <- tmp$virus_specific
names(tmp2) <- tmp$unique_index
obj$virus_specific <- tmp2

obj[[]] %>%
    group_by(expansion, virus_specific, TCR_Paired_Chains) %>%
    tally()

View(obj[[]])

SaveH5Seurat(obj, filename = "data/seurat_results_update2.h5Seurat", overwrite = TRUE)
