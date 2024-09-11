library(dplyr)
library(readr)
library(Seurat)
library(SeuratDisk)
library(rstatix)

# ranked list in Wu
df_wu <- read_csv("data/adult data from public/blood mCD8T from dataset_Wu.csv.gz")
panel_wu <- colnames(df_wu)[1:3734]
tmp <- df_wu %>%
    group_by(expansion) %>%
    summarise(across(all_of(panel_wu), mean))
od <- tmp$expansion
tmp <- t(tmp[, panel_wu])
colnames(tmp) <- od
res_fc <- tmp %>%
    as.data.frame() %>%
    mutate(logFC = expanded - non_expanded)
res_fc$gene <- rownames(res_fc)
tmp <- df_wu %>% tidyr::pivot_longer(values_to = "level", names_to = "gene", cols = all_of(panel_wu))
res_p <- tmp %>%
    group_by(gene) %>%
    t_test(level ~ expansion) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
rank_wu <- left_join(res_fc, res_p, by = "gene") %>%
    mutate(rank_metric = logFC * -log10(p.adj)) %>%
    arrange(desc(rank_metric)) %>%
    rename(p_val_adj = p.adj)
write_csv(rank_wu, "data/DE_results/adults_Wu.csv")

# ranked list in Guo
obj_guo <- LoadH5Seurat("data/adult data from public/dataset_Guo.h5Seurat")
obj_guo_memt <- subset(obj_guo,
    subset = Clone.Status != "NA" &
        sampleType %in% c("PTC", "PTH", "PTR", "PTY") & # peripheral blood
        majorCluster %in% c(
            "CD8_C2-CD28", "CD8_C3-CX3CR1", "CD8_C4-GZMK", # memory CD8T cells
            "CD8_C5-ZNF683", "CD8_C6-LAYN", "CD8_C7-SLC4A10"
        )
)
Idents(obj_guo_memt) <- "Clone.Status"
rank_guo <- FindMarkers(obj_guo_memt,
    ident.1 = "Clonal", ident.2 = "NoClonal",
    logfc.threshold = 0, min.pct = 0.5
) %>%
    mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
    arrange(desc(rank_metric))
rank_guo$gene <- rownames(rank_guo)
write_csv(rank_guo, "data/DE_results/adults_Guo.csv")

# ranked list in Zhang
obj_zhang <- LoadH5Seurat("data/adult data from public/dataset_Zhang.h5Seurat")
obj_zhang_memt <- subset(obj_zhang,
    subset = Clonal.status != "NA" &
        sampleType %in% c("PTC", "PTH", "PTR", "PTY", "PP7") & # peripheral blood cells
        majorCluster %in% c(
            "CD8_C02-GPR183", "CD8_C03-CX3CR1", "CD8_C04-GZMK", "CD8_C05-CD6",
            "CD8_C06-CD160", "CD8_C07-LAYN", "CD8_C08-SLC4A10"
        ) # exclude CD8 T cells with LEF1 expression (memory phenotype)
)
Idents(obj_zhang_memt) <- "Clonal.status"
rank_zhang <- FindMarkers(obj_zhang_memt,
    ident.1 = "Clonal", ident.2 = "NoClonal",
    logfc.threshold = 0, min.pct = 0.5
) %>%
    mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
    arrange(desc(rank_metric))
rank_zhang$gene <- rownames(rank_zhang)
write_csv(rank_zhang, "data/DE_results/adults_Zhang.csv")

# ranked list in ISAC
## exclude virus specific
obj <- LoadH5Seurat("data/seurat_results_update2.h5Seurat") %>%
    subset(subset = cell_type == "CD8T" &
        virus_specific == "non_virus_specific" &
        TCR_Paired_Chains == 1)

for (i in unique(obj$tumor_type)) {
    obj_sub <- subset(obj,
        subset = tumor_type == i
    )
    Idents(obj_sub) <- "expansion"
    rank_tmp <- FindMarkers(obj_sub,
        ident.1 = "expanded", ident.2 = "non-expanded",
        logfc.threshold = 0, min.pct = 0
    ) %>%
        mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
        arrange(desc(rank_metric))
    rank_tmp$gene <- rownames(rank_tmp)
    write_csv(rank_tmp, paste0("data/DE_results/ISAC_", i, ".csv"))
}
