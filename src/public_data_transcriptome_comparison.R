library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

data_t <- readRDS('data/adult data from public/GSE139555_tcell_integrated.rds')
meta_t <- read.delim('data/adult data from public/GSE139555_tcell_metadata.txt.gz', sep = '\t')
panel_adult <- rownames(data_t)

obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
panel_isac <- rownames(obj)

# find common panel
# library("org.Hs.eg.db")
# mat1<- select(org.Hs.eg.db, keys=panel_isac, columns=c("ENTREZID"), keytype="ALIAS") %>% filter(!is.na(ENTREZID))
# mat2<- select(org.Hs.eg.db, keys=panel_adult, columns=c("ENTREZID"), keytype="ALIAS") %>% filter(!is.na(ENTREZID))
# detach("package:org.Hs.eg.db", unload = TRUE)
# detach("package:AnnotationDbi", unload = TRUE)

# map_df <- inner_join(mat1, mat2, by='ENTREZID', suffix = c('.1', '.2')) %>% dplyr::rename(name_adult = ALIAS.2, name_isac = ALIAS.1)
# common_entrezid <- intersect(mat1$ENTREZID, mat2$ENTREZID) %>% na.omit()

common_panel <- intersect(panel_adult, panel_isac)

# construct the merged data frame

isac_vst <- GetAssayData(obj, slot = 'count') %>% as.matrix() %>% sctransform::vst()

data_merge <- rbind(t(isac_vst$y) %>% as.data.frame() %>% select(all_of(common_panel)), 
                    t(data_t[common_panel, ]) %>% as.data.frame())
meta_merge <- rbind(obj[[]] %>% select(cell_type, patient_id) %>% mutate(source = 'Blood', group = 'isac'),
                    meta_t %>% select(ident, source, patient) %>% rename(cell_type = ident, patient_id = patient) %>% mutate(group = 'adult'))

data_correct <- sva::ComBat(t(data_merge), batch = meta_merge$group) %>% t() %>% as_tibble()

df <- cbind(data_correct, meta_merge) %>% filter(source == 'Blood', cell_type %in% c('8.2-Tem', '8.3a-Trm', '8.3b-Trm', '8.3c-Trm', 'CD8T'))


# directly compare transcript -------------
ggplot(df, aes(x=group, y=LGALS1)) + geom_jitter()

tmp <- df %>% group_by(group) %>% summarise(across(all_of(common_panel), mean))
tmp <- t(tmp[,common_panel])
colnames(tmp) <- c('adult', 'isac')
res <- tmp %>% as.data.frame() %>% mutate(logFC = isac - adult)
res$gene <- rownames(res)

library(rstatix)

tmp <- df %>% tidyr::pivot_longer(values_to = 'level', names_to = 'gene', cols = all_of(common_panel))
res_tmp <- tmp %>% 
  group_by(gene) %>% 
  t_test(level ~ group) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

res <- left_join(res, res_tmp, by = 'gene')

library(EnhancedVolcano)
EnhancedVolcano(res, x='logFC', y='p.adj', lab = res$gene, 
                pCutoff = 1e-04, FCcutoff = NA, drawConnectors = TRUE, arrowheads = F, 
                selectLab = (res %>% filter(p.adj < 1e-10))$gene,
                xlim=c(-1,1),
                col= 'black',
                title = 'memory CD8 T cells in Blood',
                subtitle = 'isac vs adult') + theme_minimal()
ggsave('figures/compare with adults/volcano mCD8T in blood.pdf', width = 10, height = 10)

ggplot(df, aes(x=group, y=CCL5)) + geom_jitter()
# ggplot(df, aes(x=group, y=CCL5)) + geom_violin()
ggsave('figures/compare with adults/distribution example mCD8T in blood.pdf')

# pca
tmp <- df %>% group_by(patient_id, group) %>% summarise(across(all_of(common_panel), mean))
pca.res <- prcomp(tmp[,common_panel])
cbind(pca.res$x, tmp) %>% 
  ggplot(aes(x=PC1, y=PC2, color=group)) + geom_point() + theme_minimal()
ggsave('figures/compare with adults/pca mCD8T in blood.pdf')

# umap
library(umap)
umap.res <- umap(df[,common_panel])
cbind(umap.res$layout, df) %>% ggplot(aes(x=`1`, y=`2`, color=group)) + geom_point() + xlim(-6, 6) + ylim(-6,6) + theme_minimal()
ggsave('figures/compare with adults/umap mCD8T in blood.png')

# in adults, highly expanded vs others --------------

data_adult <- t(data_t) %>% as.data.frame()
data_adult$X = rownames(data_adult)
data_adult <- as_tibble(data_adult)

data_adult_merge <- left_join(data_adult %>% select(all_of(common_panel), X), meta_t, by='X')

df <- data_adult_merge %>% filter(source == 'Blood', ident %in% c('8.2-Tem', '8.3a-Trm', '8.3b-Trm', '8.3c-Trm'))

expanded_clones <- df %>% group_by(clonotype) %>% tally() %>% filter(n>1) %>% select(clonotype) %>% unlist() %>% na.omit()

df <- df %>% mutate(expansion = if_else(clonotype %in% expanded_clones, 'expanded', 'non_expanded'))

tmp <- df %>% group_by(expansion) %>% summarise(across(all_of(common_panel), mean))
tmp <- t(tmp[,common_panel])
colnames(tmp) <- c('expanded', 'non_expanded')
res <- tmp %>% as.data.frame() %>% mutate(logFC = expanded - non_expanded)
res$gene <- rownames(res)

library(rstatix)

tmp <- df %>% tidyr::pivot_longer(values_to = 'level', names_to = 'gene', cols = all_of(common_panel))
res_tmp <- tmp %>% 
  group_by(gene) %>% 
  t_test(level ~ expansion) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

res <- left_join(res, res_tmp, by = 'gene')

library(EnhancedVolcano)
EnhancedVolcano(res, x='logFC', y='p.adj', lab = res$gene, 
                pCutoff = 1e-04, FCcutoff = NA, drawConnectors = TRUE, arrowheads = F, 
                selectLab = (res %>% filter(p.adj < 1e-5))$gene,
                xlim=c(-1,1),
                col= 'black',
                title = 'memory CD8 T cells in Blood (adult)',
                subtitle = 'expanded vs non-expanded') + theme_minimal()
ggsave('figures/compare with adults/volcano mCD8T in blood (adult) expanded vs non-expanded.pdf', width = 10, height = 10)

# cluster the adult's memory T cells ---------

ggplot(df, aes(x=UMAP_1, y=UMAP_2, color = ident)) + geom_point(size = 1, alpha = 0.5) + theme_bw()
ggsave('figures/compare with adults/umap mCD8T in blood (adult) original label.pdf', width = 10, height = 10)

# try use seurat pipeline
ind <- meta_t$source == 'Blood' & meta_t$ident %in% c('8.2-Tem', '8.3a-Trm', '8.3b-Trm', '8.3c-Trm')
rownames(meta_t) <- meta_t$X
obj_adult <- CreateSeuratObject(counts = data_t[,ind], 
                                meta.data = meta_t[ind,])
obj_adult <- ScaleData(obj_adult)
obj_adult <- FindVariableFeatures(obj_adult)
obj_adult <- RunPCA(obj_adult, verbose = FALSE)
obj_adult <- RunUMAP(obj_adult, dims = 1:30, verbose = FALSE)
obj_adult <- FindNeighbors(obj_adult, dims = 1:30, verbose = FALSE, k.param = 25)
obj_adult <- FindClusters(obj_adult, algorithm = 3) # leiden
DimPlot(obj_adult, label = TRUE) + NoLegend()
obj_adult$subtype <- obj_adult$ident
DimPlot(obj_adult, group.by='subtype')

pbmc.markers <- FindAllMarkers(obj_adult, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(obj_adult, features = unique(top10$gene))



# TIM3
FeaturePlot(obj, 'HAVCR2')





