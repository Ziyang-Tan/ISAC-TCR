library(readr)
library(dplyr)
library(Seurat)
library(SeuratDisk)

# dataset Wu ------------------------------------
# Wu, T.D., Madireddi, S., de Almeida, P.E. et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. Nature 579, 274–278 (2020).
data_t <- readRDS('data/adult data from public/GSE139555_Wu/GSE139555_tcell_integrated.rds')
meta_t <- read.delim('data/adult data from public/GSE139555_Wu/GSE139555_tcell_metadata.txt.gz', sep = '\t')

data_Wu <- t(data_t) %>% as.data.frame()
data_Wu$X = rownames(data_Wu)
data_Wu <- as_tibble(data_Wu)
panel_Wu <- rownames(data_t)
data_Wu_merge <- left_join(data_Wu, meta_t, by='X')
df <- data_Wu_merge %>% filter(source == 'Blood', ident %in% c('8.2-Tem', '8.3a-Trm', '8.3b-Trm', '8.3c-Trm'))
expanded_clones <- df %>% group_by(clonotype) %>% tally() %>% filter(n>1) %>% select(clonotype) %>% unlist() %>% na.omit()
df <- df %>% mutate(
  expansion = if_else(clonotype %in% expanded_clones, 'expanded', 'non_expanded'),
  tumor_type = case_when(
    patient %in% c('Lung6') ~ 'Lung',
    patient %in% c('Renal1', 'Renal2', 'Renal3') ~ 'Renal'
  )
)
write_csv(df, 'data/adult data from public/blood mCD8T from dataset_Wu.csv.gz')

# dataset Guo ------------------------------------
# Guo, X., Zhang, Y., Zheng, L. et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing. Nat Med 24, 978–985 (2018). 

raw_counts <- read.delim('data/adult data from public/GSE99254_Guo/GSE99254_NSCLC.TCell.S12346.count.txt.gz', sep = '\t')
meta <- read.delim('data/adult data from public/GSE99254_Guo/meta.txt', sep = '\t')
data_scTCR <- readxl::read_excel('data/adult data from public/GSE99254_Guo/scTCR_info.xlsx')

raw_counts <- raw_counts %>% select(-geneID) %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE)

counts_seurat <- raw_counts %>% select(-symbol)
rownames(counts_seurat) <- raw_counts$symbol

meta_seurat <- meta %>% left_join(
  data_scTCR %>% select(Cell.Name, Clone.Status, Clone.Frequency) %>% rename(UniqueCell_ID = Cell.Name)
) %>% mutate(UniqueCell_ID = gsub('\\-', '\\.', UniqueCell_ID))
rownames(meta_seurat) <- meta_seurat$UniqueCell_ID

obj <- CreateSeuratObject(counts = counts_seurat,
                          meta.data = meta_seurat)
obj <- SCTransform(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE, k.param = 25)
obj <- FindClusters(obj, algorithm = 3) # leiden

SaveH5Seurat(obj,filename = 'data/adult data from public/dataset_Guo.h5Seurat', overwrite = T)

# dataset Zhang ---------------------------------
# Zhang, L., Yu, X., Zheng, L. et al. Lineage tracking reveals dynamic relationships of T cells in colorectal cancer. Nature 564, 268–272 (2018).

raw_counts <- read.delim('data/adult data from public/GSE108989_Zhang/GSE108989_CRC.TCell.S11138.count.txt.gz', sep = '\t')
meta <- read.delim('data/adult data from public/GSE108989_Zhang/meta.txt', sep = '\t')
data_scTCR <- readxl::read_excel('data/adult data from public/GSE108989_Zhang/scTCR_info.xlsx')

raw_counts <- raw_counts %>% select(-geneID) %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE)

counts_seurat <- raw_counts %>% select(-symbol)
rownames(counts_seurat) <- raw_counts$symbol

meta_seurat <- meta %>% left_join(
  data_scTCR %>% select(`Cell name`, `Clonal status`, Frequency) %>% rename(UniqueCell_ID = `Cell name`)
) %>% mutate(UniqueCell_ID = gsub('\\-', '\\.', UniqueCell_ID))
rownames(meta_seurat) <- meta_seurat$UniqueCell_ID

obj <- CreateSeuratObject(counts = counts_seurat,
                          meta.data = meta_seurat)
obj <- SCTransform(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE, k.param = 25)
obj <- FindClusters(obj, algorithm = 3) # leiden

SaveH5Seurat(obj,filename = 'data/adult data from public/dataset_Zhang.h5Seurat', overwrite = T)

