library(dplyr)
library(readr)
library(igraph)
source('src/GLIPH_related_plots.R')
source('src/clone_expansion_plots.R')

data_TRA <- read_csv('data/TZY_TRA_annotated.csv')
data_TRB <- read_csv('data/TZY_TRB_annotated.csv')

data_TRAB <- rbind(data_TRA, data_TRB) %>%
  mutate(CDR3aa_concat = paste0(cdr3a, '_', cdr3b),
         Sample_Name = paste0(sid, '_', condition),
         cluster_id = paste0(index, '_', stringr::str_extract(pattern, '^.{4}')),
         virus_specific = !is.na(CDR3a) | !is.na(CDR3b))

clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz') %>%
  mutate(clone_id = as.character(clone_id),
         patient_id = sub('_.*$', '', Sample_Name))
data_scTCR <- read_csv('data/scTCR_data_merge_paired_chain.csv.gz')
top_clones <- read_csv('data/top_clones_all_timepoints.csv') %>%
  filter(sub_name == 'CD8T')
nt2aa <- data_scTCR %>% select(CDR3_concat, CDR3aa_concat) %>% distinct()

virus_specific_list <- left_join(
  data_scTCR %>% select(CDR3aa_concat, clone_id) %>% distinct(),
  data_TRAB %>% select(CDR3aa_concat, virus_specific) %>% distinct(), by = 'CDR3aa_concat') %>%
  filter(virus_specific) %>%
  select(clone_id) %>% 
  distinct()

write_csv(virus_specific_list, 'data/virus_specific_clone_id.csv')

# network of all samples
pdf(paste0('figures/GLIPH_network/GLIPH_with_known_viral_clusters_clonecount', 3, '.pdf'))
plot.gliph(data_TRA, data_TRB, clone_exp, nt2aa, clone_thre=3)
title('all 8 samples')
dev.off()

# network per sample

for (name_i in unique(clone_exp$patient_id)){
  pdf(paste0('figures/GLIPH_network/GLIPH_with_known_viral_clusters_clonecount', 1, '_', name_i, '.pdf'))
  plot.gliph(data_TRA %>% filter(sid==name_i), 
             data_TRB %>% filter(sid==name_i), 
             clone_exp %>% filter(patient_id==name_i), 
             nt2aa, 
             clone_thre=1,
             labels = filter(top_clones, patient==name_i)$clone_id)
  title(main=name_i)
  dev.off()
}

# clonal expansion and changes on GLIPH cluster level

cluster_exp <- data_TRAB %>% 
  group_by(pattern, Sample_Name, cluster_id) %>% 
  summarise(clone_count=sum(frequency)) %>%
  rename(clone_id = cluster_id,
         CDR3_concat = pattern)
  #mutate(patient = sub('_.*$', '', Sample_Name))


g_list <- lapply(unique(data_TRAB$sid), clone_expansion_alluvium, cluster_exp, top_mod = 'large', n=10)
ggarrange(plotlist = g_list, ncol = 2, nrow = 2) %>%
  ggexport(filename = 'figures/GLIPH_cluster_changes_CD8T.pdf')


# analysis of graph
comp <- components(g)
clusters <- groups(comp)
big_clusters <- clusters[comp$csize > 3]


tmp <- data_scTCR %>% filter(clone_id %in% clusters$`5`)
tmp %>% group_by(Sample_Name, clone_id) %>% tally() %>% View()

data_TRAB %>% filter(CDR3aa_concat %in% unique(tmp$CDR3aa_concat)) %>% View()

