library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
source("src/clone_expansion_plots.R")

# load data
clone_id_map <- read_csv('data/clone_id_map.csv.gz')
data_scTCR <- read_csv('data/scTCR_data_merge.csv.gz')


# construct sample list

tmp <- tidyr::separate(data_scTCR %>% select(Sample_Name), Sample_Name, into=c('patient', NA), remove = F)
sample_list <-lapply(unique(tmp$patient), function(x){
  tmp %>% filter(patient == x) %>% select(Sample_Name) %>% distinct() %>% unlist(use.names = F) %>% gtools::mixedsort()
})
names(sample_list) <- unique(tmp$patient)

# donut chart 

for (patient in names(sample_list)){
  for (sub_name in c('CD4T', 'CD8T', 'gdT')){
    fig_dir <- file.path('figures', 'clonal expansion', patient)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    
    
    clone_exp_sub <- data_scTCR %>%
      filter(cell_type == sub_name) %>%
      # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
      group_by(CDR3_concat, Sample_Name) %>%
      summarise(clone_count = n()) %>%
      ungroup() %>%
      inner_join(clone_id_map, by='CDR3_concat')
    g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
    ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
      ggexport(filename = file.path(fig_dir, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
               width = 10, height = 20)
    #g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_gdT'))
    #ggsave(plot = g, filename = file.path(fig_dir, paste0(patient, '_top_clone_changes_', sub_name, '.pdf')))
  }
}

# public clone id

clone_exp <- data_scTCR %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')

public_clone_id <- clone_exp_sub %>% 
  tidyr::separate(Sample_Name, into = c('individual', NA), remove = F) %>% 
  select(individual, clone_id) %>% 
  distinct() %>% 
  group_by(clone_id) %>% 
  tally() %>%
  filter(n>1) %>%
  select(clone_id) %>% unlist(use.names = F)

clone_exp %>% filter(clone_id %in% public_clone_id)

