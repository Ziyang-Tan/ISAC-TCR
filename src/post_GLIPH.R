library(dplyr)
library(readr)

clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz')
sample_info <- read_csv(file = 'data/ISAC_TCR_sample_info.csv') %>%
  mutate(individual = paste0('ISAC', study_id))

raw <- read_csv(file = 'data/GLIPH_output_mix_visit.csv') 
data <- raw %>% 
  rename(GLIPH_cluster = index) %>%
  mutate(individual = sub(':.*', '', Sample)) %>%
  left_join(sample_info, by = 'individual') %>%
  mutate(log10_fisher = -log10(Fisher_score)) %>%
  mutate(log2_expansion_score = -log2(expansion_score))

highest_freq <- data %>%
  group_by(GLIPH_cluster) %>%
  summarise(highest_freq = max(Freq)) %>%
  mutate(log2_highest_freq = log2(highest_freq))


ggplot(data %>% 
         distinct(GLIPH_cluster, log10_fisher, log2_expansion_score, number_subject) %>% 
         left_join(highest_freq, by='GLIPH_cluster') %>%
         mutate(label = case_when(
           log2_highest_freq > 2 & log10_fisher > 7.5 ~ GLIPH_cluster,
           log2_highest_freq > 7 ~ GLIPH_cluster,
           log10_fisher > 20 ~ GLIPH_cluster)), 
       aes(x=log2_expansion_score, y=log10_fisher)) +
  geom_point(aes(size = number_subject)) +
  scale_radius()
  #ggrepel::geom_text_repel(aes(label=label))

ggsave('figures/GLIPH_cluster_overview_exp_score.pdf')

data %>% filter(Freq > 5) %>% select(GLIPH_cluster) %>% distinct() %>% unlist(use.names = F)









tmp <- data %>% group_by(GLIPH_cluster, individual) %>% summarise(clone_sum = sum(Freq))

tmp2 <- tmp %>% group_by(GLIPH_cluster) %>% summarise(count = sum(clone_sum))

tmp3 <- sapply(tmp2$GLIPH_cluster, function(x){
  data %>% 
    filter(GLIPH_cluster == x) %>% 
    select(individual) %>% 
    distinct() %>% 
    unlist(use.names = F) %>% 
    paste(collapse = ', ')
})
summary <- cbind(tmp2,tmp3)
