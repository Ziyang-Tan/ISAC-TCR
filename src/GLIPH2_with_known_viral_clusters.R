library(dplyr)
library(readr)
library(igraph)


data_TRA <- read_csv('data/TZY_TRA_annotated.csv')
data_TRB <- read_csv('data/TZY_TRB_annotated.csv')
GLIPH_known_specificity <- rbind(data_TRA, data_TRB) %>%
  mutate(CDR3aa_concat = paste0(cdr3a, '_', cdr3b)) %>%
  select(CDR3aa_concat, CDR3a, CDR3b) %>%
  distinct() %>%
  rename(CDR3a_known = CDR3a,
         CDR3b_known = CDR3b)

clone_exp <- read_csv(file = 'data/clone_expansion.csv.gz') %>%
  mutate(clone_id = as.character(clone_id))
data_scTCR <- read_csv('data/scTCR_data_merge_paired_chain.csv.gz')

nt2aa <- data_scTCR %>% select(CDR3_concat, CDR3aa_concat) %>% distinct() %>%
  left_join(GLIPH_known_specificity, by = "CDR3aa_concat")

data <- clone_exp %>% left_join(nt2aa, by='CDR3_concat', multiple = 'all') %>%
  filter(clone_count > 5) %>%
  arrange(desc(clone_count))

write_csv(data, 'data/clone_exp_with_known_specificity.csv')

# graph
# connection with known specificity
edge_a <- data %>% 
  #filter(!is.na(CDR3a_known), clone_count > 1) %>%
  select(clone_id, CDR3a_known) %>%
  mutate(CDR3a_known = strsplit(CDR3a_known, ":")) %>%
  tidyr::unnest(CDR3a_known) %>%
  rename(V1 = clone_id, V2 = CDR3a_known) %>%
  mutate(source = if_else(is.na(V2), NA, 'CDR3a'))

edge_b <- data %>% 
  #filter(!is.na(CDR3b_known), clone_count > 1) %>%
  select(clone_id, CDR3b_known) %>%
  mutate(CDR3b_known = strsplit(CDR3b_known, ":")) %>%
  tidyr::unnest(CDR3b_known) %>%
  rename(V1 = clone_id, V2 = CDR3b_known) %>%
  mutate(source = if_else(is.na(V2), NA, 'CDR3b'))

el <- rbind(edge_a, edge_b) %>% distinct()

g1 <- graph_from_edgelist(as.matrix(el[!is.na(el$V2),c('V1', 'V2')]), directed = F)

# connection with each other (same GLIPH1 group)

edge_a_2 <- data_TRA %>% mutate(CDR3aa_concat = paste0(cdr3a, '_', cdr3b)) %>% 
  filter(CDR3aa_concat %in% unique(data$CDR3aa_concat)) %>% 
  select(index, CDR3aa_concat) %>%
  left_join(data, by='CDR3aa_concat') %>%
  select(index, clone_id) %>%
  distinct() %>% 
  group_by(index) %>%
  filter(n()>=2) %>%
  do(data.frame(t(combn(.$clone_id, 2)))) %>% # elegant way
  mutate(source='CDR3a')

edge_b_2 <- data_TRB %>% mutate(CDR3aa_concat = paste0(cdr3a, '_', cdr3b)) %>% 
  filter(CDR3aa_concat %in% unique(data$CDR3aa_concat)) %>% 
  select(index, CDR3aa_concat) %>%
  left_join(data, by='CDR3aa_concat') %>%
  select(index, clone_id) %>%
  distinct() %>% 
  group_by(index) %>%
  filter(n()>=2) %>%
  do(data.frame(t(combn(.$clone_id, 2)))) %>% # elegant way
  mutate(source='CDR3b')

el2 <- rbind(edge_a_2, edge_b_2) %>% distinct()
g2 <- graph_from_edgelist(as.matrix(el2[,c('X1', 'X2')]), directed = F)

setdiff(data$clone_id, union(V(g1)$name, V(g2)$name))

g <- g1 + g2 + vertex(setdiff(data$clone_id, union(V(g1)$name, V(g2)$name)))
g <- simplify(g)

V_attr <- tibble(name=V(g)$name) %>%
  left_join(data %>% 
              select(clone_count, clone_id) %>% 
              group_by(clone_id) %>%
              summarise(clone_count = sum(clone_count)) %>%
              mutate(name=as.character(clone_id)),
            by = 'name') %>%
  mutate(
    size = case_when(
      name %in% c('EBV', 'CMV', 'Influenza', 'InfluenzaA') ~ 10,
      TRUE ~ log(clone_count)
    ),
    label = case_when(
      name %in% c('EBV', 'CMV', 'Influenza', 'InfluenzaA') ~ name,
      #clone_count > 10 ~ name,
      TRUE ~ NA
    ),
    color = case_when(
      name == 'EBV' ~ '#ff71ce',
      name == 'CMV' ~ '#01cdfe',
      name == 'Influenza' ~ '#05ffa1',
      name == 'InfluenzaA' ~ '#b967ff',
      TRUE ~ '#fffb96'
    )
  )

V(g)$size <- V_attr$size
V(g)$label <- V_attr$label
V(g)$color <- V_attr$color

pdf('figures/GLIPH_with_known_viral_clusters.pdf')
plot(g)
legend('topleft', pch=21, pt.cex=1,
  legend = c('expanded clone (>5)', 'EBV', 'CMV', 'Influenza', 'InfluenzaA'), 
  pt.bg = c('#fffb96', '#ff71ce', '#01cdfe', '#05ffa1', '#b967ff'))
dev.off()



