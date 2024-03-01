


# ISAC02 ---------------
# trend determination

patient <- 'ISAC02'
sub_name <- 'CD8T'
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')
g <- clone_expansion_alluvium(patient, clone_exp_sub, clone_label = T) + labs(title = paste0(patient, '_', sub_name)) + guides(fill='none', color='none')

g2 <- trend_determination_plot(patient, clone_exp_sub, top_mod='union')
# ISAC02
# increase_clone <- c(31042, 24360, 1709, 29062, 33659, 21963)
# decrease_clone <- c(4834, 19762, 32933, 20882, 2490, 26712, 11642, 32213, 26805, 
#                     30159, 12573, 11556, 30108, 20878, 27993, 20975, 34986)
# ISAC99
#increase_clone <- c(489, 2881, 8580)
#decrease_clone <- c(7517, 20443, 21102)

increase_clone <- c(48795, 37805, 2890, 45479, 53017, 33890)
decrease_clone <- c(7918, 30260, 51808, 32056, 4293, 41385, 18961, 50670, 41552, 
                    47305, 20424, 18818, 47227, 32052, 43590, 32214, 55081)



# flow by trend
df <- clone_exp_alluvium_preparation(patient, clone_exp_sub, top_mod='union')
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
  labs(title = patient)

obj_sub <- subset(obj, subset = patient_id == patient & 
                    clone_id %in% union(increase_clone, decrease_clone) &
                    cell_type == sub_name)
tmp <- obj_sub[[]] %>% left_join(df %>% select(clone_id, trend) %>% distinct(), by='clone_id') %>% select(trend)
rownames(tmp) <- obj_sub[[]]$unique_index
obj_sub$trend <- tmp
Idents(obj_sub) <- 'trend'
res <- FindMarkers(obj_sub, ident.1 = 'increase', ident.2 = 'decrease',
                   logfc.threshold = 0)

g4 <- EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
                      pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
                      selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
                      title = paste0(sub_name, ' of ', patient),
                      subtitle = 'increase vs decrease')

ggarrange(g, g3, g4, ncol = 3, common.legend=T) %>% ggexport(filename = file.path('figures', 'patient_specific_analysis', patient, paste0(patient, '_trend_RNA.pdf')), width = 20, height = 8)

# add some GO?




