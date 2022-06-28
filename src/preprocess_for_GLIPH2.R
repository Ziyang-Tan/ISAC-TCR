library(dplyr)
library(readr)
library(stringr)

# load data
data_scTCR <- read_csv('data/scTCR_data_merge_vb.csv.gz')
info_bulk <- read_csv('data/bulk_info.csv.gz') %>%
  rename(ngi_id = `NGI ID`)
sample_id_map <- read_delim('data/sample_id_map.csv', delim = ';') %>%
  tidyr::separate(SC_ID, into = c('SampleID', NA), remove = F)
data_bulk <- read_csv('data/bulk_data.csv.gz') %>%
  left_join(info_bulk, by = 'ngi_id') %>%
  rename(Bulk_ID = `User ID`) %>%
  left_join(sample_id_map, by = 'Bulk_ID')
data_HLA <- readxl::read_excel(path = '/Users/tan/Library/CloudStorage/OneDrive-KI.SE/TCR_processed_data/HLA typing ISAC.xlsx')

# prepare GLIPH_input_TCR.txt

tmp_sc <- data_scTCR %>% 
  filter(cell_type == 'CD8T') %>%
  mutate(CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant,
         TRBV = gsub('\\*.*', '', TCR_Beta_Delta_V_gene_Dominant),
         TRBJ = gsub('\\*.*', '', TCR_Beta_Delta_J_gene_Dominant),
         CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
         subject = paste0(Sample_Name, ':NA')) %>%
  select(CDR3b, TRBV, TRBJ, CDR3a, subject) %>%
  group_by(CDR3b, TRBV, TRBJ, CDR3a, subject) %>%
  tally(name = 'count')
  #filter(!grepl('TRDV', TRBV)) # exclude gdT clone type

# prepare GLIPH_input_TCR.txt 
# mix visit

tmp_sc <- data_scTCR %>% 
  filter(cell_type == 'CD8T') %>%
  tidyr::separate(Sample_Name, into = c('individual', 'visit'), remove = F) %>%
  mutate(CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant,
         TRBV = gsub('\\*.*', '', TCR_Beta_Delta_V_gene_Dominant),
         TRBJ = gsub('\\*.*', '', TCR_Beta_Delta_J_gene_Dominant),
         CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
         subject = paste0(individual, ':NA')) %>%
  select(CDR3b, TRBV, TRBJ, CDR3a, subject) %>%
  group_by(CDR3b, TRBV, TRBJ, CDR3a, subject) %>%
  tally(name = 'count')
#filter(!grepl('TRDV', TRBV)) # exclude gdT clone type




# tmp_bulk <- data_bulk %>%
#   mutate(CDR3b = aaSeqCDR3,
#          TRBV = gsub('\\*.*', '', allVHitsWithScore),
#          TRBJ = gsub('\\*.*', '', allJHitsWithScore),
#          CDR3a = NA,
#          subject = paste0(SC_ID, ':NA'),
#          count = cloneFraction) %>%
#   select(CDR3b, TRBV, TRBJ, CDR3a, subject, count)

# GLIPH_input_TCR <- rbind(tmp_sc, tmp_bulk)
GLIPH_input_TCR <- tmp_sc
write_delim(GLIPH_input_TCR, file = 'data/GLIPH_input_TCR_mix_visit.txt', delim = '\t', col_names = F)

# prepare GLIPH_input_TCR.txt
# need to do something with HLA ambiguity?

GLIPH_input_HLA <- data_HLA %>%
  filter(SampleID %in% unique(info_bulk$Sample_Name)) %>%
  right_join(sample_id_map, by = 'SampleID') %>%
  mutate(A1 = paste0('A*', sub(':$', '', str_extract(paste0(A1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         A2 = paste0('A*', sub(':$', '', str_extract(paste0(A2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         B1 = paste0('B*', sub(':$', '', str_extract(paste0(B1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         B2 = paste0('B*', sub(':$', '', str_extract(paste0(B2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         C1 = paste0('C*', sub(':$', '', str_extract(paste0(C1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         C2 = paste0('C*', sub(':$', '', str_extract(paste0(C2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DRB11 = paste0('DRB1*', sub(':$', '', str_extract(paste0(DRB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DRB12 = paste0('DRB1*', sub(':$', '', str_extract(paste0(DRB12, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DPB11 = paste0('DPB1*', sub(':$', '', str_extract(paste0(DPB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DPB12 = paste0('DPB1*', sub(':$', '', str_extract(paste0(DPB12, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DQB11 = paste0('DQB1*', sub(':$', '', str_extract(paste0(DQB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DQB12 = paste0('DQB1*', sub(':$', '', str_extract(paste0(DQB12, ':'), '(\\d*:|[[:upper:]]*:){2}')))) %>%
  select(SC_ID, A1, A2, B1, B2, C1, C2, DRB11, DRB12, DPB11, DPB12, DQB11, DQB12)
write_delim(GLIPH_input_HLA, file = 'data/GLIPH_input_HLA.txt', delim = '\t', col_names = F)  

# mix visit
GLIPH_input_HLA <- data_HLA %>%
  filter(SampleID %in% unique(info_bulk$Sample_Name)) %>%
  #right_join(sample_id_map, by = 'SampleID') %>%
  mutate(A1 = paste0('A*', sub(':$', '', str_extract(paste0(A1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         A2 = paste0('A*', sub(':$', '', str_extract(paste0(A2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         B1 = paste0('B*', sub(':$', '', str_extract(paste0(B1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         B2 = paste0('B*', sub(':$', '', str_extract(paste0(B2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         C1 = paste0('C*', sub(':$', '', str_extract(paste0(C1, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         C2 = paste0('C*', sub(':$', '', str_extract(paste0(C2, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DRB11 = paste0('DRB1*', sub(':$', '', str_extract(paste0(DRB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DRB12 = paste0('DRB1*', sub(':$', '', str_extract(paste0(DRB12, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DPB11 = paste0('DPB1*', sub(':$', '', str_extract(paste0(DPB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DPB12 = paste0('DPB1*', sub(':$', '', str_extract(paste0(DPB12, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DQB11 = paste0('DQB1*', sub(':$', '', str_extract(paste0(DQB11, ':'), '(\\d*:|[[:upper:]]*:){2}'))),
         DQB12 = paste0('DQB1*', sub(':$', '', str_extract(paste0(DQB12, ':'), '(\\d*:|[[:upper:]]*:){2}')))) %>%
  select(SampleID, A1, A2, B1, B2, C1, C2, DRB11, DRB12, DPB11, DPB12, DQB11, DQB12)
write_delim(GLIPH_input_HLA, file = 'data/GLIPH_input_HLA_mix_visit.txt', delim = '\t', col_names = F)  

# notebook -----------------------------------------------------------------------------------------------------------------------
  
tmp <- data_scTCR %>% 
  filter(cell_type == 'CD8T') %>%
  tidyr::separate(Sample_Name, into = c('individual', 'visit'), remove = F)
  
write_csv(tmp, file = 'data/scTCR_data_merge_CD8T.csv')  

library(ggseqlogo)
library(ggplot2)
library(ggpubr)

tmp <- tmp %>% mutate(CDR3b_aa_len = nchar(TCR_Beta_Delta_CDR3_Translation_Dominant))

g1 <- ggplot(tmp, aes(x=CDR3b_aa_len)) + geom_bar() + xlim(5,25) + ggtitle('our CDR3b')
g2 <- ggplot() + geom_logo((tmp %>% filter(CDR3b_aa_len == 15))$TCR_Beta_Delta_CDR3_Translation_Dominant) + theme_logo() + 
  ggtitle('logotype of aa_length = 13')  + theme(legend.position = "none")


demo_data <- read_delim('data/demo_tcr.txt',  delim = '\t', col_names = c('CDR3b', 'V', 'J', 'CDR3a', 'subject', 'count')) %>% 
  mutate(CDR3b_aa_len = nchar(CDR3b))

g3 <- ggplot(demo_data, aes(x=CDR3b_aa_len)) + geom_bar() + xlim(5,25) + ggtitle('from demo_data')
g4 <- ggplot() + geom_logo((demo_data %>% filter(CDR3b_aa_len == 15))$CDR3b) + theme_logo() + 
  ggtitle('logotype of aa_length = 15') + theme(legend.position = "none")
ggarrange(ggarrange(g1,g2, ncol = 1, nrow = 2), ggarrange(g3,g4, ncol = 1, nrow = 2))
ggsave('figures/CDR3b_distribution_comparison.pdf')

