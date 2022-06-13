library(bigrquery)
library(dplyr)
library(tibble)
library(mmR)

# there're some problem with the HLA table on bigquery. Fix it when free

sql <- "SELECT * FROM `cradle-259115.HLA.typing`"
tb <- bq_project_query('cradle-259115', sql)
HLA <- bq_table_download(tb)
sql <- "SELECT * FROM `cradle-259115.HLA.appendix`"
tb <- bq_project_query('cradle-259115', sql)
HLA_full <- bq_table_download(tb)


HLA_full <- HLA_full %>% mutate(Included_Alleles = case_when(
  is.na(Included_Alleles) ~ Reporting_Allele,
  TRUE ~ Included_Alleles
))


HLA0201 <- HLA %>% 
  filter(grepl('^02:01', A1) | grepl('^02:01',A2))
#write_csv(HLA0201, 'HLA_0201.csv')

HLA_frequency <- function(HLA_dat, allele){
  # output a frequency table showing that in how many subjects a HLA type appears in 
  # at least one of the alleles
  # HLA_dat is HLA typing table from bigQuery
  # allele is a vector of allele names
  nameList <- HLA_dat %>% select(one_of(allele)) %>% unlist(use.names = F) %>% unique()
  table <- lapply(nameList, function(x){
    tmp <- HLA_dat %>% filter(get(allele[1]) == x | get(allele[2]) == x) %>% count()
    return(tmp)
  }) %>% bind_rows() %>% add_column(HLA = nameList) %>% arrange(desc(n))
  return(table)
}

DRB <- HLA_frequency(HLA0201, c('DRB11', 'DRB12'))
DPB <- HLA_frequency(HLA0201, c('DPB11', 'DPB12'))
DQB <- HLA_frequency(HLA0201, c('DQB11', 'DQB12'))

# select the samples with A02:01:01 B07:02:01 and DRB15:01:01
HLASelect <- HLA %>% 
  filter(grepl('^02:01', A1) | grepl('^02:01',A2)) %>%
  filter(grepl('^15:01', DRB11) | grepl('^15:01', DRB12)) %>%
  filter(grepl('^07:02', B1) | grepl('^07:02', B2))
#write_csv(HLASelect, 'HLA_select.csv')


# check celiac associated DQ2 and DQ8
# DQ2: DQB1*02:01 or 02:02
# DQ8: DQB1*03:02
HLA_full_search_pattern <- function(HLA_full, pattern, locus){
  d <- HLA_full %>%
    filter(Locus == locus) %>%
    filter(sapply(.$Included_Alleles %>% strsplit('/'), function(x){
      any(grepl(pattern, x))
    }))
  return(d)
}
DQ2 <- bind_rows(HLA_full_search_pattern(HLA_full, '^02:01', 'DQB1'),
                 HLA_full_search_pattern(HLA_full, '^02:02', 'DQB1'))
DQ8 <- HLA_full_search_pattern(HLA_full, '^03:02', 'DQB1')

mm.fastwrite(bind_rows(DQ2, DQ8), 'potential celiac.xlsx')






