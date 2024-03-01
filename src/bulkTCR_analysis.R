library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

raw <- read_csv('data/bulk_data.csv.gz')
#info <- read_csv('data/bulk_info.csv.gz')

process_mixcr_raw <- function(raw){
  raw <- raw[,colSums(!is.na(raw))!=0]
  data <- raw %>% 
    mutate(V = sub('\\*.*$', '', allVHitsWithScore),
           D = sub('\\*.*$', '', allDHitsWithScore),
           J = sub('\\*.*$', '', allJHitsWithScore),
           C = sub('\\*.*$', '', allCHitsWithScore)) %>%
    rename(CDR3 = nSeqCDR3,
           CDR3aa = aaSeqCDR3) %>%
    # select(cloneCount, cloneFraction, V, D, J, C, CDR3, CDR3aa, `NGI ID`, `User ID`, Mreads, proj_id) %>%
    mutate(normCloneCount = round(cloneCount/Mreads))
  return(data)
}

count_V_fraction <- function(data) {
  res <- lapply(unique(data$`User ID`), function(cur_sample){
    sub_data <- data %>% filter(`User ID` == cur_sample)
    tmp <- sub_data %>% 
      group_by(V) %>% 
      summarise(freq = sum(cloneFraction))
    colnames(tmp) <- c('V', cur_sample)
    return(tmp)
  }) %>% Reduce(f=function(x,y) full_join(x,y,by='V'))
}

data <- process_mixcr_raw(raw)
res <- count_V_fraction(data)

rowOrder <- (res %>% rowwise() %>% mutate(mean = mean(c_across(contains('ISAC')), na.rm = T)) %>% arrange(desc(mean)))$V

res.long <- res %>%
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res), c('V')), names_to = 'sample_id', values_to = 'freq') %>%
  tidyr::separate(sample_id, into = c('subject', 'timepoint')) %>%
  mutate(timepoint = as.factor(timepoint))

res.long_sub <- res.long %>% filter(subject == 'ISAC81')
ggplot(res.long_sub, aes(x=V, y=freq, fill=timepoint)) +
  geom_bar(stat = 'identity', position="dodge") + 
  #geom_jitter(aes(group=timepoint), position=position_jitterdodge(), size=0.3, alpha=0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('ISAC81')


