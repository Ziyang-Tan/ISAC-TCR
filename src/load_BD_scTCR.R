library(dplyr)
library(stringr)
library(readr)

BD_load_VDJ <- function(x, dir_path){
  path <- Sys.glob(paste0(dir_path, '/*/*/*', x, '_VDJ_perCell.csv'))
  read_csv(path, skip = 7, show_col_types = FALSE) %>%
    mutate(at_least_one_chain = !is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant) | 
             !is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
    mutate(is_gdT = str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}") == 'TRGV') %>%
    mutate(is_gdT = case_when(
      is.na(is_gdT) ~ FALSE,
      TRUE ~ is_gdT
    )) %>%
    mutate(proj_id = x) %>%
    mutate(unique_index = paste0(x, '-', Cell_Index)) %>%
    select(-Cell_Index)
}

BD_load_sample_tag <- function(x, dir_path){
  path <- Sys.glob(paste0(dir_path, '/*/*/*', x, '_Sample_Tag_Calls.csv'))
  if (identical(path, character(0))){
    print(paste0('no sample tag file found for ', x))
  } else {
    read_csv(path, skip = 7, show_col_types = FALSE) %>%
      mutate(unique_index = paste0(x, '-', Cell_Index)) %>%
      select(-Cell_Index)
  }
}

BD_load_gene_exp <- function(x, dir_path, norm_method){
  path <- Sys.glob(paste0(dir_path, '/*/*/*', x, '_', norm_method, '_ReadsPerCell.csv'))
  read_csv(path, skip = 7, show_col_types = FALSE) %>%
    mutate(unique_index = paste0(x, '-', Cell_Index)) %>%
    select(-Cell_Index)
}