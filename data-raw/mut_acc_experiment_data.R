## code to prepare `DATASET` dataset goes here
library(tidyverse)

mut_acc_experiment_data = read_csv("data-raw/table-nb-occurences-and-mut-by-triplet-strain.csv") %>%
  rename(m_sc = `m_{s,c}`) %>%
  arrange(desc(n_c)) %>%
  mutate(context_id = n_c %>% desc %>% dense_rank())

usethis::use_data(mut_acc_experiment_data, overwrite = TRUE)


mut_acc_experiment_data_3610 = read.csv("data-raw/table-nb-occurences-and-mut-by-triplet-strain_DmutS3610.csv") %>%
  rename(m_sc = `m_.s.c.`) %>%
  filter(strain != "MMR-") %>%
  mutate(strain = gsub("DmutS3610", "MMR-", strain)) %>%
  arrange(desc(n_c)) %>%
  mutate(context_id = n_c %>% desc %>% dense_rank())

usethis::use_data(mut_acc_experiment_data_3610, overwrite = TRUE)


mut_acc_experiment_data_ts_tv_indel = read.csv("data-raw/table-nb-mut-ts-tv-indel_20230720.csv") %>%
  rename(m_sc = `m_.s.type.`) %>%
  # arrange(desc(n_c)) %>%
  # mutate(context_id = n_c %>% desc %>% dense_rank()) %>%
  mutate(context_id = mut.type %>%
           factor(levels = unique(mut.type)) %>%
           as.numeric())

usethis::use_data(mut_acc_experiment_data_ts_tv_indel, overwrite = TRUE)

mut_acc_experiment_data_ts_tv_3610 = read.csv("data-raw/table-nb-mut-ts-tv_2024_01_26.csv") %>%
  rename(m_sc = `m_.s.type.`) %>%
  # arrange(desc(n_c)) %>%
  # mutate(context_id = n_c %>% desc %>% dense_rank()) %>%
  filter(strain != "MMR-") %>%
  mutate(strain = gsub("DmutS3610", "MMR-", strain)) %>%
  mutate(context_id = mut.type %>%
           factor(levels = unique(mut.type)) %>%
           as.numeric())

usethis::use_data(mut_acc_experiment_data_ts_tv_3610, overwrite = TRUE)

input_data = read.delim(file = "data-raw/SAT_mcmc_chain_labmut_R3610MMR-3610_inputfile.csv", sep = " ")
minimal_input_data_onestrain = input_data %>% 
  filter(strain == first(strain)) %>% 
  rename(mutation_id = context_id,
         m = m_sc, 
         n = n_c,
         t = t_s) %>% 
  select(mutation_id, m, n, t) 
usethis::use_data(minimal_input_data_onestrain, overwrite = TRUE)
