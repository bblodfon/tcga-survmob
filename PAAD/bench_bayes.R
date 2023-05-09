#' Bayesian Analysis of Benchmarking Results
#' => Fit Bayesian linear mixed-effects models
library(survmob)
library(tidyverse)
library(rstanarm)

disease_code = 'PAAD'

# Reproducibility
seed = 42
set.seed(seed)

# Models tuned using Uno's C-index ----
## Benchmark data ----
mob_uno = readRDS(file = paste0(disease_code, '/bench/mob_uno.rds'))
cox = readRDS(file = paste0(disease_code, '/bench/cox.rds'))
bench_path = paste0(disease_code, '/bench/mob_uno')
dir.create(bench_path)

# merge results
res = dplyr::bind_rows(mob_uno$result, cox$result %>% select(!matches('model')))

# some checks
stopifnot(mob_uno$test_nrsmps == cox$test_nrsmps)
stopifnot(mob_uno$test_measure_ids == cox$test_measure_ids)
nrsmps = mob_uno$test_nrsmps
msr_ids = mob_uno$test_measure_ids
task_ids = unname(mlr3misc::map_chr(mob_uno$tasks, `[[`, 'id'))
print(task_ids) # for creating the omic columns

# reshape results
df = res %>%
  mutate(scores = map(boot_res, 'scores')) %>%
  select(!matches('boot_res')) %>%
  tibble::add_column(id = list(tibble(rsmp_id = 1:nrsmps))) %>%
  tidyr::unnest(cols = c(id, scores)) %>%
  dplyr::relocate(rsmp_id, .after = lrn_id) %>%
  tidyr::pivot_longer(cols = !matches('task_id|lrn_id|rsmp_id'),
    names_to = 'measure', values_to = 'value') %>%
  # One-hot encoded each omic to separate column
  mutate(
    Clinical_omic = as.integer(str_detect(task_id, 'Clinical')),
    mRNA_omic = as.integer(str_detect(task_id, 'mRNA')),
    miRNA_omic = as.integer(str_detect(task_id, 'miRNA')),
    CNA_omic = as.integer(str_detect(task_id, 'CNA')),
    Methyl_omic = as.integer(str_detect(task_id, 'Methyl'))
  )
#' `tibble` with a very general format of benchmarking results
df

# Inspect observed performance stats
# df %>%
#   filter(measure == 'harrell_c', task_id == 'Clinical') %>%
#   #filter(measure == 'harrell_c') %>%
#   group_by(lrn_id) %>%
#   summarise(average = mean(value), std = sd(value)) %>%
#   arrange(desc(average))

## Compare learners ----
for (msr_id in msr_ids) {
  message('\nMeasure: ', msr_id)
  df_sub = df %>%
    filter(measure == msr_id) %>%
    select(!matches('.*omic'))

  model = rstanarm::stan_glmer(
    data = df_sub,
    # MAIN FORMULA FOR MODEL COMPARISON
    formula = value ~ -1 + lrn_id + (1 | task_id/rsmp_id),
    chains = 4, cores = 4, iter = 5000, seed = seed
  )

  saveRDS(object = model,
    #' `mob_uno` => tuning by optimizing Uno's C-index
    #' `msr_id`  => which test measure we use to compare the models
    file = paste0(disease_code, '/bench/mob_uno/stan_glmer_', msr_id, '.rds')
  )
}

## Omic importance ----
omics = colnames(df)[grepl(pattern = 'omic', colnames(df))]
for (msr_id in msr_ids) {
  message('\nMeasure: ', msr_id)
  df_sub = df %>%
    filter(measure == msr_id) %>%
    # reduce rsmp_ids to half for faster model fit
    filter(rsmp_id %in% sample(1:nrsmps, nrsmps/2))

  # fit different model per omic
  model_list = list()
  for (omic in omics) {
    model = rstanarm::stan_glmer(
      data = df_sub,
      # MAIN FORMULA FOR PER OMIC IMPORTANCE (expected difference in
      # performance with vs without that omic)
      formula = as.formula(paste0('value ~ ', omic,
        ' + (1 | lrn_id) + (1 | task_id/rsmp_id)')),
      chains = 4, cores = 4, iter = 5000, seed = seed
    )

    model_list[[omic]] = model
  }

  saveRDS(object = model_list,
    #' `mob_uno` => tuning by optimizing Uno's C-index
    #' `msr_id`  => which test measure we used
    file = paste0(disease_code, '/bench/mob_uno/omic_fit_', msr_id, '.rds')
  )
}
