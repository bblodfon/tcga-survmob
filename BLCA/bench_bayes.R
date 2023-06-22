#' Bayesian Analysis of Benchmarking Results
#' => Fit Bayesian linear mixed-effects models
#' => Compare models and omics
#' => Run in parallel
#' => Execute: `Rscript BLCA/bench_bayes.R`
library(survmob)
library(tidyverse)
library(rstanarm)
library(future.apply)
library(progressr)

# TCGA study
disease_code = 'BLCA'

# Reproducibility
seed = 42

# Progress bars
#' Note: no output bars when using same number of workers and measures (e.g. 4)
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Parallelization
#' Note: Total maximum cores used is `workers` x `n_cores`
plan('multicore', workers = 4) # across performance metrics
n_chains = 4 # MCMC chains
n_iters = 5000 # MCMC iterations per chain
n_cores = 4 # one core per chain

# Helper functions ----

#' @description
#' Reshape the output result object from [survmob::MOBenchmark]
#'
#' @param res a `tibble` with columns `task_id`, `lrn_id` and `boot_res`
#'
#' @return `tibble` with columns `task_id`, `lrn_id`, `rsmp_id`, `measure`,
#' `value` and one column per omic/modality, which has values 1 or 0
#' depending on the presence or absence of the omic name in the `task_id` column.
reshape_mob_res = function(res) {
  # check number of test bootstrap resamplings
  nrsmps_vec = mlr3misc::map_dbl(res$boot_res, `[[`, 'test_nrsmps')
  if (length(unique(nrsmps_vec)) != 1) {
    stop('Number of test bootstrap resamplings different?!')
  }
  nrsmps = res$boot_res[[1]]$test_nrsmps

  df = res %>%
    mutate(scores = purrr::map(boot_res, 'scores')) %>%
    select(!matches('boot_res')) %>%
    tibble::add_column(id = list(tibble(rsmp_id = 1:nrsmps))) %>%
    tidyr::unnest(cols = c(id, scores)) %>%
    dplyr::relocate(rsmp_id, .after = lrn_id) %>%
    tidyr::pivot_longer(cols = !matches('task_id|lrn_id|rsmp_id'),
      names_to = 'measure', values_to = 'value') %>%
    # exclude `dcal` measure (wide range of positive values, posterior draws
    # negative values!)
    filter(measure != 'dcal')

  # Extract omic names
  omics = unique(unlist(strsplit(df$task_id, '-')))

  # Iterate over omics and create new columns
  for (omic in omics) {
    omic_column = paste0(omic, '_omic')
    df = df %>%
      mutate(!!omic_column := as.integer(str_detect(task_id, omic)))
  }

  df
}

#' @description
#' Fit bayesian models to compare multiple learners
#'
#' @param df a `data.frame`/`tibble` object with the benchmarking multi-omics
#' results.
#' Columns names must be: `task_id`, `lrn_id`, `rsmp_id`, `measure`, `value`
#' and one-hot encoded columns that have the substring `omic`.
#' See `reshape_mob_res`.
#' @param n_chains Number of Markon chains for MCMC. Default = 4
#' @param n_iters Number of iterations per chain. Default = 5000
#' @param n_cores Number of cores for MCMC parallelization. Default = 4
#' @param seed seed number for MCMC
#' @param bench_path directory to store the bayesian models
fit_bayes_model_cmp = function(df, n_chains = 4, n_iters = 5000, n_cores = 4,
  seed, bench_path) {
  message('\n### Run Bayesian models to compare learners')

  # Measures - one model per measure
  msr_ids = unique(df$measure)
  message('Measures: ', paste0(msr_ids, collapse = ', '))

  # progress bar - parallelize across measures
  pb = progressr::progressor(steps = length(msr_ids))

  y = future.apply::future_lapply(1:length(msr_ids), function(i) {
    set.seed(i)

    msr_id = msr_ids[i]
    pb(sprintf('%s', msr_id))

    df_sub = df %>%
      filter(measure == msr_id) %>%
      select(!matches('.*omic'))

    model = rstanarm::stan_glmer(
      data = df_sub,
      # MAIN FORMULA FOR MODEL COMPARISON
      formula = value ~ -1 + lrn_id + (1 | task_id/rsmp_id),
      chains = n_chains, cores = n_cores, iter = n_iters, seed = seed
    )

    saveRDS(
      object = model,
      file = paste0(bench_path, '/models_fit_', msr_id, '.rds')
    )
  }, future.seed = TRUE)
}

#' @description
#' Fit bayesian models to compare omic importance
#'
#' @param df a `data.frame`/`tibble` object with the benchmarking multi-omics
#' results.
#' Columns names must be: `task_id`, `lrn_id`, `rsmp_id`, `measure`, `value`
#' and one-hot encoded columns that have the substring `omic`.
#' See `reshape_mob_res`.
#' @param n_chains Number of Markon chains for MCMC. Default = 4
#' @param n_iters Number of iterations per chain. Default = 5000
#' @param n_cores Number of cores for MCMC parallelization. Default = 4
#' @param seed seed number for MCMC. Default = 42
#' @param nrsmps_perc Percentage of the resampling ids to keep for faster
#' bayesian model fit (number between 0 and 1).
#' Default = 0.5, so keep half of the resampling ids (the same ones for all
#' models to be fitted)
#' @param bench_path directory to store the bayesian models
fit_bayes_omic_cmp = function(df, n_chains = 4, n_iters = 5000, n_cores = 4,
  seed, nrsmps_perc = 0.5, bench_path) {
  # small check
  stopifnot(nrsmps_perc > 0, nrsmps_perc <= 1)

  message('\n### Run Bayesian models to compare omics')

  # reducing rsmp_ids for faster model fit
  set.seed(seed)
  rsmp_ids = unique(df$rsmp_id)
  nrsmps = length(rsmp_ids)
  df = df %>%
    filter(rsmp_id %in% sample(x = rsmp_ids, size = nrsmps * nrsmps_perc))
  message('Reducing to ', nrsmps * nrsmps_perc, ' rsmp_ids')

  # Measures
  msr_ids = unique(df$measure)
  message('Measures: ', paste0(msr_ids, collapse = ', '))

  # progress bar - parallelize across measures
  pb = progressr::progressor(steps = length(msr_ids))

  # omic column names (must have `omic`)
  omics = colnames(df)[grepl(pattern = 'omic', colnames(df))]

  y = future.apply::future_lapply(1:length(msr_ids), function(i) {
    set.seed(i)

    msr_id = msr_ids[i]
    pb(sprintf('%s', msr_id))

    df_sub = df %>%
      filter(measure == msr_id)

    # fit different model per omic
    model_list = list()
    for (omic in omics) {
      #message('### ', omic, ' ###')
      model = rstanarm::stan_glmer(
        data = df_sub,
        # MAIN FORMULA FOR PER OMIC IMPORTANCE (expected difference in
        # performance with vs without that omic)
        formula = as.formula(paste0('value ~ ', omic,
          ' + (1 | lrn_id) + (1 | task_id/rsmp_id)')),
        chains = n_chains, cores = n_cores, iter = n_iters, seed = seed
      )

      model_list[[omic]] = model
    }

    saveRDS(
      object = model_list,
      file = paste0(bench_path, '/omic_fit_', msr_id, '.rds')
    )
  }, future.seed = TRUE)
}

# Models tuned using Uno's C-index ----
message('\n### Uno\'s C-index ###')

## Benchmark data ----
mob_uno = readRDS(file = paste0(disease_code, '/bench/mob_uno.rds'))
bench_path_uno = paste0(disease_code, '/bench/mob_uno')
dir.create(bench_path_uno)

#' `tibble` with a very general format of benchmarking results
df_uno = reshape_mob_res(res = mob_uno$result)
df_uno

## Compare learners ----
fit_bayes_model_cmp(df = df_uno, n_chains = n_chains, n_iters = n_iters,
  n_cores = n_cores, seed = seed, bench_path = bench_path_uno)

## Omic importance ----
fit_bayes_omic_cmp(df = df_uno, n_chains = n_chains, n_iters = n_iters,
  n_cores = n_cores, seed = seed, nrsmps_perc = 0.5, bench_path = bench_path_uno)

# Models tuned using RCLL ----
message('\n### RCLL ###')

## Benchmark data ----
mob_rcll = readRDS(file = paste0(disease_code, '/bench/mob_rcll.rds'))
bench_path_rcll = paste0(disease_code, '/bench/mob_rcll')
dir.create(bench_path_rcll)

#' `tibble` with a very general format of benchmarking results
df_rcll = reshape_mob_res(res = mob_rcll$result)
df_rcll

## Compare learners ----
fit_bayes_model_cmp(df = df_rcll, n_chains = n_chains, n_iters = n_iters,
  n_cores = n_cores, seed = seed, bench_path = bench_path_rcll)

## Omic importance ----
fit_bayes_omic_cmp(df = df_rcll, n_chains = n_chains, n_iters = n_iters,
  n_cores = n_cores, seed = seed, nrsmps_perc = 0.5, bench_path = bench_path_rcll)
