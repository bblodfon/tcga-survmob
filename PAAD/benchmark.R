library(survmob)
library(mlr3verse)
library(mlr3mbo)
library(mlr3proba)
library(progressr)

# cancer TCGA study
disease_code = 'PAAD'

# create benchmark results dir
bench_path = paste0(disease_code, '/bench')
dir.create(bench_path)

# Logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Reproducibility
set.seed(42)

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Tasks, stratify by status
tasks = readRDS(file = paste0(disease_code, '/data/tasks_cindex.rds'))
for (task in tasks) task$col_roles$stratum = 'status'

# Data partition (train/test)
part = readRDS(file = paste0(disease_code, '/data/part.rds'))

# Learners
lrn_ids = c('coxph', 'coxnet', 'rsf_cindex', 'rsf_logrank',
  'rsf_maxstat', 'rsf_extratrees', 'aorsf', 'coxboost', 'xgboost_cox_early',
  'xgboost_aft_early')

# Test Measures
test_measure_ids = c('harrell_c', 'uno_c', 'ibrier_erv', 'rcll_erv', 'dcal')

# Create benchmark objects
mob_uno = MOBenchmark$new(
  tasks = tasks, part = part,
  gen_task_powerset = TRUE,
  lrn_ids = lrn_ids, nthreads_rsf = 30, nthreads_xgb = 6,
  tune_measure_id = 'uno_c', tune_nevals = 100,
  test_nrsmps = 100, test_workers = 20,
  test_measure_ids = test_measure_ids,
  tune_rsmp = rsmp('cv', folds = 5),
  quiet = FALSE, keep_models = TRUE
)

# Execute benchmark
mob_uno$run()

# Save result object
#mob_uno$drop_tasks()
saveRDS(mob_uno, file = paste0(bench_path, '/mob_uno.rds'))
mob_uno$drop_models() # reduce size by A lot :)
saveRDS(mob_uno, file = paste0(bench_path, '/mob_uno_light.rds'))