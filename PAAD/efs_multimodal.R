#' => Ensemble feature selection on all PAAD omics merged into a single
#' multi-modal task
#' => Use the C-index (OOB error) to optimize during feature selection
#' => We don't include the `aorsf` learner due to memory issues
#' => Run in parallel mode
library(mlr3verse)
library(survmob)
library(tictoc)
library(future)
library(progressr)

disease_code = 'PAAD'

# Reproducibility
set.seed(42)

# Logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3' )$set_threshold('warn')

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Data partition
part = readRDS(file = paste0(disease_code, '/data/part.rds'))

# Tasks
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))
po_featureunion = mlr3pipelines::po('featureunion')
task_multimodal = po_featureunion$train(tasks)[[1L]]
task_multimodal$id = paste0(mlr3misc::map_chr(tasks, `[[`, 'id'), collapse = '-')
task_multimodal

task_multimodal$filter(rows = part$train) # filter to train set
task_multimodal$col_roles$stratum = 'status' # Stratify

# eFS
lrn_ids = eFS$new()$supported_lrn_ids()
lrn_ids = lrn_ids[lrn_ids != 'aorsf'] # exclude aorsf

efs = eFS$new(msr_id = 'oob_error', resampling = rsmp('insample'), repeats = 1000,
  nthreads_rsf = 15, n_features = 5, subsample_ratio = 0.9, lrn_ids = lrn_ids,
  feature_fraction = 0.85)
#efs$run(task = task_multimodal)

plan('multicore', workers = 8)
efs$run_parallel(task = task_multimodal)

# Save result
file_name = paste0(disease_code, '/efs/multimodal_efs.rds')
saveRDS(efs, file = file_name)
