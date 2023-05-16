#' => Ensemble feature selection on all PAAD omics merged into a single
#' multi-modal task
#' => Use the C-index (OOB error) to optimize during feature selection
#' => Execute command: `Rscript --max-ppsize=100000 PAAD/efs_multimodal.R`
#' (aorsf fails with `protect(): protection stack overflow` for large datasets)
library(mlr3verse)
library(survmob)
library(tictoc)

disease_code = 'PAAD'

# Reproducibility
set.seed(42)

# Logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3' )$set_threshold('warn')

# Data partition
part  = readRDS(file = paste0(disease_code, '/data/part.rds'))

# Tasks
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))
po_featureunion = mlr3pipelines::po('featureunion')
task_multimodal = po_featureunion$train(tasks)[[1L]]
task_multimodal$id = paste0(mlr3misc::map_chr(tasks, `[[`, 'id'), collapse = '-')
task_multimodal

task_multimodal$filter(rows = part$train)
task_multimodal$col_roles$stratum = 'status' # Stratify

# eFS
lrn_ids = eFS$new()$supported_lrn_ids()

efs = eFS$new(msr_id = 'oob_error', resampling = rsmp('insample'), repeats = 500,
  nthreads_rsf = 30, n_features = 5, subsample_ratio = 0.9, lrn_ids = lrn_ids,
  feature_fraction = 0.8)
efs$run(task = task_multimodal)

# Save result
file_name = paste0(disease_code, '/efs/multimodal_efs.rds')
saveRDS(efs, file = file_name)
