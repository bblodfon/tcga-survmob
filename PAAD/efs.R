#' => Ensemble feature selection on PAAD omics (tasks)
#' => Use separately the C-index (OOB error) and RCLL measure to optimize
#' during feature selection
library(mlr3verse)
library(survmob)
library(tictoc)

disease_code = 'PAAD'

# Reproducibility
set.seed(42)

# Logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3' )$set_threshold('warn')

# Tasks and data partition
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))
part  = readRDS(file = paste0(disease_code, '/data/part.rds'))

# eFS variables
nthreads_rsf = 30
n_features = 5
repeats = 100

# helper function
run_efs = function(task, part, msr_id, rs, repeats, nthreads_rsf, n_features,
  lrn_ids) {
  message('Omic: ', task$id)

  task2 = task$clone()
  task2$filter(rows = part$train)
  task2$col_roles$stratum = 'status' # Stratify

  efs = eFS$new(msr_id = msr_id, resampling = rs, repeats = repeats,
    nthreads_rsf = nthreads_rsf, n_features = n_features, lrn_ids = lrn_ids)
  efs$run(task = task2)

  efs
}

# Measure: OOB (1 - C-index) ----
message('### OOB error ###')
msr_id = 'oob_error'
rs = rsmp('insample')
lrn_ids = eFS$new()$supported_lrn_ids()

data_list = list()
for (task in tasks) {
  if (task$id != 'Clinical') { # efs only for omic data
    tic()
    efs = run_efs(task = task, part = part, msr_id = msr_id, rs = rs,
      nthreads_rsf = nthreads_rsf, n_features = n_features, repeats = repeats,
      lrn_ids = lrn_ids)
    toc()

    data_list[[task$id]] = efs

    # temporary save just in case:
    file_name = paste0(disease_code, '/efs/cindex_efs_', task$id, '.rds')
    saveRDS(efs, file = file_name)
  }
}

saveRDS(data_list, file = paste0(disease_code, '/efs/cindex_efs.rds'))

# Measure: RCLL ----
message('### RCLL ###')
msr_id = 'rcll'
rs = rsmp('cv', folds = 4)
#' exclude `aorsf` for RCLL (C++ issues)
lrn_ids = lrn_ids[lrn_ids != 'aorsf']

data_list = list()
for (task in tasks) {
  if (task$id != 'Clinical') { # efs only for omic data
    tic()
    efs = run_efs(task = task, part = part, msr_id = msr_id, rs = rs,
      nthreads_rsf = nthreads_rsf, n_features = n_features, repeats = repeats,
      lrn_ids = lrn_ids)
    toc()

    data_list[[task$id]] = efs

    # temporary save just in case:
    file_name = paste0(disease_code, '/efs/rcll_efs_', task$id, '.rds')
    saveRDS(efs, file = file_name)
  }
}

saveRDS(data_list, file = paste0(disease_code, '/efs/rcll_efs.rds'))
