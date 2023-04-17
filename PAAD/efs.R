#' => Ensemble feature selection on PAAD omics (tasks)
#' => Use separately the C-index (OOB error) and RCLL measure to optimize
#' during feature selection
library(tidyverse)
library(mlr3verse)
library(survmob)
library(tictoc)

# Reproducibility
set.seed(42)

# Logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3' )$set_threshold('warn')

# Tasks and partition
tasks = readRDS(file = 'PAAD/data/tasks.rds')
part  = readRDS(file = 'PAAD/data/part.rds')

# eFS variables
nthreads_rsf = 14
n_features = 5
repeats = 1

# helper function
execute_efs = function(task, part, msr_id, rs, nthreads_rsf, n_features, repeats) {
  print(task$id)

  task2 = task$clone()
  task2$filter(rows = part$train)
  task2$col_roles$stratum = 'status'

  efs = eFS$new(msr_id = msr_id, resampling = rs, repeats = repeats,
    nthreads_rsf = nthreads_rsf, n_features = n_features)
  efs$run(task = task2)

  efs
}

# Measure: OOB (1 - C-index) ----
msr_id = 'oob_error'
rs = rsmp('insample')

data_list = list()
for (task in tasks) {
  if (task$id != 'Clinical') { # efs only for omic data
    tic()
    efs = execute_efs(task = task, part = part, msr_id = msr_id, rs = rs,
      nthreads_rsf = nthreads_rsf, n_features = n_features, repeats = repeats
    )
    toc()

    data_list[[task$id]] = efs

    # temporary save just in case:
    file_name = paste0('PAAD/efs/cindex_efs_', task$id, '.rds')
    saveRDS(efs, file = file_path)
  }
}

saveRDS(data_list, file = 'PAAD/efs/cindex_efs.rds')

# Measure: RCLL ----
msr_id = 'rcll'
rs = rsmp('cv', folds = 5)

data_list = list()
for (task in tasks) {
  if (task$id != 'Clinical') { # efs only for omic data
    efs = execute_efs(task = task, part = part, msr_id = msr_id, rs = rs,
      nthreads_rsf = nthreads_rsf, n_features = n_features, repeats = repeats
    )

    data_list[[task$id]] = efs

    # temporary save just in case:
    file_name = paste0('PAAD/efs/rcll_efs_', task$id, '.rds')
    saveRDS(efs, file = file_path)
  }
}

saveRDS(data_list, file = 'PAAD/efs/rcll_efs.rds')
