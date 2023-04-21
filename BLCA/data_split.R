library(mlr3verse)
library(mlr3proba)

set.seed(42)
disease_code = 'BLCA'
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))

# verify all tasks have the same (time, status) outcomes
for (tsk_index in 1:(length(tasks) - 1)) {
  truth1 = tasks[[tsk_index]]$data(
    cols = tasks[[tsk_index]]$target_names
  )
  truth2 = tasks[[tsk_index + 1]]$data(
    cols = tasks[[tsk_index + 1]]$target_names
  )
  # ?all.equal.data.table
  res = all.equal(truth1, truth2, check.attributes = FALSE)
  if (class(res) == 'character') {
    stop('Target values must be the same! ', res)
  }
}

#' `mlr3proba:::partition.TaskSurv` is used by default which stratifies
#' the survival task by `status`
part = partition(tasks$Clinical, ratio = 0.8, stratify = TRUE)

saveRDS(part, file = paste0(disease_code, '/data/part.rds'))
