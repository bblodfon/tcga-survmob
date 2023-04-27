#' Subset `mlr3` survival tasks to the most stable/robust features found
#' using the `eFS` ensemble feature selection algorithm results based on two
#' different measures for optimizing the feature selection (C-index vs RCLL)
library(survmob)

disease_code = 'PAAD'

# Tasks and eFS results ----
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))
efs_cindex = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
efs_rcll   = readRDS(file = paste0(disease_code, '/efs/rcll_efs.rds'))

ts = TaskSub$new(nfeats = 20, max_nfeats = 30, cutoff = 0.4, perc = 0.01)

omics = names(tasks)[!names(tasks) %in% 'Clinical']
omic_tasks = tasks[omics] # with all features (max 10000)

# Subset tasks (C-index) ----
stopifnot(names(efs_cindex) == omics) # same order
print(ts$summary_stats(eFSlist = efs_cindex, tasks = omic_tasks))
## get omic tasks with reduced features
omic_tasks_cindex = ts$subset_tasks(eFSlist = efs_cindex, tasks = omic_tasks, method = 'top_n')
tasks_cindex = c(tasks['Clinical'], omic_tasks_cindex)
saveRDS(object = tasks_cindex, file = paste0(disease_code, '/data/tasks_cindex.rds'))

# Subset tasks (RCLL) ----
stopifnot(names(efs_rcll) == omics) # same order
print(ts$summary_stats(eFSlist = efs_rcll, tasks = omic_tasks))
## get omic tasks with reduced features
omic_tasks_rcll = ts$subset_tasks(eFSlist = efs_rcll, tasks = omic_tasks, method = 'top_n')
tasks_rcll = c(tasks['Clinical'], omic_tasks_rcll)
saveRDS(object = tasks_rcll, file = paste0(disease_code, '/data/tasks_rcll.rds'))
