#' Subset the `mlr3` survival tasks/datasets used during the execution
#' of the `eFS` ensemble feature selection algorithm.
#' The most *frequently selected* features are used per omic based on some
#' manual cutoffs.
#'
#' We analyze the results from:
#' => `efs.R`: includes optimization of the C-index (OOB error) and the RCLL
#' measure during feature selection *separately* on each omic
#' => `efs_multimodal.R`: includes optimization of the C-index (OOB error) on
#' a multimodal task (all omics *combined together* into a single task)
library(survmob)
library(mlr3pipelines)
library(tidyverse)

disease_code = 'BLCA'

# Tasks and eFS results ----
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))
efs_cindex = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
efs_rcll   = readRDS(file = paste0(disease_code, '/efs/rcll_efs.rds'))
#' Only run after executing `efs_multimodal.R`
efs_multimodal = readRDS(file = paste0(disease_code, '/efs/multimodal_efs.rds'))

# top 20 features
ts = TaskSub$new(nfeats = 20, max_nfeats = 30, cutoff = 0.4, perc = 0.01)

omics = names(tasks)[!names(tasks) %in% 'Clinical']
omic_tasks = tasks[omics] # with all features (max 10000)

# Subset tasks (C-index) ----
message('Subset tasks (C-index)\n')
stopifnot(names(efs_cindex) == omics) # same order
print(ts$summary_stats(eFSlist = efs_cindex, tasks = omic_tasks))
## get omic tasks with reduced features
omic_tasks_cindex = ts$subset_tasks(eFSlist = efs_cindex, tasks = omic_tasks, method = 'top_n')
tasks_cindex = c(tasks['Clinical'], omic_tasks_cindex)
saveRDS(object = tasks_cindex, file = paste0(disease_code, '/data/tasks_cindex.rds'))

# Subset tasks (RCLL) ----
if (FALSE) {
  message('Subset tasks (RCLL)\n')
  stopifnot(names(efs_rcll) == omics) # same order
  print(ts$summary_stats(eFSlist = efs_rcll, tasks = omic_tasks))
  ## get omic tasks with reduced features
  omic_tasks_rcll = ts$subset_tasks(eFSlist = efs_rcll, tasks = omic_tasks, method = 'top_n')
  tasks_rcll = c(tasks['Clinical'], omic_tasks_rcll)
  saveRDS(object = tasks_rcll, file = paste0(disease_code, '/data/tasks_rcll.rds'))
}

# Subset multimodal task (C-index) ----
#' Only run after executing `efs_multimodal.R`
message('Subset multimodal task\n')

nfeats = sapply(tasks, function(task) {length(task$feature_names)})
message('#Features per omic:')
print(nfeats)
prop_omic = nfeats/sum(nfeats)
prop_omic_log = log(nfeats)/sum(log(nfeats))

# Proportion of omic features to the total multimodal feature set
# It's like a probability of a random chosen feature belonging to a particular omic
message('%Features per omic:')
100*prop_omic
message('%Features per omic (log scaled):')
100*prop_omic_log

# consensus ranked features (from all omics)
fs_tbl = efs_multimodal$fs_stats()$consensus

message('%Total features selected per omic:')
sel_feats = fs_tbl$feat_name
total_sel_feats = length(sel_feats)
100*
  c(Clinical = sum(str_detect(string = sel_feats, pattern = 'mRNA|miRNA|CNA|meth', negate = TRUE))/total_sel_feats,
  mRNA = sum(str_detect(string = sel_feats, pattern = 'mRNA\\.'))/total_sel_feats,
  miRNA = sum(str_detect(string = sel_feats, pattern = 'miRNA\\.'))/total_sel_feats,
  CNA = sum(str_detect(string = sel_feats, pattern = 'CNA\\.'))/total_sel_feats,
  Methyl = sum(str_detect(string = sel_feats, pattern = 'meth\\.'))/total_sel_feats
)

#' We need to *correct the selection frequency* to account for the
#' difference in number of features between the different omics.
#' We scale the multimodal selection frequency as: `freq ^ prop_omic`
#' => The lower the `prop_omic`, the higher the scaled frequency
fs_tbl = fs_tbl %>%
  mutate(freq_scaled = case_when(
    str_detect(string = feat_name, pattern =  'mRNA\\.')  ~ freq ^ prop_omic['mRNA'],
    str_detect(string = feat_name, pattern =  'miRNA\\.') ~ freq ^ prop_omic['miRNA'],
    str_detect(string = feat_name, pattern =  'CNA\\.')   ~ freq ^ prop_omic['CNA'],
    str_detect(string = feat_name, pattern =  'meth\\.')  ~ freq ^ prop_omic['Methyl'],
    TRUE ~ freq ^ prop_omic['Clinical']
  ), freq_scaled_log = case_when(
    str_detect(string = feat_name, pattern =  'mRNA\\.')  ~ freq ^ prop_omic_log['mRNA'],
    str_detect(string = feat_name, pattern =  'miRNA\\.') ~ freq ^ prop_omic_log['miRNA'],
    str_detect(string = feat_name, pattern =  'CNA\\.')   ~ freq ^ prop_omic_log['CNA'],
    str_detect(string = feat_name, pattern =  'meth\\.')  ~ freq ^ prop_omic_log['Methyl'],
    TRUE ~ freq ^ prop_omic_log['Clinical']
  ))

#' Use the 3 selection frequency values as ranking scores to select the top
#' 100 multimodal features and calculate the number of features per omic
topn = 100
feat_stat_list = list()
for (freq_val in c('freq', 'freq_scaled', 'freq_scaled_log')) {
  topn_feats = fs_tbl %>%
    arrange(desc(!!rlang::sym(freq_val))) %>%
    pull(feat_name) %>%
    (function(x) x[1:topn])

  feat_stat_list[[freq_val]] = tibble(
    freq_value = freq_val,
    topn_feature_set = list(topn_feats),
    mRNA     = length(grep(x = topn_feats, pattern = 'mRNA')),
    miRNA    = length(grep(x = topn_feats, pattern = 'miRNA')),
    CNA      = length(grep(x = topn_feats, pattern = 'CNA')),
    Methyl   = length(grep(x = topn_feats, pattern = 'meth')),
    Clinical = length(grep(x = topn_feats, pattern = 'mRNA|miRNA|CNA|meth', invert = TRUE))
  )
}
feat_stat_tbl = bind_rows(feat_stat_list)
feat_stat_tbl
#' observe the bias:
#' - `freq` => more features from omics with large p
#' - `freq_scaled` => more features from omics with small p
#' - `freq_scaled_log` => a good compromise I think :)

#' the multimodal task used in `efs_multimodal.R`
po_featureunion = mlr3pipelines::po('featureunion')
task_multimodal = po_featureunion$train(tasks)[[1L]]
task_multimodal$id = 'Multimodal-eFS'
task_multimodal

# subset multimodal task
features =
  feat_stat_tbl %>%
  filter(freq_value == 'freq_scaled_log') %>%
  pull(topn_feature_set) %>%
  `[[`(1)
task = minimize_backend(task_multimodal$clone()$select(cols = features))
task

saveRDS(object = task, file = paste0(disease_code, '/data/task_multimodal.rds'))
