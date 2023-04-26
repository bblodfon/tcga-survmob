#' => Make `mlr3` survival tasks using omics and clinical data
#' => Preprocess omics data with `ml3pipelines` and `survmob`
#' => Save preprocessed `mlr3` tasks
#'
#' Execute command: `Rscript OV/preprocess.R > OV/data/preprocess.log`
library(tidyverse)
library(survmob)
library(mlr3verse)
library(mlr3proba)
library(tictoc) # to measure the preprocessing time
set.seed(42)

disease_code = 'OV'

# Clinical ----
print('Clinical')
clin_tbl = readr::read_csv(file = paste0(disease_code, '/data/clinical.csv'))

clin_tbl2 = clin_tbl %>%
  select(-c('patient_id')) %>% # remove patient ID
  mutate(across(where(is.double), as.integer)) # integerize :)

#' convert to `TaskSurv`
clin_task = mlr3proba::as_task_surv(
  x = clin_tbl2,
  time = 'time', event = 'status', id = 'Clinical'
)
clin_task
print('Missing features:')
clin_task$missings()
print('Clinical variables:') # only 'age'
print(paste0(clin_task$feature_names, collapse = ', '))

# impute age by sampling from the empirical training data distribution
impute_hist = po('imputehist')
# center and scale numeric (not 0-1 encoded) features
po_scale = po('scale', affect_columns = selector_name(feature_names = 'age'))

pre = impute_hist %>>% po_scale

tic()
clin_task = pre$train(clin_task)[[1L]]
toc()

all(clin_task$missings() == 0) # check no missing features
print(paste0('#features remained: ', length(clin_task$feature_names)))

# number of events censoring rate
status = clin_task$truth()[,2]
print(paste0('Number of events: ', sum(status == 1)))
print(paste0('Censoring rate: ',
  floor(100 * sum(status==0)/length(status)), '%')) # 37%

# PipeOps ----
# remove features with more than 20% of NAs
po_na = po('removenas', cutoff = 0.2)
# remove features with more than 20% of zeros
po_zeros = po('removezeros', cutoff = 0.2)
#' Keep top 10000 most variant features (`NA`s are removed)
po_variance = po('filter', filter = flt('variance', na.rm = TRUE),
  filter.nfeat = 10000)
# Impute missing numeric features with xgboost regression learner
# (supports training data with missing features)
po_xgboost_impute = po('imputelearner',
  learner = lrn('regr.xgboost', nthread = 16,
    nrounds = 50, eta = 0.1, max_depth = 3),
  context_columns = selector_all(), # use all features for training
  affect_columns = selector_missing() # impute only features with missing data
)
# remove highly correlated features (>0.95 Pearson Correlation)
po_corrfeats = po('filter', filter = flt('find_correlation', method = 'pearson',
  use = 'pairwise.complete.obs'), filter.cutoff = 0.05)
# log2(x+1) transform
po_logtrans = po('logtransform')
# center and scale features
po_scale = po('scale')

# function to print statistics of an omics preprocessing pipeline
#' `pipeline` => mlr3pipelines::Graph object
print_stats = function(pipeline) {
  state = pipeline$state
  if (is.null(state)) {
    print('No state stored')
  } else {
    if (!is.null(state$remove_na)) {
      print(paste0('#features removed (NAs): ',
        state$remove_na$removed_column_num)
      )
    }
    if (!is.null(state$remove_zeros)) {
      print(paste0('#features removed (zeros): ',
        state$remove_zeros$removed_column_num)
      )
    }
    if (!is.null(state$variance)) {
      print(paste0('#features after variance filtering: ',
        length(state$variance$features))
      )
    }
    if (!is.null(state$imputelearner)) {
      print(paste0('#features imputed: ',
        length(state$imputelearner$affected_cols))
      )

    }
    if (!is.null(state$find_correlation)) {
      print(paste0('#features removed (pair-wise correlation): ',
        sum(state$find_correlation$scores <= 0.05))
      )
    }
  }
}

# Omics Preprocessing ----
## mRNA ----
print('mRNA')
mRNA = readr::read_csv(file = paste0(disease_code, '/data/mRNA.csv'), col_types = 'd')

colnames(mRNA) = make.names(colnames(mRNA)) # proper syntactic names
colnames(mRNA) = paste0('mRNA.', colnames(mRNA)) # omic prefix for feature uniqueness
mRNA_task = mlr3proba::as_task_surv(
  x = cbind(clin_task$data(cols = clin_task$target_names), mRNA),
  time = 'time', event = 'status', id = 'mRNA'
)

mRNA_ppl =
  po_na %>>%
  po_zeros %>>%
  po_variance %>>%
  po_xgboost_impute %>>%
  po_corrfeats %>>%
  po_logtrans %>>%
  po_scale

tic()
mRNA_task = mRNA_ppl$train(mRNA_task)[[1L]]
toc()
mRNA_task = survmob::minimize_backend(mRNA_task)

# check no missing features
stopifnot(all(mRNA_task$missings() == 0))
# Preprocessing stats
print_stats(mRNA_ppl)
print(paste0('#features remained: ', length(mRNA_task$feature_names)))

## CNA ----
print('CNA')
CNA = readr::read_csv(file = paste0(disease_code, '/data/CNA.csv'), col_types = 'd')

colnames(CNA) = make.names(colnames(CNA)) # proper syntactic names
colnames(CNA) = paste0('CNA.', colnames(CNA)) # omic prefix for feature uniqueness
CNA_task = mlr3proba::as_task_surv(
  x = cbind(clin_task$data(cols = clin_task$target_names), CNA),
  time = 'time', event = 'status', id = 'CNA'
)

CNA_ppl =
  po_na %>>%
  po_variance %>>%
  po_xgboost_impute %>>%
  po_scale

tic()
CNA_task = CNA_ppl$train(CNA_task)[[1L]]
toc()
CNA_task = survmob::minimize_backend(CNA_task)

# check no missing features
stopifnot(all(CNA_task$missings() == 0))
# Preprocessing stats
print_stats(CNA_ppl)
print(paste0('#features remained: ', length(CNA_task$feature_names)))

## miRNA ----
print('miRNA')
miRNA = readr::read_csv(file = paste0(disease_code, '/data/miRNA.csv'), col_types = 'd')

colnames(miRNA) = make.names(colnames(miRNA)) # proper syntactic names
colnames(miRNA) = paste0('miRNA.', colnames(miRNA)) # omic prefix for feature uniqueness
miRNA_task = mlr3proba::as_task_surv(
  x = cbind(clin_task$data(cols = clin_task$target_names), miRNA),
  time = 'time', event = 'status', id = 'miRNA'
)

miRNA_ppl =
  po_na %>>%
  po_zeros %>>%
  po_variance %>>%
  po_xgboost_impute %>>%
  po_corrfeats %>>%
  po_logtrans %>>%
  po_scale

tic()
miRNA_task = miRNA_ppl$train(miRNA_task)[[1L]]
toc()
miRNA_task = survmob::minimize_backend(miRNA_task)

# check no missing features
stopifnot(all(miRNA_task$missings() == 0))
# Preprocessing stats
print_stats(miRNA_ppl)
print(paste0('#features remained: ', length(miRNA_task$feature_names)))

## RPPA ----
print('RPPA')
RPPA = readr::read_csv(file = paste0(disease_code, '/data/RPPA.csv'), col_types = 'd')

colnames(RPPA) = make.names(colnames(RPPA)) # proper syntactic names
colnames(RPPA) = paste0('RPPA.', colnames(RPPA)) # omic prefix for feature uniqueness
RPPA_task = mlr3proba::as_task_surv(
  x = cbind(clin_task$data(cols = clin_task$target_names), RPPA),
  time = 'time', event = 'status', id = 'RPPA'
)

RPPA_ppl =
  po_na %>>%
  po_xgboost_impute %>>%
  po_corrfeats %>>%
  po_scale

tic()
RPPA_task = RPPA_ppl$train(RPPA_task)[[1L]]
toc()
RPPA_task = survmob::minimize_backend(RPPA_task)

# check no missing features
stopifnot(all(RPPA_task$missings() == 0))
# Preprocessing stats
print_stats(RPPA_ppl)
print(paste0('#features remained: ', length(RPPA_task$feature_names)))

## Methylation ----
print('Methylation')
Methyl = readr::read_csv(file = paste0(disease_code, '/data/methyl.csv'), col_types = 'd')

colnames(Methyl) = make.names(colnames(Methyl)) # proper syntactic names
colnames(Methyl) = paste0('meth.', colnames(Methyl)) # omic prefix for feature uniqueness
methyl_task = mlr3proba::as_task_surv(
  x = cbind(clin_task$data(cols = clin_task$target_names), Methyl),
  time = 'time', event = 'status', id = 'Methyl'
)

methyl_ppl =
  po_na %>>%
  po_variance %>>%
  po_xgboost_impute %>>%
  po_corrfeats %>>%
  po_logtrans %>>%
  po_scale

tic()
methyl_task = methyl_ppl$train(methyl_task)[[1L]]
toc()
methyl_task = survmob::minimize_backend(methyl_task)

# check no missing features
stopifnot(all(methyl_task$missings() == 0))
# Preprocessing stats
print_stats(methyl_ppl)
print(paste0('#features remained: ', length(methyl_task$feature_names)))

# Save tasks ----
tasks = list(
  Clinical = clin_task,
  mRNA = mRNA_task,
  miRNA = miRNA_task,
  CNA = CNA_task,
  RPPA = RPPA_task,
  Methyl = methyl_task
)
saveRDS(tasks, file = paste0(disease_code, '/data/tasks.rds'))
