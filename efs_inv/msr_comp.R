#' Compare the two measures (RCLL vs C-index) used for optimizing
#' Ensemble feature selection (eFS) algorithm in terms of produced
#' top N ranked features
library(survmob)
library(stabm)
library(tidyverse)
library(future.apply)

# set up parallel backend
plan('multicore', workers = 15)

# Reproducibility ----
set.seed(42)

# Helper functions ----
# Jaccard similarity between two vectors
jaccard = function(a, b) {
  a = unique(a)
  b = unique(b)
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' Jaccard similarity under permutation testing
#' @param a vector (the one which will be permutated)
#' @param b vector
#' @param n number of first/top elements to take from each vector
#' @param N number of permutations
pjaccard = function(a, b, n, N = 10000) {
  n = min(n, length(a), length(b))
  topn_b = b[1:n]

  jacc_hist = replicate(N, {
    perm_a = sample(x = a, size = length(a), replace = FALSE)
    topn_a = perm_a[1:n]
    jaccard(topn_a, topn_b)
  })

  mean(jacc_hist)
}

calculate_similarity = function(omic, cindex_feats, rcll_feats, topn, N) {
  #' just in case `topn` exceeds the number of total feature selected
  cindex_topn = min(topn, length(cindex_feats))
  rcll_topn   = min(topn, length(rcll_feats))

  # get topn features
  cindex_topn_feats = cindex_feats[1:cindex_topn]
  rcll_topn_feats   = rcll_feats[1:rcll_topn]

  # Jaccard similarity
  jacc = jaccard(cindex_topn_feats, rcll_topn_feats)

  # Permutation-based Jaccard similarity
  pj = pjaccard(cindex_feats, rcll_feats, n = topn, N = N)

  tibble::tibble(
    omic = omic, ntop = topn, jaccard = jacc, permut_jaccard = pj
  )
}

# get similarity table
get_sim_tbl = function(omics, cindex_efs, rcll_efs, topn_vec, N = 10000) {
  data_list = list()
  for (omic in omics) {
    # get consensus features
    cindex_feats = cindex_efs[[omic]]$fs_stats()$consensus$feat_name
    rcll_feats   = rcll_efs  [[omic]]$fs_stats()$consensus$feat_name

    message('### Omic: ', omic)
    message('# C-index features: ', length(cindex_feats))
    message('# RCLL features: ', length(rcll_feats))

    data_list[[omic]] = future_lapply(topn_vec, function(topn) {
      calculate_similarity(omic, cindex_feats, rcll_feats, topn, N)
    }, future.seed = TRUE)
  }
  dplyr::bind_rows(data_list)
}

run_analysis = function(disease_code, topn_vec) {
  message('\nStudy: ', disease_code)
  cindex_efs = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
  rcll_efs = readRDS(file = paste0(disease_code, '/efs/rcll_efs.rds'))
  stopifnot(names(rcll_efs) == names(cindex_efs))
  omics = names(cindex_efs)

  sim_tbl = get_sim_tbl(omics = omics, cindex_efs = cindex_efs,
    rcll_efs = rcll_efs, topn_vec = topn_vec)
  set1_colors = RColorBrewer::brewer.pal(n = 3, name = 'Set1')

  sim_tbl %>%
    #mutate(jaccard = zoo::rollmean(jaccard, k = 5, fill = NA)) %>%
    rename(Observed = jaccard, Expected = permut_jaccard) %>%
    pivot_longer(cols = 3:4) %>%
    mutate(name = factor(name, levels = c('Observed', 'Expected'))) %>%
    ggplot(aes(x = ntop, y = value, color = name)) +
    scale_color_brewer(palette = 'Set1', name = '') +
    geom_hline(yintercept = 0.6, linetype = 'dashed', color = set1_colors[3]) +
    geom_line() +
    facet_grid(rows = 'omic') +
    labs(x = 'Top N consensus features', y = 'Jaccard Similarity',
      title = paste0('Ensemble Feature Selection (', disease_code, ')\nC-index vs RCLL')
    ) +
    ylim(c(0,1)) +
    theme_bw(base_size = 14)
  ggsave(filename = paste0('efs_inv/', disease_code, '_msr_comp.png'),
    width = 7, height = 5, dpi = 300)
}

# Execute ----
disease_codes = c('PAAD', 'BLCA')
topn_vec = seq(2, 150)
for (disease_code in disease_codes) {
  run_analysis(disease_code, topn_vec)
}

#' simpler figure
# sim_tbl %>%
#   ggplot(aes(x = ntop, y = jaccard)) +
#   geom_line() +
#   geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'red') +
#   facet_grid(rows = 'omic') +
#   labs(x = 'Top N consensus features',
#     y = 'Jaccard Similarity',
#     title = paste0('Ensemble Feature Selection on ', disease_code, ' dataset: \n',
#       'RCLL vs C-index')
#   ) +
#   ylim(c(0,1)) +
#   theme_bw()
