#' Assess how similar are the top N features per omic modality when using
#' the eFS algorithm with each omic independently (`efs.R`) or input is the
#' combined multimodal dataset (`efs_multimodal.R`)
library(survmob)
library(tidyverse)

# Reproducibility
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

# Jaccard similarity between two vectors for varying top (first) N elements
calculate_similarity = function(omic, feats1, feats2, topn) {
  #' just in case `topn` exceeds the number of total features selected
  feats1_topn = min(topn, length(feats1))
  feats2_topn = min(topn, length(feats2))

  # get the top N features
  feats1 = feats1[1:feats1_topn]
  feats2 = feats2[1:feats2_topn]

  # Jaccard similarity
  jacc = jaccard(feats1, feats2)

  tibble::tibble(
    omic = omic, ntop = topn, jaccard = jacc
  )
}

#disease_codes = c('PAAD', 'BLCA')
disease_codes = c('PAAD')
topn_vec = seq(2, 150)

for (disease_code in disease_codes) {
  message(disease_code)

  # efs results - all optimized with OOB (1 - C-index)
  efs_single      = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
  efs_multimodal  = readRDS(file = paste0(disease_code, '/efs/multimodal_efs.rds'))
  mm_feats = efs_multimodal$fs_stats()$consensus$feat_name
  message(length(mm_feats), ' multimodal features')

  # multimodal task (top 100 multimodal features after correction)
  task_multimodal = readRDS(file = paste0(disease_code, '/data/task_multimodal.rds'))
  mm_feats100 = task_multimodal$feature_names

  omics = names(efs_single)
  omic_prefixes = omics
  omic_prefixes[which(omic_prefixes == 'Methyl')] = 'meth'

  data_list  = list()
  data_list2 = list()
  for (index in 1:length(omics)) {
    omic = omics[index]
    message('#Omic: ', omic)

    feats1 = efs_single[[index]]$fs_stats()$consensus$feat_name
    feats2 = mm_feats[str_detect(string = mm_feats, pattern = omic_prefixes[index])]

    data_list[[omic]] = lapply(topn_vec, function(topn) {
      calculate_similarity(omic, feats1, feats2, topn)
    })

    # get features from multimodal task
    feats3 = mm_feats100[str_detect(string = mm_feats100, pattern = omic_prefixes[index])]

    # subset single eFS features
    feats1 = feats1[1:length(feats3)]

    # Jaccard similarity
    jacc = jaccard(feats1, feats3)
    data_list2[[omic]] = tibble::tibble(omic = omic, jaccard = jacc)
  }
  sim_tbl = dplyr::bind_rows(data_list)
  sim_tbl2 = dplyr::bind_rows(data_list2)
  print(sim_tbl2)

  sim_tbl %>%
    ggplot(aes(x = ntop, y = jaccard)) +
    geom_hline(yintercept = 0.6, linetype = 'dashed', color = '#4DAF4A') +
    geom_line() +
    facet_grid(rows = 'omic') +
    labs(x = 'Top N consensus features', y = 'Jaccard Similarity',
      title = paste0('Ensemble Feature Selection (', disease_code, ')'),
      subtitle = 'Single omic eFS vs Multimodal eFS'
    ) +
    ylim(c(0,1)) +
    theme_bw(base_size = 14)
  ggsave(filename = paste0('efs_inv/', disease_code, '_single_vs_mm_efs.png'),
    width = 7, height = 5, dpi = 300)
}
