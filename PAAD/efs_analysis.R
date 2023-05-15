#' Analyze efs results from `efs.R`
library(survmob)
library(stabm)
library(tidyverse)

# TCGA study
disease_code = 'PAAD'

# Reproducibility
set.seed(42)

# efs results
cindex_efs = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
rcll_efs = readRDS(file = paste0(disease_code, '/efs/rcll_efs.rds'))
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))

# where to save figures
cindex_path = paste0(disease_code, '/efs/cindex')
rcll_path   = paste0(disease_code, '/efs/rcll')
dir.create(cindex_path, recursive = TRUE)
dir.create(rcll_path, recursive = TRUE)

# Plotting helpers
my_labels = c(
  'aorsf' = 'AORSF',
  'rsf_cindex' = 'RSF (C-index)',
  'rsf_extratrees' = 'RSF (ExtraTrees)',
  'rsf_logrank' = 'RSF (Logrank)',
  'rsf_maxstat' = 'RSF (Maxstat)',
  'consensus' = 'Consensus'
)
hue_colors = scale_color_hue()$palette(n = 6)
hue_colors = hue_colors[-4] # remove blue-ish color
my_colors = c(hue_colors, 'black') # `black` for consensus
names(my_colors) = names(my_labels)

# C-index eFS ----
message('### C-index eFS ###')
for (efs in cindex_efs) {
  omic = efs$task_id
  omic_name = ifelse(omic == 'Methyl', 'Methylation', omic)
  message(omic)

  #' Consensus features
  efs$ffs_plot(top_n = 10) +
    labs(title = 'Consensus Features', x = '')
  ggsave(filename = paste0(cindex_path, '/', omic, '_consensus_top10.png'),
    width = 5, height = 5, units = 'in', dpi = 300)

  #' Performance
  efs$res_plot(type = 'perf', msr_label = 'OOB (1 - C-index)', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = my_colors, breaks = names(my_labels)) +
    scale_x_discrete(labels = my_labels) # change model names on x-axis
  ggsave(filename = paste0(cindex_path, '/', omic, '_perf.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Nfeatures
  efs$res_plot(type = 'nfeat', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = my_colors, breaks = names(my_labels)) +
    scale_x_discrete(labels = my_labels)
  ggsave(filename = paste0(cindex_path, '/', omic, '_nfeats.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Stability
  p = efs$res_plot(type = 'stab', task = tasks[[omic]], title = omic_name,
    include_legend = TRUE, ylimits = c(0, 0.6))
  lrn_levels = levels(p$data$lrn_id)
  p +
    guides(fill = guide_legend(title = 'Learners')) +
    scale_fill_manual( # set same colors as prev figures but keep the desc order
      values = my_colors[lrn_levels],
      breaks = names(my_labels[lrn_levels]),
      labels = my_labels[lrn_levels]
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(filename = paste0(cindex_path, '/', omic, '_stab.png'),
    width = 6, height = 5, units = 'in', dpi = 300)
}

# RCLL eFS ----
message('### RCLL eFS ###')
for (efs in rcll_efs) {
  omic = efs$task_id
  omic_name = ifelse(omic == 'Methyl', 'Methylation', omic)
  message(omic)

  #' Consensus features
  efs$ffs_plot(top_n = 10) +
    labs(title = 'Consensus Features', x = '')
  ggsave(filename = paste0(rcll_path, '/', omic, '_consensus_top10.png'),
    width = 5, height = 5, units = 'in', dpi = 300)

  #' Performance
  efs$res_plot(type = 'perf', msr_label = 'RCLL', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = my_colors, breaks = names(my_labels)) +
    scale_x_discrete(labels = my_labels) # change model names on x-axis
  ggsave(filename = paste0(rcll_path, '/', omic, '_perf.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Nfeatures
  efs$res_plot(type = 'nfeat', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = my_colors, breaks = names(my_labels)) +
    scale_x_discrete(labels = my_labels)
  ggsave(filename = paste0(rcll_path, '/', omic, '_nfeats.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Stability
  p = efs$res_plot(type = 'stab', task = tasks[[omic]], title = omic_name,
    include_legend = TRUE, ylimits = c(0, 0.6))
  lrn_levels = levels(p$data$lrn_id)
  p +
    guides(fill = guide_legend(title = 'Learners')) +
    scale_fill_manual( # set same colors as prev figures but keep the desc order
      values = my_colors[lrn_levels],
      breaks = names(my_labels[lrn_levels]),
      labels = my_labels[lrn_levels]
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(filename = paste0(rcll_path, '/', omic, '_stab.png'),
    width = 6, height = 5, units = 'in', dpi = 300)
}
