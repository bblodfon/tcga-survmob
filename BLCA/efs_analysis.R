#' Analyze efs results from `efs.R`
library(survmob)
library(stabm)
library(tidyverse)

# TCGA study
disease_code = 'BLCA'

# Reproducibility
set.seed(42)

# efs results
cindex_efs = readRDS(file = paste0(disease_code, '/efs/cindex_efs.rds'))
rcll_efs = readRDS(file = paste0(disease_code, '/efs/rcll_efs.rds'))
tasks = readRDS(file = paste0(disease_code, '/data/tasks.rds'))

# where to save the figures
cindex_path = paste0(disease_code, '/efs/cindex')
rcll_path   = paste0(disease_code, '/efs/rcll')
dir.create(cindex_path, recursive = TRUE)
dir.create(rcll_path, recursive = TRUE)

# Plotting helpers ----
lrn_labels = c(
  'aorsf' = 'AORSF',
  'rsf_cindex' = 'RSF (C-index)',
  'rsf_extratrees' = 'RSF (ExtraTrees)',
  'rsf_logrank' = 'RSF (Logrank)',
  'rsf_maxstat' = 'RSF (Maxstat)',
  'consensus' = 'Consensus'
)
hue_colors = scale_color_hue()$palette(n = length(lrn_labels))
hue_colors = hue_colors[-4] # remove blue-ish color

# colors for the learners (including the 'consensus')
lrn_colors = c(hue_colors, 'black') # `black` for consensus
names(lrn_colors) = names(lrn_labels)

# colors for the omics
stopifnot(names(cindex_efs) == names(rcll_efs))
omic_names = names(cindex_efs)
omic_names = dplyr::if_else(
  condition = omic_names == 'Methyl',
  true = 'Methylation',
  false = omic_names
)
omic_colors = scale_color_hue()$palette(n = length(omic_names))

# C-index eFS ----
message('### C-index eFS ###')

## Per-omic plots ----
res_list = list()
stab_list = list()
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
  efs$result %>%
    mutate(score = 1 - score) %>% # OOB = 1 - Cindex
    ggplot(aes(x = lrn_id, y = score, fill = lrn_id)) +
    geom_boxplot() +
    labs(y = 'Harrell\'s C-index', x = NULL, title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = lrn_colors, breaks = names(lrn_labels)) +
    scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  ggsave(filename = paste0(cindex_path, '/', omic, '_perf.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Nfeatures
  efs$res_plot(type = 'nfeat', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = lrn_colors, breaks = names(lrn_labels)) +
    scale_x_discrete(labels = lrn_labels)
  ggsave(filename = paste0(cindex_path, '/', omic, '_nfeats.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Stability
  p = efs$res_plot(type = 'stab', task = tasks[[omic]], title = omic_name,
    include_legend = TRUE, ylimits = c(0, 0.6))
  lrn_levels = levels(p$data$lrn_id)
  p +
    guides(fill = guide_legend(title = 'Learners')) +
    scale_fill_manual( # set same colors as prev figures but keep the desc order
      values = lrn_colors[lrn_levels],
      breaks = names(lrn_labels[lrn_levels]),
      labels = lrn_labels[lrn_levels]
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(filename = paste0(cindex_path, '/', omic, '_stab.png'),
    width = 6, height = 5, units = 'in', dpi = 300)

  # aggregate results across omics
  res_list[[omic]] =
    efs$result %>%
    mutate(score = 1 - score) %>%
    add_column(omic_name)

  # nogueira stability
  stab_list[[omic]] =
    efs$stab(task = tasks[[omic]], stab_metrics = 'nogueira') %>%
    add_column(omic_name)
}
res = res_list %>% bind_rows()
stab_res = stab_list %>% bind_rows()

## Multi-omics plots ----
# Performance + Nfeatures
# order omics according to performance metric (C-index)
res = res %>%
  mutate(omic_name = forcats::fct_reorder(omic_name, score, .desc = TRUE))

omic_levels = levels(res$omic_name)
names(omic_colors) = omic_levels

res %>%
  ggplot(aes(x = lrn_id, y = score, fill = omic_name)) +
  geom_boxplot() +
  labs(y = 'Harrell\'s C-index', x = NULL) +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(cindex_path, '/all_omics_perf.png'),
  width = 5, height = 5, units = 'in', dpi = 300)

res %>%
  ggplot(aes(x = lrn_id, y = nfeatures, fill = omic_name)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = NULL) +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(cindex_path, '/all_omics_nfeats.png'),
  width = 5, height = 5, units = 'in', dpi = 300)

# Stability plot
stab_res = stab_res %>%
  filter(lrn_id != 'consensus') %>%
  mutate(omic_name = fct_reorder(omic_name, nogueira, .desc = TRUE))
omic_levels_stab = levels(stab_res$omic_name)

stab_res %>%
  ggplot(aes(x = lrn_id, y = nogueira, fill = omic_name)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = NULL, y = 'Nogueira Similarity') +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels_stab) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(cindex_path, '/all_omics_nogueira_stab.png'),
  width = 6, height = 5, units = 'in', dpi = 300)

# RCLL eFS ----
message('### RCLL eFS ###')

## Per-omic plots ----
res_list = list()
stab_list = list()
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
    scale_fill_manual(values = lrn_colors, breaks = names(lrn_labels)) +
    scale_x_discrete(labels = lrn_labels) # change model names on x-axis
  ggsave(filename = paste0(rcll_path, '/', omic, '_perf.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Nfeatures
  efs$res_plot(type = 'nfeat', title = omic_name) +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = lrn_colors, breaks = names(lrn_labels)) +
    scale_x_discrete(labels = lrn_labels)
  ggsave(filename = paste0(rcll_path, '/', omic, '_nfeats.png'),
    width = 5, height = 6, units = 'in', dpi = 300)

  #' Stability
  p = efs$res_plot(type = 'stab', task = tasks[[omic]], title = omic_name,
    include_legend = TRUE, ylimits = c(0, 0.6))
  lrn_levels = levels(p$data$lrn_id)
  p +
    guides(fill = guide_legend(title = 'Learners')) +
    scale_fill_manual( # set same colors as prev figures but keep the desc order
      values = lrn_colors[lrn_levels],
      breaks = names(lrn_labels[lrn_levels]),
      labels = lrn_labels[lrn_levels]
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(filename = paste0(rcll_path, '/', omic, '_stab.png'),
    width = 6, height = 5, units = 'in', dpi = 300)

  # aggregate results across omics
  res_list[[omic]] =
    efs$result %>%
    add_column(omic_name)

  # nogueira stability
  stab_list[[omic]] =
    efs$stab(task = tasks[[omic]], stab_metrics = 'nogueira') %>%
    add_column(omic_name)
}
res = res_list %>% bind_rows()
stab_res = stab_list %>% bind_rows()

## Multi-omics plots ----
# Performance + Nfeatures
# order omics according to performance metric (C-index)
res = res %>%
  mutate(omic_name = forcats::fct_reorder(omic_name, score, .desc = TRUE))

omic_levels = levels(res$omic_name)
names(omic_colors) = omic_levels

res %>%
  ggplot(aes(x = lrn_id, y = score, fill = omic_name)) +
  geom_boxplot() +
  labs(y = 'RCLL', x = NULL) +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(rcll_path, '/all_omics_perf.png'),
  width = 5, height = 5, units = 'in', dpi = 300)

res %>%
  ggplot(aes(x = lrn_id, y = nfeatures, fill = omic_name)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = NULL) +
  ylim(c(0, 1000)) +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(rcll_path, '/all_omics_nfeats.png'),
  width = 5, height = 5, units = 'in', dpi = 300)

# Stability plot
stab_res = stab_res %>%
  filter(lrn_id != 'consensus') %>%
  mutate(omic_name = fct_reorder(omic_name, nogueira, .desc = TRUE))
omic_levels_stab = levels(stab_res$omic_name)

stab_res %>%
  ggplot(aes(x = lrn_id, y = nogueira, fill = omic_name)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = NULL, y = 'Nogueira Similarity') +
  guides(x = guide_axis(angle = 45), fill = guide_legend(title = 'Omics')) +
  scale_fill_manual(values = omic_colors, breaks = omic_levels_stab) +
  scale_x_discrete(labels = lrn_labels) + # change model names on x-axis
  theme_bw(base_size = 14)
ggsave(filename = paste0(rcll_path, '/all_omics_nogueira_stab.png'),
  width = 6, height = 5, units = 'in', dpi = 300)
