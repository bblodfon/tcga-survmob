#' => Visualize benchmarking results using the bootstrap scores from the
#' test set
#' => Compare with Cox + Clinical baseline
#' => Results are generated from the scripts: `benchmark.R` and `benchmark_multimodal.R`
library(survmob)
library(tidyverse)
library(ggdist)

# TCGA study
disease_code = 'PAAD'

# Where to store the output figures
res_path = paste0(disease_code, '/bench/boot')
dir.create(res_path)

# Benchmark data
cox =
  readRDS(file = paste0(disease_code, '/bench/cox.rds')) %>%
  `$`(result) %>%
  select(-model) %>%
  add_column(tune = 'no')
#' Tuning with Uno's C-index vs RCLL (omic-specific eFS)
mob_uno =
  readRDS(file = paste0(disease_code, '/bench/mob_uno.rds')) %>%
  `$`(result) %>%
  add_column(tune = 'uno')
mob_rcll =
  readRDS(file = paste0(disease_code, '/bench/mob_rcll.rds')) %>%
  `$`(result) %>%
  add_column(tune = 'rcll')
#' Tuning with Uno's C-index vs RCLL (multimodal eFS)
mm_uno =
  readRDS(file = paste0(disease_code, '/bench/mm_uno.rds')) %>%
  `$`(result) %>%
  add_column(tune = 'uno')
mm_rcll =
  readRDS(file = paste0(disease_code, '/bench/mm_rcll.rds')) %>%
  `$`(result) %>%
  add_column(tune = 'rcll')

#' Add all benchmark data combined
bench_data = list(cox, mob_uno, mob_rcll, mm_uno, mm_rcll) %>% bind_rows()
bench_data

# Reproducibility
seed = 42

# Helper functions ----

#' @description
#' Plot bootstrap Confidence Intervals (CIs) from a benchmark result tibble
#'
#' @param bench_res tibble object with `task_id`, `lrn_id`, `id` (identifier of
#' the combination of `task_id` and `lrn_id`) and `value` (score)
#' @param topn How many `task_id`/`lrn_id` combos should be included in the
#' figure?
#' Default is 12.
#' The CIs are ordered by decreasing median `value`.
#' One CI per combo `id`.
#' @param probs 2 probabilities for the CIs
#' @param interval_colors 2 colors for the CIs
#' @param x_label x-axis label name (e.g. measure name)
#' @param add_cox add vertical line for the cox baseline
#' @param adjust_x how much to x-adjust the cox baseline annotation
#' @param adjust_y how much to y-adjust the cox baseline annotation
#' @param cox_size size of the cox baseline annotation
#' @param cox_angle angle of the cox baseline annotation
plot_bootCIs = function(bench_res, topn = 12, probs = c(0.8, 0.95),
  interval_colors = c('#A6CEE3', '#1F78B4'), x_label = 'Measure name',
  add_cox = TRUE, adjust_x = 0.01, adjust_y = -2.5, cox_size = 5, cox_angle = 0) {
  stopifnot(length(probs) == 2) # only 2 probabilities for the CIs
  stopifnot(length(interval_colors) == 2) # only 2 colors

  #' for matching colors with intervals in `scale_*_manual`
  probs = sort(probs, decreasing = TRUE)

  # Get top N model/dataset combos by median score
  top_ids =
    bench_res %>%
    group_by(id) %>%
    summarise(m = median(value), .groups = 'drop') %>%
    arrange(desc(m)) %>%
    slice_head(n = topn) %>%
    pull(id)

  bench_res_flt =
    bench_res %>%
    # keep only the top N task:learner ids
    filter(id %in% top_ids) %>%
    # add a more comprehensive model name
    mutate(model = case_when(
      lrn_id == 'coxph' ~ 'CoxPH',
      lrn_id == 'coxnet' ~ 'CoxNet',
      lrn_id == 'coxboost' ~ 'CoxBoost',
      lrn_id == 'rsf_cindex' ~ 'RSF (C-index)',
      lrn_id == 'rsf_extratrees' ~ 'RSF (ExtraTrees)',
      lrn_id == 'rsf_logrank' ~ 'RSF (Logrank)',
      lrn_id == 'rsf_maxstat' ~ 'RSF (Maxstat)',
      lrn_id == 'aorsf' ~ 'AORSF',
      lrn_id == 'xgboost_cox_early' ~ 'XGBoost (Cox)',
      lrn_id == 'xgboost_aft_early' ~ 'XGBoost (AFT)'
    )) %>%
    # make `id` more readable
    mutate(id = paste0(model, ' / ', task_id)) %>%
    # order model/task combos according to bootstrap median score
    mutate(
      id = factor(id),
      id = forcats::fct_reorder(id, .x = value, .desc = TRUE)
    )

  p =
    bench_res_flt %>%
    ggplot(aes(y = id, x = value)) +
    stat_pointinterval(
      aes(interval_color = after_stat(level)),
      .width = probs,
      orientation = 'horizontal',
      point_interval = 'median_qi', # median quantile interval
      point_size = 2,
      point_color = '#D55E00'
    ) +
    scale_color_manual(
      values = interval_colors,
      labels = scales::label_percent()(probs),
      aesthetics = 'interval_color',
      guide = guide_legend(
        override.aes = list(point_color = NA),
        title = 'Bootstrap Confidence Intervals'
      )
    ) +
    labs(x = x_label, y = '') +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = 'black'),
      axis.line.x.top = element_blank(),
      axis.ticks.x = element_line(color = 'black'),
      axis.ticks.x.top = element_line(color = 'gray50'),
      axis.ticks.y = element_blank(),
      legend.position = 'top'
    )

  # add cox median line + text annotation
  if (add_cox) {
    #' search for the `coxph` in the initial tibble since the CoxPH + clinical
    #' combo might not be in the top N task/learner combos
    cox_median = bench_res %>%
      filter(lrn_id == 'coxph') %>%
      summarize(m = median(value)) %>%
      pull(m)

    p = p +
      geom_vline(xintercept = cox_median, color = '#222222', lty = 2) +
      annotate(
        geom = 'text',
        x = cox_median + adjust_x,
        y = length(unique(bench_res_flt$id)) + adjust_y,
        label = 'Cox Baseline',
        size = cox_size,
        angle = cox_angle
      )
  }

  p
}

# Reshape benchmark data
bench_res =
  bench_data %>%
  mutate(scores = purrr::map(boot_res, 'scores')) %>%
  select(!matches('boot_res')) %>%
  tidyr::unnest(cols = scores) %>%
  tidyr::pivot_longer(cols = !matches('task_id|lrn_id|model|tune'),
    names_to = 'measure', values_to = 'value') %>%
  filter(measure != 'dcal') #' exclude `dcal` measure
bench_res
test_measures = unique(bench_res$measure)
test_measures

measure_labels = c(
  harrell_c = 'Harrell\'s C-index',
  uno_c = 'Uno\'s C-index',
  ibrier_erv = 'IBS (ERV)',
  rcll_erv = 'RCLL (ERV)'
)

# Plots ----
## Aggregated performance across datasets ----
#' => Doesn't make much sense to me
#' => All other models expect CoxPH have been tested in many more datasets and
#' thus it's more likely to bring performance 'down'
if (FALSE) {
  plist = list()
  for(msr_id in test_measures) {
    plist[[msr_id]] =
      bench_res %>%
      filter(measure == msr_id) %>%
      #filter(measure == msr_id, task_id == 'Multimodal-eFS') %>%
      mutate(
        lrn_id = factor(lrn_id),
        lrn_id = forcats::fct_reorder(lrn_id, .x = value, .desc = TRUE)
      ) %>%
      ggplot(aes(x = lrn_id, y = value)) +
      geom_boxplot() +
      guides(x = guide_axis(angle = 45)) +
      labs(y = msr_id, x = '') +
      theme_bw(base_size = 14)
  }
  plist$harrell_c
  plist$uno_c
  plist$ibrier_erv
  plist$rcll_erv
}

## Boostrap CIs ----
topn = 12

#' ALL model/dataset combos vs Cox + Clinical (baseline)
#' Separate the benchmark results from the different model tunings
#' Note: `multimodal-eFS` results are included

# Uno's C ----
message('### Models tuned with Uno\'s C-index')

res_uno = bench_res %>%
  filter(tune == 'uno'| tune == 'no') %>%
  mutate(id = paste0(task_id, ':', lrn_id)) %>% # add unique task/learner id
  relocate(id, .after = lrn_id)

# Harrell's C-index
msr_id = 'harrell_c'
res_flt = res_uno %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = 0.06) +
  # add random performance line
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.51, y = topn - 1, label = 'Random',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/uno_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# Uno's C-index
msr_id = 'uno_c'
res_flt = res_uno %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = 0.07) +
  # add random performance line
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.515, y = 4, label = 'Random',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/uno_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# IBS (ERV)
msr_id = 'ibrier_erv'
res_flt = res_uno %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = -0.02, adjust_y = -2, cox_size = 4, cox_angle = 90) +
  # add KM baseline
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.02, y = 4, label = 'Kaplan-Meier',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/uno_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# RCLL (ERV)
msr_id = 'rcll_erv'
res_flt = res_uno %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = -0.005, adjust_y = -2, cox_size = 4, cox_angle = 90) +
  # add KM baseline
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.005, y = 4, label = 'Kaplan-Meier',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/uno_tune_' , msr_id, '.png'), units = 'in',
  width = 9, height = 5, dpi = 300)

# RCLL ----
message('### Models tuned with RCLL')
res_rcll = bench_res %>%
  filter(tune == 'rcll'| tune == 'no') %>%
  mutate(id = paste0(task_id, ':', lrn_id)) %>% # add unique task/learner id
  relocate(id, .after = lrn_id)

# Harrell's C-index
msr_id = 'harrell_c'
res_flt = res_rcll %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = 0.06) +
  # add random performance line
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.51, y = topn - 1, label = 'Random',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/rcll_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# Uno's C-index
msr_id = 'uno_c'
res_flt = res_rcll %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = 0.07) +
  # add random performance line
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.515, y = 4, label = 'Random',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/rcll_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# IBS (ERV)
msr_id = 'ibrier_erv'
res_flt = res_rcll %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = -0.02, adjust_y = -2, cox_size = 4, cox_angle = 90) +
  # add KM baseline
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.02, y = 4, label = 'Kaplan-Meier',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/rcll_tune_' , msr_id, '.png'), units = 'in',
  width = 8, height = 5, dpi = 300)

# RCLL (ERV)
msr_id = 'rcll_erv'
res_flt = res_rcll %>% filter(measure == msr_id)

plot_bootCIs(bench_res = res_flt, x_label = measure_labels[[msr_id]],
  adjust_x = -0.005, adjust_y = -2, cox_size = 4, cox_angle = 90) +
  # add KM baseline
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.005, y = 4, label = 'Kaplan-Meier',
    color = 'black', angle = 90)
ggsave(file = paste0(res_path, '/rcll_tune_' , msr_id, '.png'), units = 'in',
  width = 9, height = 5, dpi = 300)
