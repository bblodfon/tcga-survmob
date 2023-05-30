#' => Visualize Benchmarking Results using the output from the
#' Bayesian analysis script (`bench_bayes.R`)
#' => Need quite a few memory to run this script (close to 10GB)
#' => Execute: `Rscript bench_bayes_vis.R`
library(survmob)
library(tidyverse)
library(rstanarm)
library(ggdist)

# TCGA study
disease_code = 'PAAD'

# Where to store the output images
models_path_uno  = paste0(disease_code, '/bench/mob_uno/')
models_path_rcll = paste0(disease_code, '/bench/mob_rcll/')
res_path_uno  = paste0(models_path_uno, 'results')
res_path_rcll = paste0(models_path_rcll, 'results')
res_path_cmp  = paste0(disease_code, '/bench/uno_vs_rcll/')
dir.create(res_path_uno, recursive = TRUE)
dir.create(res_path_rcll, recursive = TRUE)
dir.create(res_path_cmp, recursive = TRUE)

# Reproducibility
seed = 42
set.seed(seed)

# Helper functions ----

#' @param stan_model a Bayesian `stanreg` model
#' @param newdata a one-column `data.frame`.
#' The column's name must be `lrn_id` and it should have the names of the
#' different models used in the `formula` of the `stan_model`.
#' @param rename_models change model names to more friendly/understandable ones
#' @param seed seed number
model_post_draw = function(stan_model, newdata, rename_models = TRUE, seed = seed) {
  post_draws =
    rstanarm::posterior_epred(stan_model, newdata = newdata, re.form = NA, seed = seed) %>%
    as_tibble() %>%
    rlang::set_names(nd$lrn_id) %>%
    pivot_longer(everything(), names_to = 'model', values_to = 'posterior')

  if (rename_models) {
    # rename models
    post_draws = post_draws %>%
      mutate(model = case_when(
        model == 'coxph' ~ 'CoxPH',
        model == 'coxnet' ~ 'CoxNet',
        model == 'coxboost' ~ 'CoxBoost',
        model == 'rsf_cindex' ~ 'RSF (C-index)',
        model == 'rsf_extratrees' ~ 'RSF (ExtraTrees)',
        model == 'rsf_logrank' ~ 'RSF (Logrank)',
        model == 'rsf_maxstat' ~ 'RSF (Maxstat)',
        model == 'aorsf' ~ 'AORSF',
        model == 'xgboost_cox_early' ~ 'XGBoost (Cox)',
        model == 'xgboost_aft_early' ~ 'XGBoost (AFT)'
      )
    )
  }

  # reorder models according to median posterior
  post_draws = post_draws %>%
    mutate(
      model = factor(model),
      model = forcats::fct_reorder(model, .x = posterior, .desc = TRUE)
    )

  post_draws
}

#' @param post_draws output from `model_post_draw()`
model_perf_plot = function(post_draws, type = 'interval',
  probs = c(0.8, 0.95), x_label = 'Measure name') {
  match.arg(arg = type, choices = c('interval', 'halfeye'))

  if (length(probs) == 3) {
    my_cols = c('#A6CEE3', '#1F78B4', '#08306B')
  } else if (length(probs) == 2) {
    my_cols = c('#A6CEE3', '#1F78B4')
  } else if (length(probs) == 1) {
    my_cols = c('#1F78B4')
  } else {
    stop('`probs` should be a vector of length >= 1 and <= 3')
  }

  #' for matching colors with intervals in `scale_*_manual`
  probs = sort(probs, decreasing = TRUE)

  if (type == 'interval') {
    p =
      post_draws %>%
      ggplot(aes(y = model, x = posterior)) +
      stat_pointinterval(
        aes(interval_color = after_stat(level)),
        .width = probs,
        point_interval = 'median_qi', # median quantile interval
        point_size = 3,
        point_color = '#D55E00',
      ) +
      scale_color_manual(
        values = my_cols,
        labels = scales::label_percent()(probs),
        aesthetics = 'interval_color',
        guide = guide_legend(
          override.aes = list(point_color = NA),
          title = 'Credible Intervals'
        )
      )
  } else if (type == 'halfeye') {
    p =
      post_draws %>%
      ggplot(aes(y = model, x = posterior)) +
      stat_halfeye(
        aes(fill = after_stat(level)),
        .width = probs,
        point_size = 3,
        point_interval = 'median_qi', # median quantile interval
        point_color = '#D55E00'
      ) +
      scale_fill_manual(
        values = my_cols,
        labels = scales::label_percent()(probs),
        na.translate = FALSE,
        guide = guide_legend(
          override.aes = list(point_color = NA),
          title = 'Probability mass'
        )
      )
    }

  p = p +
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
  p
}

#' @param model_list a list of Bayesian `stanreg` models
#' @param seed seed number
omics_postdiff_draw = function(model_list, seed = seed) {
  post_list = list()
  for (index in 1:length(model_list)) {
    stan_model = model_list[[index]]
    omic = names(model_list)[[index]]
    message(omic)

    #' `1` => post draw with omic
    #' `0` => post draw without omic
    nd = as.integer(c(1,0)) %>% tibble::as_tibble_col(column_name = omic)

    post_list[[omic]] =
      rstanarm::posterior_epred(stan_model, newdata = nd, re.form = NA, seed = seed) %>%
      as_tibble() %>%
      # performance with omic (1st column) vs (minus) no omic (2nd column)
      mutate(diff = `1` - `2`) %>%
      select(diff) %>%
      rename(!!omic := diff)
  }

  dplyr::bind_cols(post_list) %>%
    # remove `_omic` from omic variable names
    tidyr::pivot_longer(cols = everything(), names_to = 'omic',
      names_pattern = '(.*)_omic', values_to = 'post_diff') %>%
    # better name for Methylation omic
    mutate(omic = case_when(
      omic == 'Methyl' ~ 'Methylation',
      TRUE ~ omic
    )) %>%
    # reorder omics according to median posterior differences
    mutate(
      omic = factor(omic),
      omic = forcats::fct_reorder(omic, .x = post_diff, .desc = TRUE)
    )
}

#' @description
#' Make a ridgeline plot using several posterior distributions and add a ROPE
#' shady area
#'
#' @param post_draws a `tibble` with 2 columns: a factor variable that has
#' **ordered groups** (e.g. model names or omic names for example) and a
#' second column that has the posterior distribution draws corresponding
#' to each group.
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param x_axis_label x-axis label (e.g. performance measure)
#' @param size practical effect size for ROPE (Region of Practical Equivalence).
#' Default is 0.02.
#' @param ROPE_center The center of the ROPE region.
#' Default is 0 which means that the posterior draws are distribution differences.
#' @param ROPE_area if `TRUE`, fill with a grey color the ROPE area and draw
#' vertical lines to demarcate the practical equivalence region.
#' @param draw_arrow add an arrow and a value atop to show the ROPE's
#' effect size (`size`).
#' @param pal Default choice is `hue` multicolored palette palette, nice for a
#' small number of groups, e.g. up to 5.
#' Hue colors will be consistent, i.e. the same across groups if this function
#' is used multiple times with the same groups but different distribution draws.
#' Another palette choice is `viridis` (nicer for 5+ groups, not consistent
#' across different group ordering).
#' You can also pass any color you may want, e.g. `lightblue` or `#A6CEE3`.
#' @param rope_stats output from `rope_stats()` function using the same inputs
#' as here: `post_draws`, `size` and `ROPE_center`.
#' @param add_right add the right ROPE probabilities per group as extra text
#' annotations
#' @param add_left add the left ROPE probabilities per group as extra text
#' annotations
#' @param x_right designate the place on the x-axis to add the right ROPE
#' probabilities.
#' If NULL, the maximum posterior value is used.
#' @param x_left designate the place on the x-axis to add the left ROPE
#' probabilities.
#' If NULL, the minimum posterior value is used.
#' @param label_right title for the right ROPE text annotation.
#' Default: **PPPD** => Posterior Probability of Positive Practical Difference
#' @param label_left title for the left ROPE text annotation.
#' Default: **PNPD** => Posterior Probability of Negative Practical Difference
ridgeline_plot = function(post_draws, title = '', subtitle = '', x_axis_label = '',
  size = 0.02, ROPE_center = 0, ROPE_area = TRUE, draw_arrow = TRUE, pal = 'hue',
  rope_stats = NULL, add_right = FALSE, add_left = FALSE,
  x_right = NULL, x_left = NULL, label_right = 'PPPD', label_left = 'PNPD'
) {
  grp_var = colnames(post_draws)[1]
  post_var = colnames(post_draws)[2]

  grps = sort(levels(post_draws[[grp_var]]))
  n_grps = length(grps)

  if (ROPE_area) {
    p = post_draws %>%
      ggplot(aes(y = .data[[grp_var]], x = .data[[post_var]],
        fill = .data[[grp_var]],
        fill_ramp = after_stat(abs(x - ROPE_center) > size)))
  } else {
    p = post_draws %>%
      ggplot(aes(y = .data[[grp_var]], x = .data[[post_var]]),
        fill = .data[[grp_var]])
  }

  is_palette = (pal == 'hue' | pal == 'viridis')
  if (is_palette) { # different color per distribution
    p = p +
      stat_slab(alpha = 0.8, height = 2, show.legend = FALSE)
  } else { #' `pal` is a single color, same for all distributions
    p = p +
      stat_slab(fill = pal, alpha = 0.8, height = 2, show.legend = FALSE)
  }

  p = p +
    scale_fill_ramp_discrete(from = 'grey80') +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_axis_label,
      y = ''
    ) +
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

  if (pal == 'hue') {
    # Consistent colors for each group level (irrelevant of ranking order)
    my_colors = scale_color_hue()$palette(n = n_grps)

    p = p +
      scale_fill_manual(
        values = my_colors,
        labels = grps,
        breaks = grps
      )
  } else if (pal == 'viridis') {
    p = p + scale_fill_viridis_d()
  }

  # add ROPE area
  if (ROPE_area) {
    p = p +
      geom_vline(xintercept = ROPE_center, lty = 2, color = '#E41A1C') +
      geom_vline(xintercept = c(ROPE_center - size, ROPE_center + size),
        lty = 2, alpha = 0.5)
    # add arrow + label for effect size
    if (draw_arrow) {
      p = p +
        annotate(
          geom = 'segment',
          x = 0, xend = size,
          y = n_grps, yend = n_grps,
          linewidth = 0.7,
          arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
        ) +
        annotate(
          geom = 'text', label = size,
          x = size/2, y = n_grps + 0.25
        )
    }
  }

  # add ROPE probabilities as extra text annotations
  if (!is.null(rope_stats)) {
    #' compute `x_right` and `x_left` if not given
    if (is.null(x_right)) {
      x_right = max(post_draws[[post_var]])
    }
    if (is.null(x_left)) {
      x_left = min(post_draws[[post_var]])
    }

    for (index in 1:nrow(rope_stats)) {
      if (add_right) {
        prob_right = rope_stats[index, 'prob_right']
        p = p +
          annotate(
            geom = 'text', label = paste0(round(prob_right*100), '%'),
            x = x_right, y = index + 0.5
          )
      }
      if (add_left) {
        prob_left  = rope_stats[index, 'prob_left']
        p = p +
          annotate(
            geom = 'text', label = paste0(round(prob_left*100), '%'),
            x = x_left, y = index + 0.5
          )
      }
    }

    # add labels to the ROPE probabilities annotations
    if (add_right) {
      p = p +
        annotate(
          geom = 'text', label = label_right, color = '#4DAF4A', fontface = 2,
          x = x_right, y = index + 1.4
        )
    }

    if (add_left) {
      p = p +
        annotate(
          geom = 'text', label = label_left, color = '#E41A1C', fontface = 2,
          x = x_left, y = index + 1.4
        )
    }
  }

  p
}

#' @description
#' Get ROPE statistics, i.e. 3 probabilities that summarize each posterior
#' distribution per group
#'
#' @return a tibble with 4 columns: the grouping variable, and 3 probabilities:
#' - `prob_within`: the integral of the posterior within the ROPE region
#' (-size, size)
#' - `prob_right`: the integral of the posterior to the right of the ROPE
#' (size, Inf)
#' - `prob_left`: the integral of the posterior to the left of the ROPE
#' (-Inf, -size)
#'
#' @param post_draws a `tibble` with 2 columns: a factor variable that has
#' **groups** (e.g. model names or omic names for example) and a
#' second column that has the posterior distribution draws corresponding to each
#' group.
#' @param size practical effect size for ROPE (Region of Practical Equivalence)
#' @param ROPE_center The center of the ROPE region.
#' Default is 0 which means that the posterior draws are distribution differences.
rope_stats = function(post_draws, size = 0.02, ROPE_center = 0) {
  grp_var = colnames(post_draws)[1]
  post_var = colnames(post_draws)[2]

  post_draws %>%
    group_by(!!rlang::sym(grp_var)) %>%
    summarise(
      prob_within = mean(abs(!!rlang::sym(post_var) - ROPE_center) <= size),
      prob_right  = mean(!!rlang::sym(post_var) > ROPE_center + size),
      prob_left   = mean(!!rlang::sym(post_var) < ROPE_center - size)
    )
}

# Uno's C ----
message('### Models tuned with Uno\'s C-index')

## Model comparison ----
message('## Model comparison')

### Load bayesian models ----
model_harrell_c  = readRDS(file = paste0(models_path_uno, 'models_fit_harrell_c.rds'))
model_uno_c      = readRDS(file = paste0(models_path_uno, 'models_fit_uno_c.rds'))
model_ibrier_erv = readRDS(file = paste0(models_path_uno, 'models_fit_ibrier_erv.rds'))
model_rcll_erv   = readRDS(file = paste0(models_path_uno, 'models_fit_rcll_erv.rds'))
models = list(model_harrell_c, model_uno_c, model_ibrier_erv, model_rcll_erv)

### Models posterior draw ----
post_draws_uno = list()
for (model in models) {
  msr_id = unique(model$data$measure)
  message('# Measure: ', msr_id)

  nd = distinct(model$data, lrn_id)
  post_draws_uno[[msr_id]] = model_post_draw(
    stan_model = model, newdata = nd, rename_models = TRUE, seed = seed
  ) %>% add_column(msr_id)
}

# clear some memory
rm(model_harrell_c, model_uno_c, model_ibrier_erv, model_rcll_erv, models)
gc()

### Plots ----
# Harrell's C-index
msr_id = unique(post_draws_uno$harrell_c$msr_id)
p = model_perf_plot(type = 'interval', post_draws = post_draws_uno[[msr_id]],
  probs = c(0.8, 0.95), x_label = 'Harrell\'s C-index')
p +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.505, y = 8, label = 'Random', color = 'black',
    angle = 90) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01))
ggsave(file = paste0(res_path_uno, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# Uno's C-index
msr_id = unique(post_draws_uno$uno_c$msr_id)
p = model_perf_plot(type = 'interval', post_draws = post_draws_uno[[msr_id]],
  probs = c(0.8, 0.95), x_label = 'Uno\'s C-index')
p +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.505, y = 8, label = 'Random', color = 'black',
    angle = 90)
ggsave(file = paste0(res_path_uno, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# IBS (ERV)
msr_id =  unique(post_draws_uno$ibrier_erv$msr_id)
p = model_perf_plot(type = 'halfeye', post_draws = post_draws_uno[[msr_id]],
  probs = c(0.9, 0.99), x_label = 'IBS (ERV)')
p +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.03, y = 8, label = 'Kaplan-Meier', color = 'black',
    angle = 90) +
  # restrict values
  xlim(c(-0.3, 0.25)) +
  # add arrows to show that distributions are more to the left
  annotate(
    geom = 'segment',
    x = -0.3, xend = -0.3,
    y = 8, yend = 8, # 8th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  ) +
  annotate(
    geom = 'segment',
    x = -0.3, xend = -0.3,
    y = 9, yend = 9, # 9th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  )
ggsave(file = paste0(res_path_uno, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# RCLL (ERV)
msr_id = unique(post_draws_uno$rcll_erv$msr_id)
p = model_perf_plot(type = 'halfeye', post_draws = post_draws_uno[[msr_id]],
  probs = c(0.9, 0.99), x_label = 'RCLL (ERV)')
p +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.007, y = 3, label = 'Kaplan-Meier', color = 'black',
    angle = 90) +
  # restrict values
  xlim(c(-0.06, 0.06)) +
  # add arrow to show that a distribution is more to the left
  annotate(
    geom = 'segment',
    x = -0.06, xend = -0.06,
    y = 9, yend = 9, # 9th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  )
ggsave(file = paste0(res_path_uno, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

## Omics ranking ----
message('## Omics ranking')

### Load bayesian models ----
#' Note: each `omic_fit` has a separate model for each omic
omic_fit_harrell_c  = readRDS(file = paste0(models_path_uno, 'omic_fit_harrell_c.rds'))
omic_fit_uno_c      = readRDS(file = paste0(models_path_uno, 'omic_fit_uno_c.rds'))
omic_fit_ibrier_erv = readRDS(file = paste0(models_path_uno, 'omic_fit_ibrier_erv.rds'))
omic_fit_rcll_erv   = readRDS(file = paste0(models_path_uno, 'omic_fit_rcll_erv.rds'))
omic_fits = list(omic_fit_harrell_c, omic_fit_uno_c, omic_fit_ibrier_erv, omic_fit_rcll_erv)

### Omics posterior diff draw ----
# Get the posterior differences (with vs without an omic) for each measure
omics_postdiff_draws_uno = list()
for(model_list in omic_fits) {
  msr_id = unique(model_list$Clinical_omic$data$measure)
  message('# Measure: ', msr_id)

  omics_postdiff_draws_uno[[msr_id]] =
    omics_postdiff_draw(model_list = model_list, seed = seed) %>%
    add_column(msr_id)
}

# clear some memory
rm(omic_fit_harrell_c, omic_fit_uno_c, omic_fit_ibrier_erv, omic_fit_rcll_erv, omic_fits)
gc()

### Plots ----
# effect size for ROPE (Region of Practical Equivalence)
size = 0.02

# Harrell's C-index
msr_id = unique(omics_postdiff_draws_uno$harrell_c$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_uno[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in Harrell\'s C-index',
  size = size, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_uno, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_uno[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_uno, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# Uno's C-index
msr_id = unique(omics_postdiff_draws_uno$uno_c$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_uno[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in Uno\'s C-index',
  size = size, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_uno, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_uno[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_uno, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# IBS (ERV)
msr_id = unique(omics_postdiff_draws_uno$ibrier_erv$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_uno[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in IBS (ERV)',
  size = size, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_uno, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_uno[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_uno, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# RCLL (ERV)
size2 = 0.01
msr_id = unique(omics_postdiff_draws_uno$rcll_erv$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_uno[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in RCLL (ERV)',
  size = size2, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_uno, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_uno[[msr_id]], size = size2)
write_csv(stats, file = paste0(res_path_uno, '/omics_ranking_ROPE_stats_' ,
  size2, '_', msr_id, '.csv'))

# RCLL ----
message('### Models tuned with RCLL')

## Model comparison ----
message('## Model comparison')

### Load bayesian models ----
model_harrell_c  = readRDS(file = paste0(models_path_rcll, 'models_fit_harrell_c.rds'))
model_uno_c      = readRDS(file = paste0(models_path_rcll, 'models_fit_uno_c.rds'))
model_ibrier_erv = readRDS(file = paste0(models_path_rcll, 'models_fit_ibrier_erv.rds'))
model_rcll_erv   = readRDS(file = paste0(models_path_rcll, 'models_fit_rcll_erv.rds'))
models = list(model_harrell_c, model_uno_c, model_ibrier_erv, model_rcll_erv)

### Models posterior draw ----
post_draws_rcll = list()
for (model in models) {
  msr_id = unique(model$data$measure)
  message('# Measure: ', msr_id)

  nd = distinct(model$data, lrn_id)
  post_draws_rcll[[msr_id]] = model_post_draw(
    stan_model = model, newdata = nd, rename_models = TRUE, seed = seed
  ) %>% add_column(msr_id)
}

# clear some memory
rm(model_harrell_c, model_uno_c, model_ibrier_erv, model_rcll_erv, models)
gc()

### Plots ----
# Harrell's C-index
msr_id = unique(post_draws_rcll$harrell_c$msr_id)
p = model_perf_plot(type = 'interval', post_draws = post_draws_rcll[[msr_id]],
  probs = c(0.8, 0.95), x_label = 'Harrell\'s C-index')
p +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.505, y = 6, label = 'Random', color = 'black',
    angle = 90) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01))
ggsave(file = paste0(res_path_rcll, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# Uno's C-index
msr_id = unique(post_draws_rcll$uno_c$msr_id)
p = model_perf_plot(type = 'interval', post_draws = post_draws_rcll[[msr_id]],
  probs = c(0.8, 0.95), x_label = 'Uno\'s C-index')
p +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.505, y = 8, label = 'Random', color = 'black',
    angle = 90)
ggsave(file = paste0(res_path_rcll, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# IBS (ERV)
msr_id = unique(post_draws_rcll$ibrier_erv$msr_id)
p = model_perf_plot(type = 'halfeye', post_draws = post_draws_rcll[[msr_id]],
  probs = c(0.9, 0.99), x_label = 'IBS (ERV)')
p +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = 0.03, y = 8, label = 'Kaplan-Meier', color = 'black',
    angle = 90) +
  # restrict values
  xlim(c(-0.32, 0.25)) +
  # add arrows to show that distributions are more to the left
  annotate(
    geom = 'segment',
    x = -0.32, xend = -0.32,
    y = 7, yend = 7, # 8th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  ) +
  annotate(
    geom = 'segment',
    x = -0.32, xend = -0.32,
    y = 8, yend = 8, # 8th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  ) +
  annotate(
    geom = 'segment',
    x = -0.32, xend = -0.32,
    y = 9, yend = 9, # 9th learner
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
  )
ggsave(file = paste0(res_path_rcll, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

# RCLL (ERV)
msr_id = unique(post_draws_rcll$rcll_erv$msr_id)
p = model_perf_plot(type = 'halfeye', post_draws = post_draws_rcll[[msr_id]],
  probs = c(0.9, 0.99), x_label = 'RCLL (ERV)')
p +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#E41A1C') +
  annotate(geom = 'text', x = -0.007, y = 3, label = 'Kaplan-Meier', color = 'black',
    angle = 90)
ggsave(file = paste0(res_path_rcll, '/model_cmp_' , msr_id, '.png'), units = 'in',
  width = 6, height = 5, dpi = 300)

## Omics ranking ----
message('## Omics ranking')

### Load bayesian models ----
#' Note: each `omic_fit` has a separate model for each omic
omic_fit_harrell_c  = readRDS(file = paste0(models_path_rcll, 'omic_fit_harrell_c.rds'))
omic_fit_uno_c      = readRDS(file = paste0(models_path_rcll, 'omic_fit_uno_c.rds'))
omic_fit_ibrier_erv = readRDS(file = paste0(models_path_rcll, 'omic_fit_ibrier_erv.rds'))
omic_fit_rcll_erv   = readRDS(file = paste0(models_path_rcll, 'omic_fit_rcll_erv.rds'))
omic_fits = list(omic_fit_harrell_c, omic_fit_uno_c, omic_fit_ibrier_erv, omic_fit_rcll_erv)

### Omics posterior diff draw ----
# Get the posterior differences (with vs without an omic) for each measure
omics_postdiff_draws_rcll = list()
for(model_list in omic_fits) {
  msr_id = unique(model_list$Clinical_omic$data$measure)
  message('# Measure: ', msr_id)

  omics_postdiff_draws_rcll[[msr_id]] =
    omics_postdiff_draw(model_list = model_list, seed = seed) %>%
    add_column(msr_id)
}

# clear some memory
rm(omic_fit_harrell_c, omic_fit_uno_c, omic_fit_ibrier_erv, omic_fit_rcll_erv, omic_fits)
gc()

### Plots ----
# effect size for ROPE (Region of Practical Equivalence)
size = 0.02

# Harrell's C-index
msr_id = unique(omics_postdiff_draws_rcll$harrell_c$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_rcll[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in Harrell\'s C-index',
  size = size, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_rcll, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_rcll[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_rcll, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# Uno's C-index
msr_id = unique(omics_postdiff_draws_rcll$uno_c$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_rcll[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in Uno\'s C-index',
  size = size, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_rcll, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_rcll[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_rcll, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# IBS (ERV)
msr_id = unique(omics_postdiff_draws_rcll$ibrier_erv$msr_id)
n_omics = nlevels(omics_postdiff_draws_rcll$ibrier_erv$omic)
ridgeline_plot(post_draws = omics_postdiff_draws_rcll[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in IBS (ERV)',
  size = size, ROPE_area = T, draw_arrow = F, pal = 'hue')  +
  # add arrow
  annotate(
    geom = 'segment',
    x = 0, xend = size,
    y = n_omics, yend = n_omics,
    linewidth = 0.7,
    arrow = arrow(ends = 'both', length = unit(0.05, 'in'))
  ) +
  annotate(
    geom = 'text', label = size,
    x = size, y = n_omics + 0.25
  )
ggsave(file = paste0(res_path_rcll, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_rcll[[msr_id]], size = size)
write_csv(stats, file = paste0(res_path_rcll, '/omics_ranking_ROPE_stats_' ,
  size, '_', msr_id, '.csv'))

# RCLL (ERV)
size2 = 0.01
msr_id = unique(omics_postdiff_draws_rcll$rcll_erv$msr_id)
ridgeline_plot(post_draws = omics_postdiff_draws_rcll[[msr_id]],
  title = 'Posterior Distribution Differences',
  subtitle = 'With vs without an omic',
  x_axis_label = 'Difference in RCLL (ERV)',
  size = size2, ROPE_area = T, draw_arrow = T, pal = 'hue')
ggsave(file = paste0(res_path_rcll, '/omics_ranking_' , msr_id, '.png'),
  units = 'in', width = 7, height = 5, dpi = 300)

stats = rope_stats(post_draws = omics_postdiff_draws_rcll[[msr_id]], size = size2)
write_csv(stats, file = paste0(res_path_rcll, '/omics_ranking_ROPE_stats_' ,
  size2, '_', msr_id, '.csv'))

# Uno's C vs RCLL ----
message('### Uno\'s C vs RCLL runing comparison')

## Get posterior differences ----
#' Compare model posterior distribution performance across test metrics, with
#' models tuned with Uno's C-index (`crank` prediction) vs tuned with RCLL
#' (`distr` prediction)
stopifnot(length(post_draws_uno) == length(post_draws_rcll))
res = list()
for(index in 1:length(post_draws_uno)) {
  # pick performance measure
  stopifnot(names(post_draws_uno)[[index]] == names(post_draws_rcll)[[index]])
  msr_id = names(post_draws_uno)[[index]]
  message('# Measure: ', msr_id)

  # get the posterior distributions of the measure across benchmarked models
  uno  = post_draws_uno[[index]]
  rcll = post_draws_rcll[[index]]
  stopifnot(nrow(uno) == nrow(rcll)) # it should be the same number of draws

  # get the model names
  stopifnot(sort(levels(uno$model)) == sort(levels(rcll$model)))
  models = levels(uno$model)

  # calculate posterior differences per model (Uno-C - RCLL)
  res[[msr_id]] = lapply(models, function(model) {
    uno_post = uno %>% filter(!!model == model) %>% pull(posterior)
    rcll_post = rcll %>% filter(!!model == model) %>% pull(posterior)
    #' SOS for interpretation => higher difference means Uno-C-derived distr is
    #' higher and since for all the metrics used higher is better, that would
    #' mean that models tuned with Uno-C should be preferred
    post_diff = uno_post - rcll_post
    tibble(model = model, post_diff = post_diff)
  }) %>%
  bind_rows() %>%
  mutate(
    model = factor(model),
    model = forcats::fct_reorder(model, .x = post_diff, .desc = TRUE)
  )
}

## Plots ----
# effect size for ROPE (Region of Practical Equivalence)
size = 0.02

for(msr_id in names(res)) {
  if (msr_id == 'harrell_c') { # Harrell's C-index
    ridgeline_plot(post_draws = res[[msr_id]],
      title = 'Posterior Distribution Differences',
      subtitle = 'Uno\'s C-index vs RCLL',
      x_axis_label = 'Difference in Harrell\'s C-index',
      ROPE_area = TRUE, draw_arrow = TRUE, size = size, pal = 'viridis')
    ggsave(file = paste0(res_path_cmp, '/' , msr_id, '.png'),
      units = 'in', width = 7, height = 5, dpi = 300)

    stats = rope_stats(post_draws = res[[msr_id]], size = size)
    write_csv(stats, file = paste0(res_path_cmp, '/ROPE_stats_' ,
      size, '_', msr_id, '.csv'))
  } else if (msr_id == 'uno_c') { # Uno's C-index
    ridgeline_plot(post_draws = res[[msr_id]],
      title = 'Posterior Distribution Differences',
      subtitle = 'Uno\'s C-index vs RCLL',
      x_axis_label = 'Difference in Uno\'s C-index',
      ROPE_area = TRUE, draw_arrow = TRUE, size = size, pal = 'viridis')
    ggsave(file = paste0(res_path_cmp, '/' , msr_id, '.png'),
      units = 'in', width = 7, height = 5, dpi = 300)

    stats = rope_stats(post_draws = res[[msr_id]], size = size)
    write_csv(stats, file = paste0(res_path_cmp, '/ROPE_stats_' ,
      size, '_', msr_id, '.csv'))
  } else if (msr_id == 'ibrier_erv') { # IBS (ERV)
    ridgeline_plot(post_draws = res[[msr_id]],
      title = 'Posterior Distribution Differences',
      subtitle = 'Uno\'s C-index vs RCLL',
      x_axis_label = 'Difference in IBS (ERV)',
      ROPE_area = TRUE, draw_arrow = FALSE, size = size, pal = 'viridis')
    ggsave(file = paste0(res_path_cmp, '/' , msr_id, '.png'),
      units = 'in', width = 7, height = 5, dpi = 300)

    stats = rope_stats(post_draws = res[[msr_id]], size = size)
    write_csv(stats, file = paste0(res_path_cmp, '/ROPE_stats_' ,
      size, '_', msr_id, '.csv'))
  } else if (msr_id == 'rcll_erv') { # RCLL (ERV)
    # change difference to: RCLL - Uno-C (only for figure)
    pds = res[[msr_id]] %>% mutate(post_diff = -post_diff)

    ridgeline_plot(post_draws = pds,
      title = 'Posterior Distribution Differences',
      subtitle = 'RCLL vs Uno\'s C-index', # changed this as well
      x_axis_label = 'Difference in RCLL (ERV)',
      ROPE_area = TRUE, draw_arrow = TRUE, size = size, pal = 'viridis') +
      xlim(c(-0.06, 0.06)) +
      # add arrows to show that distribution is more to the right
      annotate(
        geom = 'segment',
        x = 0.055, xend = 0.055,
        y = 9, yend = 9, # 9th learner
        linewidth = 0.7,
        arrow = arrow(ends = 'both', length = unit(0.1, 'in'))
      )
    ggsave(file = paste0(res_path_cmp, '/' , msr_id, '.png'),
      units = 'in', width = 7, height = 5, dpi = 300)

    stats = rope_stats(post_draws = res[[msr_id]], size = size)
    write_csv(stats, file = paste0(res_path_cmp, '/ROPE_stats_' ,
      size, '_', msr_id, '.csv'))
  }
}
