# tcga-survmob

This repository is a continuation of [paad-survival-bench](https://github.com/bblodfon/paad-survival-bench) - now the focus is to to use the developed [survmob](https://github.com/bblodfon/survmob) R library and many TCGA multimodal datasets to benchmark survival ML models and analyze the output results.

- TO ADD: link for the full downloaded datasets.

The data download, pre-processing, feature selection, benchmarking and analysis scripts can be found in each of the folders corresponding to the TCGA study names, e.g. `PAAD`, `BLCA`, `OV`, etc.
The general order of execution to get the benchmark results per cancer type is:
- `download_data.R`
- `preprocess.R`
- `data_split.R`
- `efs.R`
- `task_subset.R`
- `benchmark.R`
