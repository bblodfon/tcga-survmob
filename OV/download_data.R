#' => Download omics and clinical data
#' => Filter/Simplify omics data
#' => Find suitable omics subset for survival prediction task
#' => Save data to files
#'
#' Execute command: `Rscript OV/download_data.R` (not via Rstudio!)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)
library(tibble)
library(readr)

# Download Data ----

# check available disease codes
data('diseaseCodes', package = 'TCGAutils')
head(diseaseCodes)

# Ovarian serous cystadenocarcinoma
disease_code = 'OV'

# see all assays (no download => `dry.run = TRUE`)
curatedTCGAData::curatedTCGAData(diseaseCode = disease_code, version = '2.0.1',
  dry.run = TRUE)

# Download selected assays
my_assays = c('*RPPA*', '*Methylation_methyl27*', '*miRNASeqGene*',
  '*RNASeq2GeneNorm*', '*CNASNP*')
cancer_data = curatedTCGAData::curatedTCGAData(diseaseCode = disease_code,
  assays = my_assays, version = '2.0.1', dry.run = FALSE)
print(cancer_data)

# number of patients (590)
print(paste0('#Patients: ', nrow(colData(cancer_data))))
head(cancer_data$patientID) # patientIDs

# Data filtering ----

# check sample types
data('sampleTypes', package = 'TCGAutils')
print(sampleTypes) # Primary Solid Tumor = '01'

# what sample types do we have in the dataset across different data types?
print(TCGAutils::sampleTables(cancer_data))
## sometimes the above function doesn't work but the explicit one below works:
check_samples = function(cancer_data) {
  lapply(colnames(cancer_data), function(x) {
    scodes = TCGAutils::TCGAbarcode(x, participant = FALSE, sample = TRUE)
    scodes = substr(scodes, 1L, 2L)
    table(unname(scodes))
  })
}
print(check_samples(cancer_data))

# keep only primary tumor samples
cancer_data = TCGAutils::TCGAprimaryTumors(cancer_data)
print(check_samples(cancer_data)) # verify: (only '01')

MultiAssayExperiment::replicates(cancer_data) # no replicates

#' (SOS) hack for CNA data - if you perform more filtering in later steps,
#' `simplifyTCGA()` will fail, so do it now and subset to common patients later
cna = cancer_data[,,'OV_CNASNP-20160128']
#' Simplify features: from genomics ranges to gene symbols
cna_simplified = TCGAutils::simplifyTCGA(cna)

# keep only patients with the major `histological_type`
# serous cystadenocarcinoma => 579 patients
tbl_hist = table(cancer_data$histological_type, useNA = 'ifany')
print(tbl_hist)
major_hist_type = names(sort(tbl_hist, decreasing = TRUE))[1]
print(paste0('Major histological type: ', major_hist_type))

logical_vec = cancer_data$histological_type == major_hist_type
cancer_data = cancer_data[, logical_vec & !is.na(logical_vec), ]
print(paste0('#Patients: ', nrow(colData(cancer_data)))) # number of patients (579)

# keep only patients with 'ovary' tumor_tissue_site
tbl_tts = table(cancer_data$tumor_tissue_site, useNA = 'ifany')
print(tbl_tts)
major_tts = names(sort(tbl_tts, decreasing = TRUE))[1]
print(paste0('Most prominent tumor tissue site: ', major_tts))

logical_vec = cancer_data$tumor_tissue_site == major_tts
cancer_data = cancer_data[, logical_vec & !is.na(logical_vec), ]
print(paste0('#Patients: ', nrow(colData(cancer_data)))) # number of patients (574)

# keep only patients that have (time > 0)
death = cancer_data$days_to_death
cens  = cancer_data$days_to_last_followup
stopifnot(length(death) == length(cens))
time = ifelse(is.na(death), cens, death)
any(is.na(time)) # TRUE
cancer_data = cancer_data[, time > 0 & !is.na(time), ]
print(paste0('#Patients: ', nrow(colData(cancer_data)))) # number of patients (569)

# Decide which assays to keep
#' Criterion: as many common samples as possible + as many omics (>=4) as possible!
#' `pm` => **patient membership** 1/0 matrix - patient either has a data on a
#' particular omic (1) or not (0)
pm = do.call(
  function(...) { data.frame(..., check.names = TRUE) },
  lapply(mapToList(sampleMap(cancer_data)),
    function(minimap) {
      rownames(colData(cancer_data)) %in% minimap[['primary']] * 1L
    }
  )
)
rownames(pm) = rownames(colData(cancer_data)) # patient IDs as rownames
colnames(pm) = stringr::str_replace(string = colnames(pm), pattern = '\\.',
  replacement = '-') # change '.' back to '-' in column names
head(pm)

powcs = survmob::powerset_icounts(pm)

# checks number of omics in a combo and intersection counts
for (num in 4:5) {
  row = powcs %>%
    arrange(desc(intersect_count)) %>%
    select(!combo_name) %>%
    filter(n_omics == num) %>%
    slice(1)

  print(paste0(
    "#Omics: ", row$n_omics,
    ", Common samples: ", row$intersect_count,
    ", Omics: ", paste0(row$omic_combo[[1]], collapse = " / ")
  ))
}

# keep the 5 omics which have the largest amount of patients with all modalities
# mRNA, miRNA, CNA, Methyl, RPPA
assays_to_keep =
  powcs %>%
  arrange(desc(intersect_count)) %>%
  filter(n_omics == 5) %>%
  slice(1) %>%
  pull(omic_combo) %>%
  nth(1)

print(paste0('Assays to keep: ', paste0(assays_to_keep, collapse = ' / ')))

cancer_data = cancer_data[,,assays_to_keep]
print(cancer_data)

# keep only the common patients across assays
cancer_data = MultiAssayExperiment::intersectColumns(cancer_data)
print(paste0('#Patients: ', nrow(colData(cancer_data)))) # number of patients (221)

# check that patientIDs are in correct order across all assays (SOS)
for(i in 1:length(cancer_data)) {
  stopifnot(substr(colnames(cancer_data[[i]]), 1, 12) == cancer_data$patientID)
}

# Save Omics data ----
## mRNA ----
mrna = cancer_data[,,'OV_RNASeq2GeneNorm-20160128']
mrna_mat = t(assay(mrna))

mrna_mat[1:5,1:5] # count data
print(dim(mrna_mat))
readr::write_csv(x = as_tibble(mrna_mat), file = paste0(disease_code, '/data/mRNA.csv'))

## CNA ----
#' Subset the `cna` MultiAssayExperiment to the common patient ids
ids_cna  = rownames(colData(cna_simplified))
ids_orig = rownames(colData(cancer_data))

log_vec = ids_cna %in% ids_orig
cna_simplified = cna_simplified[, log_vec, ]

# check: same patients ids
all(rownames(colData(cna_simplified)) == rownames(colData(cancer_data)))

# get the assay
cna_mat = t(assay(cna_simplified))
cna_mat[1:5,1:5] # CNA log2 values (positive => amplifications, negative => deletions)
print(dim(cna_mat))
readr::write_csv(x = as_tibble(cna_mat), file = paste0(disease_code, '/data/CNA.csv'))

## miRNA ----
mirna = cancer_data[,,'OV_miRNASeqGene-20160128']
mirna_mat = t(assay(mirna))

mirna_mat[1:5,1:5] # counts
print(dim(mirna_mat))
readr::write_csv(x = as_tibble(mirna_mat), file = paste0(disease_code, '/data/miRNA.csv'))

## Methylation ----
meth = cancer_data[,,'OV_Methylation_methyl27-20160128']
methyl_mat = t(assay(meth))
methyl_mat[1:5,1:5] # beta values
print(dim(methyl_mat))

# Remove CpGs associated with sex chromosomes (X,Y)
methyl_metadata = rowData(experiments(meth)[[1]]) %>%
  as_tibble(rownames = 'CpGs')
methyl_metadata %>% count(Chromosome) %>% print(n = 25)
somatic_cpgs = methyl_metadata %>%
  filter(!is.na(Gene_Symbol), !Chromosome %in% c('X','Y',NA)) %>%
  pull(CpGs)
all(somatic_cpgs %in% colnames(methyl_mat)) # check

# subset CpGs
methyl_mat = methyl_mat[, somatic_cpgs]
print(dim(methyl_mat))
readr::write_csv(x = as_tibble(methyl_mat), file = paste0(disease_code, '/data/methyl.csv'))

## RPPA ----
rppa = cancer_data[,,'OV_RPPAArray-20160128']
rppa_mat = t(assay(rppa))

rppa_mat[1:5,1:5]
print(dim(rppa_mat))
readr::write_csv(x = as_tibble(rppa_mat), file = paste0(disease_code, '/data/RPPA.csv'))

## Clinical ----
# Level 4 clinical & pathological data
clinical_vars = TCGAutils::getClinicalNames(diseaseCode = disease_code)
length(clinical_vars) # 12

clin_tbl = colData(cancer_data)[, clinical_vars] %>%
  tibble::as_tibble(rownames = 'patient_id') %>%
  dplyr::rename(status = vital_status) %>%
  readr::type_convert(guess_integer = TRUE) %>%
  dplyr::mutate(time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death)) %>%
  dplyr::rename(age = years_to_birth) %>%
  dplyr::select(-c('days_to_death', 'days_to_last_followup',
    'tumor_tissue_site', 'histological_type', 'gender', # same value for all patients
    'date_of_initial_pathologic_diagnosis',
    'karnofsky_performance_score', 'residual_tumor', # many NAs
    'radiation_therapy', # almost all '0', quasi-constant variable
    'ethnicity'
  ))
clin_tbl # only 'age' remains

readr::write_csv(x = clin_tbl, file = paste0(disease_code, '/data/clinical.csv'))
