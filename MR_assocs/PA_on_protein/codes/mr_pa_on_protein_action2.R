#!/usr/bin/env Rscript

library(TwoSampleMR)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("No command-line arguments supplied. Exiting with error.\n")
    quit(save = "no", status = 1)
}

trait_name <- args[1] # e.g. acc_425
pval_threshold <- as.numeric(args[2]) # e.g. 5e-8

harmonized_files_dir <- "/xdisk/yann/arora/harmonized_pa_on_protein_action2/"
output_dir <- "/xdisk/yann/arora/mr_pa_on_protein_action2/"
log_dir <- "/xdisk/yann/arora/bash_scripts_codes_ua/logs_mr_pa_on_protein_action2/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

# The output file for this trait/threshold
thresh_string <- ifelse(pval_threshold < 1e-7, "5e-08", "5e-06")
out_file <- file.path(output_dir, paste0(trait_name, "_on_all_proteins_", thresh_string, ".csv"))

final_headers <- c(
    "exposure", "outcome",
    "Inverse_variance_weighted_nsnp", "Inverse_variance_weighted_b", "Inverse_variance_weighted_se", "Inverse_variance_weighted_pval",
    "Weighted_median_nsnp","Weighted_median_b","Weighted_median_se","Weighted_median_pval",
    "Weighted_mode_nsnp","Weighted_mode_b","Weighted_mode_se","Weighted_mode_pval",
    "Simple_mode_nsnp","Simple_mode_b","Simple_mode_se","Simple_mode_pval",
    "MR_Egger_nsnp","MR_Egger_b","MR_Egger_se","MR_Egger_pval",
    "Wald_ratio_nsnp","Wald_ratio_b","Wald_ratio_se","Wald_ratio_pval"
)
fwrite(setNames(data.table(matrix(nrow = 0, ncol = length(final_headers))), final_headers), out_file, sep=",", col.names=TRUE)

# List of all harmonized files for this trait
harmonized_files <- list.files(
    harmonized_files_dir,
    pattern = paste0("^harmonized_", trait_name, "_on_.*\\.txt(?:\\.txt)?$"),
    full.names = TRUE
)

cat("Found", length(harmonized_files), "harmonized files for trait:", trait_name, "\n")

format_output <- function(dt) {
    dt <- dt[, intersect(names(dt), final_headers), with=FALSE]
    for (col in setdiff(final_headers, names(dt))) dt[[col]] <- NA
    setcolorder(dt, final_headers)
    dt
}

run_mr_analysis <- function(harmonized_file, pval_threshold) {
    base <- basename(harmonized_file)
    split_fields <- strsplit(base, "_on_")[[1]]
    exposure_name <- sub("^harmonized_", "", split_fields[1]) # trait
    outcome_name  <- sub("\\.txt$", "", sub("\\.txt$", "", split_fields[2])) # protein

    skeleton <- as.data.table(as.list(setNames(rep(NA, length(final_headers)), final_headers)))
    skeleton$exposure <- exposure_name
    skeleton$outcome <- outcome_name

    dat <- suppressWarnings(tryCatch(fread(harmonized_file, sep="\t"), error = function(e) NULL))
    if (is.null(dat) || nrow(dat) == 0 || !("pval.exposure" %in% names(dat))) return(skeleton)
    dat_filtered <- dat[pval.exposure < pval_threshold]
    if (nrow(dat_filtered) == 0) return(skeleton)

    # Fix: explicit Wald ratio for nsnp==1
    if (nrow(dat_filtered) == 1) {
        mr_wald <- tryCatch(mr(dat_filtered), error = function(e) NULL)
        if (!is.null(mr_wald) && nrow(mr_wald) == 1 && mr_wald$method[1] == "Wald ratio") {
            skeleton$Wald_ratio_nsnp <- 1
            skeleton$Wald_ratio_b <- mr_wald$b
            skeleton$Wald_ratio_se <- mr_wald$se
            skeleton$Wald_ratio_pval <- mr_wald$pval
            return(skeleton)
        } else {
            return(skeleton)
        }
    }

    mr_results <- tryCatch(mr(dat_filtered), error = function(e) NULL)
    if (is.null(mr_results) || nrow(mr_results) == 0) return(skeleton)
    mr_results <- as.data.table(mr_results)

    mr_results[, method := gsub("Inverse variance weighted", "Inverse_variance_weighted", method)]
    mr_results[, method := gsub("Weighted median", "Weighted_median", method)]
    mr_results[, method := gsub("Weighted mode", "Weighted_mode", method)]
    mr_results[, method := gsub("Simple mode", "Simple_mode", method)]
    mr_results[, method := gsub("MR Egger", "MR_Egger", method)]
    mr_results[, exposure := exposure_name]
    mr_results[, outcome := outcome_name]

    target_methods <- c("Inverse_variance_weighted","Weighted_median","Weighted_mode","Simple_mode","MR_Egger","Wald_ratio")
    mr_results <- mr_results[method %in% target_methods]
    required_cols <- c("exposure", "outcome", "method", "nsnp", "b", "se", "pval")
    if (!all(required_cols %in% names(mr_results))) return(skeleton)

    safe_dcast <- function(dt, val) {
        dt <- dt[method %in% target_methods]
        if (!(val %in% names(dt))) return(NULL)
        wide <- tryCatch({
            dcast(dt, exposure + outcome ~ method, value.var = val)
        }, error = function(e) NULL)
        if (is.null(wide)) return(NULL)
        old_cols <- names(wide)[-(1:2)]
        new_cols <- paste0(old_cols, "_", val)
        setnames(wide, old = old_cols, new = new_cols)
        wide
    }

    nsnp_wide <- safe_dcast(mr_results, "nsnp")
    b_wide    <- safe_dcast(mr_results, "b")
    se_wide   <- safe_dcast(mr_results, "se")
    pval_wide <- safe_dcast(mr_results, "pval")
    if (all(sapply(list(nsnp_wide, b_wide, se_wide, pval_wide), is.null))) return(skeleton)
    wide_tables <- Filter(Negate(is.null), list(nsnp_wide, b_wide, se_wide, pval_wide))
    wide_res <- if (length(wide_tables)==1) wide_tables[[1]] else
        Reduce(function(x, y) merge(x, y, by = c("exposure", "outcome"), all = TRUE), wide_tables)
    if (!is.data.table(wide_res)) wide_res <- as.data.table(wide_res)
    wide_res <- format_output(wide_res)

    # Wald ratio fallback for methods with nsnp <= 1 (excluding Wald_ratio itself)
    methods <- c("Inverse_variance_weighted", "Weighted_median", "Weighted_mode", "Simple_mode", "MR_Egger")
    for (meth in methods) {
        nsnp_col  <- paste0(meth, "_nsnp")
        b_col     <- paste0(meth, "_b")
        se_col    <- paste0(meth, "_se")
        pval_col  <- paste0(meth, "_pval")
        wald_cols <- c("Wald_ratio_b", "Wald_ratio_se", "Wald_ratio_pval")
        if (all(c(nsnp_col, wald_cols) %in% names(wide_res))) {
            idx <- which(is.na(wide_res[[nsnp_col]]) | wide_res[[nsnp_col]] <= 1)
            if (length(idx) > 0) {
                wide_res[idx, (b_col) := get("Wald_ratio_b")]
                wide_res[idx, (se_col) := get("Wald_ratio_se")]
                wide_res[idx, (pval_col) := get("Wald_ratio_pval")]
            }
        }
    }
    wide_res <- format_output(wide_res)
    wide_res
}

cat("Starting MR analysis for", length(harmonized_files), "files...\n")
cat("Results will be written incrementally to CSV file:", out_file, "\n")

for (i in seq_along(harmonized_files)) {
    h_file <- harmonized_files[i]
    if (i %% 50 == 0) cat("Processed", i, "of", length(harmonized_files), "files\n")
    res <- run_mr_analysis(h_file, pval_threshold)
    res <- format_output(res)
    fwrite(res, out_file, sep=",", append=TRUE, col.names=FALSE)
}

cat("Completed MR analysis for all", length(harmonized_files), "files.\n")
cat("Final CSV results written to:", out_file, "\n")
if (file.exists(out_file)) {
    final_data <- fread(out_file, sep=",", nrows=5)
    cat("Sample output for verification:\n")
    print(final_data)
}
cat("Script completed successfully.\n")

