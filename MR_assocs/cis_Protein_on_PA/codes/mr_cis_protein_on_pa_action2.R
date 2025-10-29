
### mr of cis protein variants on pa traits

library(TwoSampleMR)
library(data.table)

# --- Parse Command Line Args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("No command-line arguments supplied. Using defaults for testing.\n")
    trait_name <- "acc_ave"
    pval_threshold <- 5e-8
} else {
    trait_name <- args[1]
    pval_threshold <- as.numeric(args[2])
}

# -----------------------------
# Prepare suffix for directories and filenames
# -----------------------------
pval_suffix <- format(pval_threshold, scientific = TRUE)

# -----------------------------
# Harmonized files directory depends on threshold
# -----------------------------
harmonized_files_dir <- paste0(
    "/xdisk/yann/arora/harmonized_cis_protein_on_pa_action2_",
    pval_suffix,
    "/"
)

#harmonized_files_dir <- "/xdisk/yann/arora/harmonized_cis_protein_on_pa_action2/"
output_dir <- "/xdisk/yann/arora/mr_cis_protein_on_pa_action2/"
out_file <- file.path(
    output_dir, paste0("all_cis_proteins_on_", trait_name, "_", pval_threshold, ".csv")
)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

trait_harmonized_files <- list.files(
    harmonized_files_dir,
    pattern = paste0("^harmonized_cis_.*_on_", trait_name, "\\.txt.*$"),
    full.names = TRUE
)

cat("Found", length(trait_harmonized_files), "harmonized files for trait:", trait_name, "\n")

# --- Output Header Template ---
final_headers <- c(
    "exposure", "outcome",
    "Inverse_variance_weighted_nsnp", "Inverse_variance_weighted_b", "Inverse_variance_weighted_se", "Inverse_variance_weighted_pval",
    "Weighted_median_nsnp", "Weighted_median_b", "Weighted_median_se", "Weighted_median_pval",
    "Weighted_mode_nsnp", "Weighted_mode_b", "Weighted_mode_se", "Weighted_mode_pval",
    "Simple_mode_nsnp", "Simple_mode_b", "Simple_mode_se", "Simple_mode_pval",
    "MR_Egger_nsnp", "MR_Egger_b", "MR_Egger_se", "MR_Egger_pval",
    "Wald_ratio_nsnp", "Wald_ratio_b", "Wald_ratio_se", "Wald_ratio_pval"
)
fwrite(setNames(data.table(matrix(nrow = 0, ncol = length(final_headers))), final_headers),
       out_file, sep = ",", col.names = TRUE)

format_output <- function(dt) {
    dt <- dt[, intersect(names(dt), final_headers), with = FALSE]
    for (col in setdiff(final_headers, names(dt))) dt[[col]] <- NA
    setcolorder(dt, final_headers)
    dt
}

run_mr_analysis <- function(harmonized_file, pval_threshold) {
    base <- basename(harmonized_file)
    split_fields <- strsplit(base, "_on_")[[1]]
    exposure_name <- sub("^harmonized_cis_", "", split_fields[1])
    exposure_name <- sub("\\.txt.*$", "", exposure_name)
    outcome_name <- sub("\\.txt.*$", "", split_fields[2])
    skeleton <- as.data.table(as.list(setNames(rep(NA, length(final_headers)), final_headers)))
    skeleton$exposure <- exposure_name
    skeleton$outcome <- outcome_name

    dat <- suppressWarnings(tryCatch(fread(harmonized_file, sep = "\t"), error = function(e) NULL))
    if (is.null(dat) || nrow(dat) == 0 || !("pval.exposure" %in% names(dat))) return(skeleton)

    # --- Robust SNP filtering ---
    keep_col <-
        if ("mr_keep.exposure" %in% names(dat)) "mr_keep.exposure"
        else if ("mr_keep" %in% names(dat)) "mr_keep"
        else NULL
    if (!is.null(keep_col)) {
        dat_filt <- dat[
            (pval.exposure < pval_threshold) &
            (get(keep_col) == TRUE) &
            (remove == FALSE)
        ]
    } else {
        dat_filt <- dat[(pval.exposure < pval_threshold) & (remove == FALSE)]
    }
    dat_filt <- unique(dat_filt, by = "SNP")

    if (nrow(dat_filt) == 0) return(skeleton)

    # --- Wald Ratio calculation for single SNP (with NAs handled!) ---
    if (nrow(dat_filt) == 1) {
        # Extract numeric values safely, allow for NA
        b_exp <- as.numeric(dat_filt$beta.exposure)
        b_out <- as.numeric(dat_filt$beta.outcome)
        se_exp <- as.numeric(dat_filt$se.exposure)
        se_out <- as.numeric(dat_filt$se.outcome)
        output <- skeleton
        output$Wald_ratio_nsnp <- 1
        # Only fill Wald stats if all present
        if (!any(is.na(c(b_exp, b_out, se_exp, se_out)))) {
            wald_b <- b_out / b_exp
            wald_se <- sqrt(
                (se_out^2 / b_exp^2) + ((b_out^2 * se_exp^2) / b_exp^4)
            )
            wald_p <- 2 * pnorm(-abs(wald_b / wald_se))
            output$Wald_ratio_b <- wald_b
            output$Wald_ratio_se <- wald_se
            output$Wald_ratio_pval <- wald_p
        }
        return(format_output(output))
    }

    # --- MR for multiple SNPs (>1) ---
    if (nrow(dat_filt) > 1) {
        mr_results <- tryCatch({
            mr(dat_filt)
        }, error = function(e) NULL)
        if (is.null(mr_results) || nrow(mr_results) == 0) return(skeleton)

        mr_results <- as.data.table(mr_results)
        mr_results[, method := gsub("Inverse variance weighted", "Inverse_variance_weighted", method)]
        mr_results[, method := gsub("Weighted median", "Weighted_median", method)]
        mr_results[, method := gsub("Weighted mode", "Weighted_mode", method)]
        mr_results[, method := gsub("Simple mode", "Simple_mode", method)]
        mr_results[, method := gsub("MR Egger", "MR_Egger", method)]
        mr_results[, exposure := exposure_name]
        mr_results[, outcome := outcome_name]

        wide <- function(col) {
            w <- tryCatch({
                dcast(mr_results, exposure + outcome ~ method, value.var = col)
            }, error = function(e) NULL)
            if (!is.null(w)) setnames(w,
                old = setdiff(names(w), c("exposure", "outcome")),
                new = paste0(setdiff(names(w), c("exposure", "outcome")), "_", col)
            )
            w
        }
        nsnp_wide <- wide("nsnp")
        b_wide <- wide("b")
        se_wide <- wide("se")
        pval_wide <- wide("pval")
        wide_tables <- Filter(Negate(is.null), list(nsnp_wide, b_wide, se_wide, pval_wide))
        if (length(wide_tables) == 0) return(skeleton)
        wide_res <- Reduce(function(x, y) merge(x, y, by = c("exposure", "outcome"), all = TRUE), wide_tables)
        wide_res <- format_output(wide_res)
        return(wide_res)
    }

    return(skeleton)
}

cat("Starting MR analysis for", length(trait_harmonized_files), "files...\n")
cat("Results will be written incrementally to CSV file:", out_file, "\n")

for (i in seq_along(trait_harmonized_files)) {
    h_file <- trait_harmonized_files[i]
    if (i %% 50 == 0 || i == length(trait_harmonized_files)) cat("Processed", i, "of", length(trait_harmonized_files), "files\n")
    res <- run_mr_analysis(h_file, pval_threshold)
    res <- format_output(res)
    fwrite(res, out_file, sep = ",", append = TRUE, col.names = FALSE)
}

cat("Completed MR analysis for all", length(trait_harmonized_files), "files.\n")
cat("Final CSV results written to:", out_file, "\n")
if (file.exists(out_file)) {
    final_data <- fread(out_file, sep = ",", nrows = 5)
    cat("Sample output for verification:\n")
    print(final_data)
}
cat("all completed successfully")

