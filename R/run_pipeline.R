#' Run the Full ImprintCapASM Pipeline for a Batch of Samples
#'
#' Convenience wrapper that executes the complete three-step pipeline:
#' (1) \code{\link{extract_bam_regions}}, (2) \code{\link{prepare_cpg_snp_input}},
#' and (3) \code{\link{ASM}} for every sample in a directory.
#' Controls and patients are always processed independently using their
#' respective \code{filter_cpgs} reference files.
#'
#' @param bam_dir      Character. Path to the directory containing \code{.bam}
#'   files (and their \code{.bai} indices).
#' @param snp_dir      Character. Path to the directory containing SNP
#'   \code{.out} files produced by your bisulfite SNP caller.
#' @param meth_dir     Character. Path to the directory containing CpG
#'   methylation files (Bismark \code{.cov} or similar).
#' @param cpg_ref_file Character. Path to the cohort-matched
#'   \code{filter_cpgs} reference \code{.xlsx} file.
#'   Use \code{inst/extdata/filter_cpgs_ctrl.xlsx} for controls and
#'   \code{inst/extdata/filter_cpgs_pat.xlsx} for patients.
#' @param output_dir   Character. Directory where all output files are written.
#'   Created automatically if it does not exist.
#' @param sample_type  Character. Either \code{"control"} or \code{"patient"}.
#'   Must match the cohort of samples in \code{bam_dir}.
#' @param bed_file     Character. Path to the BED file defining the 41 DMR
#'   regions used by \code{\link{extract_bam_regions}}.
#' @param min_depth    Integer. Minimum read depth filter passed to
#'   \code{\link{prepare_cpg_snp_input}}. Default \code{20L}.
#' @param window_bp    Integer. Window in base pairs around each CpG passed to
#'   \code{\link{prepare_cpg_snp_input}}. Default \code{60L}.
#' @param overwrite    Logical. If \code{TRUE}, re-runs extraction even when
#'   output BAM files already exist. Default \code{FALSE}.
#'
#' @return Invisibly returns a named list of ASM result objects, one per sample.
#'
#' @examples
#' if (FALSE) {
#'   ctrl_ref <- system.file("extdata", "filter_cpgs_ctrl.xlsx",
#'                           package = "ImprintCapASM")
#'   bed      <- system.file("extdata", "dmr_regions.bed",
#'                           package = "ImprintCapASM")
#'
#'   results <- run_pipeline(
#'     bam_dir      = "/data/controls/bam",
#'     snp_dir      = "/data/controls/snps",
#'     meth_dir     = "/data/controls/meth",
#'     cpg_ref_file = ctrl_ref,
#'     output_dir   = "/data/controls/output",
#'     sample_type  = "control",
#'     bed_file     = bed
#'   )
#' }
#'
#' @export
run_pipeline <- function(bam_dir,
                         snp_dir,
                         meth_dir,
                         cpg_ref_file,
                         output_dir,
                         sample_type  = c("control", "patient"),
                         bed_file,
                         min_depth    = 20L,
                         window_bp    = 60L,
                         overwrite    = FALSE) {
  
  sample_type <- match.arg(sample_type)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Collect input files
  bam_files  <- sort(list.files(bam_dir,  pattern = "_markdup\\.bam$",
                                full.names = TRUE))
  snp_files  <- sort(list.files(snp_dir,  pattern = "_all\\.SNPs\\.out$",
                                full.names = TRUE))
  meth_files <- sort(list.files(meth_dir, pattern = "\\.cov(\\.gz)?$",
                                full.names = TRUE))
  
  if (length(bam_files) == 0L)
    stop("No *_markdup.bam files found in: ", bam_dir)
  if (length(snp_files) == 0L)
    stop("No *_all.SNPs.out files found in: ", snp_dir)
  if (length(meth_files) == 0L)
    stop("No .cov/.cov.gz files found in: ", meth_dir)
  
  n <- length(bam_files)
  cat(sprintf("Found %d sample(s) [%s]\n", n, sample_type))
  
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    sample_id <- sub("_markdup\\.bam$", "", basename(bam_files[i]))
    cat(sprintf("\nSample %d/%d: %s\n", i, n, sample_id))
    
    # Step 1: Extract BAM regions
    cat("  Step 1/3  extract_bam_regions\n")
    extracted_bam <- extract_bam_regions(
      bam_file    = bam_files[i],
      bed_file    = bed_file,
      output_dir  = output_dir,
      overwrite   = overwrite,
      sample_type = sample_type
    )
    
    # Step 2: Prepare CpG-SNP input table
    cat("  Step 2/3  prepare_cpg_snp_input\n")
    cpg_snp_file <- file.path(output_dir,
                              paste0(sample_id, "_cpg_snp.xlsx"))
    
    prepare_cpg_snp_input(
      snp_file     = snp_files[i],
      meth_file    = meth_files[i],
      cpg_ref_file = cpg_ref_file,
      output_file  = cpg_snp_file,
      min_depth    = min_depth,
      window_bp    = window_bp,
      sample_type  = sample_type
    )
    
    # Step 3: ASM analysis
    cat("  Step 3/3  ASM\n")
    asm_output_file <- file.path(output_dir,
                                 paste0(sample_id, "_ASM_results.xlsx"))
    
    asm_result <- ASM(
      cpg_snp_file     = cpg_snp_file,
      sam_file         = extracted_bam,
      filter_cpgs_file = cpg_ref_file,
      output_file      = asm_output_file,
      sample_type      = sample_type
    )
    
    results[[i]] <- asm_result
    names(results)[i] <- sample_id
  }
  
  cat(sprintf("\nPipeline complete. Results for %d sample(s) written to: %s\n",
              n, output_dir))
  
  invisible(results)
}
