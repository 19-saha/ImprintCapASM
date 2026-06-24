#' Run the Full ImprintCapASM Pipeline
#'
#' Convenience wrapper that executes the complete three-step pipeline:
#' (1) \code{\link{prepare_cpg_snp_input}}, (2) \code{\link{extract_bam_regions}},
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
#'   Use the file bundled in \code{inst/extdata/filter_cpgs_ctrl.xlsx}
#'   for controls and \code{inst/extdata/filter_cpgs_pat.xlsx} for patients.
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
#' @param verbose      Logical. If \code{TRUE}, progress messages are printed
#'   to the console via \code{message()}. Default \code{TRUE}.
#'   Set to \code{FALSE} or wrap the call in \code{suppressMessages()} to
#'   silence all output.
#'
#' @return Invisibly returns a named list of ASM result objects, one per
#'   sample. Each element is the list returned by \code{\link{ASM}},
#'   containing \code{$asm}, \code{$snp_cpg}, and \code{$meth_summary}.
#'
#' @examples
#' \donttest{
#' if (nchar(Sys.which("samtools")) > 0L) {
#'   extdata    <- system.file("extdata", package = "ImprintCapASM")
#'   cpg_ref    <- file.path(extdata, "example_filter_cpgs.xlsx")
#'   bed        <- file.path(extdata, "example_cpg_snp.bed")
#'   output_dir <- tempdir()
#'
#'   results <- run_pipeline(
#'     bam_dir      = extdata,
#'     snp_dir      = extdata,
#'     meth_dir     = extdata,
#'     cpg_ref_file = cpg_ref,
#'     output_dir   = output_dir,
#'     sample_type  = "control",
#'     bed_file     = bed
#'   )
#'   head(results[[1]]$asm)
#' }
#' }
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
                         overwrite    = FALSE,
                         verbose      = TRUE) {
  
  sample_type <- match.arg(sample_type)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Collect input files
  bam_files  <- sort(list.files(bam_dir,  pattern = "_markdup\\.bam$",
                                full.names = TRUE))
  snp_files  <- sort(list.files(snp_dir,  pattern = "_all\\.SNPs\\.out$",
                                full.names = TRUE))
  meth_files <- sort(list.files(meth_dir, pattern = "\\.cov(\\.gz)?$|\\.txt$",
                                full.names = TRUE))
  
  if (length(bam_files) == 0L)
    stop("No *_markdup.bam files found in: ", bam_dir)
  if (length(snp_files) == 0L)
    stop("No *_all.SNPs.out files found in: ", snp_dir)
  if (length(meth_files) == 0L)
    stop("No .cov/.cov.gz/.txt methylation files found in: ", meth_dir)
  
  n <- length(bam_files)
  if (verbose)
    message(sprintf("Found %d sample(s) [%s]", n, sample_type))
  
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    sample_id <- sub("_markdup\\.bam$", "", basename(bam_files[i]))
    if (verbose)
      message(sprintf("Sample %d/%d: %s", i, n, sample_id))
    
    # Step 1 ── prepare CpG-SNP input table
    if (verbose) message("  Step 1/3  prepare_cpg_snp_input")
    cpg_snp_file <- file.path(output_dir,
                              paste0(sample_id, "_cpg_snp.xlsx"))
    prepare_cpg_snp_input(
      snp_file     = snp_files[i],
      meth_file    = meth_files[i],
      cpg_ref_file = cpg_ref_file,
      output_file  = cpg_snp_file,
      min_depth    = min_depth,
      window_bp    = window_bp,
      sample_type  = sample_type,
      verbose      = verbose
    )
    
    # Step 2 ── extract BAM regions
    if (verbose) message("  Step 2/3  extract_bam_regions")
    extracted_bam <- extract_bam_regions(
      bam_file    = bam_files[i],
      bed_file    = bed_file,
      output_dir  = output_dir,
      overwrite   = overwrite,
      sample_type = sample_type,
      verbose     = verbose
      
    )
    
    # Step 3 ── ASM analysis
    if (verbose) message("  Step 3/3  ASM")
    asm_output_file <- file.path(output_dir,
                                 paste0(sample_id, "_ASM_results.xlsx"))
    results[[i]] <- ASM(
      cpg_snp_file     = cpg_snp_file,
      sam_file         = extracted_bam,
      filter_cpgs_file = cpg_ref_file,
      output_file      = asm_output_file,
      sample_type      = sample_type
    )
    names(results)[i] <- sample_id
  }
  
  # ── 3. Done ────────────────────────────────────────────────────────────────
  if (verbose)
    message(sprintf(
      "\u2713 Pipeline complete. Results for %d sample(s) written to: %s",
      n, output_dir
    ))
  
  invisible(results)
}