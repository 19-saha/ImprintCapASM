#' Extract BAM Regions Around SNP Windows
#'
#' Calls \code{samtools view} and \code{samtools sort} to subset a
#' whole-genome bisulfite BAM file to the genomic windows defined in a BED
#' file, producing a smaller indexed BAM suitable for input to
#' \code{\link{ASM}}.
#'
#' @param bam_file   Character. Path to the input \code{_markdup.bam} file
#'   (a \code{.bai} index must be present alongside it).
#' @param bed_file   Character. Path to the BED file produced by
#'   \code{\link{prepare_cpg_snp_input}}.
#' @param output_dir Character. Directory for output files. Defaults to
#'   \code{dirname(bam_file)}.
#' @param overwrite  Logical. If \code{FALSE} (default), skips extraction
#'   when the output BAM already exists.
#' @param sample_type Character. Either \code{"control"} or
#'   \code{"patient"}.
#' @param verbose    Logical. If \code{TRUE} (default), progress messages
#'   are written via \code{message()}. Set to \code{FALSE} or wrap the call
#'   in \code{suppressMessages()} to silence all output.
#'
#' @return The path to the output BAM file (returned invisibly).
#'
#' @details
#' Requires \code{samtools} to be available on the system \code{PATH}.
#' The function calls \code{samtools view -b}, \code{samtools sort}, and
#' \code{samtools index} in sequence. If no \code{.bai} index is found
#' alongside \code{bam_file}, \code{\link[Rsamtools]{indexBam}} is called
#' automatically.
#'
#' @examples
#' \donttest{
#' if (nchar(Sys.which("samtools")) > 0L) {
#'   extdata <- system.file("extdata", package = "ImprintCapASM")
#'
#'   out_bam <- extract_bam_regions(
#'     bam_file    = file.path(extdata, "example_markdup.bam"),
#'     bed_file    = file.path(extdata, "example_cpg_snp.bed"),
#'     output_dir  = tempdir(),
#'     sample_type = "control"
#'   )
#'   file.exists(out_bam)
#' }
#' }
#'
#' @importFrom Rsamtools indexBam
#' @importFrom data.table fread
#' @export
extract_bam_regions <- function(bam_file,
                                bed_file,
                                output_dir  = dirname(bam_file),
                                overwrite   = FALSE,
                                sample_type = c("control", "patient"),
                                verbose     = TRUE) {
  
  sample_type <- match.arg(sample_type)
  # Check samtools is available
  if (nchar(Sys.which("samtools")) == 0L)
    stop("samtools is not available on PATH. Please install samtools and ensure it is accessible.")
  
  sample_id   <- sub("_markdup\\.bam$", "", basename(bam_file))
  out_bam     <- file.path(output_dir,
                           paste0(sample_type, "_", sample_id, "_wide.bam"))
  
  if (file.exists(out_bam) && !overwrite) {
    if (verbose) message(" BAM already exists, skipping extraction: ", out_bam)
    return(invisible(out_bam))
  }  
  
  if (!file.exists(bam_file)) stop("BAM not found: ", bam_file)
  if (!file.exists(bed_file)) stop("BED not found: ", bed_file)
  
  if (!file.exists(paste0(bam_file, ".bai"))) {
    if (verbose) message("  Index not found, indexing BAM...")
    indexBam(bam_file)
  }
  
  bed     <- fread(bed_file, header = FALSE,
                   col.names = c("chr", "start", "end"))
  if (verbose)
    message("  Extracting ", nrow(bed), " target regions from: ",
            basename(bam_file), " [", sample_type, "]")
  
  regions  <- paste(paste0(bed$chr, ":", bed$start + 1L, "-", bed$end),
                    collapse = " ")
  tmp_bam  <- file.path(output_dir,
                        paste0(sample_type, "_", sample_id, "_wide_unsorted.bam"))
  
  cmd_view <- paste("samtools view -b", shQuote(bam_file), regions,
                    ">", shQuote(tmp_bam))
  ret <- system(cmd_view)
  if (ret != 0) stop("samtools view failed for: ", bam_file)
  
  cmd_sort <- paste("samtools sort -o", shQuote(out_bam), shQuote(tmp_bam))
  ret <- system(cmd_sort)
  if (ret != 0) stop("samtools sort failed for: ", tmp_bam)
  
  cmd_index <- paste("samtools index", shQuote(out_bam))
  ret <- system(cmd_index)
  if (ret != 0) stop("samtools index failed for: ", out_bam)
  
  file.remove(tmp_bam)
  if (verbose) message("\u2713 Written: ", out_bam)
  return(invisible(out_bam))
}
