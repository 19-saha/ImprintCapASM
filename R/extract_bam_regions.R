#' Extract BAM Regions Around SNP Windows
#'
#' @param bam_file    Path to the input markdup BAM file
#' @param bed_file    Path to the BED file from prepare_cpg_snp_input()
#' @param output_dir  Output directory (defaults to dirname of bam_file)
#' @param overwrite   If FALSE, skip if output BAM already exists
#' @param sample_type Either "control" or "patient"
#'
#' @return Path to the output BAM file (invisibly)
#' 
#' @importFrom Rsamtools indexBam
#' @importFrom data.table fread

#' @export
extract_bam_regions <- function(bam_file,
                                bed_file,
                                output_dir  = dirname(bam_file),
                                overwrite   = FALSE,
                                sample_type = c("control", "patient")) {
 
  sample_type <- match.arg(sample_type)
  sample_id   <- sub("_markdup\\.bam$", "", basename(bam_file))
  out_bam     <- file.path(output_dir,
                           paste0(sample_type, "_", sample_id, "_wide.bam"))
  
  if (file.exists(out_bam) && !overwrite) {
    cat("  BAM already exists, skipping extraction:", out_bam, "\n")
    return(invisible(out_bam))
  }
  
  if (!file.exists(bam_file))  stop("BAM not found: ", bam_file)
  if (!file.exists(bed_file))  stop("BED not found: ", bed_file)
  if (!file.exists(paste0(bam_file, ".bai"))) {
    cat("  Index not found, indexing BAM...\n")
    indexBam(bam_file)
  }
  
  bed <- fread(bed_file, header = FALSE,
               col.names = c("chr", "start", "end"))
  
  cat("  Extracting", nrow(bed), "target regions from:", basename(bam_file),
      "[", sample_type, "]\n")
  
  regions  <- paste(paste0(bed$chr, ":", bed$start + 1L, "-", bed$end),
                    collapse = " ")
  tmp_bam  <- file.path(output_dir,
                        paste0(sample_type, "_", sample_id, "_wide_unsorted.bam"))
  
  cmd_view  <- paste("samtools view -b", shQuote(bam_file), regions,
                     ">", shQuote(tmp_bam))
  ret <- system(cmd_view)
  if (ret != 0) stop("samtools view failed for: ", bam_file)
  
  cmd_sort  <- paste("samtools sort -o", shQuote(out_bam), shQuote(tmp_bam))
  ret <- system(cmd_sort)
  if (ret != 0) stop("samtools sort failed for: ", tmp_bam)
  
  cmd_index <- paste("samtools index", shQuote(out_bam))
  ret <- system(cmd_index)
  if (ret != 0) stop("samtools index failed for: ", out_bam)
  
  file.remove(tmp_bam)
  cat("Written:", out_bam, "\n")
  return(invisible(out_bam))
}