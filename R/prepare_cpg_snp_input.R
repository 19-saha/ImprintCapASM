utils::globalVariables(c(
  "chr", "5_location", "3_location", "DMR", "start", "end",
  "AD", "total_depth", "ref_depth", "alt_depth", "GT", "pos",
  "REF", "ALT", "DP", "i.start", "w_meth", "w_cov", "c_meth",
  "c_cov", "meth_frac", "total_cov", "total_meth", "context",
  "cpg_pos", "snp_pos", "win_start", "win_end", "x.cpg_pos",
  "x.DMR", "i.snp_pos", "i.DMR", "i.sample_id", "i.REF", "i.ALT",
  "i.GT", "i.AD", "i.ref_depth", "i.alt_depth", "i.total_depth",
  "ref_DMR"
))

#' Prepare CpG-SNP Input Table
#'
#' Loads a bisulfite SNP \code{.out} file and a CpG methylation file,
#' intersects heterozygous SNPs with the 41-DMR CpG panel, and returns
#' a per-sample CpG-SNP pair table ready for input to \code{\link{ASM}}.
#'
#' @param snp_file     Character. Path to a VCF-like file ending in
#'   \code{_all.SNPs.out}.
#' @param meth_file    Character. Path to the gemBS/Bismark CpG methylation
#'   table (\code{.cov} or \code{.txt}).
#' @param cpg_ref_file Character. Path to the CpG panel reference Excel file
#'   (\code{filter_cpgs_ctrl.xlsx} or \code{filter_cpgs_pat.xlsx}).
#' @param output_file  Character or \code{NULL}. Path for the output
#'   \code{.xlsx}. Auto-named as
#'   \code{cpg_snps_CG_<sample_type>_<sample_id>.xlsx} if \code{NULL}.
#' @param min_depth    Integer. Minimum total SNP read depth. Default
#'   \code{20L}.
#' @param window_bp    Integer. Window in base pairs around each SNP used
#'   to search for CpG positions. Default \code{60L}.
#' @param sample_type  Character. Either \code{"control"} or
#'   \code{"patient"}.
#' @param verbose      Logical. If \code{TRUE} (default), progress messages
#'   are written via \code{message()}. Set to \code{FALSE} or wrap the call
#'   in \code{suppressMessages()} to silence all output.
#'
#' @return A \code{data.table} of CpG-SNP pairs (returned invisibly). The
#'   \code{.xlsx} and \code{.bed} files are written as side effects.
#'
#' @examples
#' \donttest{
#'   extdata    <- system.file("extdata", package = "ImprintCapASM")
#'   cpg_snp_tmp <- tempfile(fileext = ".xlsx")
#'
#'   result <- prepare_cpg_snp_input(
#'     snp_file     = file.path(extdata, "example_snp.out"),
#'     meth_file    = file.path(extdata, "example_cgmeth.txt"),
#'     cpg_ref_file = file.path(extdata, "example_filter_cpgs.xlsx"),
#'     output_file  = cpg_snp_tmp,
#'     sample_type  = "control"
#'   )
#'   head(result)
#' }
#'
#' @importFrom data.table as.data.table data.table setkey foverlaps fread setnames := fifelse
#' @importFrom readxl read_xlsx
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom utils head write.table
#' @importFrom writexl write_xlsx
#' @export
prepare_cpg_snp_input <- function(snp_file,
                                  meth_file,
                                  cpg_ref_file,
                                  output_file  = NULL,
                                  min_depth    = 20L,
                                  window_bp    = 60L,
                                  sample_type  = c("control", "patient"),
                                  verbose      = TRUE) {
  
  sample_type <- match.arg(sample_type)
  sample_id   <- sub("_all\\.SNPs\\.out$", "", basename(snp_file))
  if (verbose) message("Processing sample: ", sample_id, " [", sample_type, "]")
  
  if (verbose) message("  Loading CpG reference: ", cpg_ref_file)
  cpg_ref <- as.data.table(read_xlsx(cpg_ref_file))
  
  if ("5_location" %in% names(cpg_ref)) {
    cpg_intervals <- cpg_ref[, .(
      chr   = chr,
      start = `5_location` + 1L - window_bp,
      end   = `3_location` + 1L + window_bp,
      DMR   = DMR
    )]
  } else {
    cpg_intervals <- cpg_ref[, .(
      chr   = chr,
      start = start - window_bp,
      end   = end   + window_bp,
      DMR   = DMR
    )]
  }
  setkey(cpg_intervals, chr, start, end)
  
  if (verbose) message("  Loading SNP file: ", snp_file)
  vr      <- read.vcfR(snp_file, verbose = FALSE)
  gt_data <- extract.gt(vr, element = "GT", return.alleles = FALSE)
  ad_data <- extract.gt(vr, element = "AD", return.alleles = FALSE)
  dp_data <- extract.gt(vr, element = "DP", return.alleles = FALSE)
  
  snps <- data.table(
    chr = vr@fix[, "CHROM"],
    pos = as.integer(vr@fix[, "POS"]),
    REF = vr@fix[, "REF"],
    ALT = vr@fix[, "ALT"],
    GT  = if (!is.null(gt_data)) gt_data[, 1] else NA_character_,
    AD  = if (!is.null(ad_data)) ad_data[, 1] else NA_character_,
    DP  = if (!is.null(dp_data)) as.integer(dp_data[, 1]) else NA_integer_
  )
  snps[, `:=`(
    ref_depth = as.integer(sapply(strsplit(AD, ","), `[`, 1)),
    alt_depth = as.integer(sapply(strsplit(AD, ","), `[`, 2))
  )]
  snps[, total_depth := ref_depth + alt_depth]
  snps <- snps[GT == "0/1" & !is.na(total_depth) & total_depth >= min_depth]
  if (verbose) message("  Heterozygous SNPs (depth >= ", min_depth, "): ", nrow(snps))
  
  snp_intervals <- snps[, .(
    chr, start = pos, end = pos,
    REF, ALT, GT, AD, DP,
    ref_depth, alt_depth, total_depth
  )]
  setkey(snp_intervals, chr, start, end)
  
  cpg_with_snps <- foverlaps(
    snp_intervals, cpg_intervals,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "start", "end"),
    type = "any", nomatch = NULL
  )
  cpg_with_snps[, `:=`(
    snp_pos = i.start,
    start   = NULL, end     = NULL,
    i.start = NULL, i.end   = NULL
  )]
  cpg_with_snps[, sample_id := sample_id]
  cpg_with_snps <- unique(cpg_with_snps, by = c("chr", "snp_pos", "DMR", "sample_id"))
  if (verbose) message("  Unique SNP-DMR pairs: ", nrow(cpg_with_snps))
  
  if (verbose) message("  Loading methylation file: ", meth_file)
  meth <- fread(meth_file)
  setnames(meth,
           old = c("#CHROM", "POS", "CONTEXT",
                   "Watson-METH", "Watson-COVERAGE", "Watson-QUAL",
                   "Crick-METH",  "Crick-COVERAGE",  "Crick-QUAL"),
           new = c("chr", "pos", "context",
                   "w_meth", "w_cov", "w_qual",
                   "c_meth", "c_cov", "c_qual")
  )
  meth[, `:=`(
    w_meth = as.integer(w_meth), w_cov = as.integer(w_cov),
    c_meth = as.integer(c_meth), c_cov = as.integer(c_cov)
  )]
  meth[, `:=`(
    total_meth = fifelse(!is.na(w_meth) & !is.na(c_meth), w_meth + c_meth, NA_integer_),
    total_cov  = fifelse(!is.na(w_cov)  & !is.na(c_cov),  w_cov  + c_cov,  NA_integer_)
  )]
  meth[, meth_frac := fifelse(total_cov > 0, total_meth / total_cov, NA_real_)]
  meth <- meth[context == "CG", .(chr, pos, context, total_meth, total_cov, meth_frac)]
  setkey(meth, chr, pos)
  
  if ("5_location" %in% names(cpg_ref)) {
    cpg_ref_positions <- cpg_ref[, .(chr, cpg_pos = `5_location` + 1L, DMR)]
  } else {
    cpg_ref_positions <- cpg_ref[, .(chr, cpg_pos = start + 1L, DMR)]
  }
  setkey(cpg_ref_positions, chr, cpg_pos)
  
  cpg_with_snps[, `:=`(win_start = snp_pos - window_bp,
                       win_end   = snp_pos + window_bp)]
  
  result <- cpg_ref_positions[cpg_with_snps,
                              on = .(chr, cpg_pos >= win_start, cpg_pos <= win_end),
                              nomatch = NULL, allow.cartesian = TRUE,
                              .(chr,
                                pos       = x.cpg_pos,
                                ref_DMR   = x.DMR,
                                snp_pos   = i.snp_pos,
                                DMR       = i.DMR,
                                sample_id = i.sample_id,
                                REF       = i.REF,
                                ALT       = i.ALT,
                                GT        = i.GT,
                                AD        = i.AD,
                                ref_depth = i.ref_depth,
                                alt_depth = i.alt_depth,
                                total_depth = i.total_depth)
  ]
  result <- result[ref_DMR == DMR]
  result[, ref_DMR := NULL]
  if (verbose) message("  SNP-CpG pairs found: ", nrow(result))
  
  result <- meth[result, on = .(chr, pos), nomatch = NULL]
  result <- result[pos != snp_pos]
  result <- unique(result)
  cpg_with_snps[, c("win_start", "win_end") := NULL]
  if (verbose) message("  Final SNP+CpG rows: ", nrow(result))
  
  bed <- unique(result[, .(
    chr,
    start = snp_pos - window_bp - 1L,
    end   = snp_pos + window_bp + 1L
  )])
  bed_file <- sub("\\.xlsx$", ".bed",
                  ifelse(is.null(output_file),
                         paste0("cpg_snps_CG_", sample_type, "_", sample_id, ".xlsx"),
                         output_file))
  write.table(bed, file = bed_file, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (verbose) message("  BED file written: ", bed_file)
  
  if (is.null(output_file))
    output_file <- paste0("cpg_snps_CG_", sample_type, "_", sample_id, ".xlsx")
  writexl::write_xlsx(as.data.frame(result), output_file)
  if (verbose) message("\u2713 Excel written: ", output_file)
  
  return(invisible(result))
}
