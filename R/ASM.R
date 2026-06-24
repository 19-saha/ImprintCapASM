utils::globalVariables(c(
  "allele_type", "n_methylated", "n_unmethylated", "ref_allele",
  "alt_allele", "chr", "DMR", "snp_pos", "cpg_pos", "allele_order",
  ".N", ".SD", ":=", ".", "mean_methylation", "sd_methylation",
  "MSDvar", "SDvar", "Mvar", "Category", "LOWvar", "5_location",
  "Allele_Type", "Methylation_Fraction", "Padded_Sequence",
  "marked_sequence", "N", "padding_amount", "min_read_start_per_snp"
))

#' Run Allele-Specific Methylation Analysis
#'
#' Performs SNP-phased allele-specific methylation (ASM) analysis across
#' the canonical human imprinted differentially methylated regions (DMRs).
#' Reads are assigned to REF or ALT alleles based on bisulfite-aware SNP
#' detection, and per-allele CpG methylation fractions are computed.
#'
#' @param cpg_snp_file    Character. Path to the Excel output from
#'   \code{\link{prepare_cpg_snp_input}}.
#' @param sam_file        Character. Path to the extracted wide BAM file
#'   produced by \code{\link{extract_bam_regions}}.
#' @param filter_cpgs_file Character. Path to the CpG filter/variability
#'   reference Excel file (\code{filter_cpgs_ctrl.xlsx} or
#'   \code{filter_cpgs_pat.xlsx}).
#' @param output_file     Character or \code{NULL}. Path for the ASM output
#'   \code{.xlsx}. Auto-named as \code{asm_<sample_type>_<sample_id>.xlsx}
#'   if \code{NULL}.
#' @param sample_type     Character. Either \code{"control"} or
#'   \code{"patient"}.
#' @param verbose         Logical. If \code{TRUE} (default), progress
#'   messages are written via \code{message()}. Set to \code{FALSE} or
#'   wrap the call in \code{suppressMessages()} to silence all output.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{asm}{Full read-level ASM results as a \code{data.table}.}
#'   \item{snp_cpg}{Per SNP x CpG summary table.}
#'   \item{meth_summary}{Per CpG allele methylation summary table.}
#' }
#'
#' @examples
#' \dontrun{
#'   extdata     <- system.file("extdata", package = "ImprintCapASM")
#'   cpg_snp_tmp <- tempfile(fileext = ".xlsx")
#'
#'   prepare_cpg_snp_input(
#'     snp_file     = file.path(extdata, "example_snp.out"),
#'     meth_file    = file.path(extdata, "example_cgmeth.txt"),
#'     cpg_ref_file = file.path(extdata, "example_filter_cpgs.xlsx"),
#'     output_file  = cpg_snp_tmp,
#'     sample_type  = "control"
#'   )
#'
#'   result <- ASM(
#'     cpg_snp_file     = cpg_snp_tmp,
#'     sam_file         = file.path(extdata, "example_markdup.bam"),
#'     filter_cpgs_file = file.path(extdata, "example_filter_cpgs.xlsx"),
#'     output_file      = tempfile(fileext = ".xlsx"),
#'     sample_type      = "control"
#'   )
#'   head(result$snp_cpg)
#' }
#'
#' @importFrom data.table as.data.table data.table rbindlist dcast melt setorder uniqueN := .N .SD fcase
#' @importFrom readxl read_xlsx
#' @importFrom writexl write_xlsx
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_shape_manual facet_wrap scale_color_manual scale_x_continuous scale_y_continuous labs theme_bw theme element_text element_rect
#' @importFrom grDevices pdf dev.off
#' @importFrom stats sd
#' @importFrom utils write.table capture.output
#' @export
ASM <- function(cpg_snp_file,
                sam_file,
                filter_cpgs_file,
                output_file  = NULL,
                sample_type  = c("control", "patient"),
                verbose      = TRUE) {
  
  sample_type <- match.arg(sample_type)
  
  
  # Bisulfite-aware allele assignment
  
  assign_allele_bisulfite <- function(snp_base, ref_al, alt_al,
                                      strand_label, flag, md_tag) {
    allele_type     <- NA_character_
    na_allele       <- NA_character_
    assignment_note <- NA_character_
    snp_type <- paste0(ref_al, alt_al)
    b        <- toupper(snp_base)
    fwd      <- (strand_label == "Forward")
    rev      <- (strand_label == "Reverse")
    
    if (snp_type == "CA") {
      if (fwd) {
        if      (b == "C") { allele_type <- "ALT"; assignment_note <- "CA;fwd;C meth" }
        else if (b == "T") { allele_type <- "REF"; assignment_note <- "CA;fwd;C>T=REF" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "CA;fwd;A=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "C") { allele_type <- "REF"; assignment_note <- "CA;rev;C=REF RC" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "CA;rev;A=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "AC") {
      if (fwd) {
        if      (b == "A") { allele_type <- "REF"; assignment_note <- "AC;fwd;A=REF" }
        else if (b == "C") { allele_type <- "REF"; assignment_note <- "AC;fwd;C meth" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "AC;fwd;C>T=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "A") { allele_type <- "REF"; assignment_note <- "AC;rev;A=REF RC" }
        else if (b == "C") { allele_type <- "ALT"; assignment_note <- "AC;rev;C=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "CT") {
      if (rev) {
        if      (b == "C") { allele_type <- "REF"; assignment_note <- "CT;rev;C=REF RC" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "CT;rev;T=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "TC") {
      if (rev) {
        if      (b == "T") { allele_type <- "REF"; assignment_note <- "TC;rev;T=REF RC" }
        else if (b == "C") { allele_type <- "ALT"; assignment_note <- "TC;rev;C=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "GA") {
      if (fwd) {
        if      (b == "G") { allele_type <- "REF"; assignment_note <- "GA;fwd;G=REF" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "GA;fwd;A=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "AG") {
      if (fwd) {
        if      (b == "A") { allele_type <- "REF"; assignment_note <- "AG;fwd;A=REF" }
        else if (b == "G") { allele_type <- "ALT"; assignment_note <- "AG;fwd;G=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "TG") {
      if (fwd) {
        if      (b == "T") { allele_type <- "REF"; assignment_note <- "TG;fwd;T=REF" }
        else if (b == "G") { allele_type <- "ALT"; assignment_note <- "TG;fwd;G=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "T") { allele_type <- "REF"; assignment_note <- "TG;rev;T=REF RC" }
        else if (b == "G") { allele_type <- "ALT"; assignment_note <- "TG;rev;G=C meth" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "TG;rev;C>T=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "GT") {
      if (fwd) {
        if      (b == "G") { allele_type <- "REF"; assignment_note <- "GT;fwd;G=REF" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "GT;fwd;T=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "G") { allele_type <- "REF"; assignment_note <- "GT;rev;G=C meth" }
        else if (b == "A") { allele_type <- "REF"; assignment_note <- "GT;rev;C>T=REF RC" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "GT;rev;T=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "AT") {
      if (fwd) {
        if      (b == "A") { allele_type <- "REF"; assignment_note <- "AT;fwd;A=REF" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "AT;fwd;T=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "A") { allele_type <- "REF"; assignment_note <- "AT;rev;A=REF RC" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "AT;rev;T=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "TA") {
      if (fwd) {
        if      (b == "T") { allele_type <- "REF"; assignment_note <- "TA;fwd;T=REF" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "TA;fwd;A=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "T") { allele_type <- "REF"; assignment_note <- "TA;rev;T=REF RC" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "TA;rev;A=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "CG") {
      if (fwd) {
        if      (b == "C") { allele_type <- "REF"; assignment_note <- "CG;fwd;C=C meth" }
        else if (b == "T") { allele_type <- "REF"; assignment_note <- "CG;fwd;C>T=REF" }
        else if (b == "G") { allele_type <- "ALT"; assignment_note <- "CG;fwd;G=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "C") { allele_type <- "REF"; assignment_note <- "CG;rev;C=REF RC" }
        else if (b == "G") { allele_type <- "ALT"; assignment_note <- "CG;rev;G=C meth" }
        else if (b == "A") { allele_type <- "ALT"; assignment_note <- "CG;rev;C>T=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else if (snp_type == "GC") {
      if (fwd) {
        if      (b == "G") { allele_type <- "REF"; assignment_note <- "GC;fwd;G=REF" }
        else if (b == "C") { allele_type <- "ALT"; assignment_note <- "GC;fwd;C=C meth" }
        else if (b == "T") { allele_type <- "ALT"; assignment_note <- "GC;fwd;C>T=ALT" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      } else if (rev) {
        if      (b == "G") { allele_type <- "REF"; assignment_note <- "GC;rev;G=C meth" }
        else if (b == "A") { allele_type <- "REF"; assignment_note <- "GC;rev;C>T=REF RC" }
        else if (b == "C") { allele_type <- "ALT"; assignment_note <- "GC;rev;C=ALT RC" }
        else               { na_allele   <- "?";   assignment_note <- "?" }
      }
    } else {
      if      (b == toupper(ref_al)) { allele_type <- "REF"; assignment_note <- paste0(snp_type, ";", b, "=REF") }
      else if (b == toupper(alt_al)) { allele_type <- "ALT"; assignment_note <- paste0(snp_type, ";", b, "=ALT") }
      else                           { na_allele   <- "?";   assignment_note <- paste0(snp_type, "?") }
    }
    list(allele_type = allele_type, na_allele = na_allele,
         assignment_note = assignment_note)
  }
  
  
  # SNP/CpG summary table

  make_snp_cpg_table <- function(dt) {
    result <- dt[, {
      ref_idx   <- !is.na(allele_type) & allele_type == "REF"
      alt_idx   <- !is.na(allele_type) & allele_type == "ALT"
      ref_total <- sum(ref_idx, na.rm = TRUE)
      alt_total <- sum(alt_idx, na.rm = TRUE)
      ref_meth  <- sum(n_methylated[ref_idx],   na.rm = TRUE)
      alt_meth  <- sum(n_methylated[alt_idx],   na.rm = TRUE)
      ref_base  <- ref_allele[!is.na(ref_allele)][1L]
      alt_base  <- alt_allele[!is.na(alt_allele)][1L]
      list(
        sample_id     = sample_id[1],
        chr           = chr[1],
        DMR           = DMR[1],
        ref_allele    = ref_base,
        alt_allele    = alt_base,
        REF_m         = ref_meth,
        REF_um        = sum(n_unmethylated[ref_idx], na.rm = TRUE),
        ALT_m         = alt_meth,
        ALT_um        = sum(n_unmethylated[alt_idx], na.rm = TRUE),
        REF_tot       = ref_total,
        ALT_tot       = alt_total,
        MI            = (ref_meth + alt_meth) / (ref_total + alt_total),
        REF_f         = as.numeric(if (!is.na(ref_total) && ref_total > 0) round(ref_meth / ref_total, 3) else NA_real_),
        ALT_f         = as.numeric(if (!is.na(alt_total) && alt_total > 0) round(alt_meth / alt_total, 3) else NA_real_),
        ref_alt_ratio = as.numeric(if (!is.na(ref_total) && !is.na(alt_total) && alt_total > 0) round(ref_total / alt_total, 3) else NA_real_)
      )
    },
    by = .(snp_pos, cpg_pos)
    ][order(DMR, snp_pos, cpg_pos)]
    result
  }
  

  # Methylation summary table
  
  make_meth_summary <- function(dt) {
    dt2 <- dt[!allele_type %in% c("NA", "?") & !is.na(allele_type)]
    result <- dt2[,
                  .(total_reads  = .N,
                    methylated   = sum(n_methylated,   na.rm = TRUE),
                    unmethylated = sum(n_unmethylated, na.rm = TRUE),
                    meth_frac    = round(sum(n_methylated, na.rm = TRUE) / .N, 3)),
                  by = .(sample_id, snp_pos, cpg_pos, DMR, allele_type)
    ]
    result[, allele_order := fcase(
      allele_type == "REF", 1L,
      allele_type == "ALT", 2L,
      default = 3L
    )]
    setorder(result, snp_pos, cpg_pos, DMR, allele_order)
    result[, allele_order := NULL]
    result
  }
  
  
  # Line plot PDF per DMR
  make_lineplot_pdf <- function(snp_cpg_dt, pdf_file) {
    if (verbose) message("Generating DMR line plots: ", pdf_file)
    plot_data <- melt(
      snp_cpg_dt,
      id.vars       = c("cpg_pos", "snp_pos", "DMR", "sample_id",
                        "mean_methylation", "sd_methylation", "Category"),
      measure.vars  = c("REF_f", "ALT_f"),
      variable.name = "Allele_Type",
      value.name    = "Methylation_Fraction"
    )
    plot_data[, Allele_Type := ifelse(Allele_Type == "REF_f",
                                      "REF Allele", "ALT Allele")]
    dmr_list <- unique(plot_data$DMR)
    pdf(pdf_file, width = 14, height = 8)
    for (dmr in dmr_list) {
      plot_data_dmr <- plot_data[DMR == dmr & !is.na(Methylation_Fraction)]
      if (nrow(plot_data_dmr) == 0) next
      unique_cpgs <- sort(unique(plot_data_dmr$cpg_pos))
      cpg_index   <- data.table(cpg_pos   = unique_cpgs,
                                cpg_index = seq_along(unique_cpgs))
      plot_data_dmr <- merge(plot_data_dmr, cpg_index, by = "cpg_pos")
      sid <- plot_data_dmr$sample_id[1]
      p <- ggplot(plot_data_dmr,
                  aes(x = cpg_index, y = Methylation_Fraction,
                      color = Allele_Type, group = Allele_Type)) +
        geom_line(linewidth = 0.8, alpha = 0.7) +
        geom_point(aes(shape = Category), size = 2.5, alpha = 0.9) +
        scale_shape_manual(
          values   = c("MSDvar" = 17, "SDvar" = 15, "Mvar" = 18, "LOWvar" = 19),
          na.value = 1,
          name     = "CpG Variability"
        ) +
        facet_wrap(~ snp_pos, scales = "free_x", ncol = 3) +
        scale_color_manual(values = c("REF Allele" = "#FF0000",
                                      "ALT Allele" = "#0000FF")) +
        scale_x_continuous(breaks = cpg_index$cpg_index,
                           labels = cpg_index$cpg_pos,
                           expand = c(0.02, 0)) +
        scale_y_continuous(limits = c(0, 1),
                           breaks = seq(0, 1, by = 0.25),
                           labels = c("0.00","0.25","0.50","0.75","1.00")) +
        labs(title    = paste0(dmr, " \u2014 ", sid, " [", sample_type, "]"),
             subtitle = paste0("SNPs: ", uniqueN(plot_data_dmr$snp_pos)),
             x = "CpG Position", y = "Methylation Fraction", color = "Allele") +
        theme_bw(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", size = 14),
              plot.subtitle    = element_text(size = 10, colour = "grey40"),
              strip.background = element_rect(fill = "#f0f0f0"),
              strip.text       = element_text(face = "bold", size = 10),
              axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
              legend.position  = "top")
      print(p)
      if (verbose) message("  Plotted DMR:", dmr)
    }
    dev.off()
    if (verbose) message("Saved: ", pdf_file)
  }
  

  if (verbose) message("Loading CpG/SNP reference file:", cpg_snp_file)
  if (!file.exists(cpg_snp_file)) stop("CpG/SNP file not found: ", cpg_snp_file)
  cpg_snp_data <- as.data.table(read_xlsx(cpg_snp_file))
  if (verbose) message("  Rows:", nrow(cpg_snp_data))
  
  if (verbose) message("Loading CpG filter file:", filter_cpgs_file)
  if (!file.exists(filter_cpgs_file)) stop("CpG file not found: ", filter_cpgs_file)
  processed_data <- as.data.table(read_xlsx(filter_cpgs_file))
  if (verbose) message("  Rows:", nrow(processed_data))
  
  if (verbose) message("Building CpG variation lookup...")
  cpg_ref_raw  <- as.data.table(read_xlsx(filter_cpgs_file))
  sample_cols  <- grep("^(Control|Patient)_", names(cpg_ref_raw), value = TRUE)
  if (verbose) message("  Sample columns detected:", length(sample_cols))
  
  cpg_ref_raw[, mean_methylation := round(rowMeans(.SD, na.rm = TRUE), 2),
              .SDcols = sample_cols]
  cpg_ref_raw[, sd_methylation   := round(apply(.SD, 1, sd, na.rm = TRUE), 2),
              .SDcols = sample_cols]
  cpg_ref_raw[, `:=`(
    SDvar  =  sd_methylation > 5,
    Mvar   = !(mean_methylation >= 40 & mean_methylation <= 60),
    LOWvar =  (sd_methylation < 5 & mean_methylation >= 40 & mean_methylation <= 60)
  )]
  cpg_ref_raw[, MSDvar := SDvar & Mvar]
  cpg_ref_raw[, Category := fcase(
    MSDvar, "MSDvar", SDvar, "SDvar", Mvar, "Mvar", LOWvar, "LOWvar",
    default = NA_character_
  )]
  cpg_ref_raw[, Category := factor(Category,
                                   levels = c("MSDvar","SDvar","Mvar","LOWvar"))]
  cpg_ref_raw[, c("SDvar","Mvar","LOWvar","MSDvar") := NULL]
  
  cpg_var_lookup <- cpg_ref_raw[, .(
    cpg_pos          = `5_location` + 1L,
    DMR,
    mean_methylation,
    sd_methylation,
    Category
  )]
  if (verbose) message("  Variation category counts in reference:")
  print(cpg_var_lookup[!is.na(Category), .N, by = Category][order(Category)])
  if (verbose) message()
  
  if (verbose) message("Loading BAM file:", sam_file)
  if (!file.exists(sam_file)) stop("BAM file not found: ", sam_file)
  
  bam_scan <- scanBam(
    sam_file,
    param = ScanBamParam(
      what = c("qname", "flag", "rname", "pos", "seq"),
      tag  = c("MD", "YD")
    )
  )[[1]]
  
  sam_dt <- data.table(
    read_id    = bam_scan$qname,
    flag       = as.integer(bam_scan$flag),
    read_chr   = as.character(bam_scan$rname),
    read_start = as.numeric(bam_scan$pos),
    read_seq   = as.character(bam_scan$seq),
    md_tag     = sapply(bam_scan$tag$MD,
                        function(x) if (is.null(x)) NA_character_ else as.character(x)),
    yd_tag     = sapply(bam_scan$tag$YD,
                        function(x) if (is.null(x)) NA_character_ else as.character(x))
  )
  sam_dt <- sam_dt[flag %in% c(99L, 147L, 83L, 163L)]
  if (verbose) message("  Total reads loaded:", nrow(sam_dt))
  
  sample_id  <- sub("_all_wide$", "",
                    tools::file_path_sans_ext(basename(sam_file)))
  
  target_cpgs <- data.table(
    c_pos = processed_data$`5_location` + 1L,
    g_pos = processed_data$`5_location` + 2L
  )
  
  ref_chrs    <- unique(cpg_snp_data$chr)
  result_list <- vector("list", nrow(sam_dt) * 5L)
  list_idx    <- 1L
  
  if (verbose) message("Processing BAM reads...")
  for (i in seq_len(nrow(sam_dt))) {
    if (verbose && i %% 50000 == 0) message("  Processed ", i, " reads")
    
    flag       <- sam_dt$flag[i]
    read_chr   <- sam_dt$read_chr[i]
    read_start <- sam_dt$read_start[i]
    read_seq   <- sam_dt$read_seq[i]
    read_id    <- sam_dt$read_id[i]
    md_tag     <- sam_dt$md_tag[i]
    yd_tag     <- sam_dt$yd_tag[i]
    read_len   <- nchar(read_seq)
    read_end   <- read_start + read_len - 1L
    
    if (!(read_chr %in% ref_chrs)) next
    
    if (flag %in% c(99L, 147L)) {
      strand_label <- "Forward"
      flag_context <- paste0("Forward ", flag, "/147")
    } else if (flag %in% c(83L, 163L)) {
      strand_label <- "Reverse"
      flag_context <- paste0("Reverse 83/", flag)
    }
    
    combined_tags   <- paste0("MD:Z:", md_tag, ";YD:Z:", yd_tag)
    snps_this_chr   <- unique(cpg_snp_data[chr == read_chr, .(snp_pos, DMR)])
    
    for (snp_idx in seq_len(nrow(snps_this_chr))) {
      current_snp_pos <- snps_this_chr$snp_pos[snp_idx]
      current_dmr     <- snps_this_chr$DMR[snp_idx]
      if (current_snp_pos < read_start || current_snp_pos > read_end) next
      
      ref_match <- cpg_snp_data[chr == read_chr &
                                  snp_pos == current_snp_pos &
                                  DMR == current_dmr]
      if (nrow(ref_match) == 0) next
      
      ref_al     <- ref_match$REF[1]
      alt_al     <- ref_match$ALT[1]
      snp_offset <- current_snp_pos - read_start + 1L
      if (snp_offset < 1L || snp_offset > read_len) next
      
      snp_base   <- substr(read_seq, snp_offset, snp_offset)
      assignment <- assign_allele_bisulfite(snp_base, ref_al, alt_al,
                                            strand_label, flag, md_tag)
      if (is.na(assignment$allele_type) ||
          assignment$allele_type == "NA" ||
          !is.na(assignment$na_allele)) next
      
      snp_marked <- paste0(
        substr(read_seq, 1L, snp_offset - 1L),
        "[", snp_base, "]",
        substr(read_seq, snp_offset + 1L, read_len)
      )
      
      cpg_window_start <- current_snp_pos - 60L
      cpg_window_end   <- current_snp_pos + 60L
      
      for (cpg_idx in seq_len(nrow(target_cpgs))) {
        c_pos_genome <- target_cpgs$c_pos[cpg_idx]
        g_pos_genome <- target_cpgs$g_pos[cpg_idx]
        if (c_pos_genome < cpg_window_start || c_pos_genome > cpg_window_end) next
        if (c_pos_genome < read_start || g_pos_genome > read_end) next
        
        c_offset <- c_pos_genome - read_start + 1L
        g_offset <- g_pos_genome - read_start + 1L
        dinuc    <- substr(read_seq, c_offset, g_offset)
        
        if (strand_label == "Forward") {
          if      (dinuc == "CG") { is_methylated <- TRUE  }
          else if (dinuc == "TG") { is_methylated <- FALSE }
          else                    { next }
        } else if (strand_label == "Reverse") {
          if      (dinuc == "CG") { is_methylated <- TRUE  }
          else if (dinuc == "CA") { is_methylated <- FALSE }
          else                    { next }
        }
        
        result_list[[list_idx]] <- data.table(
          sample_id       = sample_id,
          sample_type     = sample_type,
          id              = read_id,
          read_sequence   = read_seq,
          read_start      = read_start,
          flag            = flag,
          flag_context    = flag_context,
          combined_tags   = combined_tags,
          strand          = strand_label,
          chr             = read_chr,
          DMR             = current_dmr,
          snp_pos         = current_snp_pos,
          cpg_pos         = c_pos_genome,
          allele_type     = assignment$allele_type,
          ref_allele      = ref_al,
          alt_allele      = alt_al,
          assignment_note = assignment$assignment_note,
          n_methylated    = as.integer(is_methylated),
          n_unmethylated  = as.integer(!is_methylated),
          meth_frac       = as.numeric(is_methylated),
          marked_sequence = snp_marked
        )
        list_idx <- list_idx + 1L
      }
    }
  }
  
  if (verbose) message("Finished processing BAM file")
  
  if (list_idx == 1L) { message("No observations found"); return(NULL) }
  
  final_results <- rbindlist(result_list[1:(list_idx - 1L)])
  final_results <- final_results[cpg_pos != snp_pos]
  
  if (verbose) message("Filtering CpG positions: removing those with REF < 20 or ALT < 20 reads...")
  cpg_counts      <- final_results[allele_type %in% c("REF","ALT"),
                                   .(n = .N), by = .(cpg_pos, snp_pos, allele_type)]
  cpg_counts_wide <- dcast(cpg_counts, cpg_pos + snp_pos ~ allele_type,
                           value.var = "n", fill = 0L)
  valid_cpg_snp <- cpg_counts_wide[
    get("REF") >= 20L &
      get("ALT") >= 20L &
      (get("REF") / get("ALT")) >= 0.5 &
      (get("REF") / get("ALT")) <= 1.5,
    .(cpg_pos, snp_pos)
  ]
  
  n_before <- uniqueN(final_results$cpg_pos)
  final_results <- final_results[valid_cpg_snp, on = .(cpg_pos, snp_pos), nomatch = 0L]
  n_after  <- uniqueN(final_results$cpg_pos)
  if (verbose) message("  CpG positions before filter:", n_before)
  if (verbose) message("  CpG positions after filter: ", n_after)
  if (verbose) message("  CpG positions removed:      ", n_before - n_after)
  
  if (verbose) message("Computing per-SNP alignment padding...")
  snp_groups    <- final_results[,
                                 .(min_read_start_per_snp = min(read_start, na.rm = TRUE)),
                                 by = .(DMR, snp_pos)]
  final_results <- merge(final_results, snp_groups, by = c("DMR","snp_pos"),
                         all.x = TRUE)
  final_results[, padding_amount := ifelse(
    is.na(min_read_start_per_snp) | is.na(read_start), 0L,
    as.integer(read_start - min_read_start_per_snp)
  )]
  final_results[, Padded_Sequence := paste0(strrep("*", padding_amount),
                                            marked_sequence)]
  final_results[, c("marked_sequence","min_read_start_per_snp",
                    "padding_amount") := NULL]
  final_results[, allele_order := fcase(
    allele_type == "REF", 1L, allele_type == "ALT", 2L, default = 3L
  )]
  setorder(final_results, DMR, read_start, snp_pos, cpg_pos, allele_order)
  final_results[, allele_order := NULL]
  final_results <- unique(
    final_results,
    by = c("sample_id", "id", "chr", "DMR", "snp_pos", "cpg_pos", "allele_type")
  )
  
  if (verbose) message("Building SNP/CpG summary table...")
  snp_cpg_table <- make_snp_cpg_table(final_results)
  if (verbose) message("Building methylation summary table...")
  meth_summary_table <- make_meth_summary(final_results)
  
  if (verbose) message("SUMMARY STATISTICS")
  if (verbose) message("  Total rows:",           nrow(final_results))
  if (verbose) message("  Unique read IDs:",      uniqueN(final_results$id))
  if (verbose) message("  Unique DMRs:",          uniqueN(final_results$DMR))
  if (verbose) message("  Unique CpG positions:", uniqueN(final_results$cpg_pos))
  if (verbose) message("  REF reads:",            sum(final_results$allele_type == "REF", na.rm = TRUE))
  if (verbose) message("  ALT reads:",            sum(final_results$allele_type == "ALT", na.rm = TRUE))
  if (verbose) message("DMR distribution:")
  print(final_results[, .N, by = DMR][order(-N)])
  
  if (verbose) message()
  
  if (verbose) message("Joining CpG variation info into output tables...")
  join_var <- function(dt) {
    merge(dt, cpg_var_lookup, by = c("cpg_pos","DMR"), all.x = TRUE, sort = FALSE)
  }
  final_results      <- join_var(final_results)
  snp_cpg_table      <- join_var(snp_cpg_table)
  meth_summary_table <- join_var(meth_summary_table)
  if (verbose) message("  Done.")
  
  if (is.null(output_file))
    output_file <- paste0("asm_", sample_type, "_", sample_id, ".xlsx")
  snp_cpg_file  <- paste0("snp_cpg_",      sample_type, "_", sample_id, ".xlsx")
  meth_sum_file <- paste0("meth_summary_", sample_type, "_", sample_id, ".xlsx")
  pdf_file      <- paste0("dmr_plots_",    sample_type, "_", sample_id, ".pdf")
  
  if (verbose) message("Writing ASM table:           ", output_file)
  write_xlsx(final_results,      output_file)
  if (verbose) message("Writing SNP/CpG summary:     ", snp_cpg_file)
  write_xlsx(snp_cpg_table,      snp_cpg_file)
  if (verbose) message("Writing methylation summary: ", meth_sum_file)
  write_xlsx(meth_summary_table, meth_sum_file)
  
  make_lineplot_pdf(snp_cpg_table, pdf_file)
  
  if (verbose) message("Done!")
  return(list(asm          = final_results,
              snp_cpg      = snp_cpg_table,
              meth_summary = meth_summary_table))
}