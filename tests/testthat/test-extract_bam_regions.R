library(ImprintCapASM)

test_that("rejects invalid sample_type", {
  expect_error(
    extract_bam_regions(
      bam_file    = "dummy.bam",
      bed_file    = "dummy.bed",
      sample_type = "CONTROL"
    ),
    regexp = "should be one of"
  )
})

test_that("errors when samtools is missing or bam_file not found", {
  expect_error(
    extract_bam_regions(
      bam_file    = "nonexistent.bam",
      bed_file    = "dummy.bed",
      sample_type = "control"
    ),
    regexp = "samtools|not found|does not exist|cannot open"
  )
})
