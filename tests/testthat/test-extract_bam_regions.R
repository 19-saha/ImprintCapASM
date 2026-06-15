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

test_that("rejects missing bam_file", {
  expect_error(
    extract_bam_regions(
      bam_file    = "nonexistent.bam",
      bed_file    = "dummy.bed",
      sample_type = "control"
    )
  )
})
