test_that("rejects invalid sample_type", {
  expect_error(
    ASM(
      cpg_snp_file     = "dummy.xlsx",
      sam_file         = "dummy.bam",
      filter_cpgs_file = "dummy.xlsx",
      sample_type      = "xyz"
    ),
    regexp = "should be one of"
  )
})

test_that("rejects missing cpg_snp_file", {
  expect_error(
    ASM(
      cpg_snp_file     = "nonexistent.xlsx",
      sam_file         = "dummy.bam",
      filter_cpgs_file = "dummy.xlsx",
      sample_type      = "patient"
    )
  )
})
