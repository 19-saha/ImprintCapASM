library(ImprintCapASM)
test_that("rejects invalid sample_type", {
  expect_error(
    prepare_cpg_snp_input(
      snp_file     = "dummy.vcf",
      meth_file    = "dummy.txt",
      cpg_ref_file = "dummy.xlsx",
      sample_type  = "banana"
    ),
    regexp = "should be one of"
  )
})

test_that("rejects missing snp_file", {
  expect_error(
    prepare_cpg_snp_input(
      snp_file     = "nonexistent.vcf",
      meth_file    = "dummy.txt",
      cpg_ref_file = "dummy.xlsx",
      sample_type  = "control"
        ),
    regexp = "not found|does not exist|cannot open"
  )
})
