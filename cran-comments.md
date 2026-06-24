## Resubmission (v0.1.1)

Changes since initial submission:

* Added `run_pipeline()` — the main convenience wrapper executing the full
  three-step pipeline (prepare → extract → ASM) across all samples in a
  directory.
* Fixed non-ASCII characters in `R/run_pipeline.R`.
* Fixed incorrect argument names in `ASM()` call within `run_pipeline()`.
* Added `Sys.which("samtools")` availability guard in `extract_bam_regions()`
  before calling `system()`, with a clear error message when samtools is absent.
* Switched `\donttest{}` to `\dontrun{}` in `ASM()` and `run_pipeline()`
  examples — these require BAM files and samtools and are not suitable for
  automated example checking.
* Fixed `DESCRIPTION` field to contain complete sentences ending with a full
  stop.
* Added `SystemRequirements: samtools (>= 1.10)` to `DESCRIPTION`.
* Updated `test-extract_bam_regions.R` regexp to match the new samtools
  error message.


## R CMD check results

0 errors | 0 warnings | 4 notes
