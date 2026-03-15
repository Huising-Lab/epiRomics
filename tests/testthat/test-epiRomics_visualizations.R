library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

histone_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "histone"]
enhancers <- epiRomics::epiRomics_enhancers(epiRomics_dB,
  histone_marks[1], histone_marks[2])
enhanceosome <- epiRomics::epiRomics_enhanceosome(enhancers, epiRomics_dB)

# --- call_accessible_regions ---
# (plot_super_enhancers, plot_signal_enrichment, plot_tf_chromatin_state,
#  plot_tf_chromatin_celltype are not in stable — see experimental fork)

test_that("call_accessible_regions returns logical vector", {
  skip_if_no_extdata()
  bw_file <- system.file("extdata", "BigWigs",
    "Kaestner_Beta.atac.signal.subset.chr1.bigwig", package = "epiRomics")
  skip_if(!file.exists(bw_file), "BigWig file not available")

  test_regions <- enhanceosome@annotations[1:min(50, length(enhanceosome@annotations))]
  result <- call_accessible_regions(bw_file, test_regions, z_threshold = 2.0)

  expect_true(is.logical(result))
  expect_equal(length(result), length(test_regions))
})

test_that("call_accessible_regions with strict threshold returns fewer TRUE", {
  skip_if_no_extdata()
  bw_file <- system.file("extdata", "BigWigs",
    "Kaestner_Beta.atac.signal.subset.chr1.bigwig", package = "epiRomics")
  skip_if(!file.exists(bw_file), "BigWig file not available")

  test_regions <- enhanceosome@annotations[1:min(50, length(enhanceosome@annotations))]
  result_loose <- call_accessible_regions(bw_file, test_regions, z_threshold = 1.0)
  result_strict <- call_accessible_regions(bw_file, test_regions, z_threshold = 3.0)

  expect_true(sum(result_strict) <= sum(result_loose))
})

test_that("call_accessible_regions errors on invalid inputs", {
  skip_if_no_extdata()
  expect_error(call_accessible_regions("nonexistent.bw",
    enhanceosome@annotations[1:5]))
  bw_file <- system.file("extdata", "BigWigs",
    "Kaestner_Beta.atac.signal.subset.chr1.bigwig", package = "epiRomics")
  skip_if(!file.exists(bw_file), "BigWig file not available")
  expect_error(call_accessible_regions(bw_file, "not_a_granges"))
})
