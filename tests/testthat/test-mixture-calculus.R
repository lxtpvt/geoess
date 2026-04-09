ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "mix_mvnorm_calculus.R"))
}

test_that("geoess_Ipi_from_mix matches precision for single MVN", {
  Sigma <- diag(c(2, 3))
  mix <- structure(
    list(
      family = "mvnorm",
      weights = 1,
      mu = list(c(0, 0)),
      Sigma = list(Sigma),
      prec = list(diag(c(1 / 2, 1 / 3)))
    ),
    class = "geoess_mix"
  )
  Ipi <- geoess_Ipi_from_mix(mix, c(0, 0))
  expect_equal(Ipi, diag(c(1 / 2, 1 / 3)), tolerance = 1e-8)
})
