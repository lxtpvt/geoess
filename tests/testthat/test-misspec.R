ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "ess_misspec.R"))
}

test_that("geoess_godambe reduces to curvature when J = I", {
  I <- diag(c(2, 3))
  J <- I
  Ipi <- diag(c(1, 1))
  res_g <- geoess_godambe(I, J, Ipi)
  res_c <- geoess_curvature(I, Ipi)
  expect_equal(res_g$ess, res_c$ess, tolerance = 1e-8)
})
