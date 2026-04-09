ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "mix_beta.R"))
}

test_that("mixbeta creates a betaMix with normalized weights", {
  mix <- mixbeta(c(0.6, 2, 5), c(0.4, 5, 2), param = "ab")
  expect_true(inherits(mix, "betaMix"))
  expect_equal(sum(mix["w", ]), 1, tolerance = 1e-12)
})

test_that("dmix and pmix return valid values", {
  mix <- mixbeta(c(1, 2, 5), param = "ab")
  x <- c(0.2, 0.5)
  d <- dmix(mix, x)
  p <- pmix(mix, x)
  expect_true(all(d >= 0))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("qmix inverts pmix", {
  mix <- mixbeta(c(0.5, 2, 5), c(0.5, 5, 2), param = "ab")
  p <- c(0.1, 0.5, 0.9)
  q <- qmix(mix, p)
  p2 <- pmix(mix, q)
  expect_true(max(abs(p2 - p)) < 1e-6)
})

test_that("rmix samples within [0,1]", {
  mix <- mixbeta(c(0.5, 2, 5), c(0.5, 5, 2), param = "ab")
  samp <- rmix(mix, 200)
  expect_true(all(samp >= 0 & samp <= 1))
})
