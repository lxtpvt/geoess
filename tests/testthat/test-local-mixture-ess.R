ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "mix_beta.R"))
  source(file.path(pkg_root, "R", "mix_mvnorm_calculus.R"))
  source(file.path(pkg_root, "R", "local_mixture_ess.R"))
}

test_that("geoess_local_mixture applies exact weight scaling in 1D normal mixture", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.7, 0.3),
      mu = list(c(1), c(2)),
      Sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(1, 1, 1))
    ),
    class = "geoess_mix"
  )

  res <- geoess_local_mixture(
    mix = mix,
    theta_star = c(1, 2),
    I1 = matrix(1, 1, 1),
    scale_by_weight = TRUE
  )
  tab <- res$results

  expect_equal(tab$ess_unscaled, c(1, 1), tolerance = 1e-12)
  expect_equal(tab$ess_scaled_weight, c(0.7, 0.3), tolerance = 1e-12)

  # Unified pointwise summary reports attribution and full-mixture ESS at
  # every evaluated reference point.
  expect_equal(res$pointwise$ess_attr_weight, c(1, 1), tolerance = 1e-12)
  expect_equal(res$pointwise$ess_attr_resp, c(1, 1), tolerance = 1e-12)
  expect_equal(res$pointwise$ess_attr_hard_resp, c(0.7, 0.7), tolerance = 1e-12)
})

test_that("geoess_local_mixture computes responsibility scaling", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.7, 0.3),
      mu = list(c(1), c(2)),
      Sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(1, 1, 1))
    ),
    class = "geoess_mix"
  )

  res <- geoess_local_mixture(
    mix = mix,
    theta_star = c(1, 1),
    I1 = matrix(1, 1, 1),
    scale_by_weight = FALSE,
    scale_by_responsibility = TRUE
  )
  tab <- res$results

  expect_true(tab$responsibility[1] > tab$responsibility[2])
  expect_true(tab$responsibility[1] > 0.75)
  expect_true(tab$responsibility[2] < 0.25)
  expect_equal(tab$ess_scaled_resp, tab$responsibility * tab$ess_unscaled, tolerance = 1e-12)
})

test_that("existing mixture ESS default remains unchanged", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.25, 0.75),
      mu = list(c(0, 0), c(0, 0)),
      Sigma = list(diag(2), diag(c(2, 2)))
    ),
    class = "geoess_mix"
  )
  I1 <- diag(c(2, 2))

  old_style <- geoess_ess(I1 = I1, Ipi = mix, method = "trace")
  with_none <- geoess_ess(I1 = I1, Ipi = mix, method = "trace",
                          mixture_attribution = "none")

  expect_equal(old_style$ess_components, c(0.5, 0.25), tolerance = 1e-12)
  expect_equal(old_style$ess_weighted, 0.3125, tolerance = 1e-12)
  expect_equal(with_none$ess, old_style$ess, tolerance = 1e-12)
})

test_that("geoess_ess supports opt-in full-mixture local ESS", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.7, 0.3),
      mu = list(c(1), c(2)),
      Sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(1, 1, 1))
    ),
    class = "geoess_mix"
  )
  I1 <- matrix(1, 1, 1)
  theta_star <- 1.5

  expected <- geoess_trace(I1, geoess_Ipi_from_mix(mix, theta_star), symmetrize = TRUE)
  res_mix <- geoess_ess(
    I1 = I1,
    Ipi = mix,
    method = "trace",
    mixture_curvature = "mix",
    theta_star = theta_star
  )

  expect_equal(res_mix$ess, expected, tolerance = 1e-12)
  expect_equal(res_mix$ess_mixture, expected, tolerance = 1e-12)
  expect_equal(res_mix$mixture_curvature, "mix")
})

test_that("geoess_ess supports opt-in responsibility attribution", {
  mix <- structure(
    list(
      weights = c(0.5, 0.5),
      mu = list(c(0), c(4)),
      Sigma = list(matrix(1, 1, 1), matrix(4, 1, 1))
    ),
    class = "geoess_mix"
  )
  I1 <- matrix(1, 1, 1)

  res_default <- geoess_ess(I1 = I1, Ipi = mix, method = "trace")
  res_resp <- geoess_ess(
    I1 = I1,
    Ipi = mix,
    method = "trace",
    mixture_attribution = "responsibility",
    theta_star = 0
  )

  expect_equal(res_default$ess, 0.625, tolerance = 1e-12)
  expect_gt(res_resp$ess, res_default$ess)
})

test_that("geoess_local_mixture pointwise summary supports unified reporting at each theta", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.7, 0.3),
      mu = list(c(1), c(2)),
      Sigma = list(matrix(1, 1, 1), matrix(4, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(0.25, 1, 1))
    ),
    class = "geoess_mix"
  )

  res <- geoess_local_mixture(
    mix = mix,
    theta_star = 1,
    component = "given",
    component_id = 1,
    I1 = matrix(1, 1, 1),
    scale_by_weight = TRUE,
    scale_by_responsibility = FALSE
  )

  expect_equal(nrow(res$pointwise), 1)
  expect_equal(res$pointwise$ess_attr_weight, 0.775, tolerance = 1e-12)
  expect_equal(res$pointwise$selected_attribution, "weight")
  expect_equal(res$pointwise$ess_attributed, res$pointwise$ess_attr_weight, tolerance = 1e-12)
  expect_true(is.finite(res$pointwise$ess_mixture))

  contrib <- res$contributions
  expect_equal(nrow(contrib), 2)
  expect_equal(contrib$ess_contrib_weight, c(0.7, 0.075), tolerance = 1e-12)
})

test_that("geoess_ess supports hard responsibility attribution", {
  mix <- structure(
    list(
      weights = c(0.5, 0.5),
      mu = list(c(0), c(4)),
      Sigma = list(matrix(1, 1, 1), matrix(4, 1, 1))
    ),
    class = "geoess_mix"
  )
  I1 <- matrix(1, 1, 1)

  res_resp <- geoess_ess(
    I1 = I1,
    Ipi = mix,
    method = "trace",
    mixture_attribution = "responsibility",
    theta_star = 0
  )
  res_hard <- geoess_ess(
    I1 = I1,
    Ipi = mix,
    method = "trace",
    mixture_attribution = "hard_resp",
    theta_star = 0
  )

  expect_equal(res_hard$ess, 0.5, tolerance = 1e-12)
  expect_lt(res_hard$ess, res_resp$ess)
})

test_that("hard responsibility averages weighted ESS across ties", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.5, 0.5),
      mu = list(c(-1), c(1)),
      Sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(1, 1, 1))
    ),
    class = "geoess_mix"
  )

  res <- geoess_local_mixture(
    mix = mix,
    theta_star = 0,
    component = "given",
    component_id = 1,
    I1 = matrix(1, 1, 1),
    scale_by_weight = TRUE,
    scale_by_responsibility = FALSE
  )

  # At theta = 0 the responsibilities tie (0.5, 0.5), so hard attribution
  # uses the mean of weighted component ESS: mean(c(0.5 * 1, 0.5 * 1)) = 0.5.
  expect_equal(res$pointwise$ess_attr_hard_resp, 0.5, tolerance = 1e-12)
})
