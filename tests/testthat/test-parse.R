# Tests for formula parsing

test_that("parse_fixed handles intercept-only", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- bayesreml:::.parse_fixed(y ~ 1, dat)
  expect_equal(info$response, "y")
  expect_true(info$intercept)
  expect_length(info$terms, 0)
})

test_that("parse_fixed handles no-intercept with factor", {
  dat <- data.frame(y = rnorm(20), env = factor(rep(1:4, 5)))
  info <- bayesreml:::.parse_fixed(y ~ 0 + env, dat)
  expect_false(info$intercept)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "factor")
  expect_equal(info$terms[[1]]$var, "env")
  expect_equal(info$terms[[1]]$n_levels, 4)
})

test_that("parse_fixed handles continuous covariate", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- bayesreml:::.parse_fixed(y ~ x, dat)
  expect_true(info$intercept)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "continuous")
  expect_equal(info$terms[[1]]$var, "x")
})

test_that("parse_fixed handles factor interaction", {
  dat <- data.frame(y = rnorm(20),
                    A = factor(rep(1:2, 10)),
                    B = factor(rep(1:5, 4)))
  info <- bayesreml:::.parse_fixed(y ~ A:B, dat)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "factor_interaction")
})

test_that("parse_fixed handles I() expressions", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- bayesreml:::.parse_fixed(y ~ I(x^2), dat)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "expression")
})

test_that("parse_fixed errors on missing response", {
  dat <- data.frame(x = rnorm(10))
  expect_error(
    bayesreml:::.parse_fixed(y ~ x, dat),
    "not found"
  )
})

test_that("parse_formula handles simple random", {
  dat <- data.frame(geno = factor(rep(1:5, 4)))
  terms <- bayesreml:::.parse_formula(~ geno, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "simple")
  expect_equal(terms[[1]]$var, "geno")
  expect_equal(terms[[1]]$var_n, 5)
})

test_that("parse_formula handles crossed random", {
  dat <- data.frame(geno = factor(rep(1:5, 4)),
                    env  = factor(rep(1:4, each = 5)))
  terms <- bayesreml:::.parse_formula(~ geno + env, dat)
  expect_length(terms, 2)
  expect_equal(terms[[1]]$type, "simple")
  expect_equal(terms[[2]]$type, "simple")
})

test_that("parse_formula handles nested (colon)", {
  dat <- data.frame(block = factor(rep(1:3, 4)),
                    rep   = factor(rep(1:2, 6)))
  terms <- bayesreml:::.parse_formula(~ block:rep, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "nested")
})

test_that("parse_formula handles vm()", {
  dat <- data.frame(geno = factor(1:5))
  terms <- bayesreml:::.parse_formula(~ vm(geno, Gmat), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "vm")
  expect_equal(terms[[1]]$var, "geno")
  expect_equal(terms[[1]]$mat, "Gmat")
})

test_that("parse_formula handles at(env):geno", {
  dat <- data.frame(geno = factor(rep(1:5, 4)),
                    env  = factor(rep(1:4, each = 5)))
  terms <- bayesreml:::.parse_formula(~ at(env):geno, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_simple")
})

test_that("parse_formula handles us(env):id(geno)", {
  dat <- data.frame(geno = factor(rep(1:5, 4)),
                    env  = factor(rep(1:4, each = 5)))
  terms <- bayesreml:::.parse_formula(~ us(env):id(geno), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "us_gxe")
})

test_that("parse_formula handles fa(env,2):id(geno)", {
  dat <- data.frame(geno = factor(rep(1:5, 4)),
                    env  = factor(rep(1:4, each = 5)))
  terms <- bayesreml:::.parse_formula(~ fa(env, 2):id(geno), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "fa_gxe")
  expect_equal(terms[[1]]$k, 2L)
})

test_that("parse_formula errors on fa(env, 0)", {
  dat <- data.frame(geno = factor(rep(1:5, 4)),
                    env  = factor(rep(1:4, each = 5)))
  expect_error(
    bayesreml:::.parse_formula(~ fa(env, 0):id(geno), dat),
    "k >= 1"
  )
})

test_that("parse_formula handles ar1(row):id(col)", {
  dat <- data.frame(row = factor(rep(1:5, 3)),
                    col = factor(rep(1:3, each = 5)))
  terms <- bayesreml:::.parse_formula(~ ar1(row):id(col), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "ar1_spatial")
})

test_that("parse_formula handles spl(x)", {
  dat <- data.frame(x = rnorm(20))
  terms <- bayesreml:::.parse_formula(~ spl(x), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "spline")
})

test_that("parse_formula handles at(env):units (rcov)", {
  dat <- data.frame(env = factor(rep(1:4, each = 5)))
  terms <- bayesreml:::.parse_formula(~ at(env):units, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_units")
})

test_that("parse_formula handles units", {
  terms <- bayesreml:::.parse_formula(~ units, data.frame(x = 1))
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "units")
})

test_that("parse_formula handles dsum(~units|env)", {
  dat <- data.frame(env = factor(rep(1:3, 5)))
  terms <- bayesreml:::.parse_formula(~ dsum(~units | env), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_units")
  expect_equal(terms[[1]]$var, "env")
})

test_that("parse_formula handles three-way combo", {
  dat <- data.frame(a = factor(1:2), b = factor(1:2), c = factor(1:2))
  terms <- bayesreml:::.parse_formula(~ a:b:c, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "combo")
  expect_length(terms[[1]]$vars, 3)
})

test_that("resolve_family works", {
  expect_equal(bayesreml:::.resolve_family("gaussian", NULL)$link, "identity")
  expect_equal(bayesreml:::.resolve_family("binomial", NULL)$link, "logit")
  expect_equal(bayesreml:::.resolve_family("binomial", "probit")$link, "probit")
  expect_equal(bayesreml:::.resolve_family("poisson", NULL)$link, "log")
  expect_error(bayesreml:::.resolve_family("foo", NULL), "Unsupported")
})
