# NW2024 equations (13)--(16), with the existing output-gap data specification.
# Run: Rscript NKPC_1.R [directory containing china_cpi.csv and china_gdp.csv]

nkpc_quarter_adjust <- function(x, quarters) {
  # The intercept and quarterly dummies remove the mean and fixed seasonality.
  design <- model.matrix(~ factor(quarters))
  as.numeric(lm.fit(design, as.numeric(x))$residuals)
}

nkpc_cycle <- function(x, quarters) {
  cycle <- mFilter::hpfilter(log(x), freq = 1600)$cycle
  nkpc_quarter_adjust(cycle, quarters)
}

read_nkpc_inputs <- function(data_dir = ".") {
  # Construct the price index from unadjusted monthly CPI ratios.
  cpi <- read.csv(file.path(data_dir, "china_cpi.csv"), check.names = FALSE)
  gdp <- read.csv(file.path(data_dir, "china_gdp.csv"), check.names = FALSE)
  regions <- as.character(cpi$Region)
  months <- as.matrix(cpi[, -1L])
  n_quarters <- ncol(months) %/% 3L
  stopifnot(ncol(months) == 144L, !anyDuplicated(regions),
            setequal(regions, gdp$Region))
  quarterly_cpi <- matrix(NA_real_, nrow(months), n_quarters)
  for (j in seq_len(n_quarters)) {
    quarterly_cpi[, j] <- apply(months[, 3L * (j - 1L) + 1:3, drop = FALSE], 1L,
      function(x) {
        value <- 1
        for (month in x) value <- value * month / 100
        value * 100
      })
  }
  price <- matrix(1, nrow(months), n_quarters)
  for (j in seq.int(n_quarters - 1L, 1L)) {
    price[, j] <- price[, j + 1L] * (quarterly_cpi[, j] / 100)
  }
  # Convert year-to-date nominal GDP to quarterly flows before CPI deflation.
  # Columns descend Q4, Q3, Q2, Q1 within each year; Q1 is already a flow.
  nominal <- as.matrix(gdp[match(regions, gdp$Region), 1L + seq_len(n_quarters)])
  non_q1 <- which(rep(4:1, length.out = n_quarters) != 1L)
  nominal[, non_q1] <- nominal[, non_q1, drop = FALSE] -
    nominal[, non_q1 + 1L, drop = FALSE]
  real_gdp <- t(nominal / price)[rev(seq_len(n_quarters)), , drop = FALSE]
  quarters <- paste0(rep(2012:2023, each = 4L), "Q", rep(1:4, 12L))
  quarter_of_year <- substr(quarters, 6L, 6L)
  gap <- vapply(seq_along(regions), function(i) {
    nkpc_cycle(real_gdp[, i], quarter_of_year)
  }, numeric(n_quarters))
  inflation <- t(quarterly_cpi - 100)[rev(seq_len(n_quarters)), , drop = FALSE]
  inflation <- apply(inflation, 2L, nkpc_quarter_adjust, quarters = quarter_of_year)
  dimnames(inflation) <- list(quarters, regions)
  dimnames(gap) <- list(quarters, regions)
  list(y = inflation, s = gap, quarters = quarters, regions = regions)
}

nkpc_factors <- function(moment, r_max = 8L, relative_tol = 1e-10, count = NULL) {
  stopifnot(r_max >= 1L, r_max == as.integer(r_max),
            relative_tol > 0, relative_tol < 1)
  eig <- eigen((moment + t(moment)) / 2, symmetric = TRUE)
  values <- eig$values                       # descending, as returned by R
  negative_allowance <- 64 * nrow(moment) * .Machine$double.eps * max(abs(values))
  cutoff <- max(relative_tol * max(values[1L], 0), negative_allowance)
  if (any(values < -negative_allowance)) {
    stop("Factor moment matrix is not positive semidefinite")
  }
  values[values < 0] <- 0                    # roundoff-sized negative values only
  candidates <- seq_len(min(r_max, length(values) - 1L))
  ratios <- rep(NA_real_, length(candidates))
  eligible <- values[candidates + 1L] > cutoff
  ratios[eligible] <- values[candidates[eligible]] / values[candidates[eligible] + 1L]
  if (is.null(count)) {
    if (!any(eligible)) stop("No admissible eigenvalue ratio; specify a factor count")
    # ER searches positive factor counts; exact ties select the smallest count.
    count <- candidates[which.max(ratios)]
  }
  stopifnot(length(count) == 1L, is.finite(count), count == as.integer(count),
            count >= 0L, count <= sum(values > cutoff))
  factors <- sqrt(nrow(moment)) * eig$vectors[, seq_len(count), drop = FALSE]
  list(count = as.integer(count), factors = factors, eigenvalues = eig$values,
       cutoff = cutoff, negative_allowance = negative_allowance,
       numerical_rank = sum(values > cutoff),
       criterion = data.frame(factors = candidates, ratio = ratios))
}

nkpc_projection <- function(factors) {
  if (!ncol(factors)) return(diag(nrow(factors)))
  decomposition <- qr(factors)
  if (decomposition$rank != ncol(factors)) stop("Factor window has lost column rank")
  basis <- qr.Q(decomposition)
  diag(nrow(factors)) - tcrossprod(basis)
}

nkpc_iv <- function(y, X, Z, projection) {
  nobs <- sum(vapply(y, length, integer(1L)))
  A <- Reduce(`+`, Map(function(z, x) crossprod(z, projection %*% x), Z, X)) / nobs
  B <- Reduce(`+`, lapply(Z, function(z) crossprod(z, projection %*% z))) / nobs
  g <- Reduce(`+`, Map(function(z, yi) crossprod(z, projection %*% yi), Z, y)) / nobs
  bread <- crossprod(A, solve(B, A))
  theta <- solve(bread, crossprod(A, solve(B, g)))
  rownames(theta) <- colnames(X[[1L]])
  list(theta = theta, A = A, B = B, g = g, nobs = nobs,
       bread = bread, condition = kappa(bread, exact = TRUE))
}

nkpc_covariance <- function(step, Z, residuals1, projection) {
  periods <- nrow(projection)
  sigma2 <- vapply(residuals1, function(u) {
    as.numeric(crossprod(u, projection %*% u)) / periods
  }, numeric(1L))
  Omega <- Reduce(`+`, Map(function(z, variance) {
    variance * crossprod(z, projection %*% z)
  }, Z, sigma2)) / step$nobs
  H <- solve(step$bread, t(step$A) %*% solve(step$B))
  Sigma <- H %*% Omega %*% t(H)
  vcov <- (Sigma + t(Sigma)) / (2 * step$nobs)
  list(vcov = vcov, Sigma = Sigma, Omega = Omega, sigma2 = sigma2)
}

Jhat <- function(residuals, Z, projection, Omega) {
  # Equation (17) quadratic form. Calibration for this B-weighted fit is separate.
  nobs <- sum(vapply(residuals, length, integer(1L)))
  score <- Reduce(`+`, Map(function(z, u) crossprod(z, projection %*% u), Z, residuals))
  as.numeric(crossprod(score, solve(Omega, score)) / nobs)
}

estimate_nkpc <- function(input, max_lag = 5L, r_max_g = 8L, r_max_f = 8L,
                          relative_tol = 1e-10, m_g = NULL, m_f = NULL) {
  N <- ncol(input$s)
  T_raw <- nrow(input$s)
  stopifnot(max_lag >= 1L, max_lag == as.integer(max_lag), T_raw > max_lag + 1L)
  target <- seq.int(max_lag + 1L, T_raw - 1L)
  T_eff <- length(target)
  offsets <- c(0L, seq_len(max_lag), -1L)
  G <- nkpc_factors(tcrossprod(input$s) / (N * T_raw),
                    r_max_g, relative_tol, m_g)
  # Equation (13): shift the variable AND factors, then project on that window.
  windows <- lapply(offsets, function(lag) target - lag)
  M_G <- lapply(windows, function(index) nkpc_projection(G$factors[index, , drop = FALSE]))
  Z <- lapply(seq_len(N), function(i) {
    z <- vapply(seq_along(windows), function(j) {
      as.vector(M_G[[j]] %*% input$s[windows[[j]], i])
    }, numeric(T_eff))
    colnames(z) <- c("current", paste0("lag", seq_len(max_lag)), "lead1")
    z
  })
  y <- lapply(seq_len(N), function(i) input$y[target, i])
  X <- lapply(seq_len(N), function(i) {
    # Equation (8) retains the forcing variable before the G projection.
    cbind(pi_lag = input$y[target - 1L, i],
          pi_lead = input$y[target + 1L, i], y_gap = input$s[target, i])
  })
  first <- nkpc_iv(y, X, Z, diag(T_eff))
  residuals1 <- Map(function(yi, x) as.vector(yi - x %*% first$theta), y, X)
  U1 <- do.call(cbind, residuals1)
  F <- nkpc_factors(tcrossprod(U1) / (N * T_eff), r_max_f, relative_tol, m_f)
  M_F <- nkpc_projection(F$factors)
  second <- nkpc_iv(y, X, Z, M_F)
  covariance <- nkpc_covariance(second, Z, residuals1, M_F)
  coefficient <- as.numeric(second$theta)
  se <- sqrt(diag(covariance$vcov))
  statistic <- coefficient / se
  coefficients <- data.frame(term = rownames(second$theta), estimate = coefficient,
    std.error = se, statistic = statistic, p.value = 2 * pnorm(-abs(statistic)),
    conf.low = coefficient - qnorm(0.975) * se,
    conf.high = coefficient + qnorm(0.975) * se, row.names = NULL)
  residuals2 <- Map(function(yi, x) as.vector(yi - x %*% second$theta), y, X)
  list(coefficients = coefficients, first_step = first, second_step = second,
       covariance = covariance, factors = list(G = G, F = F),
       sample = list(regions = input$regions, quarters = input$quarters[target],
                     N = N, T_raw = T_raw, T_eff = T_eff, nobs = N * T_eff),
       inputs = list(y = y, X = X, Z = Z, M_G = M_G, M_F = M_F,
                     residuals1 = residuals1, residuals2 = residuals2),
       overidentification = list(statistic = Jhat(residuals2, Z, M_F, covariance$Omega),
                                  df = ncol(Z[[1L]]) - ncol(X[[1L]]), p.value = NA_real_))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  nkpc_result <- estimate_nkpc(read_nkpc_inputs(if (length(args)) args[[1L]] else "."))
  print(nkpc_result$coefficients, row.names = FALSE, digits = 6)
  cat("Factors: G =", nkpc_result$factors$G$count, "; F =", nkpc_result$factors$F$count,
      "\nObservations:", nkpc_result$sample$nobs, "\n")
}
