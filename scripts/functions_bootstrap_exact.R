# ==============================================================================
# functions_bootstrap_exact.R
#
# Deterministic exact nonparametric-bootstrap distribution for the sample
# median. No external R packages are required.
#
# Numerical note: terms in the even-N branch are never pruned by an
# epsilon-based magnitude check. A term this small can still be multiplied
# by an enormous binomial coefficient (e.g. choose(80, 40)) and contribute
# materially to the total probability mass, so a term is skipped only when
# it evaluates to exactly zero (or non-positive) after floating-point
# arithmetic. This keeps the total probability mass at 1 and matches
# independently derived exact-bootstrap reference values. Floating-point
# support values are never converted to factors or characters, to avoid
# silent precision loss.
#
# Inferential target: an exact empirical nonparametric percentile-bootstrap
# confidence interval for the cohort median.
# ==============================================================================

exact_bootstrap_median_distribution <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]

  n <- length(v)

  if (n == 0L) {
    return(
      data.frame(
        value = NA_real_,
        probability = 1,
        stringsAsFactors = FALSE
      )
    )
  }

  values <- sort(unique(v))

  counts <- vapply(
    values,
    function(z) sum(v == z),
    integer(1)
  )

  probs_emp <- counts / n

  F <- pmin(
    1,
    pmax(
      0,
      cumsum(probs_emp)
    )
  )

  out_value <- numeric()
  out_prob <- numeric()

  add_mass <- function(value, probability) {
    if (!is.finite(value) ||
        !is.finite(probability)) {
      stop(
        "Non-finite value/probability in exact bootstrap distribution."
      )
    }

    # Permit only negligible floating-point negative noise.
    if (probability < -1e-12) {
      stop(
        "Negative exact-bootstrap probability beyond rounding tolerance: ",
        format(probability, digits = 17)
      )
    }

    probability <- max(
      0,
      as.numeric(probability)
    )

    if (probability > 0) {
      out_value <<- c(
        out_value,
        as.numeric(value)
      )
      out_prob <<- c(
        out_prob,
        probability
      )
    }

    invisible(NULL)
  }

  if ((n %% 2L) == 1L) {
    # Odd n = 2r - 1.  Median = X_(r).
    r <- (n + 1L) %/% 2L

    cdf_at_support <- vapply(
      F,
      function(fj) {
        if (fj <= 0) return(0)
        if (fj >= 1) return(1)

        # P[ Binomial(n, fj) >= r ].
        stats::pbinom(
          q = r - 1L,
          size = n,
          prob = fj,
          lower.tail = FALSE
        )
      },
      numeric(1)
    )

    masses <- c(
      cdf_at_support[[1L]],
      diff(cdf_at_support)
    )

    for (ii in seq_along(values)) {
      add_mass(
        values[[ii]],
        masses[[ii]]
      )
    }

  } else {
    # Even n = 2r. R median = (X_(r) + X_(r+1))/2.
    #
    # Exact support consists of:
    #   (i)  X_(r) = X_(r+1) = v_i
    #   (ii) X_(r) = v_i < v_j = X_(r+1)
    #
    # Each mass is computed directly from the empirical category
    # probabilities.  No floating-point support grouping is performed.

    r <- n %/% 2L
    choose_nr <- choose(n, r)
    s <- length(values)

    for (ii in seq_len(s)) {
      q_left <- if (ii == 1L) {
        0
      } else {
        F[[ii - 1L]]
      }

      p_i <- probs_emp[[ii]]

      # ------------------------------------------------------------------------
      # Same central value:
      #   X_(r) = X_(r+1) = v_i
      #
      # Let A = # observations strictly below v_i.
      # Need A <= r-1 and, conditional on A=a, at least r+1-a observations
      # equal v_i.
      # ------------------------------------------------------------------------

      if (q_left >= 1) {
        q_cond_i <- 0
      } else {
        q_cond_i <- p_i / (1 - q_left)
      }

      q_cond_i <- min(
        1,
        max(
          0,
          q_cond_i
        )
      )

      a <- 0:(r - 1L)

      p_a <- stats::dbinom(
        x = a,
        size = n,
        prob = q_left
      )

      # B >= r+1-a  <=>  B > r-a.
      p_b_tail <- stats::pbinom(
        q = r - a,
        size = n - a,
        prob = q_cond_i,
        lower.tail = FALSE
      )

      p_same <- sum(
        p_a * p_b_tail
      )

      add_mass(
        values[[ii]],
        p_same
      )

      # ------------------------------------------------------------------------
      # Distinct central values:
      #   X_(r) = v_i < v_j = X_(r+1)
      #
      # Exactly r bootstrap observations must be <= v_i and exactly r must be
      # >= v_j, with at least one observation equal to each boundary value.
      # Hence no bootstrap observation lies strictly between v_i and v_j.
      # ------------------------------------------------------------------------

      left_term <-
        F[[ii]]^r -
        q_left^r

      # IMPORTANT: do NOT compare left_term with machine epsilon.
      # It is multiplied by choose(n, r), so a value far below epsilon can
      # contribute non-negligible probability mass.
      if (left_term > 0 &&
          ii < s) {

        for (jj in (ii + 1L):s) {
          F_j_minus_1 <- F[[jj - 1L]]
          F_j <- F[[jj]]

          right_term <-
            (1 - F_j_minus_1)^r -
            (1 - F_j)^r

          p_pair <-
            choose_nr *
            left_term *
            right_term

          add_mass(
            (
              values[[ii]] +
                values[[jj]]
            ) / 2,
            p_pair
          )
        }
      }
    }
  }

  if (length(out_value) == 0L) {
    stop(
      "Exact-bootstrap median distribution is empty."
    )
  }

  # Sort numerically and DO NOT aggregate by a floating-point grouping key.
  ord <- order(
    out_value,
    method = "radix"
  )

  dist <- data.frame(
    value = out_value[ord],
    probability = out_prob[ord],
    stringsAsFactors = FALSE
  )

  total_probability <- sum(
    dist$probability
  )

  if (!is.finite(total_probability) ||
      total_probability <= 0) {
    stop(
      "Invalid total probability in exact-bootstrap median distribution."
    )
  }

  if (abs(total_probability - 1) > 1e-10) {
    stop(
      "Exact-bootstrap probability mass does not sum to 1 within tolerance: ",
      format(total_probability, digits = 17)
    )
  }

  # Normalize only floating-point summation noise.
  dist$probability <-
    dist$probability /
    total_probability

  rownames(dist) <- NULL
  dist
}


exact_bootstrap_distribution_quantile <- function(
    distribution,
    probability) {

  if (!is.numeric(probability) ||
      length(probability) != 1L ||
      !is.finite(probability) ||
      probability < 0 ||
      probability > 1) {
    stop(
      "Invalid exact-bootstrap quantile probability."
    )
  }

  if (nrow(distribution) == 0L) {
    return(NA_real_)
  }

  if (anyNA(distribution$value) ||
      anyNA(distribution$probability)) {
    return(NA_real_)
  }

  if (any(diff(distribution$value) < 0)) {
    stop(
      "Exact-bootstrap distribution is not numerically sorted."
    )
  }

  cdf <- cumsum(
    distribution$probability
  )

  # Exact inverse ECDF:
  # inf{x : F(x) >= probability}.
  idx <- which(
    cdf >= probability
  )[[1L]]

  if (is.na(idx)) {
    idx <- nrow(distribution)
  }

  as.numeric(
    distribution$value[[idx]]
  )
}


exact_bootstrap_median_ci <- function(
    v,
    level = 0.95,
    min_n = 1L) {

  v <- as.numeric(v)
  v <- v[is.finite(v)]

  if (length(v) < min_n) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  if (!is.numeric(level) ||
      length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "Invalid confidence level."
    )
  }

  dist <-
    exact_bootstrap_median_distribution(
      v
    )

  alpha <- 1 - level

  c(
    lower =
      exact_bootstrap_distribution_quantile(
        dist,
        alpha / 2
      ),
    upper =
      exact_bootstrap_distribution_quantile(
        dist,
        1 - alpha / 2
      )
  )
}


exact_bootstrap_independent_median_difference_distribution <- function(
    x,
    y) {

  x <- as.numeric(x)
  y <- as.numeric(y)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) == 0L ||
      length(y) == 0L) {
    return(
      data.frame(
        value = NA_real_,
        probability = 1,
        stringsAsFactors = FALSE
      )
    )
  }

  dx <-
    exact_bootstrap_median_distribution(
      x
    )

  dy <-
    exact_bootstrap_median_distribution(
      y
    )

  n_pairs <-
    nrow(dx) *
    nrow(dy)

  if (n_pairs > 5e6) {
    stop(
      "Exact independent median-difference distribution would require ",
      n_pairs,
      " support pairs; exceeds safety limit."
    )
  }

  values <- as.vector(
    outer(
      dx$value,
      dy$value,
      FUN = "-"
    )
  )

  probs <- as.vector(
    outer(
      dx$probability,
      dy$probability,
      FUN = "*"
    )
  )

  ord <- order(
    values,
    method = "radix"
  )

  dist <- data.frame(
    value = values[ord],
    probability = probs[ord],
    stringsAsFactors = FALSE
  )

  total_probability <- sum(
    dist$probability
  )

  if (!is.finite(total_probability) ||
      total_probability <= 0) {
    stop(
      "Invalid total probability in exact median-difference distribution."
    )
  }

  if (abs(total_probability - 1) > 1e-10) {
    stop(
      "Exact median-difference probability mass does not sum to 1 ",
      "within tolerance: ",
      format(total_probability, digits = 17)
    )
  }

  dist$probability <-
    dist$probability /
    total_probability

  rownames(dist) <- NULL
  dist
}


exact_bootstrap_independent_median_difference_ci <- function(
    x,
    y,
    level = 0.95) {

  x <- as.numeric(x)
  y <- as.numeric(y)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) == 0L ||
      length(y) == 0L) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  dist <-
    exact_bootstrap_independent_median_difference_distribution(
      x,
      y
    )

  alpha <- 1 - level

  c(
    lower =
      exact_bootstrap_distribution_quantile(
        dist,
        alpha / 2
      ),
    upper =
      exact_bootstrap_distribution_quantile(
        dist,
        1 - alpha / 2
      )
  )
}


# ==============================================================================
# Independent mathematical self-tests
# ==============================================================================

bootstrap_mass_at <- function(
    dist,
    target,
    tolerance = 1e-14) {

  sum(
    dist$probability[
      abs(
        dist$value -
          target
      ) <= tolerance
    ]
  )
}


validate_exact_bootstrap_helpers <- function(
    tolerance = 1e-12) {

  # ---------------------------------------------------------------------------
  # Odd n=3, distinct values.
  # Full exhaustive bootstrap has 3^3 = 27 resamples:
  # median probabilities 7/27, 13/27, 7/27.
  # ---------------------------------------------------------------------------

  d1 <-
    exact_bootstrap_median_distribution(
      c(0, 1, 2)
    )

  expected_value_1 <- c(
    0, 1, 2
  )

  expected_prob_1 <- c(
    7, 13, 7
  ) / 27

  for (ii in seq_along(expected_value_1)) {
    observed <-
      bootstrap_mass_at(
        d1,
        expected_value_1[[ii]]
      )

    if (abs(
      observed -
        expected_prob_1[[ii]]
    ) > tolerance) {
      stop(
        "Exact-bootstrap odd-n self-test FAILED."
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Even n=4, distinct values.
  # Full exhaustive bootstrap has 4^4 = 256 resamples.
  # ---------------------------------------------------------------------------

  d2 <-
    exact_bootstrap_median_distribution(
      c(0, 1, 2, 3)
    )

  expected_value_2 <- c(
    0, 0.5, 1, 1.5, 2, 2.5, 3
  )

  expected_prob_2 <- c(
    13, 30, 55, 60, 55, 30, 13
  ) / 256

  for (ii in seq_along(expected_value_2)) {
    observed <-
      bootstrap_mass_at(
        d2,
        expected_value_2[[ii]]
      )

    if (abs(
      observed -
        expected_prob_2[[ii]]
    ) > tolerance) {
      stop(
        "Exact-bootstrap even-n self-test FAILED."
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Even n=4 with ties.
  # ---------------------------------------------------------------------------

  d3 <-
    exact_bootstrap_median_distribution(
      c(0, 0, 1, 2)
    )

  expected_value_3 <- c(
    0, 0.5, 1, 1.5, 2
  )

  expected_prob_3 <- c(
    80, 72, 61, 30, 13
  ) / 256

  for (ii in seq_along(expected_value_3)) {
    observed <-
      bootstrap_mass_at(
        d3,
        expected_value_3[[ii]]
      )

    if (abs(
      observed -
        expected_prob_3[[ii]]
    ) > tolerance) {
      stop(
        "Exact-bootstrap tied even-n self-test FAILED."
      )
    }
  }

  # Odd n=3 with ties.
  d4 <-
    exact_bootstrap_median_distribution(
      c(0, 0, 1)
    )

  if (abs(
    bootstrap_mass_at(
      d4,
      0
    ) -
      20 / 27
  ) > tolerance ||
      abs(
        bootstrap_mass_at(
          d4,
          1
        ) -
          7 / 27
      ) > tolerance) {
    stop(
      "Exact-bootstrap tied odd-n self-test FAILED."
    )
  }

  # ---------------------------------------------------------------------------
  # Large even-N tail-mass regression: for N=80, choose(80,40) is huge, so
  # terms below machine epsilon must not be pruned before multiplication
  # (see the numerical note in the file header).
  # ---------------------------------------------------------------------------

  d_large <-
    exact_bootstrap_median_distribution(
      seq_len(80)
    )

  large_mass <- sum(
    d_large$probability
  )

  if (abs(
    large_mass - 1
  ) > tolerance) {
    stop(
      "Exact-bootstrap large-even-N probability-mass self-test FAILED: ",
      format(large_mass, digits = 17)
    )
  }

  # Constant-vector identity.
  ci_const <-
    exact_bootstrap_median_ci(
      rep(
        0.25,
        7
      )
    )

  if (max(
    abs(
      ci_const -
        0.25
    )
  ) > tolerance) {
    stop(
      "Exact-bootstrap constant-vector self-test FAILED."
    )
  }

  # Independent difference of constants.
  ci_diff <-
    exact_bootstrap_independent_median_difference_ci(
      rep(0.75, 5),
      rep(0.25, 7)
    )

  if (max(
    abs(
      ci_diff -
        0.50
    )
  ) > tolerance) {
    stop(
      "Exact-bootstrap median-difference self-test FAILED."
    )
  }

  TRUE
}


EXACT_BOOTSTRAP_HELPERS_VALIDATED <-
  validate_exact_bootstrap_helpers()
