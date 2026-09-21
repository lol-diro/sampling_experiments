# ==============================================================================
# functions_sampling.R
#
# Shared, deterministic functions for exhaustive sector downsampling.
#
# Design principles
#   - Patient is the inferential unit.
#   - All subsets are enumerated exactly; no Monte Carlo is required here.
#   - "Unconstrained/random" means uniformly sampling among all k-sector subsets.
#   - Spatial strategy functions use sector ORDER/RANK only; they do not invent
#     physical distances.
#   - Reference terminology is ubiquitous / non-ubiquitous / private.
# ==============================================================================

# ------------------------------ basic checks ---------------------------------

.assert_integer_scalar <- function(x, name, min_value = 1L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      abs(x - round(x)) > .Machine$double.eps^0.5 ||
      x < min_value) {
    stop(name, " must be one integer >= ", min_value, ".")
  }
  as.integer(round(x))
}

# ---------------------------- subset generators ------------------------------

enumerate_all_subsets <- function(n, k) {
  n <- .assert_integer_scalar(n, "n")
  k <- .assert_integer_scalar(k, "k")
  if (k > n) stop("k cannot exceed n.")
  utils::combn(seq_len(n), k)
}

# Grid-based dispersed downsampling on ORDERED sectors.
#
# k = 1:
#   choose the central rank(s); for even n there are two equivalent centres.
#
# k >= 2:
#   require coverage of both ends (rank 1 and rank n), then choose the subset(s)
#   minimizing squared deviation from k equally spaced ideal rank positions.
#
# Multiple exactly equivalent solutions are retained, never broken by sample ID.
grid_subsets_rank <- function(n, k, tolerance = 1e-12) {
  n <- .assert_integer_scalar(n, "n")
  k <- .assert_integer_scalar(k, "k")
  if (k > n) stop("k cannot exceed n.")

  if (k == 1L) {
    target <- (n + 1) / 2
    ranks <- seq_len(n)
    loss <- abs(ranks - target)
    keep <- which(abs(loss - min(loss)) <= tolerance)
    return(matrix(as.integer(keep), nrow = 1L))
  }

  cmb <- enumerate_all_subsets(n, k)
  keep_endpoints <- cmb[1L, ] == 1L & cmb[k, ] == n
  cmb <- cmb[, keep_endpoints, drop = FALSE]

  if (ncol(cmb) == 0L) {
    stop("Internal grid error: no subset covers both endpoints.")
  }

  targets <- seq(1, n, length.out = k)
  loss <- vapply(seq_len(ncol(cmb)), function(j) {
    sum((as.numeric(cmb[, j]) - targets)^2)
  }, numeric(1))

  keep <- which(abs(loss - min(loss)) <= tolerance)
  cmb[, keep, drop = FALSE]
}

# Clustered/adjacent sampling on ORDERED sectors.
# Every possible contiguous k-sector window is retained.
clustered_subsets_rank <- function(n, k) {
  n <- .assert_integer_scalar(n, "n")
  k <- .assert_integer_scalar(k, "k")
  if (k > n) stop("k cannot exceed n.")

  starts <- seq_len(n - k + 1L)
  out <- vapply(starts, function(s) seq.int(s, length.out = k), integer(k))

  # vapply simplifies k=1 to a vector; restore the common k x n_windows shape.
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = k)
  }
  storage.mode(out) <- "integer"
  out
}

# --------------------------- incidence construction ---------------------------

# Build one logical event x sector matrix for one patient/event set.
#
# patient_reference must contain one row per patient-specific event and:
#   event_key, is_ubiquitous, is_nonubiquitous, is_private
#
# patient_presence may contain all event presences; rows not present in
# patient_reference are ignored. Thus the same function can be used later for
# the protein-altering sensitivity simply by passing a filtered reference.
build_patient_incidence <- function(patient_presence,
                                    patient_reference,
                                    sample_ids) {
  required_presence <- c("sample_id", "event_key")
  required_reference <- c(
    "event_key", "is_ubiquitous", "is_nonubiquitous", "is_private"
  )

  if (!all(required_presence %in% names(patient_presence))) {
    stop("patient_presence lacks: ",
         paste(setdiff(required_presence, names(patient_presence)), collapse = ", "))
  }
  if (!all(required_reference %in% names(patient_reference))) {
    stop("patient_reference lacks: ",
         paste(setdiff(required_reference, names(patient_reference)), collapse = ", "))
  }

  sample_ids <- as.character(sample_ids)
  if (length(sample_ids) < 1L || anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids) > 0L) {
    stop("sample_ids must be unique, non-missing, non-empty strings.")
  }

  ref <- as.data.frame(patient_reference, stringsAsFactors = FALSE)
  if (anyNA(ref$event_key) || anyDuplicated(ref$event_key) > 0L) {
    stop("patient_reference event_key must be unique and non-missing.")
  }

  event_keys <- ref$event_key
  n_events <- length(event_keys)
  n_samples <- length(sample_ids)

  mat <- matrix(FALSE, nrow = n_events, ncol = n_samples)
  colnames(mat) <- sample_ids

  if (n_events > 0L) {
    pp <- as.data.frame(patient_presence, stringsAsFactors = FALSE)
    pp <- pp[pp$event_key %in% event_keys, required_presence, drop = FALSE]

    if (nrow(pp) > 0L) {
      if (anyDuplicated(pp[, c("sample_id", "event_key"), drop = FALSE]) > 0L) {
        stop("Duplicate sample_id x event_key presence detected.")
      }

      rr <- match(pp$event_key, event_keys)
      cc <- match(as.character(pp$sample_id), sample_ids)

      if (anyNA(rr)) stop("Internal event-key matching failure.")
      if (anyNA(cc)) {
        stop("At least one event presence refers to a sample outside sample_ids.")
      }
      mat[cbind(rr, cc)] <- TRUE
    }

    # Every reference event must be observed in at least one available sector.
    observed_counts <- rowSums(mat)
    if (any(observed_counts < 1L)) {
      stop("At least one reference event has no presence in the incidence matrix.")
    }

    # The reference classes must agree with the reconstructed incidence matrix.
    ref_ubiq <- observed_counts == n_samples
    ref_nonubiq <- observed_counts > 0L & observed_counts < n_samples
    ref_private <- observed_counts == 1L

    if (!identical(as.logical(ref$is_ubiquitous), as.logical(ref_ubiq)) ||
        !identical(as.logical(ref$is_nonubiquitous), as.logical(ref_nonubiq)) ||
        !identical(as.logical(ref$is_private), as.logical(ref_private))) {
      stop("Reference class labels disagree with reconstructed event presences.")
    }
  }

  pairwise_jaccard <- pairwise_jaccard_distance_matrix(mat)
  full_jaccard <- mean_jaccard_from_distance_matrix(
    pairwise_jaccard, seq_len(n_samples)
  )

  list(
    matrix = mat,
    reference = ref,
    event_keys = event_keys,
    sample_ids = sample_ids,
    n_events = n_events,
    n_samples = n_samples,
    pairwise_jaccard = pairwise_jaccard,
    full_jaccard_ith = full_jaccard
  )
}

# ------------------------------ Jaccard ITH ----------------------------------

# Pairwise sector Jaccard-distance matrix:
#   1 - |A intersect B| / |A union B|
#
# If both sector event sets are empty, distance is defined as 0, preserving the
# legacy convention. This situation does not occur in the validated HCC data.
pairwise_jaccard_distance_matrix <- function(binary_matrix) {
  if (is.null(dim(binary_matrix))) {
    stop("binary_matrix must be a matrix.")
  }

  n <- ncol(binary_matrix)
  d <- matrix(0, nrow = n, ncol = n)

  if (n < 2L) return(d)

  for (i in seq_len(n - 1L)) {
    a <- binary_matrix[, i]
    for (j in seq.int(i + 1L, n)) {
      b <- binary_matrix[, j]
      union_n <- sum(a | b)
      dij <- if (union_n == 0L) {
        0
      } else {
        1 - sum(a & b) / union_n
      }
      d[i, j] <- dij
      d[j, i] <- dij
    }
  }

  colnames(d) <- colnames(binary_matrix)
  rownames(d) <- colnames(binary_matrix)
  d
}

mean_jaccard_from_distance_matrix <- function(distance_matrix, subset_idx) {
  subset_idx <- as.integer(subset_idx)
  if (length(subset_idx) < 2L) return(NA_real_)

  z <- distance_matrix[subset_idx, subset_idx, drop = FALSE]
  mean(z[upper.tri(z)])
}

# ----------------------------- subset metrics --------------------------------

compute_subset_metrics <- function(patient_data, subset_idx) {
  mat <- patient_data$matrix
  ref <- patient_data$reference

  subset_idx <- as.integer(subset_idx)
  if (length(subset_idx) < 1L ||
      anyNA(subset_idx) ||
      any(subset_idx < 1L) ||
      any(subset_idx > ncol(mat)) ||
      anyDuplicated(subset_idx) > 0L) {
    stop("subset_idx must contain unique valid sector-column indices.")
  }

  k <- length(subset_idx)
  sub_counts <- if (nrow(mat) == 0L) {
    integer(0)
  } else if (k == 1L) {
    as.integer(mat[, subset_idx])
  } else {
    rowSums(mat[, subset_idx, drop = FALSE])
  }

  nonubiq <- as.logical(ref$is_nonubiquitous)
  private <- as.logical(ref$is_private)

  n_events_reference <- length(sub_counts)
  n_nonubiq_reference <- sum(nonubiq)
  n_private_reference <- sum(private)

  # 1) Primary endpoint: detection of full-reference non-ubiquitous events.
  n_nonubiq_detected <- if (n_nonubiq_reference == 0L) {
    0L
  } else {
    sum(sub_counts[nonubiq] > 0L)
  }

  recall_nonubiquitous_detection <- if (n_nonubiq_reference == 0L) {
    NA_real_
  } else {
    n_nonubiq_detected / n_nonubiq_reference
  }

  # 2) Key secondary: event is both detected and remains non-ubiquitous within
  #    the selected sectors.
  n_nonubiq_correctly_classified <- if (n_nonubiq_reference == 0L) {
    0L
  } else {
    sum(sub_counts[nonubiq] > 0L & sub_counts[nonubiq] < k)
  }

  recall_heterogeneity_classification <- if (n_nonubiq_reference == 0L) {
    NA_real_
  } else {
    n_nonubiq_correctly_classified / n_nonubiq_reference
  }

  # Full-reference non-ubiquitous events that are detected in every selected
  # sector and therefore appear ubiquitous under the reduced sampling design.
  n_apparent_ubiquity <- if (n_nonubiq_reference == 0L) {
    0L
  } else {
    sum(sub_counts[nonubiq] == k)
  }

  apparent_ubiquity_error <- if (n_nonubiq_reference == 0L) {
    NA_real_
  } else {
    n_apparent_ubiquity / n_nonubiq_reference
  }

  conditional_apparent_ubiquity_rate <- if (n_nonubiq_detected == 0L) {
    NA_real_
  } else {
    n_apparent_ubiquity / n_nonubiq_detected
  }

  # Identity required by the three-state decomposition:
  # detected = correctly heterogeneous + apparently ubiquitous.
  if (n_nonubiq_reference > 0L) {
    lhs <- recall_nonubiquitous_detection
    rhs <- recall_heterogeneity_classification + apparent_ubiquity_error
    if (abs(lhs - rhs) > 1e-12) {
      stop("Detection/classification decomposition identity failed.")
    }
  }

  # 3) Full-reference private-event detection.
  recall_private <- if (n_private_reference == 0L) {
    NA_real_
  } else {
    sum(sub_counts[private] > 0L) / n_private_reference
  }

  # 4) Non-ubiquitous fraction among events observed in the subset.
  n_observed_subset <- sum(sub_counts > 0L)
  n_nonubiquitous_subset <- sum(sub_counts > 0L & sub_counts < k)

  subset_nonubiquitous_fraction <- if (n_observed_subset == 0L) {
    NA_real_
  } else {
    n_nonubiquitous_subset / n_observed_subset
  }

  full_nonubiquitous_fraction <- if (n_events_reference == 0L) {
    NA_real_
  } else {
    n_nonubiq_reference / n_events_reference
  }

  nonubiquitous_fraction_abs_error <-
    if (is.na(subset_nonubiquitous_fraction) ||
        is.na(full_nonubiquitous_fraction)) {
      NA_real_
    } else {
      abs(subset_nonubiquitous_fraction - full_nonubiquitous_fraction)
    }

  # 5) Mean pairwise Jaccard ITH and its absolute error.
  subset_jaccard_ith <- mean_jaccard_from_distance_matrix(
    patient_data$pairwise_jaccard, subset_idx
  )
  full_jaccard_ith <- patient_data$full_jaccard_ith

  jaccard_ith_abs_error <-
    if (is.na(subset_jaccard_ith) || is.na(full_jaccard_ith)) {
      NA_real_
    } else {
      abs(subset_jaccard_ith - full_jaccard_ith)
    }

  data.frame(
    k = k,
    n_events_reference = n_events_reference,
    n_nonubiquitous_reference = n_nonubiq_reference,
    n_private_reference = n_private_reference,

    n_nonubiquitous_detected = n_nonubiq_detected,
    n_nonubiquitous_correctly_classified = n_nonubiq_correctly_classified,
    n_apparent_ubiquity = n_apparent_ubiquity,

    recall_nonubiquitous_detection = recall_nonubiquitous_detection,
    recall_heterogeneity_classification = recall_heterogeneity_classification,
    apparent_ubiquity_error = apparent_ubiquity_error,
    conditional_apparent_ubiquity_rate = conditional_apparent_ubiquity_rate,
    recall_private = recall_private,

    full_nonubiquitous_fraction = full_nonubiquitous_fraction,
    subset_nonubiquitous_fraction = subset_nonubiquitous_fraction,
    nonubiquitous_fraction_abs_error = nonubiquitous_fraction_abs_error,

    full_jaccard_ith = full_jaccard_ith,
    subset_jaccard_ith = subset_jaccard_ith,
    jaccard_ith_abs_error = jaccard_ith_abs_error,

    stringsAsFactors = FALSE
  )
}

compute_all_subset_metrics <- function(patient_data, k) {
  k <- .assert_integer_scalar(k, "k")
  n <- patient_data$n_samples
  if (k > n) stop("k cannot exceed patient's number of sectors.")

  cmb <- enumerate_all_subsets(n, k)

  out <- lapply(seq_len(ncol(cmb)), function(j) {
    idx <- cmb[, j]
    z <- compute_subset_metrics(patient_data, idx)
    z$subset_id <- j
    z$subset_indices <- paste(idx, collapse = ",")
    z$subset_samples <- paste(patient_data$sample_ids[idx], collapse = ",")
    z
  })

  ans <- do.call(rbind, out)
  ans <- ans[, c(
    "subset_id", "subset_indices", "subset_samples",
    setdiff(names(ans), c("subset_id", "subset_indices", "subset_samples"))
  )]
  rownames(ans) <- NULL
  ans
}
