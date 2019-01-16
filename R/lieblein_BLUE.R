liebP <- function (n, m, i, t) {

  # symbol translation from the Lieblein paper to the R function dhyper
  # dhyper m = i
  # dhyper x = t
  # dhyper n = n - i
  # dhyper k = m
  dhyper(x = t, m = i, n = n - i, k = m)
}

blueCoeff <- function (m) {


  if (m == 2) {
    value <- list(a = c(0.916373, 0.083627),
                  b = c(-0.721348, 0.721348))
  } else if (m == 3) {
    value <- list(a = c(0.656320, 0.255714, 0.087966),
                  b = c(-0.630541, 0.255816, 0.374725))
  } else if (m == 4) {
    value <- list(a = c(0.510998, 0.263943, 0.153680, 0.071380),
                  b = c(-0.558619, 0.085903, 0.223919, 0.248797))
  } else if (m == 5) {
    value <- list(a = c(0.418934, 0.246282, 0.167609, 0.108824, 0.058350),
                  b = c(-0.503127, 0.006534, 0.130455, 0.181656, 0.184483))
  } else if (m == 6) {
    value <- list(a = c(0.355450, 0.225488, 0.165620, 0.121054, 0.083522, 0.048867),
                  b = c(-0.459273, -0.035992, 0.073199, 0.126724, 0.149534, 0.145807))
  } else if (m == 7) {
    value <- list(a = c(0.309008, 0.206260, 0.158590, 0.123223, 0.093747,
                        0.067331, 0.041841),
                  b = c(-0.423700, -0.060698, 0.036192, 0.087339, 0.114868,
                        0.125859, 0.120141))
  } else if (m == 8) {
    value <- list(a = c(0.273535, 0.189428, 0.150200, 0.121174,
                        0.097142, 0.075904, 0.056132, 0.036485),
                  b = c(-0.394187, -0.075767, 0.011124, 0.058928,
                        0.087162, 0.102728, 0.108074, 0.101936))
  } else if (m == 9) {
    value <- list(a = c(0.245539, 0.174882, 0.141789, 0.117357, 0.097218,
                        0.079569, 0.063400, 0.047957, 0.032291),
                  b = c(-0.369242, -0.085203, -0.006486, 0.037977, 0.065574,
                        0.082654, 0.091965, 0.094369, 0.088391))
  } else if (m == 10) {
    value <- list(a = c(0.222867, 0.162308, 0.133845, 0.112868, 0.095636,
                        0.080618, 0.066988, 0.054193, 0.041748, 0.028929),
                  b = c(-0.347830, -0.091158, -0.019210, 0.022179, 0.048671,
                        0.066064, 0.077021, 0.082771, 0.083552, 0.077940))
  } else if (m == 11) {
    value <- list(a = c(0.204123, 0.151384, 0.126522, 0.108226, 0.093234, 0.080222,
                        0.068485, 0.057578, 0.047159, 0.036886, 0.026180),
                  b = c(-0.329210, -0.094869, -0.028604, 0.010032, 0.035284, 0.052464,
                        0.064071, 0.071381, 0.074977, 0.074830, 0.069644))
  } else if (m == 12) {
    value <- list(a = c(0.188361, 0.141833, 0.119838, 0.103673, 0.090455, 0.079018,
                        0.068747, 0.059266, 0.050303, 0.041628, 0.032984, 0.023894),
                  b = c(-0.312840, -0.097086, -0.035655, 0.000534, 0.024548,
                        0.041278, 0.053053, 0.061112, 0.066122, 0.068357,
                        0.067671, 0.062906))
  } else if (m == 13) {
    value <- list(a = c(0.174916, 0.133422, 0.113759, 0.099323, 0.087540,
                        0.077368, 0.068264, 0.059900, 0.052047,
                        0.044528, 0.037177, 0.029790, 0.021965),
                  b = c(-0.298313, -0.098284, -0.041013, -0.006997,
                        0.015836, 0.032014, 0.043710, 0.052101, 0.057862,
                        0.061355, 0.062699, 0.061699, 0.057330))
  } else if (m == 14) {
    value <- list(a = c(0.163309, 0.125966, 0.108230, 0.095223, 0.084619,
                        0.075484, 0.067331, 0.059866, 0.052891, 0.046260,
                        0.039847, 0.033526, 0.027131, 0.020317),
                  b = c(-0.285316, -0.098775, -0.045120, -0.013039,
                        0.008690, 0.024282, 0.035768, 0.044262, 0.050418,
                        0.054624, 0.057083, 0.057829, 0.056652, 0.052642))
  } else if (m == 15) {
    value <- list(a = c(0.153184, 0.119134, 0.103196, 0.091384, 0.081767,
                        0.073495, 0.066128, 0.059401, 0.053140, 0.047217,
                        0.041529, 0.035984, 0.030484, 0.024887, 0.018894),
                  b = c(-0.273606, -0.098768, -0.048285, -0.017934, 0.002773,
                        0.017779, 0.028988, 0.037452, 0.043798, 0.048415,
                        0.051534, 0.053267, 0.053603, 0.052334, 0.048648))
  } else if (m == 16) {
    value <- list(a = c(0.144271, 0.113346, 0.098600, 0.087801,
                        0.079021, 0.071476, 0.064771, 0.058660,
                        0.052989, 0.047646, 0.042539, 0.037597,
                        0.032748, 0.027911, 0.022969, 0.017653),
                  b = c(-0.262990, -0.098406, -0.050731, -0.021933,
                        -0.002167, 0.012270, 0.023168, 0.031528,
                        0.037939, 0.042787, 0.046308, 0.048646,
                        0.049860, 0.049912, 0.048602, 0.045207))
  } else {
    stop('must have integer n between 2 and 16 inclusively for blueCoeff')
  }

  value
}

# #' @export
liebCoeff <- function (n, m = NULL, blue_if_possible = TRUE) {

  if (n < 2) {
    stop('libCoeff must have n >= 2')
  } else if (n >= 2 && n <= 16 && blue_if_possible) {
    value <- blueCoeff(m = n)
  } else {

    if (is.null(m)) {

      m <- min(n - 1, 16)
      blue_coeffs <- blueCoeff(m = m)
      a_prime <- matrix(nrow = m, ncol = n)
      b_prime <- matrix(nrow = m, ncol = n)
      for (i in seq_len(n)) {

        a_prime[, i] <- blue_coeffs$a*
          (seq_len(m)/i)*
          liebP(n = n, m = m, i = i, t = seq_len(m))
        b_prime[, i] <- blue_coeffs$b*
          (seq_len(m)/i)*
          liebP(n = n, m = m, i = i, t = seq_len(m))
      }
      value <- list(a = colSums(a_prime),
                    b = colSums(b_prime))
    } else if (!is.null(m) && m >= 2 && m <= 16 && m < n) {

      blue_coeffs <- blueCoeff(m = m)
      a_prime <- matrix(nrow = m, ncol = n)
      b_prime <- matrix(nrow = m, ncol = n)
      for (i in seq_len(n)) {

        a_prime[, i] <- blue_coeffs$a*
          (seq_len(m)/i)*
          liebP(n = n, m = m, i = i, t = seq_len(m))
        b_prime[, i] <- blue_coeffs$b*
          (seq_len(m)/i)*
          liebP(n = n, m = m, i = i, t = seq_len(m))
      }
      value <- list(a = colSums(a_prime),
                    b = colSums(b_prime))
    } else {
      stop('non NULL m in liebCoeff must be an integer between 2 and 16 inclusively and less than n')
    }
  }

  value
}

# #' @export
liebBLUEParamEst <- function (x) {

  xn <- sort(x, decreasing = FALSE)

  blue_coeffs <- liebCoeff(length(x))

  mu <- sum(blue_coeffs$a*xn)
  sigma <- sum(blue_coeffs$b*xn)

  list(par = c(mu, sigma),
       y = x)
}

partSeries <- function (x, n_parts, n_buffer) {

  n_col <- ceiling(length(x)/n_parts)
  if (n_col - 2*n_buffer - 1 <= 0) {

    stop("must have n_col - 2*n_buffer - 1 > 0\n")
  }
  n_total <- n_parts*n_col
  new_x <- c(x, rep(x = -Inf, times = n_total - length(x)))
  new_x <- matrix(data = new_x, nrow = n_parts, ncol = n_col, byrow = TRUE)

  partitioned_series <- double(n_parts)
  for (i in 1:n_parts) {

    if (i == 1) {

      partitioned_series[i] <- max(new_x[i, 1:(n_col - n_buffer)])
    } else if (i == n_parts) {

      partitioned_series[i] <- max(new_x[i, (n_buffer + 1):n_col])
    } else {

      partitioned_series[i] <- max(new_x[i, (n_buffer + 1):(n_col - n_buffer)])
    }
  }

  partitioned_series
}

# #' @export
liebBLUE <- function (x, n_parts, target_n_parts = n_parts[1], probs = 0.5704) {

  partitioned_series <- partSeries(x = x, n_parts = n_parts[1], n_buffer = 0)

  lieb_params <- liebBLUEParamEst(x = partitioned_series)

  results <- matrix(nrow = length(probs), ncol = length(target_n_parts))
  colnames(results) <- as.character(target_n_parts)
  rownames(results) <- as.character(probs)

  for (i in seq_along(probs)) {

    for (j in seq_along(target_n_parts)) {

      results[i, j] <- lieb_params$par[2]*log(target_n_parts[j]) -
        lieb_params$par[2]*log(-log(probs[i])) + lieb_params$par[1]
    }
  }

  if (length(probs) == 1 && length(target_n_parts) > 1) {

    results <- results[1, ]
  } else if (length(probs) > 1 && length(target_n_parts) == 1) {

    results <- results[, 1]
  } else if (length(probs) == 1 && length(target_n_parts) == 1) {

    results <- results[1, 1]
  }

  results
}
