#' Rank Energy Change Point Test via Permutation (REp)
#'
#' Detects a single change point in a multivariate sequence using a rank-based
#' energy statistic. Ranks are obtained via optimal transport to a Halton grid,
#' making the method robust to heavy-tailed distributions. The p-value is
#' computed via permutation.
#'
#' @param x A numeric matrix of observations with dimension \code{n x d},
#'   where rows represent observations in sequence.
#' @param n.perm A positive integer specifying the number of permutation
#'   replications. Defaults to \code{999}.
#' @param beta A positive numeric scalar specifying the power applied to
#'   the pairwise distance matrix. Defaults to \code{1}.
#' @param min_size A positive integer specifying the minimum slice size.
#'   If \code{NULL} (default), it is set automatically: \code{2} when
#'   \code{n < 30} and \code{10} otherwise.
#'
#' @return An object of class \code{"rcptest"} with the following components:
#'   \describe{
#'     \item{statistic}{The observed test statistic.}
#'     \item{p.value}{P-value computed via permutation.}
#'     \item{location}{Estimated change point location.}
#'     \item{method}{Character string \code{"REp"}.}
#'     \item{n}{Number of observations.}
#'     \item{n.perm}{Number of permutation replications.}
#'     \item{beta}{Power applied to the pairwise distance matrix.}
#'   }
#'
#' @seealso \code{\link{rank_energy_eig}}
#'
#' @export
rank_energy_perm <- function(x, n.perm = 999, beta = 1, min_size = NULL){
  n <- nrow(x); d <- ncol(x)

  if (is.null(min_size)){
    min_size <- if (n < 30) 2 else 10
  }

  # optimal transport ranks
  D_rank  <- compute_ranks(x,beta = beta)

  # observed statistic
  stats_obs <- split_re(s_ = 1, e_ = n, D_ = D_rank, min_size_ = min_size)
  stat      <- stats_obs[2]
  location  <- stats_obs[1]

  # permutation distribution
  stat0 <- sapply(1:n.perm, function(i) perm_once(D_rank, n, min_size))
  p.value <- (sum(stat0 >= stat) + 1) / (n.perm + 1)

  return(structure(
    list(
      statistic = stat,
      p.value   = p.value,
      location  = location,
      method    = "REp",
      n         = n,
      n.perm    = n.perm,
      beta      = beta
    ),
    class = "recp"
  ))
}
#' Rank Energy Change Point Test via Limiting Distribution (REs)
#'
#' Detects a single change point in a multivariate sequence using a rank-based
#' energy statistic. Ranks are obtained via optimal transport to a Halton grid,
#' making the method robust to heavy-tailed distributions. The p-value is
#' computed via asymptotic eigenvalue approximation using Brownian bridge
#' simulations.
#'
#' @param x A numeric matrix of observations with dimension \code{n x d},
#'   where rows represent observations in sequence.
#' @param n.perm A positive integer specifying the number of simulation
#'   replications for the asymptotic approximation. Defaults to \code{999}.
#' @param beta A positive numeric scalar specifying the power applied to
#'   the pairwise distance matrix. Defaults to \code{1}.
#' @param min_size A positive integer specifying the minimum slice size.
#'   If \code{NULL} (default), it is set automatically: \code{2} when
#'   \code{n < 30} and \code{10} otherwise.
#' @param eig_num A positive integer specifying the number of eigenvalues
#'   used in the asymptotic approximation. If \code{NULL} (default), it is
#'   set automatically: \code{50} when \code{n >= 50} and \code{n} otherwise.
#'
#' @return An object of class \code{"rcptest"} with the following components:
#'   \describe{
#'     \item{statistic}{The observed test statistic.}
#'     \item{p.value}{P-value computed via asymptotic eigenvalue approximation.}
#'     \item{location}{Estimated change point location.}
#'     \item{method}{Character string \code{"REs"}.}
#'     \item{n}{Number of observations.}
#'     \item{n.perm}{Number of simulation replications.}
#'     \item{beta}{Power applied to the pairwise distance matrix.}
#'     \item{eig_num}{Number of eigenvalues used in the approximation.}
#'   }
#'
#' @seealso \code{\link{rank_energy_perm}}
#'
#' @export
rank_energy_eig <- function(x, n.perm = 999, beta = 1,
                            min_size = NULL, eig_num = NULL){
  n <- nrow(x); d <- ncol(x)

  if (is.null(min_size)){
    min_size <- if (n < 30) 2 else 10
  }
  if (is.null(eig_num)){
    eig_num <- if (n >= 50) 50 else n
  }

  # optimal transport ranks
  D_rank  <- compute_ranks(x,beta = beta)

  # observed statistic
  stats_obs <- split_re(s_ = 1, e_ = n, D_ = D_rank, min_size_ = min_size)
  stat      <- stats_obs[2]
  location  <- stats_obs[1]

  # eigenvalue approximation
  eig_Hn <- eigencompute(D_rank, n)
  eigs   <- eig_Hn[order(abs(eig_Hn), decreasing = TRUE)[1:eig_num]]

  time_ind <- seq(0, 1, length = 1000)
  ttt      <- time_ind * (1 - time_ind)

  stat0 <- sapply(1:n.perm, function(i){
    rr  <- sapply(1:eig_num, function(j){
      sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = 999)^2
    })
    Yt <- colSums(abs(eigs) * t(ttt - rr))
    return(max(abs(Yt)))
  })

  p.value <- (sum(stat0 >= stat) + 1) / (n.perm + 1)

  return(structure(
    list(
      statistic = stat,
      p.value   = p.value,
      location  = location,
      method    = "REs",
      n         = n,
      n.perm    = n.perm,
      beta      = beta,
      eig_num   = eig_num
    ),
    class = "recp"
  ))
}

#' Print method for recp objects
#'
#' @param x An object of class \code{"recp"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
print.recp <- function(x, ...){
  cat("\n\tRank Energy Change Point Test\n\n")
  cat("Method:", x$method, "\n")
  cat("Statistic:", round(x$statistic, 4), "\n")
  cat("Estimated change point location:", x$location, "\n")
  cat("P-value:", round(x$p.value, 4), "\n")
  invisible(x)
}

#' Centered Kernel Matrix Eigendecomposition
#'
#' Computes the eigenvalues of the centered kernel matrix \eqn{H_n}, defined
#' as \eqn{H_n = C D C^\top / n}, where \eqn{C = I_n - \mathbf{1}\mathbf{1}^\top / n}
#' is the centering matrix and \eqn{D} is a pairwise distance matrix.
#'
#' @param D A numeric \code{n x n} pairwise distance matrix.
#' @param n A positive integer giving the number of observations.
#'
#' @return A vector returned by \code{eigen()}, containing only eigenvalues
#'   (computed with \code{only.values = TRUE}).
#'
#' @keywords internal
eigencompute <- function(D,n){
  cc <- diag(n)-1/n*matrix(rep(1,n),ncol=1)%*%t(matrix(rep(1,n),ncol=1))
  #centering the kernel matrix D#
  Hn <- cc%*%D%*%t(cc)/n
  eig_Hn <- eigen(Hn,only.values = T)
  return(eig_Hn$values)
}

#' Compute Optimal Transport Ranks
#'
#' Maps observations to multivariate ranks via optimal transport to a Halton
#' grid, and returns the powered pairwise distance matrix of the ranked points.
#'
#' @param x A numeric matrix of observations with dimension \code{n x d},
#'   where rows represent observations in sequence.
#' @param beta A positive numeric scalar specifying the power applied to
#'   the pairwise distance matrix.
#'
#' @return A numeric \code{n x n} matrix of powered pairwise distances among
#'   the rank-transformed observations.
#'
#' @keywords internal
compute_ranks <- function(x, beta){
  n <- nrow(x); d <- ncol(x)
  grid     <- matrix(randtoolbox::halton(n, dim = d), ncol = d)
  distmat  <- Rfast::dista(xnew = x, x = grid)^2
  solution <- transport::transport(rep(1/n, n), rep(1/n, n), p = 2, method = 'networkflow', costm = distmat)
  co_rank <- matrix(grid[solution[, 2],], ncol = d)
  return((Rfast::Dist(co_rank))^beta)
}

#' Single Permutation Statistic
#'
#' Computes the energy statistic for a single random permutation of the
#' distance matrix, used in the permutation distribution of REp.
#'
#' @param D A numeric \code{n x n} pairwise distance matrix.
#' @param n A positive integer giving the number of observations.
#' @param min_size A positive integer specifying the minimum slice size.
#'
#' @return A numeric scalar giving the permuted test statistic.
#'
#' @keywords internal
perm_once <- function(D, n, min_size){
  shuffle <- sample(1:n)
  return(split_re(s_ = 1, e_ = n, D_ = D[shuffle, shuffle], min_size_ = min_size)[2])
}
