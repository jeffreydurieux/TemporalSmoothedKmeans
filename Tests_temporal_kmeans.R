# Mon 9 12:01:34 2025 ------------------------------
# temporal smoothing k-means

## ---------------------------------------------------------------
## Temporal k-means with centroid smoothing + Hungarian matching
## ---------------------------------------------------------------

temporal_kmeans <- function(X_list, K, eta = 1, max_iter = 50, tol = 1e-4,
                            use_hungarian = TRUE, verbose = TRUE) {
  # X_list: list of length T, each element: matrix n_t x p
  # K:      number of clusters
  # eta:    temporal smoothness penalty
  # max_iter: maximum number of EM-like iterations
  # tol:    convergence threshold on centroid change
  # use_hungarian: logical, whether to enforce label continuity
  # verbose: print objective per iteration
  
  if (use_hungarian) {
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' is required for Hungarian matching. Install via install.packages('clue').")
    }
  }
  
  Tn <- length(X_list)
  if (Tn < 1) stop("X_list must have at least one time point.")
  
  p <- ncol(X_list[[1]])
  if (any(sapply(X_list, ncol) != p)) {
    stop("All time points must have same number of columns (features).")
  }
  
  # Helper: build 1D Laplacian matrix for penalty sum_{t=2}^T (m_t - m_{t-1})^2
  build_L <- function(Tn) {
    if (Tn == 1) return(matrix(0, 1, 1))
    L <- matrix(0, Tn, Tn)
    for (t in 1:Tn) {
      if (t == 1) {
        L[t, t] <- 1
        L[t, t + 1] <- -1
      } else if (t == Tn) {
        L[t, t] <- 1
        L[t, t - 1] <- -1
      } else {
        L[t, t] <- 2
        L[t, t - 1] <- -1
        L[t, t + 1] <- -1
      }
    }
    L
  }
  L <- build_L(Tn)
  
  ## -----------------------------------------------------------
  ## Initialization
  ## -----------------------------------------------------------
  
  # Pooled k-means for rough starting centers
  X_all <- do.call(rbind, X_list)
  km_init <- stats::kmeans(X_all, centers = K, nstart = 5)
  init_centers <- km_init$centers  # K x p
  
  # Initialize cluster assignments at each time by nearest pooled center
  z_list <- vector("list", Tn)          # cluster labels per time
  M_list <- vector("list", Tn)          # centroid matrices per time (K x p)
  
  for (t in 1:Tn) {
    X_t <- X_list[[t]]                 # n_t x p
    n_t <- nrow(X_t)
    dists <- as.matrix(dist(rbind(init_centers, X_t)))[1:K, (K+1):(K+n_t)]
    # assign to nearest initial centroid
    z_t <- apply(dists, 2, which.min)
    z_list[[t]] <- z_t
    
    # compute initial centroids at time t
    M_t <- matrix(NA_real_, K, p)
    for (k in 1:K) {
      idx <- which(z_t == k)
      if (length(idx) > 0) {
        M_t[k, ] <- colMeans(X_t[idx, , drop = FALSE])
      } else {
        # if empty cluster, fall back to global init center
        M_t[k, ] <- init_centers[k, ]
      }
    }
    M_list[[t]] <- M_t
  }
  
  ## -----------------------------------------------------------
  ## Objective function (for monitoring)
  ## -----------------------------------------------------------
  
  compute_objective <- function(X_list, z_list, M_list, eta, L) {
    Tn <- length(X_list)
    K <- nrow(M_list[[1]])
    
    # Within-cluster sum of squares
    wcss <- 0
    for (t in 1:Tn) {
      X_t <- X_list[[t]]
      M_t <- M_list[[t]]
      z_t <- z_list[[t]]
      for (k in 1:K) {
        idx <- which(z_t == k)
        if (length(idx) > 0) {
          diff <- X_t[idx, , drop = FALSE] - 
            matrix(M_t[k, ], nrow = length(idx),
                   ncol = ncol(X_t), byrow = TRUE)
          wcss <- wcss + sum(diff^2)
        }
      }
    }
    
    # Temporal smoothness penalty
    penalty <- 0
    if (Tn > 1) {
      for (t in 2:Tn) {
        diff <- M_list[[t]] - M_list[[t - 1]]
        penalty <- penalty + sum(diff^2)
      }
    }
    
    wcss + eta * penalty
  }
  
  ## -----------------------------------------------------------
  ## Hungarian matching helper
  ## -----------------------------------------------------------
  
  align_labels_hungarian <- function(M_list, z_list) {
    Tn <- length(M_list)
    K  <- nrow(M_list[[1]])
    
    if (Tn <= 1) return(list(M_list = M_list, z_list = z_list))
    
    for (t in 2:Tn) {
      M_prev <- M_list[[t - 1]]  # K x p (reference)
      M_curr <- M_list[[t]]      # K x p (to be permuted)
      
      # Cost matrix: rows = prev clusters, cols = current clusters
      cost <- matrix(0, K, K)
      for (j in 1:K) {
        for (i in 1:K) {
          diff <- M_prev[j, ] - M_curr[i, ]
          cost[j, i] <- sum(diff^2)
        }
      }
      
      # Hungarian matching: perm[j] = index of current cluster assigned to prev cluster j
      perm <- as.vector(clue::solve_LSAP(cost))
      
      # Reorder current centroids so that row j corresponds to same "identity" as in t-1
      M_curr_new <- M_curr[perm, , drop = FALSE]
      
      # Relabel current assignments:
      # old labels are in {1,...,K}; each old label i maps to new label j with perm[j] == i
      map_old_to_new <- match(1:K, perm)  # for each old i, gives new label j
      z_curr_old <- z_list[[t]]
      z_curr_new <- map_old_to_new[z_curr_old]
      
      M_list[[t]] <- M_curr_new
      z_list[[t]] <- z_curr_new
    }
    
    list(M_list = M_list, z_list = z_list)
  }
  
  ## -----------------------------------------------------------
  ## Main loop: assignment + centroid smoothing + Hungarian match
  ## -----------------------------------------------------------
  
  obj_values <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    # 1) Assignment step: given M_list, update z_list
    for (t in 1:Tn) {
      X_t <- X_list[[t]]         # n_t x p
      M_t <- M_list[[t]]         # K x p
      n_t <- nrow(X_t)
      
      # compute squared distances to centroids at time t
      d2 <- matrix(0, n_t, K)
      for (k in 1:K) {
        diff <- X_t - matrix(M_t[k, ], nrow = n_t, ncol = p, byrow = TRUE)
        d2[, k] <- rowSums(diff^2)
      }
      z_list[[t]] <- max.col(-d2)  # nearest centroid
    }
    
    # 2) Centroid smoothing step: given z_list, update M_list
    M_list_old <- M_list
    
    for (k in 1:K) {
      # For each cluster k, build n_k_t and xbar_k_t over time
      n_k_t <- numeric(Tn)
      xbar_k_t <- matrix(0, Tn, p)
      
      for (t in 1:Tn) {
        X_t <- X_list[[t]]
        z_t <- z_list[[t]]
        idx <- which(z_t == k)
        if (length(idx) > 0) {
          n_k_t[t] <- length(idx)
          xbar_k_t[t, ] <- colMeans(X_t[idx, , drop = FALSE])
        } else {
          n_k_t[t] <- 0
          xbar_k_t[t, ] <- 0
        }
      }
      
      # For each feature j, solve (A + eta*L) m = A * xbar
      A <- diag(n_k_t, nrow = Tn, ncol = Tn)
      
      for (j in 1:p) {
        rhs <- n_k_t * xbar_k_t[, j]
        if (Tn == 1) {
          m_hat <- if (n_k_t[1] > 0) rhs / n_k_t[1] else 0
        } else {
          M_mat <- A + eta * L
          m_hat <- solve(M_mat, rhs)
        }
        
        # Put m_hat into M_list
        for (t in 1:Tn) {
          if (is.null(M_list[[t]])) {
            M_list[[t]] <- matrix(0, K, p)
          }
          M_list[[t]][k, j] <- m_hat[t]
        }
      }
    }
    
    # 3) Hungarian matching to enforce label continuity
    if (use_hungarian && K > 1 && Tn > 1) {
      aligned <- align_labels_hungarian(M_list, z_list)
      M_list <- aligned$M_list
      z_list <- aligned$z_list
    }
    
    # 4) Check convergence (after smoothing + matching)
    max_change <- 0
    for (t in 1:Tn) {
      diff <- M_list[[t]] - M_list_old[[t]]
      max_change <- max(max_change, sqrt(sum(diff^2)))
    }
    
    obj_values[iter] <- compute_objective(X_list, z_list, M_list, eta, L)
    
    if (verbose) {
      cat(sprintf("Iter %3d: objective = %.4f, max centroid change = %.6f\n",
                  iter, obj_values[iter], max_change))
    }
    
    if (max_change < tol) {
      if (verbose) cat("Converged.\n")
      obj_values <- obj_values[1:iter]
      break
    }
  }
  
  list(
    centers = M_list,        # list of K x p centroid matrices for each time
    cluster = z_list,        # list of integer vectors for each time
    objective = obj_values,
    K = K,
    eta = eta,
    use_hungarian = use_hungarian
  )
}

## ---------------------------------------------------------------
## Minimal example: drifting 2D clusters over 5 time points
## ---------------------------------------------------------------

set.seed(123)

simulate_temporal_data <- function(Tn = 5, n_per_cluster = 100, K = 2) {
  X_list <- vector("list", Tn)
  for (t in 1:Tn) {
    # cluster 1 center drifts right, cluster 2 drifts up
    c1 <- c(0 + 0.7 * (t - 1), 0)
    c2 <- c(3, 3 + 0.7 * (t - 1))
    
    x1 <- sweep(matrix(rnorm(n_per_cluster * 2), ncol = 2), 2, c1, "+")
    x2 <- sweep(matrix(rnorm(n_per_cluster * 2), ncol = 2), 2, c2, "+")
    X_list$X[[t]] <- rbind(x1, x2)
    X_list$Centroids[[t]] <- rbind(c1, c2)
  }
  X_list
}

# Simulate data
X_list <- simulate_temporal_data(Tn = 5, n_per_cluster = 50, K = 2)

# Fit temporal k-means with Hungarian label alignment
Xlist_dat <- X_list$X
fit <- temporal_kmeans(Xlist_dat, K = 2, eta = 1, max_iter = 50, tol = 1e-3,
                       use_hungarian = TRUE)

# Quick visualization for the last time point
t_plot <- length(X_list)
X_t <- X_list[[t_plot]]
z_t <- fit$cluster[[t_plot]]
M_t <- fit$centers[[t_plot]]

plot(X_t, col = z_t, pch = 19, main = paste("Time", t_plot))
points(M_t, pch = 4, cex = 2, lwd = 2)





## ---------------------------------------------------------------
## Plot centroid trajectories over time
## ---------------------------------------------------------------
## Assumes output from temporal_kmeans()
## centers: list of length T, each element K x p
## ---------------------------------------------------------------

plot_centroid_trajectories <- function(fit, mode = c("2d", "time_series"),
                                       feature_indices = c(1, 2),
                                       main = NULL) {
  mode <- match.arg(mode)
  
  M_list <- fit$centers
  Tn <- length(M_list)
  K  <- nrow(M_list[[1]])
  p  <- ncol(M_list[[1]])
  
  if (is.null(main)) main <- "Centroid trajectories"
  
  ## -----------------------------------------------------------
  ## Mode 1: 2D trajectories in feature space (p must be >= 2)
  ## -----------------------------------------------------------
  if (mode == "2d") {
    if (p < 2) stop("Need at least 2 features for mode = '2d'.")
    
    f1 <- feature_indices[1]
    f2 <- feature_indices[2]
    
    # Collect all centroid positions for range
    all_pts <- do.call(rbind, lapply(M_list, function(M) M[, c(f1, f2), drop = FALSE]))
    x_range <- range(all_pts[, 1])
    y_range <- range(all_pts[, 2])
    
    plot(NA, xlim = x_range, ylim = y_range,
         xlab = paste("Feature", f1), ylab = paste("Feature", f2),
         main = main)
    
    cols <- 1:K  # basic colours per cluster
    
    for (k in 1:K) {
      # Extract trajectory of cluster k over time
      traj <- sapply(M_list, function(M) M[k, c(f1, f2)])
      traj <- t(traj)      # T x 2
      lines(traj[, 1], traj[, 2], col = cols[k], lwd = 2)
      points(traj[, 1], traj[, 2], col = cols[k], pch = 19)
      
      # mark start and end
      points(traj[1, 1],  traj[1, 2],  col = cols[k], pch = 15, cex = 1.5)
      points(traj[Tn, 1], traj[Tn, 2], col = cols[k], pch = 17, cex = 1.5)
    }
    
    legend("topleft", legend = paste("Cluster", 1:K),
           col = cols, lwd = 2, pch = 19, bty = "n")
  }
  
  ## -----------------------------------------------------------
  ## Mode 2: time-series of centroid coordinates
  ## -----------------------------------------------------------
  if (mode == "time_series") {
    time <- 1:Tn
    
    # If p > 1, facet by feature
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    par(mfrow = c(p, 1), mar = c(4, 4, 3, 1))
    
    cols <- 1:K
    
    for (j in 1:p) {
      mat <- sapply(M_list, function(M) M[, j])  # K x T
      mat <- t(mat)  # T x K
      
      matplot(time, mat, type = "l", lty = 1, lwd = 2,
              col = cols, xlab = "Time index",
              ylab = paste("Feature", j),
              main = paste(main, "- feature", j))
      points(rep(time, each = K),
             as.vector(mat),
             col = rep(cols, each = Tn),
             pch = 19, cex = 0.7)
      
      legend("topleft", legend = paste("Cluster", 1:K),
             col = cols, lwd = 2, pch = 19, bty = "n")
    }
  }
}


## 1) 2D trajectories in feature space (for 2D data)
plot_centroid_trajectories(fit, mode = "2d",
                           feature_indices = c(1, 2),
                           main = "Centroid paths in (x1, x2)")

## 2) Time series of centroid coordinates (general p)
plot_centroid_trajectories(fit, mode = "time_series",
                           main = "Centroid coordinates over time")


Xlist_dat <- X_list$X
fit <- temporal_kmeans(Xlist_dat, K = 2, eta = 10, max_iter = 50, tol = 1e-3,
                       use_hungarian = TRUE)



library(mclust)

cl <- c(rep(1,50), rep(2,50))
adjustedRandIndex(cl, fit$cluster[[1]])
adjustedRandIndex(cl, fit$cluster[[2]])
adjustedRandIndex(cl, fit$cluster[[3]])
adjustedRandIndex(cl, fit$cluster[[4]])
adjustedRandIndex(cl, fit$cluster[[5]])



psych::factor.congruence(X_list$Centroids[[1]], fit$centers[[1]])
psych::factor.congruence(X_list$Centroids[[2]], fit$centers[[2]])
psych::factor.congruence(X_list$Centroids[[3]], fit$centers[[3]])
psych::factor.congruence(X_list$Centroids[[4]], fit$centers[[4]])
psych::factor.congruence(X_list$Centroids[[5]], fit$centers[[5]])


fit$centers[[1]]
fit$centers[[2]]
fit$centers[[3]]
fit$centers[[4]]
fit$centers[[5]]
