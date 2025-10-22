plot_ecdf_with_sims <- function(m = 200, alpha = 0.05, nsim = 1e3, 
                                dist = "norm", dist_params = list(mean = NULL, sd = 1)) {
  theta_m <- sin(seq(m))
  
  # Helper to draw sample from specified distribution
  generate_sample <- function(size, mean_vec, distribution, params) {
    if (distribution == "norm") {
      mean <- if (is.null(params$mean)) mean_vec else params$mean
      sd <- if (is.null(params$sd)) 1 else params$sd
      return(rnorm(size, mean, sd))
    } else if (distribution == "unif") {
      minv <- if (is.null(params$min)) mean_vec else params$min
      maxv <- if (is.null(params$max)) mean_vec + 1 else params$max
      return(runif(size, minv, maxv))
    } else if (distribution == "pois") {
      lambda <- if (is.null(params$lambda)) abs(mean_vec) else params$lambda
      return(rpois(size, lambda))
    } else {
      stop("Unsupported distribution")
    }
  }
  
  Xs <- matrix(nrow = nsim, ncol = m)
  for(i in 1:nsim){
    x <- generate_sample(m, theta_m, dist, dist_params)
    x_sorted <- sort(x)
    Xs[i,] <- x_sorted
  }
  
  # take a new one (main sample)
  x <- generate_sample(m, theta_m, dist, dist_params)
  
  # ECDF and evaluation
  ecdf_fn <- ecdf(x)
  x_sorted <- sort(x)
  ecdf_x_sorted <- ecdf_fn(x_sorted)
  
  eps <- sqrt(-log(alpha/4)*2/m)
  
  u_bound <- pmin(ecdf_x_sorted + eps, 1)
  l_bound <- pmax(ecdf_x_sorted - eps, 0)
  
  # Plot all the simulated ECDFs (step lines), then the main ECDF with its bounds
  plot(x_sorted, ecdf_x_sorted, type="n", ylim=c(0,1), 
       main="ECDF with Confidence Bounds and Simulations",
       ylab="ECDF", xlab="x")
  
  # Overlay all simulated ECDFs in light gray
  for(i in 1:nsim){
    x_sim <- Xs[i,]
    ecdf_sim <- ecdf(x_sim)
    lines(x_sorted, ecdf_sim(x_sorted), type="s", col=rgb(0.7,0.7,0.7,0.15))
  }
  
  # Plot actual ECDF of the new sample
  lines(x_sorted, ecdf_x_sorted, type="s", lwd=2, col="black")
  
  # Plot bounds
  lines(x_sorted, u_bound, col="red", lty=2, lwd=2)
  lines(x_sorted, l_bound, col="blue", lty=2, lwd=2)
}

plot_ecdf_with_sims(dist = "unif", alpha = 0.5)