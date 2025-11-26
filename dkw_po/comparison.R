
# Compare our Makarov quantile method and save simulation results
set.seed(123)

source("generate.R")
source("bounds.R")
library(RIQITE)
library(MASS)
library(ggplot2)
plot.dir <- "../fig/"
compare_makarov_quantile <- function(params) {

  result <- do.call(get_bounds, params)
  bounds <- result$bounds
  data <- result$data
  obs_Y1 <- sort(result$obs_Y1)
  obs_Y0 <- sort(result$obs_Y0)
  # naive bounds
  naive_u <- max(obs_Y1) - min(obs_Y0)
  naive_l <- min(obs_Y1) - max(obs_Y0)
  # Compute Makarov bounds for all t values
  diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)

  quants_check <- seq(0,1,by = 0.25)
  quants_true <- quantile(data$tau, quants_check)
  conf_int_idx <- lapply(quants_check, \(z) which(diff_bounds[,2] >= z))
  # just do lower bound since that is what riqite does well in 
  conf_int_us <- lapply(conf_int_idx, \(z)  min(diff_bounds[z, 3]))
  conf_riqite <- ci_quantile(data$Z, data$Y, quants_check*nrow(data), 
                             nperm = 1e4, alpha = params$alpha) # because two sided
  for(i in length(conf_int_us)){
    if(abs(conf_int_us[[i]] - naive_l) < 0.25) {
      conf_int_us[[i]] <- -Inf
    }
  }

  cover_us <- quants_true > conf_int_us
  cover_riqite <- quants_true > conf_riqite$lower
  width_diff <- unlist(conf_int_us) - conf_riqite$lower

  return(
    data.frame(
      m = params$m,
      p = params$p,
      tau = params$tau,
      rho = params$rho,
      tau = params$tau,
      tight_bounds = params$tight_bounds,
      width_us = unlist(conf_int_us),
      width_riqite = conf_riqite$lower,
      cover_us = cover_us,
      cover_riqite = cover_riqite,
      quants_check = quants_check
    )
  )
}

set.seed(123) # for reproducibility
NSIM <- 20 # or set to desired number of simulations
m_values <- c(100, 500,1e3,2e3, 1e4)
sim_results <- data.frame()
for(rho in seq(0,1, length = 4)){
  for(tight_bounds in c(FALSE, TRUE)){
    for(tau in c(1, 2, 3)){
      for (size in m_values) {
        params <- list(
          m = size,
          p = 0.5,
          alpha = 0.1,
          tau = tau,
          rho = rho,
          tight_bounds = tight_bounds
        )
        
        for (sim in 1:NSIM) {
          print(sim)
          sim_results <- rbind(sim_results, compare_makarov_quantile(params))
        }
      }
    }
  }
  
}



# Optionally, save results to a file
write.csv(sim_results, file = "makarov_quantile_sim_results.csv", row.names = FALSE)



plot_quantile_width <- function(df, 
                               xvar = c("rho", "tau", "m"), 
                               tight_bounds_bool = FALSE) {
  xvar <- match.arg(xvar)
  if(!tight_bounds_bool){
    df <- df |> 
      dplyr::filter(tight_bounds == tight_bounds_bool,
                    m >= 1000)
    
  }else{
    df <- df |> dplyr::filter(tight_bounds == tight_bounds_bool)
  }
  df$width[is.infinite(df$width)] <- max(subset(sim_results, 
                                                !is.infinite(width))$width)
  
  # Ensure categorical if not numeric
  if (xvar %in% c("rho", "tau")) df[[xvar]] <- as.factor(round(df[[xvar]], 3))
  agg <- aggregate(width ~ quants_check + get(xvar), df, mean, na.rm = TRUE)
  xcol <- names(agg)[2]
  names(agg)[2] <- "xval"
  ggplot(agg, aes(x = xval, y = width, color = factor(quants_check), group = quants_check)) +
    geom_line(aes(linetype = factor(quants_check)), size = 1, alpha = 0.7) +
    geom_point() +
    labs(x = xvar, 
         y = "Magnitude Better Lower Bound", 
         color = "For Quantile", 
         linetype = NULL,
         title = sprintf("Comparing Better Lower Bound by %s", xvar)) +
    theme_minimal()
}

# figure 2
png(paste0(plot.dir, "comparison-by-size-study.png"))
plot_quantile_width(sim_results, xvar = "m", tight_bounds = FALSE)
dev.off()
png(paste0(plot.dir, "comparison-by-rho-study.png"))
plot_quantile_width(sim_results, xvar = "rho", tight_bounds = FALSE)
dev.off()
png(paste0(plot.dir, "comparison-by-tau-study.png"))
plot_quantile_width(sim_results, xvar = "tau", tight_bounds = FALSE)
dev.off()


# figure 4

png(paste0(plot.dir, "comparison-by-size-study-tight.png"))
plot_quantile_width(sim_results, xvar = "m", tight_bounds = TRUE)
dev.off()
png(paste0(plot.dir, "comparison-by-rho-study-tight.png"))
plot_quantile_width(sim_results, xvar = "rho", tight_bounds = TRUE)
dev.off()
png(paste0(plot.dir, "comparison-by-tau-study-tight.png"))
plot_quantile_width(sim_results, xvar = "tau", tight_bounds = TRUE)
dev.off()


