# Johnson SB PDFs - rescale all theta distributions using LNG scaling factor
# Output summary includes: mean, median, sd, skewness, kurtosis, Q1, Q3
# Updated to print both original and uniformly scaled (LNG) summary stats
# Updated: 2025-12-16

library(moments)
library(ggplot2)
library(reshape2)

# ---------------- Johnson SB PDF / helper functions ----------------
djsb <- function(x, gamma, delta) {
  z <- log(x / (1 - x))
  dens <- (delta / sqrt(2*pi)) * (1 / (x*(1-x))) *
    exp(-0.5 * (gamma + delta*z)^2)
  dens[x <= 0 | x >= 1] <- NA
  dens
}

# Logitnormal moments via integration (Johnson SB mapping)
logitnorm_mean <- function(mu, sigma) {
  integrate(function(z) plogis(z) * dnorm(z, mu, sigma),
            -Inf, Inf)$value
}
logitnorm_var <- function(mu, sigma) {
  m <- logitnorm_mean(mu, sigma)
  integrate(function(z) (plogis(z)-m)^2 * dnorm(z, mu, sigma),
            -Inf, Inf)$value
}

# ---------------- Parameter sets ----------------
param_sets <- list(
  list(gamma=-0.2,    delta=0.8742, label="Coal"),
  list(gamma=-1.0637, delta=1.0935, label="LNG"),
  list(gamma=1.2444,  delta=0.9692, label="Offshore wind"),
  list(gamma=1,       delta=1.4,    label="Solar PV")
)

labels <- sapply(param_sets, `[[`, "label")
cols <- c("#8B3A3A", "#607B8B", "#979461", "#CD5733")

# ---------------- Original means (for diagnostics) ----------------
orig_means <- sapply(param_sets, function(ps) {
  mu    <- -ps$gamma / ps$delta
  sigma <- 1 / ps$delta
  logitnorm_mean(mu, sigma)
})
names(orig_means) <- labels
cat("Original means (per technology):\n")
print(round(orig_means, 6))

# ---------------- Uniform scaling: use LNG scaling factor ----------------
# User requested uniform scaling by LNG scaling factor = 0.1436254
scaling_factor <- 0.1436254
scaling_factors <- rep(scaling_factor, length(param_sets))
names(scaling_factors) <- labels
cat("Using uniform scaling factor (LNG):", scaling_factor, "\n\n")

# ---------------- Sample and compute statistics ----------------
set.seed(1)
n <- 10000

theta_orig_list <- vector("list", length(param_sets))
theta_scaled_list <- vector("list", length(param_sets))

# Prepare containers for summary stats
summary_stats <- data.frame(
  Technology = labels,
  OrigMean = NA, ScaledMean = NA,
  OrigMedian = NA, ScaledMedian = NA,
  OrigSD = NA, ScaledSD = NA,
  OrigSkew = NA, ScaledSkew = NA,
  OrigKurtosis = NA, ScaledKurtosis = NA,
  OrigQ1 = NA, ScaledQ1 = NA,
  OrigQ3 = NA, ScaledQ3 = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(param_sets)) {
  ps <- param_sets[[i]]
  mu    <- -ps$gamma / ps$delta
  sigma <- 1 / ps$delta
  
  # Sample original theta (logitnormal via plogis of normal)
  theta_orig <- plogis(rnorm(n, mu, sigma))
  
  # Apply uniform LNG-based scaling to every technology
  theta_scaled <- scaling_factor * theta_orig
  # Clamp to [0, 1]
  theta_scaled <- pmax(pmin(theta_scaled, 1), 0)
  
  theta_orig_list[[i]] <- theta_orig
  theta_scaled_list[[i]] <- theta_scaled
  
  # Compute summary stats for original
  summary_stats$OrigMean[i]      <- mean(theta_orig)
  summary_stats$OrigMedian[i]    <- median(theta_orig)
  summary_stats$OrigSD[i]        <- sd(theta_orig)
  summary_stats$OrigSkew[i]      <- skewness(theta_orig)
  summary_stats$OrigKurtosis[i]  <- kurtosis(theta_orig)
  summary_stats$OrigQ1[i]        <- quantile(theta_orig, 0.25)
  summary_stats$OrigQ3[i]        <- quantile(theta_orig, 0.75)
  
  # Compute summary stats for scaled
  summary_stats$ScaledMean[i]      <- mean(theta_scaled)
  summary_stats$ScaledMedian[i]    <- median(theta_scaled)
  summary_stats$ScaledSD[i]        <- sd(theta_scaled)
  summary_stats$ScaledSkew[i]      <- skewness(theta_scaled)
  summary_stats$ScaledKurtosis[i]  <- kurtosis(theta_scaled)
  summary_stats$ScaledQ1[i]        <- quantile(theta_scaled, 0.25)
  summary_stats$ScaledQ3[i]        <- quantile(theta_scaled, 0.75)
}

# Round for printing
summary_print <- summary_stats
num_cols <- setdiff(names(summary_print), "Technology")
summary_print[num_cols] <- round(summary_print[num_cols], 6)

cat("Summary statistics (original vs uniformly scaled by LNG factor):\n")
print(summary_print)

# ---------------- Plot: overlayed densities (full range and zoom) ----------------
plot_overlayed_theta <- function(theta_list, means, medians, cols, labels, xlim = c(0,1),
                                 main_title, add_median=FALSE, add_mean=FALSE, draw_legend=FALSE) {
  plot(NULL, xlim = xlim, ylim = c(0, 20),
       xlab = expression(theta), ylab = "Density", main = main_title, cex.lab=1.15, bty="n")
  grid(nx=10, ny=8, col="#efefef", lty=1)
  for (i in seq_along(theta_list)) {
    dens <- density(theta_list[[i]], from=xlim[1], to=xlim[2])
    polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0,length(dens$y))),
            col = adjustcolor(cols[i], 0.25), border = NA)
    lines(dens, col = cols[i], lwd = 2.2)
    if(add_mean)   abline(v = means[i], col = cols[i], lwd=1.8, lty=1)
    if(add_median) abline(v = medians[i], col = cols[i], lwd=1.8, lty=2)
  }
  if(draw_legend) legend("topright", legend = labels, col = cols, lwd=3, bty="n", cex=0.9)
}

# Compute scaled means/medians vectors for plotting purposes
theta_means_scaled <- sapply(theta_scaled_list, mean)
theta_medians_scaled <- sapply(theta_scaled_list, median)

par(mfrow = c(1,2), mar = c(4, 4, 3, 1) + 0.1)

# Full range 0-1
plot_overlayed_theta(theta_scaled_list, theta_means_scaled, theta_medians_scaled,
                     cols, labels, xlim = c(0, 1),
                     main_title = "Rescaled Disruption Distributions (Uniform LNG Scaling)",
                     add_mean = FALSE, add_median = FALSE, draw_legend = FALSE)

# Zoomed-in (close-up)
plot_overlayed_theta(theta_scaled_list, theta_means_scaled, theta_medians_scaled,
                     cols, labels, xlim = c(0, 0.45),
                     main_title = "Zoom (0 - 0.45)",
                     add_mean = FALSE, add_median = FALSE, draw_legend = TRUE)

par(mfrow = c(1,1))

# ---------------- ggplot version (faceted) ----------------
# Combine for ggplot (original vs uniformly scaled)
theta_df <- data.frame(
  value = c(unlist(theta_orig_list), unlist(theta_scaled_list)),
  tech = rep(rep(labels, each=n), 2),
  version = rep(c("Original", "Scaled (LNG)"), each = n * length(labels))
)

ggplot(theta_df, aes(x = value, color = version, fill = version)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~ tech, nrow = 2) +
  scale_color_manual(values = c("Original" = "grey40", "Scaled (LNG)" = "black")) +
  scale_fill_manual(values = c("Original" = "grey60", "Scaled (LNG)" = "black")) +
  labs(title = "Theta Distributions: Original vs. Uniformly Scaled (LNG factor)",
       x = expression(theta), y = "Density", color="Version", fill="Version") +
  theme_minimal()