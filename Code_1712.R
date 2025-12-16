# Updated, streamlined version of Code_1712.R
# - Consolidated libraries
# - Fixed run_scenario signature and return values
# - Added mean_gen_by_tech_year and sim_realized_gen to return
# - Consolidated plotting helpers, robust color mapping
# - Compact export + zip

library(readxl)
library(nloptr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SuppDists)
library(scales)

# ---------------- USER/PARAMETERS ----------------
input_path   <- "C:/Users/Lunac/OneDrive/Dokumente/UCPH/Master Thesis/Data/model_data.xlsx"
sheet_name   <- "R_technologies"
CARBON_FEES  <- c(48.24, 0)
SCEN_LABELS  <- c("With Carbon Fee", "No Carbon Fee")
SCEN_COLORS  <- c("With Carbon Fee" = "#3A4E48", "No Carbon Fee" = "#BEB0A7")
T            <- 10
nsim         <- 30
set.seed(123)

# ---------------- LOAD DATA ----------------
df <- read_excel(input_path, sheet = sheet_name) %>%
  rename(
    tech = tech, I = I_i, c = c_i,
    starting_cap = starting_cap, capfac = cap_factor,
    e = e_i, mean_theta = mean_theta, eta = depreciation
  )
techs <- df$tech; N <- nrow(df)
cap_scale <- 1e3; gen_scale <- 1e3

d_raw <- 246344038.11; s <- 0.95; sd_raw <- s * d_raw
d <- sd_raw / gen_scale # demand in GWh/year

I_vec <- as.numeric(df$I); c_vec <- as.numeric(df$c)
capfac_vec <- as.numeric(df$capfac); eta_vec <- as.numeric(df$eta)
x0_raw <- as.numeric(df$starting_cap); x0 <- x0_raw / cap_scale
emissions <- as.numeric(df$e)
# Johnson SB params
gamma_vec <- df$gamma; delta_vec <- df$delta; xi_vec <- df$xi; lambda_vec <- df$lambda

# Proportional theta scaling (kept as before)
orig_means <- as.numeric(df$mean_theta)
scaling_factor <- 0.143616257
scaled_mean_theta <- orig_means * scaling_factor
scaling_factors <- rep(scaling_factor, length(orig_means))

theta_scaling_table <- data.frame(Technology = techs, OrigMeanTheta = orig_means, ScalingFactor = scaling_factors, ScaledMeanTheta = scaled_mean_theta)
write.csv(theta_scaling_table, "theta_scaling_summary.csv", row.names = FALSE)

# ---------------- SCENARIO RUNNER ----------------
run_scenario <- function(carbon_fee, T, nsim,
                         gamma_vec, delta_vec, xi_vec, lambda_vec,
                         scaling_factors, capfac_vec, I_vec, c_vec, eta_vec, x0, emissions, d, techs, scaled_mean_theta) {
  N <- length(techs)
  nvars <- 2 * N * T
  idx_invest <- function(i, t) (i-1)*T + t
  idx_gen    <- function(i, t) N*T + (i-1)*T + t
  I_scaled <- I_vec * 1e6

  compute_capacity_path <- function(v){
    X <- matrix(0, nrow=N, ncol=T+1); X[,1] <- x0
    for(t in 1:T) for(i in 1:N) X[i, t+1] <- (1 - eta_vec[i]) * X[i, t] + v[idx_invest(i, t)]
    X
  }

  objective_fn <- function(v){
    X <- compute_capacity_path(v); total <- 0
    for(t in 1:T){
      disc <- 1 / ((1 + 0.03)^(t - 1))
      for(i in 1:N){
        inv_gw <- v[idx_invest(i,t)]; gen_gwh <- v[idx_gen(i,t)]
        gen_mwh <- gen_gwh * gen_scale
        total <- total + disc * ( I_scaled[i] * inv_gw^2 + c_vec[i] * gen_mwh + carbon_fee * emissions[i] * gen_mwh )
      }
    }
    as.numeric(total)
  }

  eval_g_ineq <- function(v){
    X <- compute_capacity_path(v)
    constr <- c()
    for(t in 1:T){
      sumg_gwh <- sum(sapply(1:N, function(i) v[idx_gen(i,t)]))
      constr <- c(constr, d - sumg_gwh)
    }
    for(i in 1:N){
      for(t in 1:T){
        g_it_gwh <- v[idx_gen(i,t)]; cap_gw <- X[i, t]
        avail_gwh <- (1 - scaled_mean_theta[i]) * capfac_vec[i] * cap_gw * 8760
        constr <- c(constr, 1/1000000 * (g_it_gwh - avail_gwh))
      }
    }
    for(i in 1:N) for(t in 1:(T+1)) constr <- c(constr, - X[i, t])
    constr
  }

  # initial guess and solve
  invest_init <- rep(0.6, N * T); gen_init <- rep(d / N, N * T)
  res <- nloptr(x0 = c(invest_init, gen_init), eval_f = objective_fn, lb = rep(0, nvars), ub = c(rep(20, N*T), rep(d, N*T)), eval_g_ineq = eval_g_ineq,
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-6, ftol_rel = 1e-8, maxeval = 100000))
  v_opt <- res$solution

  invest_opt <- matrix(v_opt[1:(N*T)], nrow = N, byrow = TRUE)
  gen_opt    <- matrix(v_opt[(N*T + 1):(2*N*T)], nrow = N, byrow = TRUE)
  X_opt <- compute_capacity_path(v_opt)

  det_cost <- sapply(1:T, function(t) sum(sapply(1:N, function(i){ inv_gw <- invest_opt[i,t]; gen_gwh <- gen_opt[i,t]; gen_mwh <- gen_gwh * gen_scale; I_scaled[i] * inv_gw^2 + c_vec[i] * gen_mwh + carbon_fee * emissions[i] * gen_mwh })))
  det_emis <- sapply(1:T, function(t) sum(sapply(1:N, function(i) emissions[i]*gen_opt[i,t]*gen_scale)))

  # Simulation
  sim_realized_gen <- array(0, dim = c(nsim, N, T))
  sim_cost_per_year <- matrix(0, nsim, T)
  sim_emis_year <- matrix(0, nsim, T)
  sim_total_prod <- matrix(0, nsim, T)
  sim_shortfall_flag <- matrix(0, nsim, T)
  sim_shortfall_amount <- matrix(0, nsim, T)

  for(sim in 1:nsim){
    theta_mat <- matrix(0, nrow=N, ncol=T)
    for(i in 1:N){ mu <- -gamma_vec[i]/delta_vec[i]; sigma <- 1/delta_vec[i]; theta_raw <- plogis(rnorm(T, mu, sigma)); theta_mat[i,] <- pmin(pmax(scaling_factors[i] * theta_raw, 0), 1) }
    for(t in 1:T){
      cap_start <- X_opt[, t]; total_y <- 0; emis_y <- 0
      for(i in 1:N){
        gen_gwh <- (1 - theta_mat[i, t]) * capfac_vec[i] * 8760 * cap_start[i]
        sim_realized_gen[sim, i, t] <- gen_gwh
        total_y <- total_y + gen_gwh
        emis_y <- emis_y + emissions[i] * gen_gwh
      }
      sim_total_prod[sim, t] <- total_y
      sim_shortfall_flag[sim, t] <- as.integer(total_y < d)
      sim_shortfall_amount[sim, t] <- max(d - total_y, 0)

      # cost: convert gen GWh -> MWh for c_vec and carbon_fee
      total_cost <- sum(sapply(1:N, function(i){ inv_gw <- invest_opt[i, t]; gen_gwh_sim <- sim_realized_gen[sim, i, t]; gen_mwh_sim <- gen_gwh_sim * gen_scale; I_scaled[i] * inv_gw^2 + c_vec[i] * gen_mwh_sim + carbon_fee * emissions[i] * gen_mwh_sim }))
      sim_cost_per_year[sim, t] <- total_cost; sim_emis_year[sim, t] <- emis_y
    }
  }

  mean_sim_cost <- colMeans(sim_cost_per_year)
  cost_p5 <- apply(sim_cost_per_year, 2, quantile, probs = 0.05)
  cost_p95 <- apply(sim_cost_per_year, 2, quantile, probs = 0.95)
  mean_sim_emis <- colMeans(sim_emis_year)
  prob_shortfall_by_year <- colMeans(sim_shortfall_flag)
  mean_shortfall_by_year <- colMeans(sim_shortfall_amount)
  nonzero_shortfall_mean <- sapply(1:T, function(tt){ vals <- sim_shortfall_amount[,tt]; nnz <- vals[vals>0]; if(length(nnz)>0) mean(nnz) else 0 })

  # mean per-tech generation matrix (N x T) in GWh
  mean_gen_by_tech_year <- matrix(NA_real_, nrow = N, ncol = T)
  for(i in 1:N) mean_gen_by_tech_year[i, ] <- colMeans(sim_realized_gen[, i, , drop = FALSE])
  rownames(mean_gen_by_tech_year) <- techs; colnames(mean_gen_by_tech_year) <- paste0("Year", 1:T)

  # utilization table
  max_possible_by_tech_year <- sapply(1:N, function(i) capfac_vec[i] * 8760 * X_opt[i, 1:T])
  mean_gen_by_tech_year_tmp <- mean_gen_by_tech_year
  share_used <- mean_gen_by_tech_year_tmp / max_possible_by_tech_year
  util_df <- data.frame(tech = rep(techs, each = T), year = rep(1:T, times = N), share = as.vector(t(share_used)))

  list(
    det_cost = det_cost, det_emis = det_emis,
    mean_sim_cost = mean_sim_cost, cost_p5 = cost_p5, cost_p95 = cost_p95,
    mean_sim_emis = mean_sim_emis, prob_shortfall = prob_shortfall_by_year,
    mean_shortfall = mean_shortfall_by_year, mean_shortfall_if = nonzero_shortfall_mean,
    util_df = util_df, cap_end = X_opt[, T+1], techs = techs,
    invest_opt = invest_opt, X_opt = X_opt, gen_opt = gen_opt,
    mean_theta_per_tech = scaled_mean_theta,
    mean_sim_total_prod = colMeans(sim_total_prod), sim_cost_per_year = sim_cost_per_year,
    sim_total_prod = sim_total_prod, scaling_factors = scaling_factors,
    mean_gen_by_tech_year = mean_gen_by_tech_year, sim_realized_gen = sim_realized_gen
  )
}

# Run scenarios
results <- list()
results[[SCEN_LABELS[1]]] <- run_scenario(CARBON_FEES[1], T, nsim, gamma_vec, delta_vec, xi_vec, lambda_vec, scaling_factors, capfac_vec, I_vec, c_vec, eta_vec, x0, emissions, d, techs, scaled_mean_theta)
results[[SCEN_LABELS[2]]] <- run_scenario(CARBON_FEES[2], T, nsim, gamma_vec, delta_vec, xi_vec, lambda_vec, scaling_factors, capfac_vec, I_vec, c_vec, eta_vec, x0, emissions, d, techs, scaled_mean_theta)

# Helper: safe color mapping
safe_color_map <- function(vec_levels, master_colors) {
  levels <- unique(as.character(vec_levels))
  cmap <- master_colors[names(master_colors) %in% levels]
  if(length(cmap) != length(levels)) cmap <- setNames(unname(master_colors)[seq_along(levels)], levels)
  cmap
}

# Example: per-technology variable cost decomposition + export
make_var_decomp_and_export <- function(results, out_dir = "outputs"){
  dir.create(out_dir, showWarnings = FALSE)
  dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
  gen_scale <- 1e3
  var_list <- list()
  for(sc in names(results)){
    r <- results[[sc]]
    mg <- as.matrix(r$mean_gen_by_tech_year) # N x T
    if(nrow(mg) != length(r$techs)) mg <- t(mg)
    mean_gen_mwh <- mg * gen_scale
    var_cost_mat <- sweep(mean_gen_mwh, 1, c_vec, `*`)
    emis_cost_mat <- sweep(mean_gen_mwh, 1, emissions, `*`) * (ifelse(sc==SCEN_LABELS[1], CARBON_FEES[1], CARBON_FEES[2]))
    df <- as.data.frame(t(var_cost_mat)); colnames(df) <- r$techs; df$Year <- 1:nrow(df)
    df_long <- pivot_longer(df, cols = all_of(r$techs), names_to = "Tech", values_to = "VarCost_USD")
    df_long$Scenario <- sc
    write.csv(df_long, file.path(out_dir, "tables", paste0("per_tech_varcost_long_", gsub(" ","_",sc), ".csv")), row.names = FALSE)
    var_list[[sc]] <- df_long
  }
  var_all <- bind_rows(var_list)
  p <- ggplot(var_all, aes(x = Year, y = VarCost_USD/1e6, fill = Tech)) + geom_area(position = "stack") + facet_wrap(~Scenario, nrow = 2) + theme_minimal()
  ggsave(file.path(out_dir, "plots", "per_tech_varcost_stack.png"), p, width=10, height=6)
  invisible(list(var_all = var_all, plot = p))
}

# Run export for variable decomposition
var_results <- make_var_decomp_and_export(results)

# Compact export: save key csv/rds and zip
compact_export_zip <- function(results, out_dir = "outputs"){
  dir.create(out_dir, showWarnings = FALSE)
  saveRDS(results, file = file.path(out_dir, "results.rds"))
  # save comp_df_all and others if present
  if(exists("comp_df_all")) write.csv(comp_df_all, file.path(out_dir, "comp_df_all.csv"), row.names = FALSE)
  if(exists("decomp_df")) write.csv(decomp_df, file.path(out_dir, "decomp_df.csv"), row.names = FALSE)
  zipfile <- file.path(out_dir, "exports.zip")
  files <- list.files(out_dir, full.names = TRUE)
  if(requireNamespace("zip", quietly=TRUE)) zip::zip(zipfile, files = files) else try(utils::zip(zipfile, files))
  zipfile
}

zipfile <- compact_export_zip(results)
message("Done. Zip: ", zipfile)
