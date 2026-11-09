# ---- Load data generating models -------
source("models/utils.R")
source("models/lorenz.R")

# ----- Simulate data ---------
## ----- Lorenz --------
lrz_init <- list(x_0 = 10, y_0 = 10, z_0 = 10)
lrz_params <- list(sigma = 10, rho = 28, beta = 8/3)
dat_lrz <- run_mod(lorenz_mod,
                       c(
                         lrz_init, # initial condition
                         lrz_params # params
                       ),
                       timestep = .05,
                       duration = 100)

dat_lrz_info <- list(
  dat = dat_lrz,
  params = lrz_params,
  print_func = \(pars) {
    c(
      sprintf("dx/dt = %.4g(y - x)", pars$sigma),
      sprintf("dy/dt = x(%.4g - z) - y", pars$rho),
      sprintf("dz/dt = xy - %.4gz", pars$beta)
    )
  }
)
write_rds(dat_lrz_info, "data/lrz_dat.rds")
plot_time_series(dat_lrz_info$dat)
plot_phase_space(dat_lrz_info$dat)

