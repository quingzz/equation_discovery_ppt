source("00_depenencies.R")
source("models/utils.R")
source("ed_methods/utils.R") # load utils for equation discovery methods
source("ed_methods/GPoM.R")


# ====== Load data =========
lrz_dat <- read_rds("data/lrz_dat.rds")

# ======= Lorenz system =========
# try obtaining the model fr Lorenz system data naively
n_var <- 3 # number of variables
n_degree <- 2 # max degree for the polynomials
pmin <- 6 # min number of terms for the derived system
pmax <- 11 # max number of terms for the derived system

# compute the number of equations to be checked at each iteration
models_to_check(n_eq = n_var, n_terms=d2pMax(n_var, n_degree), n_pmin=pmin, n_pmax=pmax)

## Interpolate data ---------
# Original stepsize = .05 -> interpolate such that stepsize= .01
lrz_dat_interp <- lrz_dat

# get the variables
interp_vars <- c("x", "y", "z")
# get the time
t_in <- lrz_dat_interp$dat$t
t_new <- seq(min(t_in), max(t_in), by = .001)

# interpolate data
input_dat_interp <- sapply(interp_vars, \(var) {
  spline(t_in, lrz_dat_interp$dat[, var], xout = t_new, method = "fmm")$y
})
colnames(input_dat_interp) <- interp_vars
lrz_dat_interp$dat <- data.frame(t = t_new, input_dat_interp)

## Generate candidate models ---------
gpom_interp_out <- gPoMo(lrz_dat_interp$dat[,c("x", "y", "z")],
                  tin = lrz_dat_interp$dat[, "t"],
                  dMax = n_degree, nS = c(1,1,1), show = 1,
                  nPmin = 6, nPmax = 11, verbose = TRUE, IstepMin = 500, IstepMax = 5000)
gpom_interp_out_df <- gpom_out_to_df(gpom_interp_out, lrz_dat_interp,
                              n_var = n_var, n_degree = n_degree,
                              varnames = c("x", "y", "z"))

### cache output -------
# write_rds(gpom_interp_out_df, "outputs/lrz_w_interp_gpom.rds")

### load output for reruns -------
# gpom_interp_out_df <- read_rds("outputs/lrz_w_interp_gpom.rds")

# Figures for the ppt =====
## Candidate models ---------
gpom_interp_plots <- gpom_out_plots(gpom_interp_out_df, lrz_dat_interp)
gpom_interp_plots_xz <- gpom_out_plots(gpom_interp_out_df, lrz_dat_interp, vars = c("x", "z"))

dim(gpom_interp_out_df)

# collect plots into panels of multiple plots
gpom_interp_panels <- gpom_out_panels(gpom_interp_out_df, gpom_interp_plots,
                               max_plots_per_panel = 15,
                               exclude_class = 0, # exclude diverging models
                               ncol = 5, plot_ratio = 6/12,
                               axes_fontsize = 6,
                               title_fontsize = 8,
                               subtitle_fontsize = 7,
                               caption_fontsize = 5.5)
# save panels
iwalk(gpom_interp_panels, function(plt, i) {
  ggsave(
    filename = file.path("figs/GPoM/", sprintf("lrz_w_interp_%03d.png", as.integer(i))),
    plot = plt, width = 12, height = 6, dpi = 200
  )
})

## Good fits -----------
good_mods_id <- c(5, 15, 18, 19, 30, 33, 37, 38, 39)
good_fit_panels <- gpom_out_panels(
                gpom_interp_out_df[good_mods_id,], gpom_interp_plots[good_mods_id],
                max_plots_per_panel = 15,
                exclude_class = 0, # exclude diverging models
                ncol = 5, plot_ratio = 6/12,
                axes_fontsize = 8,
                title_fontsize = 9,
                subtitle_fontsize = 8,
                caption_fontsize = 6.25)
iwalk(good_fit_panels, function(plt, i) {
  ggsave(
    filename = file.path("figs/GPoM/", sprintf("good_fits_%03d.png", as.integer(i))),
    plot = plt, width = 12, height = 6, dpi = 200
  )
})

good_fits_panels_xz <- gpom_out_panels(
  gpom_interp_out_df[good_mods_id,], gpom_interp_plots_xz[good_mods_id],
  max_plots_per_panel = 10,
  ncol = 5, plot_ratio = 6/12,
  axes_fontsize = 8,
  title_fontsize = 9,
  subtitle_fontsize = 8,
  caption_fontsize = 6.25)

iwalk(good_fits_panels_xz, function(plt, i) {
  ggsave(
    filename = file.path("figs/GPoM/", sprintf("good_fits_xz_%03d.png", as.integer(i))),
    plot = plt, width = 12, height = 6, dpi = 200
  )
})

## Candidate models for 6 terms -----------
candidate_6_terms <- gpom_interp_out_df %>% mutate(index = 1:n()) %>%  filter(n_terms == 6) %>% pull(index)
candidate_panels <- gpom_out_panels(
  gpom_interp_out_df[candidate_6_terms,], gpom_interp_plots[candidate_6_terms],
  max_plots_per_panel = 10,
  ncol = 5, plot_ratio = 6/12,
  axes_fontsize = 8,
  title_fontsize = 9,
  subtitle_fontsize = 8,
  caption_fontsize = 6.25)
iwalk(candidate_panels, function(plt, i) {
  ggsave(
    filename = file.path("figs/GPoM/", sprintf("candidates_6terms_%03d.png", as.integer(i))),
    plot = plt, width = 12, height = 6, dpi = 200
  )
})

## Rejected models ---------
gpom_interp_out_df %>%
  filter(okMod == 0, n_terms == 6)

rej_example <- gpom_interp_out_df %>%
  filter(id == 10)
rej_example$equation %>% message()


