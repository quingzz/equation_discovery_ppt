source("00_depenencies.R")
source("models/utils.R")
source("ed_methods/SINDy.R") # load SINDy alg
source("ed_methods/utils.R") # load utils for equation discovery methods

# load simulated datasets
lrz_dat <- read_rds("data/lrz_dat.rds")

# ========= High level SINDy process =======
sindy_steps <- function(input_ls, max_deg = 2, nDeriv = 1,
                        iter=200, threshold=.15, ...){

  # Generate data to fit -----------
  # first build the library of terms
  terms_lib <- input_ls$dat %>% build_polynomial_lib(max_deg=max_deg, ...)
  # then estimate the derivative
  dat_deriv <- input_ls$dat %>% compute_deriv(nDeriv=nDeriv, ...)
  # get the data to fit the model
  fit_dat <- terms_lib %>%
    right_join(
      dat_deriv, by = join_by(t == new_t)
    ) %>%
    relocate(t, .before=1)

  # get the indices for terms
  nterms <- 2:ncol(terms_lib)
  message(paste0("Terms used: ", paste0(colnames(fit_dat)[nterms], collapse = ", ")))

  # Derive the equation for each var-----------
  # apply sindy alg
  dat_deriv %>%
    select(-new_t) %>%
    colnames() %>%
    map_dfr(\(var){
      out <- sindy_alg(fit_dat[, nterms], fit_dat[,var], iter=iter, threshold = threshold)
      tibble(
        variable = var,
        coefs = list(out$coefs),
        formula = out$formula
      )
    }) %>%
    bind_cols(
      tibble(true_formula = input_ls$print_func(input_ls$params))
    )
}


# ------- Lorenz -------
## Interpolate data ------
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

## Derive model ---------
interp_out <- sindy_steps(lrz_dat_interp, max_deg = 2, threshold = .1, vars=c("x", "y", "z"))

## Cache model -------
write_rds(interp_out, "outputs/lrz_w_interp_sindy.rds")
# interp_out <- read_rds("outputs/lrz_w_interp_sindy.rds")

## Reconstruct model -------
K9 <- interp_out$coefs %>% bind_cols() %>% as.matrix()  # x y z xy xz yz x2 y2 z2
K <- matrix(0, nrow = 10, ncol = 3)
K[c(7, 4, 2, 9, 8, 5, 10, 6, 3), ] <- K9

reconstr_lrz <- GPoM::numicano(
    3, 2,
    KL=K,
    onestep = lrz_dat_interp$dat[2, "t"] - lrz_dat_interp$dat[1, "t"],
    Istep = nrow(lrz_dat_interp$dat),
    v0 = as.numeric(lrz_dat_interp$dat[1,c("x", "y", "z")]),
    method = "rk4"
  )
reconstr_lrz <- reconstr_lrz$reconstr
colnames(reconstr_lrz) <- c("t", "x", "y", "z")
plot_time_series(as.data.frame(reconstr_lrz))



# Figures for ppt =========
## Time series -------
lrz_dat_interp$dat %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50) %>%
  ggplot() +
    geom_line(aes(x = t, y = pop, color = state)) +
    facet_wrap(~state, ncol = 1, scales = "free_y") +
    guides(color = "none") +
    labs(y = "") +
    theme_bw() +
    scale_color_manual(
      values = c("x" = "#F8766D", y = "#00BA38", z = "#619CFF")
    )
ggsave("figs/equation_discovery_ppt/lrz_time_series.png",
       width = 8, height = 5, dpi = 300)

## Phase space --------
lrz_dat_interp$dat %>%
  plot_phase_space()

## Derivative ---------
lrz_deriv <- lrz_dat_interp$dat %>% compute_deriv(nDeriv=1)
lrz_deriv %>%
  rename(
    t = new_t,
    `dx/dt` = x_deriv1,
    `dy/dt` = y_deriv1,
    `dz/dt` = z_deriv1
  ) %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50) %>%
  ggplot() +
  geom_line(aes(x = t, y = pop, color = state)) +
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  guides(color = "none") +
  labs(y = "") +
  theme_bw() +
  scale_color_manual(
    values = c("dx/dt" = "#F8766D", "dy/dt" = "#00BA38", "dz/dt" = "#619CFF")
  )
ggsave("figs/equation_discovery_ppt/lrz_deriv.png",
       width = 8, height = 5, dpi = 300)

## Deriv simplified ----------
lrz_deriv %>%
  rename(
    t = new_t,
    `dx/dt` = x_deriv1,
    `dy/dt` = y_deriv1,
    `dz/dt` = z_deriv1
  ) %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50) %>%
  ggplot() +
  geom_path(aes(x = pop, y = t, color = state)) +
  facet_wrap(~state, ncol = 3, scales = "free_y") +
  guides(color = "none") +
  labs(y = "", x = "") +
  theme_bw() +
  scale_color_manual(
    values = c("dx/dt" = "#F8766D", "dy/dt" = "#00BA38", "dz/dt" = "#619CFF")
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank()
  )

ggsave("figs/equation_discovery_ppt/lrz_deriv_verticle.png",
       width = 5, height = 10, dpi = 300)

## Libraries ----------
lrz_libs <- lrz_dat$dat %>% build_polynomial_lib(max_deg=2)
lrz_libs %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50) %>%
  ggplot() +
  geom_line(aes(x = t, y = pop), color = "blueviolet") +
  facet_wrap(~state, ncol = 2, scales = "free_y") +
  guides(color = "none") +
  labs(y = "") +
  theme_bw()

ggsave("figs/equation_discovery_ppt/lrz_libs.png",
       width = 6, height = 8, dpi = 300)

## Reconstructed ----------
lrz_reconstr_plot_dat <- as.data.frame(reconstr_lrz) %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50)

lrz_plot_dat <- lrz_dat_interp$dat %>%
  pivot_longer(cols = -"t", values_to = "pop", names_to = "state") %>%
  filter(t < 50)

## Reconstructed phase space ----------
ggplot() +
  geom_path(aes(x = x, y = y, color = "derived system"), data = as.data.frame(reconstr_lrz) %>% filter(t < 50)) +
  geom_path(aes(x = x, y = y, color = "true system"), linetype = "dashed", data = lrz_dat_interp$dat %>% filter(t < 50)) +
  scale_color_manual(
    values = c(
      "derived system" = "cornflowerblue", "true system" = "red"
    )
  ) +
  theme_bw()
ggsave("figs/equation_discovery_ppt/lrz_reconstructed_phase_space.png", width = 9, height = 5, dpi = 300)

## Reconstructed time series ----------
ggplot() +
  geom_line(aes(x = t, y = pop, color = "derived system"), data = lrz_reconstr_plot_dat) +
  geom_line(aes(x = t, y = pop, color = "true system"), linetype = "dashed", data = lrz_plot_dat) +
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c(
      "derived system" = "cornflowerblue", "true system" = "red"
    )
  ) +
  labs(y = "") +
  theme_bw()
ggsave("figs/equation_discovery_ppt/lrz_reconstructed_time_series.png", width = 10, height = 5, dpi = 300)


