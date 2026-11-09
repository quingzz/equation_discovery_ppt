library(GPoM)
library(patchwork)
library(future)
library(furrr)

# Helper function to compute no. of models to check at every iteration
# n_eq number of equations
# n_terms maximum number of terms
# n_pmin min number of parameters for the model
# n_pmax max number of parameters for the model
models_to_check <- function(n_eq=3, n_terms=10, n_pmin=1, n_pmax=14){
  # compute expected number of models to check every runs
  df <- do.call(expand_grid, lapply(1:n_eq, function(i) 1:n_terms))
  total <- rowSums(df)
  nrow(df[total >= n_pmin & total <= n_pmax, ])
}

# Helper function to generate an output data.frame which is easier to work w/
gpom_out_to_df <- function(gpom_out, input_dat, n_var = 3, n_degree = 2,
                           varnames = c("x", "y", "z"), workers=NULL, mod_type=NULL){

  # set up for parallelisation
  n_models <- length(gpom_out$okMod)
  if (is.null(workers)) {
    workers <- max(1, availableCores() - 1)
  }
  workers <- min(workers, n_models)
  old_plan <- plan(multisession, workers = workers)
  on.exit(plan(old_plan), add = TRUE)
  opts <- furrr_options(packages = "GPoM", seed = TRUE)


  out <- tibble(
    id    = seq_along(gpom_out$okMod),
    mat   = map(id, ~ gpom_out$models[[paste0("mToTest", .x)]]),
    okMod = gpom_out$okMod,
    # okMod_label = okMod_labels[as.character(okMod)],
    n_terms  = map_int(mat, ~ sum(.x != 0)),
    # print out equations
    equation = map_chr(mat, ~ suppressMessages(
      visuEq(.x, nVar = n_var, dMax = n_degree, substit = varnames, approx=TRUE)
    ) %>%
      unlist() %>% paste(collapse = "\n")
    )
  )

  if(!is.null(mod_type)){
    out <- out %>%
      filter(okMod == mod_type)
  }

  out %>%
    mutate(
      # reconstruct time series with current derived equation
      out_df = future_map(mat, \(KL){
        out <- numicano(
          n_var, n_degree, KL=KL,
          onestep = input_dat$dat[2, "t"] - input_dat$dat[1, "t"],
          Istep = nrow(input_dat$dat),
          v0 = as.numeric(input_dat$dat[1,varnames]),
          method = "rk4"
        )

        colnames(out$reconstr) <- c("time", varnames)
        out$reconstr
      }, .options = opts)
    )
}

# Helper function to generate a list of plots for output models
gpom_out_plots <- function(gpom_out_df, input_dat, vars = c("x", "y")){
  pmap(
    list(id = gpom_out_df$id, okMod = gpom_out_df$okMod, n_terms = gpom_out_df$n_terms,
         eq = gpom_out_df$equation, sim_df = gpom_out_df$out_df),
    \(id, okMod, n_terms, eq, sim_df) {
      ggplot()+
        geom_path(aes(x = .data[[vars[1]]], y = .data[[vars[2]]], color="input"), data = input_dat$dat) +
        geom_path(aes(x = .data[[vars[1]]], y = .data[[vars[2]]], color="reconstructed"), data = sim_df) +
        scale_color_manual(values = c(
          "reconstructed" = "cornflowerblue",
          "input" = "gray"
        )) +
        labs(
          title = paste0("Model ", id),
          subtitle = paste0("classification: ", okMod, " no. terms: ", n_terms),
          caption = eq,
          color = NULL
        ) +
        theme(plot.caption = element_text(hjust = 0))
    })
}

# Helper function to collect model plots into panels of multiple plots
gpom_out_panels <- function(gpom_out_df, gpom_plots,
                          max_plots_per_panel = 20,
                          exclude_class = NULL, # filter out models that match the classification
                          ncol = 4,
                          plot_ratio = 6/12,
                          axes_fontsize = 6,
                          title_fontsize = 7,
                          subtitle_fontsize = 6,
                          caption_fontsize = 5.5){
  filt_mods <- if(is.null(exclude_class)) 1:nrow(gpom_out_df) else
    which(gpom_out_df$okMod != exclude_class)

  gpom_plots <- gpom_plots[filt_mods]

  gpom_plots %>%
    map(\(p) p + theme(aspect.ratio = plot_ratio)) %>%
    split(ceiling(seq_along(gpom_plots) / max_plots_per_panel)) %>%
    map(\(chunk) {
      wrap_plots(chunk, ncol = ncol, guides = "collect") &
        theme(
          plot.title = element_text(size = title_fontsize),
          plot.subtitle = element_text(size = subtitle_fontsize),
          plot.caption = element_text(hjust = 0, size = caption_fontsize),
          axis.title = element_text(size = axes_fontsize),
          axis.text = element_text(size = axes_fontsize),
        )
    })
}

