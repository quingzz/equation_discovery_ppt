run_mod <- function(mod, pars, duration=100, timestep=0.05){
  # --- initialize simulation time ----
  times <- seq(0, duration, timestep)

  sys <- dust_system_create(mod, pars)
  dust_system_set_state_initial(sys)
  out <- dust_system_simulate(sys, times)
  out <- dust_unpack_state(sys, out)

  out %>%
    as.data.frame() %>%
    mutate(
      t = times
    )
}


plot_time_series <- function(dat,
                             linewidth=1, color = "cornflowerblue",
                             ncol=NA, nrow=NA, scales = "fixed", ...){
  ncol <- if (!is.na(ncol)) ncol else 1
  nrow <- if (!is.na(nrow)) nrow else ncol(dat) - 1

  dat %>%
    pivot_longer(
      cols = -"t",
      names_to = "state",
      values_to = "pop"
    ) %>%
    ggplot(aes(t, pop)) +
      geom_line(linewidth = linewidth, color = color) +
      facet_wrap(~state, scales=scales, nrow=nrow, ncol=ncol)
}

plot_phase_space <- function(dat, line = TRUE, size = 0.8, dims = NA, ...){

  dims <- if(all(is.na(dims))) colnames(dat)[colnames(dat) != "t"] else dims

  if(length(dims) == 2){
    p <- dat %>%
      ggplot(aes(x = .data[[dims[1]]], y = .data[[dims[2]]], colour = t))

    if(line){
      p <- p + geom_path(linewidth = size, ...)
    } else {
      p <- p + geom_point(size = size, ...)
    }
    p + labs(x = dims[1], y = dims[2])

  } else if(length(dims) == 3){
    p <- dat %>%
      plot_ly(
        x = ~.data[[dims[1]]],
        y = ~.data[[dims[2]]],
        z = ~.data[[dims[3]]],
        color = ~t
      )

    if(line){
      p <- p %>% add_paths(line = list(width = size * 4), ...)
    } else {
      p <- p %>% add_markers(marker = list(size = size), ...)
    }
    p %>% layout(
      scene = list(
        xaxis = list(title = dims[1]),
        yaxis = list(title = dims[2]),
        zaxis = list(title = dims[3])
      )
    )

  } else {
    stop("Require 2 or 3 variables for phase space plot")
  }
}
