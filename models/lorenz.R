source("00_depenencies.R")

lorenz_mod <- odin2::odin({
  deriv(x) <-  sigma*(y - x)
  deriv(y) <-  x*(rho - z) - y
  deriv(z) <-  x*y - beta*z

  initial(x) <- x_0
  initial(y) <- y_0
  initial(z) <- z_0

  # user specified init vals
  x_0 <- parameter(10)
  y_0 <- parameter(10)
  z_0 <- parameter(10)

  # user specified params
  sigma <- parameter(10)
  rho <- parameter(28)
  beta <- parameter(8/3)
})
