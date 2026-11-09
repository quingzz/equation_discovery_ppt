source("00_depenencies.R")

# ====== Generate library of terms ==========
# generate library of polynomial terms from degree 2 to max_degree
build_polynomial_lib <- function(dat, max_deg=2, const=FALSE, vars=NULL){
  vars <- if(!is.null(vars)) vars else colnames(dat)[colnames(dat) != "t"]

  # generate library of terms from 2:max_deg
  term_list <- map(2:max_deg, \(deg) {
    terms <- combn(vars, deg)
    # add x^deg terms
    terms %>%
      cbind(rep(vars, deg) %>% matrix(nrow = deg, byrow = TRUE))
  }) %>%
    map(~ asplit(.x, MARGIN=2)) %>% # convert from list of matrix of terms (at each degree) to list of terms
    flatten()

  # new columns (i.e. polynomial terms) for input data
  map(term_list, \(terms){
    # function to generate the name for the column
    # e.g. c("x", "y", "z") -> "xyz" and c("x", "y", "y") -> c("x", "y2")
    gen_term <- function(terms) {
      tab <- table(terms)
      paste0(names(tab), ifelse(tab == 1, "", tab), collapse = "")
    }

    # define name of the current term/column
    term_name <- gen_term(terms)
    setNames(list(reduce(dat[terms], `*`)), term_name) %>% as.data.frame(check.names = FALSE)
  }) %>%
    c(dat %>%
        select(all_of(c("t", vars))), .) %>%
    bind_cols()
}

# ======== Differentiate data ==============
compute_deriv <- function(dat, vars=NULL, nDeriv = 1){
  vars <- if(!is.null(vars)) c("t", vars) else colnames(dat)
  dat %>%
    select(all_of(vars)) %>%
    pivot_longer(-t, names_to = "var", values_to = "value") %>%
    group_by(var) %>%
    reframe({
      # use GPoM::drvSucc for derivative here
      out <- GPoM::drvSucc(t, value, nDeriv = nDeriv)
      as_tibble(out$seriesDeriv[, 2:(nDeriv + 1), drop = FALSE],
                .name_repair = ~ paste0("deriv", 1:nDeriv)) %>%
        mutate(new_t = out$tout, .before = 1)
    }) %>%
    pivot_longer(starts_with("deriv"), names_to = "order", values_to = "value") %>%
    unite("var", var, order) %>%
    pivot_wider(names_from = "var", values_from = "value")
}


# ======== Custom Integrator ===============
