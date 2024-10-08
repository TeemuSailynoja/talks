library(calculus)
library(SBC)

prior_predictive_generator_years <- function(n_years){
  # Some prior predictive samples produce negative populations.
  # We reject those, but out of curiosity I'll add a counter to keep
  # track of the number of rejected samples.
  valid_populations <- FALSE
  n_tries <- 0
  while (!valid_populations) {
    n_tries <- n_tries + 1
    # Draw parameter values from the prior:
    alpha <- qnorm(runif(1, pnorm(0, 1, .5), 1), 1, .5)
    beta <- qnorm(runif(1, pnorm(0, .05, .05), 1), .05, .05)
    gamma <- qnorm(runif(1, pnorm(0, 1, .5), 1), 1, .5)
    delta <- qnorm(runif(1, pnorm(0, .05, .05), 1), .05, .05)

    l0 <- rlnorm(1, log(10), 1)
    h0 <- rlnorm(1, log(10), 1)

    # Solve the ODE system to get populations for years 2-21.
    pops <- ode(f = function(h, l) {
      c((alpha - beta * l)*h, (-gamma + delta * h)*l)
    },
    var = c(h = h0, l = l0),
    times = seq(1,n_years, .1))
    if (all(pops > 0)) {
      valid_populations <- TRUE
    }
  }

  # With valid populations, draw observation noise from priors and
  # simulate observations.
  sigma_h <- rlnorm(1, -1, 1)
  sigma_l <- rlnorm(1, -1, 1)
  hare_pelts <- rlnorm(n_years, log(pops[seq(1,201,10),1]), sigma_h)
  lynx_pelts <- rlnorm(n_years, log(pops[seq(1,201,10),2]), sigma_l)

  # Return the ground truth as `variables` and observations as `generated`
  list(
    variables = list(
      theta = c(alpha, beta, gamma, delta),
      z_init = c(h0, l0),
      sigma = c(sigma_h, sigma_l),
      loglik = sum(dlnorm(hare_pelts, log(pops[seq(1,201,10), 1]), sigma_h, log = T)) +
        sum(dlnorm(lynx_pelts, log(pops[seq(1,201,10), 2]), sigma_l, log = T))
    ),
    generated = list(
      N = n_years - 1,
      ts = 2:n_years,
      y_init = c(hare_pelts[1], lynx_pelts[1]),
      y = matrix(c(hare_pelts[-1], lynx_pelts[-1]), ncol = 2),
      n_tries = n_tries
    )
  )
}

pp_sets <- function(n_sets, n_years = 21) {
  prior_predictive_generator <- SBC_generator_function(
    prior_predictive_generator_years, n_years = n_years)
  sets <- generate_datasets(
    prior_predictive_generator, n_sims = n_sets)
  1:n_sets |>
    lapply(\(id) {
      data.frame(rbind(
        sets$generated[[id]]$y_init,
        sets$generated[[id]]$y
      ))   |>
        dplyr::rename(Hare = X1, Lynx = X2) |>
        cbind(data.frame(gen_id = paste("Sim", sprintf("%02d", id)), Year = 1900:1920))
    })
}
