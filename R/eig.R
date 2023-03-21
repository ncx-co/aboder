#' Batch Nested Monte Carlo Expected Information Gain
#' 
#' Calculate EIG for a batch of designs.
#' 
#' @param designs a list of design dataframes
#' @param model a fitted brms object
#' @param N int, how many outer loop Monte Carlo samples
#' @param M int, how many inner loop Monte Carlo samples
#' @param cores int, how many cores and chains to run MCMC with
#' @param update_model bool, update model with new samples
#' @param ... optional parameters passed to `update_model`
#' 
#' @return expected information gain (float)
#' @export 
batch_nmc_eig <- function(designs, model, N = 1000, M = 1000, cores = 1, update_model = TRUE, ...){ 
  if(update_model) model <- update_model(model, N, M, cores, ...)
  eigs <- furrr::future_map_dbl(
    designs,
    ~ nmc_eig(
        design = .x,
        model = model,
        N = N,
        M = M,
        cores = cores,
        update_model = FALSE
    ),
    .options = furrr::furrr_options(seed = TRUE)
  )
  eigs
}

#' Nested Monte Carlo Expected Information Gain
#' 
#' Esimate the expected information gain (EIG) for a given design (data frame with potential independent variables to collect)
#' 
#' #' Uses nested Monte Carlo estimation of the EIG. 
#' See equation 4 in https://arxiv.org/pdf/1903.05480.pdf or
#' equation 9 & 10 in https://arxiv.org/pdf/1108.4146.pdf
#' 
#' @param design a dataframe of independent variable values
#' @param model a fitted brms object
#' @param N int, how many outer loop Monte Carlo samples
#' @param M int, how many inner loop Monte Carlo samples
#' @param cores int, how many cores and chains to run MCMC with
#' @param update_model bool, update model with new samples
#' @param ... optional parameters passed to `update_model`
#' 
#' @return expected information gain of design
#' @export
nmc_eig <- function(design, model, N = 1000, M = 1000, cores = 1, update_model = TRUE, ...){
  if(update_model) model <- update_model(model, N, M, cores)
  # predict N new y for new x value
  # we are treating each y_i as if it were observed in the next step of the experiment
  P <- nrow(design)
  design <- predict_design(design, model, N)

  # generate likelihoods with all draws
  all_lik <- brms::log_lik(model, newdata = design)
  # extract likelihood from parameters that produced the y values
  # these are the numerator values in eq. 4
  # p(y| theta, d), likelihood of Y with theta that produced value
  lp_N <- purrr::map_dbl(
    1:N,
    ~{
      sum(all_lik[.x, seq(.x, N * P, by = N)])
    }
  )
  lq_N <- inner_monte_carlo(all_lik, N = N, P = P) # take monte carlo expectation of p(y| d), inner monte carlo, likelihood of Y with prior distribution of thetas

  eig_N <- lp_N - lq_N

  # outer monte carlo expectation of KL divergence
  mean(eig_N)
}

update_model <- function(model, N, M, cores = 1, recompile = FALSE, ...){
  # updates model with proper number of posterior draws for use in nmc_eig
  NM <- N + M # N samples for outer expectation, M for inner expectation
  
  # need to update the model to produce more draws for the Monte Carlo estimates
  if(model$algorithm == "sampling"){
    # fit with mcmc
    model <- stats::update(
      model, 
      iter = as.integer(ceiling(NM * 2 / cores)),  # double iter for warmup, divide by number of cores/chains
      chains = cores, 
      cores = cores, 
      recompile = recompile,
      ...
    )
  } else{
    # fit with vb
    model <- stats::update(model, output_samples = NM, recompile = recompile, ...)
  }
  model
}

predict_design <- function(newdata, model, N){
  # add predicted draws to dataframe
  assertthat::assert_that(!"y_" %in% colnames(newdata), msg="Column `y_` already exists and cannot overwrite.")
  response <- all.vars(model$formula$formula)[[1]]
  y_N <- brms::posterior_predict(model, newdata = newdata)[1:N, , drop = FALSE] # don't drop to vector if only one column
  newdata |>
    dplyr::mutate(
      y_ = split(t(y_N), 1:ncol(y_N))
    ) |>
    tidyr::unnest(y_) |> # there is only one x value so we need to unnest to expand the design variables
    dplyr::mutate(!!response := y_) |>
    dplyr::select(-y_)
}

inner_monte_carlo <- function(log_likelihoods_NM, N, P){
  # input matrix of log likelihoods M draws x N points
  # return the log expected likelihood log(E(p(y)))
  # foreach column calc log(mean(exp(log_likelihood)))
  # to get log expected likelihood
  purrr::map_dbl(
    1:N,
    ~ logsumexp(rowSums(log_likelihoods_NM[, seq(.x, N * P, by = N), drop = FALSE])) - log(nrow(log_likelihoods_NM))
  )
}

logsumexp <- function(x){
  x_ <- max(x)
  x_ + log(sum(exp(x - x_)))
}
