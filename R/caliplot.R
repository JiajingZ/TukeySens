#' @title Calibrate sensitivity parameters based on pre-treatment variables
#'
#' @description Calibrate the magnitude of sensitivity parameters to the amount of variation in the treatment 
#'              assignment \code{T} that is explained by outcome \code{Y(t)}, above and beyond what is accounted
#'              for by observed pre-treatment variables \code{X}.
#' @param x a \code{tibble} or data frame with observed pre-treatment variables
#' @param trt a vector with binary treatment indicators
#' @param y a vector with outcomes
#' @param gamma_seq a vector with chosen values for sensitivity parameter \eqn{\gamma_t} 
#' @section Details: \code{caliplot} returns a plot of sensitivity parameter \eqn{\gamma_t} vs partial coefficient of 
#'          variation from outcome \code{Y}, \eqn{\rho^2_{Y|X}}. For comparison, the largest four partial coefficients 
#'          of variation from covariates are also plotted if the number of observed covariates is larger or equal to four;
#'          otherwise, all the partial coefficients of variation from covariates are plotted.
#' @export
#'
#' @examples
#' # Observed Pre-treatment Variables 
#' x = NHANES %>% select(-one_of("trt_dbp", "ave_dbp"))
#'
#' # Treatment 
#' trt = NHANES %>% select(trt_dbp)
#'
#' # Outcomes 
#' y = NHANES %>% select(ave_dbp)
#'
#' # Sensitivity Parameter Sequence 
#' gamma = seq(0.01, 0.1, by = 0.001)
#'
#' # plot 
#' caliplot(x, trt, y, gamma)


caliplot <- function( x, trt, y, gamma_seq){
  # scale the observed Covariates #
  x_scaled = x %>% mutate_all(funs(scale(.)))
  names(trt) = "trt"

  # estimated m(x) #
  propensity_fit = glm(trt ~., data = cbind(x_scaled,trt), family = binomial(logit))
  lp <- propensity_fit$linear.predictors

  # rho_{X_j|{X_{-j}}
  max_change <- sapply(colnames(x), function(nm){
    glm_fit <- glm(as.formula(paste0("trt ~ .-", nm)), data = cbind(x,trt), family = binomial(logit))
    lp2 <- glm_fit$linear.predictors
    Rw_full <- var(lp) / (var(lp) + pi^2/3)
    Rw_sub <-  var(lp2) / (var(lp2) + pi^2/3)
    (Rw_full - Rw_sub) / (1 - Rw_sub)
  })

  max_change = sort(max_change, decreasing = TRUE)

  # rho_{Y|X} #
  ## the imput "gamma.tilde" of this function should be the coefficient of scaled Y
  ## gamma.tilde = gamma*sigma(Y)
  compute_rho <- function(gamma.tilde){                                    
    reduced <- (pi^2/3) / (var(lp) + pi^2/3)  ## 1-rho_x^2
    full <- (pi^2/3) / (var(lp) + gamma.tilde^2 + pi^2/3) ## 1-rho_{x,y}^2
    (reduced - full) / reduced ## eq 40
  }

  Index_covariates <- function(max_change){
    if (length(max_change) < 4)
      index = names(max_change)
    else
      index = names(max_change)[1:4]
    return(index)
  }

  # plot #
  toplot <- max_change[Index_covariates(max_change)]
  tibble(gamma = gamma_seq, rho2 = sapply(gamma_seq*sd(as.matrix(y)), compute_rho)) %>%
    ggplot2::ggplot() + geom_line(aes(x = gamma, y=rho2), size = 1.5) +
    geom_hline(data = tibble(hline = toplot, Predictor = Index_covariates(max_change)),
               aes(yintercept = toplot, col=Predictor), linetype="dashed", size = 1.2) +
    xlab(expression(gamma[t])) + ylab(substitute(rho[m]^2, list(m = "Y|X"))) +
    guides(color = guide_legend(title = substitute(rho[m]^2, list(m = bquote(X[j]~"|"~X[-j]))),
                              title.hjust = 0.5, legend.title = element_text(size = 14))) +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16))
}


