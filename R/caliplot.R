#' @title Calibrate sensitivity parameters based on pre-treatment variables
#'
#' @description create a plot....
#' @param x a \code{tibble} or data frame of observed pre-treatment covariates
#' @param trt a vector of binary treatment indicators
#' @param y a vector of outcomes
#' @param gamma_seq a vector of values for sensitivity parameter to be investiaged
#' @return \code{caliplot} returns a plot comparing the effect of outcome with the effect of the most significant
#'         four observed explanatory variables in the assumed treatment model for sensitivity analysis.
#' @export
#'
#' @examples
#' # load data #
#' library(foreign)
#' NHANES <- read.dta(system.file("extdata", "NHANES3hbp_dbp.dta",
#'           package = "TukeySensitivity")) %>% as_tibble
#' NHANES <- NHANES %>% dplyr::filter(ave_dbp > 30)
#' x_train_joint <- NHANES %>% select(-one_of(c("ave_dbp", "ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))
#'
#' # Observed Dependent Variables #
#' x = x_train_joint %>% select(-one_of("trt_dbp"))
#'
#' # Treatment #
#' trt = x_train_joint %>% select(trt_dbp)
#'
#' # Outcomes #
#' y_train_joint = NHANES %>% select(ave_dbp)
#'
#' # Sensitivity Parameter Sequence #
#' gamma = seq(0.01, 0.05, by=0.001)
#'
#' # plot #
#' caliplot(x,trt,y_train_joint,gamma)



caliplot <- function( x, trt, y, gamma_seq)
{
  # scale the observed Covariates #
  x_scaled = x %>% mutate_all(funs(scale(.)))
  names(trt)="trt"

  # estimated m(x) #
  propensity_fit = glm(trt ~., data=cbind(x_scaled,trt), family=binomial(logit))
  lp <- propensity_fit$linear.predictors

  # rho_{X_j|{X_{-j}}
  max_change <- sapply(colnames(x), function(nm) {
    glm_fit <- glm(as.formula(paste0("trt ~ .-", nm)), data=cbind(x,trt), family=binomial(logit))
    lp2 <- glm_fit$linear.predictors
    Rw_full <- var(lp) / (var(lp) + pi^2/3)
    Rw_sub <-  var(lp2) / (var(lp2) + pi^2/3)
    (Rw_full - Rw_sub)/(1-Rw_sub)
  })

  max_change=sort(max_change, decreasing=TRUE)

  # rho_{Y|X} #
  compute_rho <- function(gamma.tilde) ## the imput "gamma.tilde" of this function should be the coefficient of scaled Y
  {                                    ## gamma.tilde = gamma*sigma(Y)
    reduced <- (pi^2/3) / (var(lp) + pi^2/3)  ## 1-rho_x^2
    full <- (pi^2/3) / (var(lp) + gamma.tilde^2 + pi^2/3) ## 1-rho_{x,y}^2
    (reduced - full) / reduced ## eq 40
  }

  Index_covariates = function(max_change)
  {
    if (length(max_change) < 4)
      index = names(max_change)
    else
      index = names(max_change)[1:4]
    return(index)
  }

  # plot #
  toplot <- max_change[Index_covariates(max_change)]
  tibble(gamma=gamma_seq, rho2=sapply(gamma_seq*sd(as.matrix(y_train_joint)), compute_rho)) %>%
    ggplot2::ggplot() + geom_line(aes(x=gamma, y=rho2), size=1.5) +
    geom_hline(data=tibble(hline=toplot, Predictor=Index_covariates(max_change)),
               aes(yintercept=toplot, col=Predictor), linetype="dashed", size=1.2) +
    xlab(expression(gamma[t])) + ylab(substitute(rho[m]^2, list(m="Y|X"))) +
    guides(color=guide_legend(title=substitute(rho[m]^2, list(m=bquote(X[j]~"|"~X[-j]))),
                              title.hjust=0.5, legend.title=element_text(size=14))) +
    theme(axis.title=element_text(size=20), axis.text=element_text(size=16))

}


