#' @title Visualize ATE Estimates by Heatmap
#'
#' @description Visualize Average Treatment Effect (ATE) estimates for a grid of sensitivity parameters. "NS" is used to denote 
#'              "not significant", meaning that the 95\% posterior credible interval of the ATE contains 0.
#' @usage heatmap_ate(tukeyFit)
#' 
#' @param tukeyFit an object of class \code{tukeyFit} as returned by \code{fit_outcome}
#' 
#' @section Details: 
#'          The Average Treatment Effect is defined as:
#'          \deqn{\tau^{ATE} := E[Y(1) - Y(0)] = E[Y(1)] - E(Y(0))}
#'          For each t, the complete-data distribution for each potential outcome can be written as a mixture of 
#'          the distribution of observed and missing outcomes:
#'          \deqn{f(Y(t) \mid X) = f(T = t \mid X)f_t^{obs}(Y(t) \mid T = t, X) + 
#'          f(T = 1-t \mid X)f_t^{mis}(Y(t) \mid T = 1 - t, X).}
#'          The logistic selection with mixtures of exponential families (logistic-mEF models) have been considered, 
#'          the marginal selection functions in each arm are specified as logistic in the potential outcomes, and the 
#'          observed data is modeled with a mixture of exponential family distributions, which can be identified using
#'          flexible nonparametric or machine learning method. 
#'          Specifically, the treatment assignment model is posited as:
#'          \deqn{f(T=1 \mid Y(t),X) = \text{logit}^{-1}\{ \alpha_t(X)+\gamma_t's_t(Y(t)) \},}
#'          where \eqn{\text{logit}^{-1}(x) = (1 + exp(-x))^{-1}}. This specification has sensitivity parameters
#'          \eqn{\gamma = (\gamma_0, \gamma_1)}, which describe how treatment assignment depends marginally on
#'          each potential outcome, and a parmeter \eqn{\alpha_t(X)} in each arm that is identified by the observed 
#'          data once \eqn{gamma_t} is specified. \cr
#'          Under these settings, the missing outcome distribution can be infered as a tilt of the observed outcome 
#'          distribution.
#'          
#' @export
#'
#' @examples
#' # Observed data in treatment group
#' NHANES_trt <- NHANES %>% dplyr::filter(trt_dbp == 1)
#' x_trt <- NHANES_trt %>% select(-one_of("trt_dbp", "ave_dbp"))
#' y_trt <- NHANES_trt %>% select(ave_dbp)
#'
#' # Observed data in control group
#' NHANES_ctrl <- NHANES %>% dplyr::filter(trt_dbp == 0)
#' x_ctrl <- NHANES_ctrl %>% select(-one_of("trt_dbp", "ave_dbp"))
#' y_ctrl <- NHANES_ctrl %>% select(ave_dbp)
#' 
#' # Ribbon Plot of QTE
#' outcome_model <- fit_outcome(x_trt, y_trt, x_ctrl, y_ctrl, largest_gamma=0.05, joint=FALSE)
#' plot(outcome_model, type="ate")

ate_plot = function(tukeyFit){
   
   stopifnot(class(tukeyFit) == "tukeyFit")
   te_mat <- tukeyFit$te_mat
   ate_mean <- tukeyFit$ate_mean
   ns_elements <- tukeyFit$ns_elements
   
   superheat::superheat(te_mat,
                        heat.pal = rev(c("#d7191c","#fdae61","#ffffbf","#abd9e9", "#2c7bb6")),
                        heat.lim = round(c(-max(abs(ate_mean)), max(abs(ate_mean)))),
                        title = "Heat Map for ATE",
                        title.alignment = "center",
                        column.title = "gamma0",
                        row.title = "gamma1",
                        legend.breaks = round(seq(from = -max(abs(ate_mean)), to = max(abs(ate_mean)), length.out = 5)),
                        legend.vspace = 0.05,
                        left.label.size = .11,
                        column.title.size = 4,
                        row.title.size = 4,
                        X.text.size = 6,
                        X.text = ns_elements)
   
}
