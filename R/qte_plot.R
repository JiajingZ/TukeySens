#' @title Visualize QTE Estimates by Ribbon Plot
#'
#' @description Visualize the 95\% posterior credible band of the quantile treatment effects for
#'        selected settings of the sensitivity parameters.
#' @usage ribbon_qte(x_trt, y_trt, x_ctrl, y_ctrl, gamma_select, joint = FALSE)
#' @param tukeyFit \code{tibble} or data frame with observed pre-treatment variables for the treatment group
#' @param gamma_select a matrix with selected settings for sensitivity parameters in rows
#' 
#' @section Details: 
#'          For analysis details, please see \code{\link{heatmap_qte}} 
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
#' gamma_select = rbind(c(0, 0), c(-0.03, 0.05), c(0.05, -0.02), c(-0.01, 0.01))
#' plot(outcome_model, type="qte", gamma_select=gamma_select)


qte_plot <- function(tukeyFit, gamma_select=NULL){

  stopifnot(class(tukeyFit)=="tukeyFit")
  
  ## First and last
  if(is.null(gamma_select)) {
    gamma_select <- cbind(c(max(tukeyFit$gamma0), min(tukeyFit$gamma0)),
                          c(min(tukeyFit$gamma1), max(tukeyFit$gamma1)))
  }
  
  
  # estimated result for T=1 #
  mu_trt_obs <- tukeyFit$outcome_means$mu_trt_obs
  mu_trt_test <- tukeyFit$outcome_means$mu_trt_test
  sig_trt_obs <- tukeyFit$outcome_sds$sig_trt_obs
  
  # estimated results in T=0 #
  mu_ctrl_obs <- tukeyFit$outcome_means$mu_ctrl_obs
  mu_ctrl_test <- tukeyFit$outcome_means$mu_ctrl_test
  sig_ctrl_obs <- tukeyFit$outcome_sds$sig_ctrl_obs
  
  ## Correct mean of missing potential outcomes ##
  ## dim = n(draws) * n(obs) * n(gamma)
  mu_ctrl_test_corrected_select <- mu_ctrl_test %o% rep(1, nrow(gamma_select)) +
    matrix(sig_ctrl_obs^2, nrow=nrow(mu_ctrl_test), ncol=ncol(mu_ctrl_test)) %o% gamma_select[, 1]
  
  mu_trt_test_corrected_select <- mu_trt_test %o% rep(1, nrow(gamma_select)) -
    matrix(sig_trt_obs^2, nrow=nrow(mu_trt_test), ncol=ncol(mu_trt_test)) %o% gamma_select[, 2]
  
  
  nsample <- nrow(mu_trt_obs)
  n <- ncol(mu_ctrl_test) + ncol(mu_ctrl_obs)
  probs <- seq(0.05, 0.95, by = 0.05)
  
  get_quantile <- function(probs, cumfun, xmin_init, xmax_init, prec = 1e-4){
    binary_search <- function(prob, xmin, xmax){
      xmid <- (xmin + xmax)/2
      prob_mid <- cumfun(xmid)
      
      if(abs(prob - prob_mid) < prec)
        return(xmid)
      else if (prob_mid < prob)
        return(binary_search(prob, xmid, xmax))
      else 
        return(binary_search(prob, xmin, xmid))
    }
    
    ord <- order(probs)
    qmin  <- xmin_init
    sapply(sort(probs), function(prob){
      cur <- binary_search(prob, qmin, xmax_init)
      qmin <- cur
    })[ord]
  }
  
  

  
  max_mu <- max(mu_ctrl_obs, mu_ctrl_test_corrected_select, mu_trt_obs, mu_trt_test_corrected_select)
  min_mu <- min(mu_ctrl_obs, mu_ctrl_test_corrected_select, mu_trt_obs, mu_trt_test_corrected_select)
  max_sigma <- max(sig_ctrl_obs, sig_trt_obs)
  xmin <- qnorm(min(probs), mean=min_mu, sd=max_sigma)
  xmax <- qnorm(max(probs), mean=max_mu, sd=max_sigma)
  
  qte <- array(dim=c(nsample, nrow(gamma_select), length(probs)))
  pb <- progress_bar$new(format = " computing quantiles [:bar] :percent eta: :eta ", total = nsample*nrow(gamma_select), clear = FALSE, width=60)
  for(i in 1:nsample)
  {
    for(k in 1:nrow(gamma_select))
    {
      ## complete distribution of Y(0) ##
      Y0comp_cum <- Vectorize(function(y)
        1/n*sum(c(pnorm(y, mu_ctrl_obs[i,], sig_ctrl_obs[i]), 
                  pnorm(y, mu_ctrl_test_corrected_select[i, ,k], sig_ctrl_obs[i]))))
      ## complete distribution of Y(1) ##
      Y1comp_cum <- Vectorize(function(y)
        1/n*sum(c(pnorm(y, mu_trt_obs[i,], sig_trt_obs[i]), 
                  pnorm(y, mu_trt_test_corrected_select[i, ,k], sig_trt_obs[i]))))
      tm <- Sys.time()

      q0 <- get_quantile(probs, function(y) Y0comp_cum(y), xmin, xmax)
      q1 <- get_quantile(probs, function(y) Y1comp_cum(y), xmin, xmax)
      
      qte[i, k, ] <- q1 - q0
      pb$tick()
    }
  }
  pb$terminate()
  
  q025 <- apply(qte, c(2, 3), function(x) quantile(x, 0.025))
  q975 <- apply(qte, c(2, 3), function(x) quantile(x, 0.975))
  
  SensitivityParams <- factor(rep(apply(gamma_select, 1, function(x) paste(x, collapse=", ")), 
                                  length(probs)))
  
  gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  ribbon_plot <- 
    tibble(q025=as.vector(q025), q975=as.vector(q975), ## vectorize by column
           quantile = rep(probs, each=nrow(gamma_select)),  ## repeat each element of qtiles n(gamma_select)=4 times
           SensitivityParams=SensitivityParams) %>%
    ggplot() +
    geom_ribbon(aes(x=quantile, ymin=q025, ymax=q975, fill=SensitivityParams), alpha=0.5) + 
    xlab("Quantile q") +
    ylab(substitute(paste(Q[q], a, Q[q], b), list(a="(Y(1)) - ", b="(Y(0))"))) +
    ggtitle("Quantile Treatment Effect") + 
    geom_hline(aes(yintercept=0), linetype="dashed", col="black", size=0.5) +
    guides(color=guide_legend(override.aes=list(size=1.5, alpha=1), title.hjust=0.5, legend.title=element_text(size=14))) +
    scale_fill_manual(values=rev(gg_color_hue(4)), name=substitute(paste(a, gamma[0], b, gamma[1], c), list(a="(", b=", ", c=")")))  + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=8), 
          legend.text=element_text(size=8), legend.position="right")
  
  return(ribbon_plot)
  
}












