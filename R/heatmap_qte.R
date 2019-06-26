#' @title Visualize QTE estimates by Heatmap
#'
#' @description Visualize  Quantile Treatment Effect (QTE) estimates for a grid of sensitivity parameters. "NS" is used to denote 
#'              "not significant", meaning that the 95\% posterior credible interval of the ATE contains 0.
#' @usage heatmap_ate(x_trt, y_trt, x_ctrl, y_ctrl, 
#'             largest_effect, gamma_length = 11)
#' @param x_trt a \code{tibble} or data frame that contains observed pre-treatment variables of the treatment group
#' @param y_trt a vector of outcomes of the treatment group
#' @param x_ctrl a \code{tibble} or data frame that contains observed pre-treatment variables of the control group
#' @param y_ctrl a vector of outcomes of the control group
#' @param p a number of probablity between 0 and 1.
#' @param largest_effect the largest magnitude of sensitivity parameter to be investigated
#' @param gamma_length desired length of sensitivity parameter sequence, which needs to be an odd integer
#' @param joint logical. If TURE, the mean surface and residual variance will be estimated jointly for both treatment 
#'              groups; if FALSE (default), the mean surface and residual variance will be estimated independently for
#'              each treatment group.
#' @section Details: 
#'          The treatment assignment model is specified as:
#'          \deqn{f(T=1 \mid Y(t),X) = \text{logit}^{-1}\{ \alpha_t(X)+\gamma_t's_t(Y(t)) \}}
#'          
#' @export
#'
#' @examples
#' ## Observed data in treatment group ##
#' NHANES_trt <- NHANES %>% dplyr::filter(trt_dbp == 1)
#' x_trt <- NHANES_trt %>% select(-one_of("trt_dbp", "ave_dbp"))
#' y_trt <- NHANES_trt %>% select(ave_dbp)
#'
#' ## Observed data in control group ##
#' NHANES_ctrl <- NHANES %>% dplyr::filter(trt_dbp == 0)
#' x_ctrl <- NHANES_ctrl %>% select(-one_of("trt_dbp", "ave_dbp"))
#' y_ctrl <- NHANES_ctrl %>% select(ave_dbp)
#' 
#' ## ATE Heatmap ##
#' heatmap_qte(x_trt, y_trt, x_ctrl, y_ctrl, p = 0.63, largest_effect = 0.05)
#' heatmap_qte(x_trt, y_trt, x_ctrl, y_ctrl, p = 0.63, largest_effect = 0.05, joint = TRUE)




heatmap_qte = function(x_trt, y_trt, x_ctrl, y_ctrl, p, largest_effect, gamma_length = 11, joint = FALSE){
  if(joint){
    ## prepare joint dataset ##
    x_joint = rbind(x_trt, x_ctrl)
    trt = c(rep(1, nrow(x_trt)), rep(0, nrow(x_ctrl)))
    x_train_joint = cbind(x_joint, trt)  ## trt as a predictor when joint = T
    x_test_joint = cbind(x_joint, 1-trt)
    y_train_joint = rbind(y_trt, y_ctrl)
    
    ## jointly fit the model using BART ##
    joint_bart_fit <- BART::wbart(as.matrix(x_train_joint), as.matrix(y_train_joint), x.test=as.matrix(x_test_joint))
    
    mu_ctrl_obs_joint <- joint_bart_fit$yhat.train[ , trt == 0] ## Y.hat(Y(0)|T=0), observed
    mu_ctrl_test_joint <- joint_bart_fit$yhat.test[ , trt == 1] ## Y.hat(Y(1)|T=0), missing
    mu_trt_obs_joint <- joint_bart_fit$yhat.train[ , trt == 1]  ## Y.hat(Y(1)|T=1), observed
    mu_trt_test_joint <- joint_bart_fit$yhat.test[ , trt == 0] ## Y.hat(Y(0)|T=1), missing
    
    sigma_joint <- joint_bart_fit$sigma[101:1100]
    
    nctrl <- ncol(mu_ctrl_obs_joint)
    ntreat <- ncol(mu_ctrl_test_joint)
    q <- rbeta(1000, ntreat + 1, nctrl + 1)
    
    gamma_0 <- c(seq(from = -largest_effect, to = 0, length.out = (gamma_length + 1) / 2),
                 seq(from= 0, to = largest_effect, length.out = (gamma_length + 1) / 2)[-1])
    gamma_1 <- gamma_0
    gamma_grid <- expand.grid(gamma_0, gamma_1)
    
    tt_joint <- mu_trt_obs_joint %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test_joint %o% rep(1, nrow(gamma_grid)) -
      matrix(sigma_joint^2, nrow = length(sigma_joint), ncol = ncol(mu_trt_obs_joint)) %o% gamma_grid[, 1]
    tc_joint <- mu_trt_test_joint %o% rep(1, nrow(gamma_grid)) - mu_ctrl_obs_joint %o% rep(1, nrow(gamma_grid)) -
      matrix(sigma_joint^2, nrow = length(sigma_joint), ncol = ncol(mu_trt_test_joint)) %o% gamma_grid[, 2]
    
    qtt_joint <- apply(tt_joint, c(1, 3), quantile, p)
    qtc_joint <- apply(tc_joint, c(1, 3), quantile, p)
    
    qte_joint <- q * qtt_joint + (1-q) * qtc_joint
    
    qte_mean_joint <- matrix(apply(qte_joint, 2, mean), byrow = T,
                           nrow = length(gamma_0), ncol = length(gamma_0))
    probs_joint1 <- matrix(apply(qte_joint, 2, function(x) mean(x < 0)), byrow = T,
                           nrow = length(gamma_0), ncol = length(gamma_0))
    probs_joint2 <- matrix(apply(qte_joint, 2, function(x) mean(x > 0)), byrow = T,
                           nrow = length(gamma_0), ncol = length(gamma_0))
    
    ns_elements_bool <- (probs_joint1 > 0.025) & (probs_joint2 > 0.025)
    ns_elements <- matrix("", nrow = nrow(ns_elements_bool), ncol = ncol(ns_elements_bool))
    ns_elements[ns_elements_bool] <- "NS"
    
    colnames(qte_mean_joint) <- gamma_0
    rownames(qte_mean_joint) <- gamma_1
    
    superheat::superheat(qte_mean_joint,
                         heat.pal = c("blue", "light blue", "white", "orange", "red"),
                         heat.lim = round(c(-max(abs(ate_mean_joint)), max(abs(ate_mean_joint)))),
                         title = paste("Heat Map for ", as.character(p*100),"% QTE", sep=""),
                         title.alignment = "center",
                         column.title = "gamma0", row.title = "gamma1",
                         legend.breaks = round(seq(from = -max(abs(qte_joint)), to = max(abs(qte_joint)), length.out = 5)),
                         legend.vspace = 0.05,
                         left.label.size = .11,
                         column.title.size = 4,
                         row.title.size = 4,
                         X.text.size = 6,
                         X.text = ns_elements)
  }
  
  if (!joint){ 
    # Observed X and y in T=1 #
    x_train_trt = x_trt
    y_train_trt = y_trt
    
    # Observed X and y in T=0 #
    x_train_ctrl = x_ctrl
    y_train_ctrl = y_ctrl
    
    # X for missing outcomes #
    x_test_trt = x_ctrl
    x_test_ctrl = x_trt
    
    # Fitting Bart for T=1 #
    trt_bart_fit <- BART::wbart(as.matrix(x_train_trt), as.matrix(y_train_trt), as.matrix(x_test_trt))
    # estimated result for T=1 #
    mu_trt_obs <- trt_bart_fit$yhat.train  ## 1000 draws * n_trt obs
    mu_trt_test <- trt_bart_fit$yhat.test  ## 1000 draws * n_ctrl obs
    sig_trt_obs <- trt_bart_fit$sigma[101:1100]
    
    
    # Fitting Bart in T=0 #
    ctrl_bart_fit <- BART::wbart(as.matrix(x_train_ctrl), as.matrix(y_train_ctrl), as.matrix(x_test_ctrl))
    # estimated results in T=0 #
    mu_ctrl_obs <- ctrl_bart_fit$yhat.train  ## 1000 draws * n_ctrl obs
    mu_ctrl_test <- ctrl_bart_fit$yhat.test  ## 1000 draws * n_trt obs
    sig_ctrl_obs <- ctrl_bart_fit$sigma[101:1100]
    
    # sample weight from posterior of f(T) #
    q <- rbeta(1000, nrow(x_train_trt) + 1, nrow(x_train_ctrl) + 1)
    
    # gamma.grid #
    gamma_0 <- c(seq(from = -largest_effect, to = 0, length.out = (gamma_length + 1) / 2),
                 seq(from= 0, to = largest_effect, length.out = (gamma_length + 1) / 2)[-1])
    gamma_1 <- gamma_0
    gamma_grid <- expand.grid(gamma_0, gamma_1)
    
    tt <- mu_trt_obs %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test %o% rep(1, nrow(gamma_grid)) -
      matrix(sig_trt_obs^2, nrow = length(sig_trt_obs), ncol = ncol(mu_trt_obs)) %o% gamma_grid[, 1]
    tc <- mu_trt_test %o% rep(1, nrow(gamma_grid)) -
      matrix(sig_ctrl_obs^2, nrow = length(sig_ctrl_obs), ncol = ncol(mu_trt_test)) %o% gamma_grid[, 2] -
      mu_ctrl_obs %o% rep(1, nrow(gamma_grid))
    # E(Y(1)-Y(0)|T=1) #
    qtt <- apply(tt, c(1, 3), quantile, p)
    # E(Y(1)-Y(0)|T=0) #
    qtc <- apply(tc, c(1, 3), quantile, p)
    # average treatment effect E(Y(1)-Y(0)) #
    qte <- q * qtt + (1-q) * qtc ## ndpost*(n_gamma0*n_gamma1)
    
    qte_mean_mat <- matrix(apply(qte, 2, mean), byrow = T,
                     nrow = length(gamma_1), ncol = length(gamma_0))
    colnames(qte_mean_mat) <- gamma_0
    rownames(qte_mean_mat) <- gamma_1
    
    prob_mat1 <- matrix(apply(qte, 2, function(x) mean(x < 0)), byrow = T,
                        nrow = length(gamma_1), ncol = length(gamma_0)) ## the prob(ate_mean < 0) in different case of gamma
    prob_mat2 <- matrix(apply(qte, 2, function(x) mean(x > 0)), byrow = T,
                        nrow = length(gamma_1), ncol = length(gamma_0)) ## the prob(ate_mean > 0) in different case of gamma
    
    ns_elements_bool <- (prob_mat1 > 0.025) & (prob_mat2 > 0.025)
    ns_elements_indep <- matrix("", nrow = nrow(ns_elements_bool), ncol = ncol(ns_elements_bool))
    ns_elements_indep[ns_elements_bool] <- "NS"
    
    superheat::superheat(qte_mean_mat,
                         heat.pal = c("blue", "light blue", "white", "orange", "red"),
                         heat.lim = round(c(-max(abs(qte)), max(abs(qte)))),
                         title =paste("Heat Map for ", as.character(p*100), "% QTE", sep=""),
                         title.alignment = "center",
                         column.title = "gamma0",
                         row.title = "gamma1",
                         legend.breaks = round(seq(from = -max(abs(qte)), to = max(abs(qte)), length.out = 5)),
                         legend.vspace = 0.05,
                         left.label.size = .11,
                         column.title.size = 4,
                         row.title.size = 4,
                         X.text.size = 6,
                         X.text = ns_elements_indep)
  }
}
