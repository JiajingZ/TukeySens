#' @title Visualize ATE estimates by Contour Plot
#'
#' @description Visualize the same level of Average Treatment Effect (ATE)
#'              for different choices of sensitivity parameters by contour lines. 
#' @usage contour_ate(x_trt, y_trt, x_ctrl, y_ctrl, 
#'             largest_effect, gamma_length = 11, ...)
#' @param x_trt a \code{tibble} or data frame that contains observed pre-treatment variables of the treatment group
#' @param y_trt a vector of outcomes of the treatment group
#' @param x_ctrl a \code{tibble} or data frame that contains observed pre-treatment variables of the control group
#' @param y_ctrl a vector of outcomes of the control group
#' @param largest_effect the largest magnitude of sensitivity parameter to be investigated
#' @param gamma_length desired length of sensitivity parameter sequence, which needs to be an odd integer
#' @param joint logical. If TURE, the mean surface and residual variance will be estimated jointly for both treatment 
#'              groups; if FALSE (default), the mean surface and residual variance will be estimated independently for
#'              each treatment group.
#' @param ... additional argument to be passed to the function \code{BART::wbart}
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
#' ## ATE Contour Plot ##
#' contour_ate(x_trt, y_trt, x_ctrl, y_ctrl, largest_effect = 0.05)
#' contour_ate(x_trt, y_trt, x_ctrl, y_ctrl, largest_effect = 0.05, joint = TRUE)


contour_ate = function(x_trt, y_trt, x_ctrl, y_ctrl, largest_effect, gamma_length = 11, joint = FALSE, ...){
  if(joint){
    ## prepare joint dataset ##
    x_joint = rbind(x_trt, x_ctrl)
    trt = c(rep(1, nrow(x_trt)), rep(0, nrow(x_ctrl)))
    x_train_joint = cbind(x_joint, trt)  ## trt as a predictor when joint = T
    x_test_joint = cbind(x_joint, 1-trt)
    y_train_joint = rbind(y_trt, y_ctrl)
    
    ## jointly fit the model using BART ##
    joint_bart_fit <- BART::wbart(as.matrix(x_train_joint), as.matrix(y_train_joint), x.test=as.matrix(x_test_joint), ...)
    
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
    colnames(gamma_grid) = c("gamma_0", "gamma_1")
    
    att_joint <- mu_trt_obs_joint %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test_joint %o% rep(1, nrow(gamma_grid)) -
      matrix(sigma_joint^2, nrow = length(sigma_joint), ncol = ncol(mu_trt_obs_joint)) %o% gamma_grid[, 1]
    atc_joint <- mu_trt_test_joint %o% rep(1, nrow(gamma_grid)) - mu_ctrl_obs_joint %o% rep(1, nrow(gamma_grid)) -
      matrix(sigma_joint^2, nrow = length(sigma_joint), ncol = ncol(mu_trt_test_joint)) %o% gamma_grid[, 2]
    
    att_mean_joint <- apply(att_joint, c(1, 3), mean)
    atc_mean_joint <- apply(atc_joint, c(1, 3), mean)
    
    ate_mean_joint <- q * att_mean_joint + (1-q) * atc_mean_joint
    ate = apply(ate_mean_joint, 2, mean)
    ate_contour = cbind(gamma_grid, ate)
    
    contour_plot = ggplot() +
      theme_bw() +
      xlab(expression(gamma[0])) +
      ylab(expression(gamma[1])) +
      ggtitle("Contour Plot") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_contour(data = ate_contour, 
                   aes(x = gamma_0, y = gamma_1, z = ate, color = ..level..),
                   binwidth = 1)+
      scale_color_continuous(name = "ATE")
    
    directlabels::direct.label(contour_plot,"bottom.pieces")
  }
  if(!joint)
  {
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
    colnames(gamma_grid) = c("gamma_0", "gamma_1")
    
    att <- mu_trt_obs %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test %o% rep(1, nrow(gamma_grid)) -
      matrix(sig_trt_obs^2, nrow = length(sig_trt_obs), ncol = ncol(mu_trt_obs)) %o% gamma_grid[, 1]
    atc <- mu_trt_test %o% rep(1, nrow(gamma_grid)) -
      matrix(sig_ctrl_obs^2, nrow = length(sig_ctrl_obs), ncol = ncol(mu_trt_test)) %o% gamma_grid[, 2] -
      mu_ctrl_obs %o% rep(1, nrow(gamma_grid))
    # E(Y(1)-Y(0)|T=1) #
    att_mean <- apply(att, c(1, 3), mean)
    # E(Y(1)-Y(0)|T=0) #
    atc_mean <- apply(atc, c(1, 3), mean)
    # average treatment effect E(Y(1)-Y(0)) #
    ate_mean <- q * att_mean + (1-q) * atc_mean ## ndpost*(n_gamma0*n_gamma1)
    ate = apply(ate_mean, 2, mean)
    ate_contour = cbind(gamma_grid, ate)
    
    contour_plot = ggplot() +
      theme_bw() +
      xlab(expression(gamma[0])) +
      ylab(expression(gamma[1])) +
      ggtitle("Contour Plot") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_contour(data = ate_contour, 
                   aes(x = gamma_0, y = gamma_1, z = ate, color = ..level..),
                   binwidth = 1)+
      scale_color_continuous(name = "ATE")
    
    directlabels::direct.label(contour_plot,"bottom.pieces")
  }
}


























