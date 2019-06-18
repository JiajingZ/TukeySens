#' @title Presenting the analysis result of ATE by superheat
#'
#' @description create a plot....
#' @usage ATE.heatmap(x_trt, y_trt, x_ctrl, y_ctrl, largest_effect, gamma_length = 11)
#' @param x_trt a \code{tibble} or data frame that constains observed explanatory variables in the treatment group
#' @param y_trt a vector of observed outcomes in the treatment group
#' @param x_ctrl a \code{tibble} or data frame that contains observed explanatory variables in the control group
#' @param y_ctrl a vector of observed outcomes in the control group
#' @param largest_effect the largest scale of sensitivity parameter to be investigated
#' @param gamma_length desired length of sensitivity parameter sequence, which needs to be an odd integer
#' @param joint logical. If TURE, ...., if FALSE (default), ......
#' @return \code{ATE.heatmap} returns a superheat graph presenting the result of ATE with various
#'         choices of sensitivity parameter along with 95% confidence interval
#' @export
#'
#' @examples
#' # load data #
#' library(foreign)
#' NHANES <- read.dta(system.file("extdata", "NHANES3hbp_dbp.dta",
#'           package = "TukeySensitivity")) %>% as_tibble
#' NHANES <- NHANES %>% dplyr::filter(ave_dbp > 30)
#' trt_df <- NHANES %>% dplyr::filter(trt_dbp == 1)  ## data frame for treated
#' ctrl_df <- NHANES %>% dplyr::filter(trt_dbp == 0)  ## data frame for control
#'
#' # X and y in T=1 #
#' x.trt <- trt_df %>% select(-one_of(c("ave_dbp", "trt_dbp", "ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))
#' y.trt <- trt_df %>% select(ave_dbp)
#'
#' # X and y in T=0 #
#' x.ctrl <- ctrl_df %>% select(-one_of(c("ave_dbp", "trt_dbp", "ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))
#' y.ctrl <- ctrl_df %>% select(ave_dbp)
#'
#' gamma.max = 0.05  # as the coefficient of Y
#' ATE.heatmap(x.trt,y.trt,x.trl,y.ctrl,gamma.max)
#' ATE.heatmap(x.trt,y.trt,x.ctrl,y.ctrl, gamma.max, joint = TRUE)

ATE.heatmap = function(x_trt, y_trt, x_ctrl, y_ctrl, largest_effect, gamma_length = 11, joint = FALSE){

if(joint)
{
  x_joint = rbind(x_trt,x_ctrl)
  trt = c(rep(1,nrow(x_trt)), rep(0,nrow(x_ctrl)))
  x_train_joint = cbind(x_joint, trt) ## trt as a predictor when joint=T
  x_test_joint = cbind(x_joint, 1-trt)

  y_train_joint = rbind(y_trt,y_ctrl)

  joint_bart_fit <- BART::wbart(as.matrix(x_train_joint), as.matrix(y_train_joint), x.test=as.matrix(x_test_joint))

  mu_ctrl_obs_joint <- joint_bart_fit$yhat.train[, trt == 0] ## Y.hat(Y(0)|T=0), observed
  mu_ctrl_test_joint <- joint_bart_fit$yhat.test[, trt == 1] ## Y.hat(Y(1)|T=0), missing
  mu_trt_obs_joint <- joint_bart_fit$yhat.train[, trt == 1]  ## Y.hat(Y(1)|T=1), observed
  mu_trt_test_joint <- joint_bart_fit$yhat.test[, trt == 0] ## Y.hat(Y(0)|T=1), missing

  sigma_joint <- joint_bart_fit$sigma[101:1100]

  nctrl <- ncol(mu_ctrl_obs_joint)
  ntreat <- ncol(mu_ctrl_test_joint)
  q <- rbeta(1000, ntreat + 1, nctrl + 1)

  gamma_0 <- c(seq(from=-largest_effect, to=0, length.out=(gamma_length+1)/2),
               seq(from=0, to=largest_effect, length.out=(gamma_length+1)/2)[-1])
  gamma_1 <- gamma_0
  gamma_grid <- expand.grid(gamma_0, gamma_1)

  att_joint <- mu_trt_obs_joint %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test_joint %o% rep(1, nrow(gamma_grid)) -
    matrix(sigma_joint^2, nrow=length(sigma_joint), ncol=ncol(mu_trt_obs_joint)) %o% gamma_grid[, 1]
  atc_joint <- mu_trt_test_joint  %o% rep(1, nrow(gamma_grid))- mu_ctrl_obs_joint %o% rep(1, nrow(gamma_grid)) +
    matrix(sigma_joint^2, nrow=length(sigma_joint), ncol=ncol(mu_trt_test_joint)) %o% gamma_grid[, 2]

  att_mean_joint <- apply(att_joint, c(1, 3), mean)
  atc_mean_joint <- apply(atc_joint, c(1, 3), mean)

  ate_mean_joint <- q * att_mean_joint + (1-q) * atc_mean_joint

  te_mat_joint <- matrix(apply(ate_mean_joint, 2, mean), byrow=T,
                         nrow=length(gamma_0), ncol=length(gamma_0))
  probs_joint1 <- matrix(apply(ate_mean_joint, 2, function(x) mean(x < 0) ), byrow=T,
                         nrow=length(gamma_0), ncol=length(gamma_0))
  probs_joint2 <- matrix(apply(ate_mean_joint, 2, function(x) mean(x > 0) ), byrow=T,
                         nrow=length(gamma_0), ncol=length(gamma_0))

  ns_elements_bool <- (probs_joint1 > 0.025) & (probs_joint2 > 0.025)
  ns_elements <- matrix("", nrow=nrow(ns_elements_bool), ncol=ncol(ns_elements_bool))
  ns_elements[ns_elements_bool] <- "NS"

  colnames(te_mat_joint) <- gamma_0
  rownames(te_mat_joint) <- gamma_1

  superheat::superheat(te_mat_joint,
            heat.pal = c("blue", "light blue", "white", "orange", "red"),
            heat.lim = round(c(-max(abs(ate_mean_joint)), max(abs(ate_mean_joint)))),
            title = "Heat Map for ATE",
            title.alignment = "center",
            column.title=" ", row.title=" ",
            legend.breaks=round(seq(from=-max(abs(ate_mean_joint)), to=max(abs(ate_mean_joint)), length.out=5)),
            legend.vspace=0.05,
            left.label.size=.11,
            column.title.size = 8,
            row.title.size = 8,
            X.text.size=6,
            X.text=ns_elements)
}

if (!joint)
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
   gamma_0 <- c(seq(from=-largest_effect, to=0, length.out=(gamma_length+1)/2),
                seq(from=0, to=largest_effect, length.out=(gamma_length+1)/2)[-1])
   gamma_1 <- gamma_0
   gamma_grid <- expand.grid(gamma_0, gamma_1)

   att <- mu_trt_obs %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test %o% rep(1, nrow(gamma_grid)) -
          matrix(sig_trt_obs^2, nrow=length(sig_trt_obs), ncol=ncol(mu_trt_obs)) %o% gamma_grid[, 1]
   atc <- mu_trt_test  %o% rep(1, nrow(gamma_grid)) +
          matrix(sig_ctrl_obs^2, nrow=length(sig_ctrl_obs), ncol=ncol(mu_trt_test)) %o% gamma_grid[, 2] -
          mu_ctrl_obs %o% rep(1, nrow(gamma_grid))
   # E(Y(1)-Y(0)|T=1) #
   att_mean <- apply(att, c(1, 3), mean)
   # E(Y(1)-Y(0)|T=0) #
   atc_mean <- apply(atc, c(1, 3), mean)
   # average treatment effect E(Y(1)-Y(0)) #
   ate_mean <- q * att_mean + (1-q) * atc_mean ## ndpost*(n_gamma0*n_gamma1)

   te_mat <- matrix(apply(ate_mean, 2, mean), byrow=T,
                    nrow=length(gamma_1), ncol=length(gamma_0))
   colnames(te_mat) <- gamma_0
   rownames(te_mat) <- gamma_1

   prob_mat1 <- matrix(apply(ate_mean, 2, function(x) mean(x < 0) ), byrow = T,
                       nrow=length(gamma_1), ncol=length(gamma_0)) ## the prob(ate_mean < 0) in different case of gamma
   prob_mat2 <- matrix(apply(ate_mean, 2, function(x) mean(x > 0) ), byrow=T,
                       nrow=length(gamma_1), ncol=length(gamma_0)) ## the prob(ate_mean > 0) in different case of gamma

   ns_elements_bool <- (prob_mat1 > 0.025) & (prob_mat2 > 0.025)
   ns_elements_indep <- matrix("", nrow=nrow(ns_elements_bool), ncol=ncol(ns_elements_bool))
   ns_elements_indep[ns_elements_bool] <- "NS"

   superheat::superheat(te_mat,
             heat.pal = c("blue", "light blue", "white", "orange", "red"),
             heat.lim = round(c(-max(abs(ate_mean)), max(abs(ate_mean)))),
             title = "Heat Map for ATE",
             title.alignment = "center",
             column.title=" ",
             row.title=" ",
             legend.breaks=round(seq(from=-max(abs(ate_mean)), to=max(abs(ate_mean)), length.out=5)),
             legend.vspace=0.05,
             left.label.size=.11,
             column.title.size = 8,
             row.title.size = 8,
             X.text.size=6,
             X.text=ns_elements_indep)
   }

}
