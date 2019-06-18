# load data #
library(foreign)
NHANES <- read.dta(system.file("extdata", "NHANES3hbp_dbp.dta",
                               package = "TukeySensitivity")) %>% as_tibble
NHANES <- NHANES %>% dplyr::filter(ave_dbp > 30)
trt_df <- NHANES %>% dplyr::filter(trt_dbp == 1)  ## data frame for treated
ctrl_df <- NHANES %>% dplyr::filter(trt_dbp == 0)  ## data frame for control

# X and y in T=1 #
x.trt <- trt_df %>% select(-one_of(c("ave_dbp", "trt_dbp", "ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))
y.trt <- trt_df %>% select(ave_dbp)

# X and y in T=0 #
x.ctrl <- ctrl_df %>% select(-one_of(c("ave_dbp", "trt_dbp", "ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))
y.ctrl <- ctrl_df %>% select(ave_dbp)

###########################################################################################
x_trt = x.trt
y_trt = y.trt
x_ctrl = x.ctrl
y_ctrl = y.ctrl

################################################################################################

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

largest_effect = 0.05 ###########
gamma_length = 11  #################
gamma_0 <- c(seq(from=-largest_effect, to=0, length.out=(gamma_length+1)/2),
             seq(from=0, to=largest_effect, length.out=(gamma_length+1)/2)[-1])
gamma_1 <- gamma_0
gamma_grid <- expand.grid(gamma_0, gamma_1)
colnames(gamma_grid) = c("gamma_0","gamma_1")

att_joint <- mu_trt_obs_joint %o% rep(1, nrow(gamma_grid)) -  mu_ctrl_test_joint %o% rep(1, nrow(gamma_grid)) -
  matrix(sigma_joint^2, nrow=length(sigma_joint), ncol=ncol(mu_trt_obs_joint)) %o% gamma_grid[, 1]
atc_joint <- mu_trt_test_joint  %o% rep(1, nrow(gamma_grid))- mu_ctrl_obs_joint %o% rep(1, nrow(gamma_grid)) +
  matrix(sigma_joint^2, nrow=length(sigma_joint), ncol=ncol(mu_trt_test_joint)) %o% gamma_grid[, 2]

att_mean_joint <- apply(att_joint, c(1, 3), mean)
atc_mean_joint <- apply(atc_joint, c(1, 3), mean)

ate_mean_joint <- q * att_mean_joint + (1-q) * atc_mean_joint
ate = apply(ate_mean_joint,2,mean)
ate_contour = cbind(gamma_grid,ate)


ggplot(ate_contour,aes(gamma_0,gamma_1,z=ate)) +
  geom_raster(aes(fill=ate)) +
  scale_fill_gradientn(colours=c("blue","light blue","white","orange","red")) +
  geom_contour(color="grey") +
  labs(title = "ATE Contour", x = expression(gamma[0]),
       y = expression(gamma[1]), fill = "ATE") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(ate_contour,aes(gamma_0,gamma_1,z=ate)) +
  geom_contour(aes(colour = stat(level)),size=1.2) +
  scale_color_gradientn(colours=c("blue","light blue","white","orange","red"),name="ATE") +
  labs(title = "ATE Contour", x = expression(gamma[0]),
       y = expression(gamma[1])) +
  theme(plot.title = element_text(hjust = 0.5))









