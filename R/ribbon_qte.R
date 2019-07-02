#' @title Visualize QTE estimates by Ribbon Plot
#'
#' @description ....
#' @usage ribbon_qte(x_trt, y_trt, x_ctrl, y_ctrl, gamma_select, joint = FALSE)
#' @param x_trt a \code{tibble} or data frame that contains observed pre-treatment variables of the treatment group
#' @param y_trt a vector of outcomes of the treatment group
#' @param x_ctrl a \code{tibble} or data frame that contains observed pre-treatment variables of the control group
#' @param y_ctrl a vector of outcomes of the control group
#' @param gamma_select a matrix that contains choices of sensitivity parameters in rows
#' 
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
#' ## Ribbon Plot of QTE ##
#' gamma_select = rbind(c(0, 0), c(-0.03, 0.05), c(0.05, -0.02), c(-0.01, 0.01))
#' ribbon_qte(x_trt, y_trt, x_ctrl, y_ctrl, gamma_select)
#' ribbon_qte(x_trt, y_trt, x_ctrl, y_ctrl, gamma_select, joint = TRUE)


ribbon_qte = function(x_trt, y_trt, x_ctrl, y_ctrl, gamma_select, joint = FALSE){
  
  if(joint){
    ## joint dataset ##
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
    
    ## Correct mean of missing potential outcomes ##
    ## dim = n(draws) * n(obs) * n(gamma)
    mu_ctrl_test_joint_corrected_select <- mu_ctrl_test_joint %o% rep(1, nrow(gamma_select)) -
      matrix(sigma_joint^2, nrow = nrow(mu_ctrl_test_joint), ncol=ncol(mu_ctrl_test_joint)) %o% gamma_select[, 1]
    
    mu_trt_test_joint_corrected_select <- mu_trt_test_joint %o% rep(1, nrow(gamma_select)) -
      matrix(sigma_joint^2, nrow = nrow(mu_trt_test_joint), ncol=ncol(mu_trt_test_joint)) %o% gamma_select[, 2]
    
    
    nsample = nrow(mu_trt_obs_joint)
    n = nrow(y_trt) + nrow(y_ctrl)
    probs = seq(0.05, 0.95, by = 0.05)
    
    ## Function to Calculate Quantiles ##
    get_quantile = function(porbs, cumfun, xmin_init, xmax_init, prec = 1e-5){
      binary_search = function(prob, xmin, xmax){
        xmid = (xmin + xmax)/2
        prob_mid = cumfun(xmid)
        
        if(abs(prob - prob_mid) < prec)
          return(xmid)
        else if (prob_mid < prob)
          return(binary_search(prob, xmid, xmax))
        else 
          return(binary_search(prob, xmin, xmid))
      }
      
      ord = order(probs)
      qmin  = xmin_init
      sapply(sort(probs), function(prob){
        cur = binary_search(prob, qmin, xmax_init)
        qmin = cur
      })[ord]
    }
    
    ## QTE Array: n(draws)*n(gamma)*n(probs) ##
    qte = array(dim=c(nsample, nrow(gamma_select), length(probs)))
    for(i in 1:nsample)
    {
      print(i)
      for(k in 1:nrow(gamma_select))
      {
        print(k)
        ## complete distribution of Y(0) ##
        Y0comp_cum = Vectorize(function(y)
          1/n*sum(c(pnorm(y, mu_ctrl_obs_joint[i,], sigma_joint[i]), 
                    pnorm(y, mu_ctrl_test_joint_corrected_select[i, ,k], sigma_joint[i]))))
        ## complete distribution of Y(1) ##
        Y1comp_cum = Vectorize(function(y)
          1/n*sum(c(pnorm(y, mu_trt_obs_joint[i,], sigma_joint[i]), 
                    pnorm(y, mu_trt_test_joint_corrected_select[i, ,k], sigma_joint[i]))))
        
        q0 = get_quantile(probs, function(y) Y0comp_cum(y), -1000, 1000)
        
        q1 = get_quantile(probs, function(y) Y1comp_cum(y), -1000, 1000)
        
        qte[i, k, ] = q1 - q0
      }
    }
    
    q025 = apply(qte, c(2, 3), function(x) quantile(x, 0.025))
    q975 = apply(qte, c(2, 3), function(x) quantile(x, 0.975))
    
    SensitivityParams <- factor(rep(apply(gamma_select, 1, function(x) paste(x, collapse=", ")), 
                                    length(probs)))
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    ribbon_plot = 
      tibble(q025=as.vector(q025), q975=as.vector(q975), ## vectorize by column
           quantile = rep(probs, each=nrow(gamma_select)),  ## repeat each element of qtiles n(gamma_select)= 4 times
           SensitivityParams=SensitivityParams) %>%
      ggplot() +
      geom_ribbon(aes(x=quantile, ymin=q025, ymax=q975, fill=SensitivityParams), alpha=0.5) + 
      xlab("Quantile q") +
      ylab(substitute(paste(Q[q], a, Q[q], b), list(a="(Y(1)) - ", b="(Y(0))"))) +
      ggtitle("Quantile Treatment Effect for Pooled Model") + 
      geom_hline(aes(yintercept=0), linetype="dashed", col="black", size=0.5) +
      guides(color=guide_legend(override.aes=list(size=1.5, alpha=1), title.hjust=0.5, legend.title=element_text(size=14))) +
      scale_fill_manual(values=rev(gg_color_hue(4)), name=substitute(paste(a, gamma[0], b, gamma[1], c), list(a="(", b=", ", c=")")))  + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.title=element_text(size=12), axis.text=element_text(size=10), legend.title=element_text(size=8), 
            legend.text=element_text(size=8), legend.position="right")
    }
  
  if(!joint){
  
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
    
  
    ## Correct mean of missing potential outcomes ##
    ## dim = n(draws) * n(obs) * n(gamma)
    mu_ctrl_test_corrected_select <- mu_ctrl_test %o% rep(1, nrow(gamma_select)) -
      matrix(sig_ctrl_obs^2, nrow=length(sig_ctrl_obs), ncol=ncol(mu_trt_obs)) %o% gamma_select[, 1]
    
    mu_trt_test_corrected_select <- mu_trt_test %o% rep(1, nrow(gamma_select)) -
      matrix(sig_trt_obs^2, nrow=length(sig_trt_obs), ncol=ncol(mu_trt_test)) %o% gamma_select[, 2]
    
    
    nsample = nrow(mu_trt_obs)
    n = nrow(y_trt) + nrow(y_ctrl)
    probs = seq(0.05, 0.95, by = 0.05)
  
    get_quantile = function(porbs, cumfun, xmin_init, xmax_init, prec = 1e-5){
      binary_search = function(prob, xmin, xmax){
        xmid = (xmin + xmax)/2
        prob_mid = cumfun(xmid)
        
        if(abs(prob - prob_mid) < prec)
          return(xmid)
        else if (prob_mid < prob)
          return(binary_search(prob, xmid, xmax))
        else 
          return(binary_search(prob, xmin, xmid))
      }
      
      ord = order(probs)
      qmin  = xmin_init
      sapply(sort(probs), function(prob){
        cur = binary_search(prob, qmin, xmax_init)
        qmin = cur
      })[ord]
    }
  
    qte = array(dim=c(nsample, nrow(gamma_select), length(probs)))
    for(i in 1:nsample)
    {
      print(i)
      for(k in 1:nrow(gamma_select))
      {
        print(k)
        ## complete distribution of Y(0) ##
        Y0comp_cum = Vectorize(function(y)
          1/n*sum(c(pnorm(y, mu_ctrl_obs[i,], sig_ctrl_obs[i]), 
                    pnorm(y, mu_ctrl_test_corrected_select[i, ,k], sig_ctrl_obs[i]))))
        ## complete distribution of Y(1) ##
        Y1comp_cum = Vectorize(function(y)
          1/n*sum(c(pnorm(y, mu_trt_obs[i,], sig_trt_obs[i]), 
                    pnorm(y, mu_trt_test_corrected_select[i, ,k], sig_trt_obs[i]))))
        
        q0 = get_quantile(probs, function(y) Y0comp_cum(y), -1000, 1000)
        
        q1 = get_quantile(probs, function(y) Y1comp_cum(y), -1000, 1000)
        
        qte[i, k, ] = q1 - q0
      }
    }
    
    q025 = apply(qte, c(2, 3), function(x) quantile(x, 0.025))
    q975 = apply(qte, c(2, 3), function(x) quantile(x, 0.975))
    
    SensitivityParams <- factor(rep(apply(gamma_select, 1, function(x) paste(x, collapse=", ")), 
                                    length(probs)))
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
  
    ribbon_plot = 
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
  }
  return(ribbon_plot)
  
}












