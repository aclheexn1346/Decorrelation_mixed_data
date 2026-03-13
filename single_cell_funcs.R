run_bcd <- function(
    X, 
    block_size,
    zeropos_list,
    block_idx = NULL,
    lambda1=1, 
    lambda2=1,
    tol=1e-7,
    maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  thetahat <- diag(n)
  num_blocks <- length(zeropos_list)
  
  if(!is.null(block_idx)){
    N <- length(block_idx)
    iid_idx <- sapply(sim_data$block_idx, FUN = sample, size=1)
    X_iid <- X[iid_idx, ]
  }else{
    # if no block structure, treat x as i.i.d. samples as approximation
    X_iid <- X
    # X_iid <- X[seq(1, n, by = block_size),]
  }
  omega2_hat_iid <- estimate_omega_square(X_iid)
  old_hbeta <- matrix(0, p, p)
  old_htheta <- matrix(0, n, n)
  diff_beta <- 1
  diff_theta <- 1
  iter <- 0
  loss_his <- vector(length = maxIter)
  
  # correct_order = NULL
  # if (!is.null(block_idx)){
  #   correct_order <- match(rownames(X), names(unlist(block_idx)))
  #   # paste0('[INFO] Iter: ', iter, "\n"))
  #   # cat(rownames(X), names(unlist(block_idx)) )  
  # }
  
  while((diff_beta > tol || diff_theta > tol) && iter < maxIter){
    bhat <- estimate_b(
      n = n, p = p, X = X,
      theta_hat = thetahat,
      lambda= rep(lambda1, p)
    )
    S <- get_sample_cov(X, sqrt(omega2_hat_iid), bhat)
    thetahat <- estimate_theta(
      S = S, p = p, 
      lambda2 = lambda2,
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx
      # correct_order=correct_order
    )
    # err_beta <- norm(estimands$b - bhat, type = '2')^2 / p^2 
    # err_theta <- norm(estimands$theta - thetahat, type = '2')^2 / n^2
    curloss <- eval_loss_func(
      X = X, 
      block_idx = block_idx,
      homega = sqrt(omega2_hat_iid),
      hbeta = bhat, 
      htheta = thetahat,
      S = S, 
      lam1 = lambda1,
      lam2 = lambda2,
      iter = iter
    )
    cat(paste0('[INFO] Iter: ', iter, "\n"))
    cat(paste0('[INFO] Loss: ', round(curloss, 7), "\n"))
    # cat(paste0('[INFO] err_beta: ', round(err_beta,7), "\n"))
    # cat(paste0('[INFO] err_theta: ', round(err_theta,7), "\n"))
    cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
    cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
    cat(paste0("---------------------------------------------","\n"))
    diff_beta <- norm(bhat - old_hbeta, type = 'f')^2 / p^2
    diff_theta <- norm(thetahat - old_htheta, type = 'f')^2 / n^2
    old_hbeta <- bhat
    old_htheta <- thetahat
    iter <- iter + 1
    loss_his[iter] <- curloss
    if(iter == 1){
      cat("[INFO] Saving estimates after one iteration. \n")
      bhat_1iter <- bhat
      thetahat_1iter <- thetahat
    }
    Sys.sleep(0.01)
    flush.console()
  }
  cat(paste0('[INFO] Iter: ', iter, "\n"))
  cat(paste0('[INFO] Loss: ', round(curloss, 7), "\n"))
  # cat(paste0('[INFO] err_beta: ', round(err_beta,7), "\n"))
  # cat(paste0('[INFO] err_theta: ', round(err_theta,7), "\n"))
  cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
  cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
  cat(paste0("---------------------------------------------","\n"))
  # cat(paste0('[INFO] diff_beta: ', round(diff_beta, 7), "\n"))
  # cat(paste0('[INFO] diff_theta: ', round(diff_theta, 7), "\n"))
  dimnames(bhat) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  return(list(
    bhat=bhat, 
    thetahat=thetahat,
    bhat_1iter=bhat_1iter,
    thetahat_1iter=thetahat_1iter,
    omegahat = omega2_hat_iid,
    losses=loss_his[loss_his!=0]
  ))
}

get_sample_cov <- function(
    X,
    h_omega, 
    h_beta
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  rho_root <- 1 / h_omega
  test.temp <- X%*%(diag(rho_root) - sweep(h_beta, MARGIN = 2, FUN = '*', STATS = rho_root))
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res), n, n)
  S.scale <- S/p
  return(S.scale)
}

estimate_b <- function(
    n, p, X, 
    theta_hat=diag(n),
    lambda
){
  choles <- chol(theta_hat) 
  est_B <- matrix(0, p, p) 
  XX <- cbind(0, X)
  B_temp <- foreach(
    i=2:p,
    .combine = "cbind",
    .packages = "glmnet"
  ) %dopar% {
    lars.temp <- glmnet(
      x = as.matrix(choles)%*%XX[,1:i],
      y = as.matrix(choles)%*%XX[,i+1],
      alpha = 1,
      intercept = F,
      lambda = lambda[i]/(2*n),
      standardize = F,
      family = "gaussian",
      thresh = 1e-10
    )
    c(as.numeric(lars.temp$beta)[-1], rep(0, p - i + 1))
  }
  est_B[,2:p] <- B_temp
  return(est_B)
}


estimate_theta <- function(
    S,
    p,
    lambda2, # this should be chosen by grid search
    block_size,
    block_idx=NULL,
    zeropos_list=NULL,
    # correct_order=NULL,
    seed=1
){
  n <- dim(S)[1]
  set.seed(seed)
  if(!is.null(zeropos_list)){
    cat('[INFO]  Using the input zero list to form blocks...', '\n')
    num_blocks = length(zeropos_list)  
  }else{
    cat('[INFO]  Using the input block_idx to form blocks...', '\n')
    num_blocks = length(block_idx)  
  }
  sig.blocks <- vector(mode = "list", length = num_blocks)
  if(!is.null(block_idx)){
    for(i in 1:num_blocks){
      if(length(block_idx[[i]]) == 1){
        sig.blocks[[i]] = 1
        next
      }
      zeros = zeropos_list[[i]]
      if(is.null(zeros) || dim(zeros)[1] == 0){
        zeros = NULL
      }
      if(length(block_idx[[i]]) < p){
        lam2 = 1
      }else{
        lam2 = lambda2
      }
      temp_sig <- glasso(s = S[block_idx[[i]], block_idx[[i]]],
                         nobs = p,
                         rho = lam2/p,
                         zero = zeros,
                         penalize.diagonal = T,
                         trace = F)
      # sig.blocks[[i]] <- cov2cor(temp_sig$w)
      sig.blocks[[i]] <- temp_sig$w
    }
  }else{
    for(i in 1:num_blocks){
      zeros = zeropos_list[[i]]
      if(length(as.numeric(zeros)) == 0 | is.null(zeros))
        zeros = NULL
      temp_sig <- glasso(s = S[(1+block_size*(i-1)) : min(block_size*i,n), 
                               (1+block_size*(i-1)) : min(block_size*i,n)],
                         nobs = p,
                         rho = lambda2/p,
                         zero = zeros,
                         penalize.diagonal = T)
      # sig.blocks[[i]] <- cov2cor(temp_sig$w)
      sig.blocks[[i]] <- temp_sig$w
    }    
  }
  sig.est <- as.matrix(do.call(bdiag, sig.blocks))
  # if(!is.null(correct_order)){
  #   sig.est <- sig.est[correct_order, correct_order]    
  # }
  theta_est <- round(solve(sig.est),5)
  return(theta_est)
}



estimate_omega_square <- function(X){
  p <- dim(X)[2]
  n <- dim(X)[1]
  res = numeric(length = p)
  res[1] = 1
  # if(p > n){
  #   for(i in 2:p){
  #     x = as.matrix(cbind(0,X[,1:(i-1)]))
  #     y = as.matrix(X[,i])  
  #     model = lm(y~x-1)
  #     res[i] = sum(model$residuals^2) / (n-1)
  #   }
  # }else{
  for(i in 2:p){
    x = as.matrix(cbind(0,X[,1:(i-1)]))
    y = as.matrix(X[,i])
    lambda <- sqrt(log(p)/n)
    model = glmnet(x = x,
                   y = y,
                   alpha = 1,
                   lambda = lambda,
                   intercept = F,
                   standardize = F,
                   family = "gaussian",
                   thresh = 1e-10) 
    res[i] <- norm(y - x%*%model$beta, '2')^2/(2*n) + lambda*norm(model$beta, '1')
    # }
  }
  
  return(res)
}


eval_loss_func <- function(X, block_idx, homega, hbeta, htheta, S, lam1, lam2, iter){
  #' This loss differs from the negative log likelihood by the omega term
  n <- dim(X)[1]
  p <- dim(X)[2]
  hphi <- sweep(hbeta, 2, homega^2, '/')
  
  if(!is.null(block_idx)){
    theta_term <- -p*sum(sapply(block_idx, function(x){sum(log(eigen(as.matrix(htheta[x,x]))$values))}))  
  }else{
    theta_term <- -p*sum(log(eigen(htheta)$values))
  }
  cat(paste0('[INFO] theta_term: ' ,theta_term, '\n'))
  if(is.infinite(theta_term)){
    warning(paste0('[INFO] Iter: ', iter, '. theta_term is infinite!!!\n'))
  }
  trace_term <- p*sum(diag(S%*%htheta))
  if(is.infinite(trace_term)){
    warning(paste0('[INFO] Iter: ', iter, '. trace_term is infinite!!!\n'))
  }
  cat(paste0('[INFO] trace_term: ' ,trace_term, '\n'))
  lam1_term <- lam1*sum(abs(hphi)) 
  if(is.infinite(lam1_term)){
    warning('lam1_term is infinite!!!\n')
  }
  cat(paste0('[INFO] lam1_term: ' ,lam1_term, '\n'))
  lam2_term <- lam2*sum(abs(htheta))
  if(is.infinite(lam2_term)){
    warning('lam2_term is infinite!!!\n')
  }
  cat(paste0('[INFO] lam2_term: ' ,lam2_term, '\n'))
  result <- theta_term + trace_term + lam1_term + lam2_term
  return(result)  
}

# check this one
# Might be the log likelihood function I need.  Can set omgmle to 1 i hope

BIC_dag <- function(X, block_idx=NULL, bmle, omgmle, theta){
  n <- dim(X)[1]
  p <- dim(X)[2]
  s0 <- sum(abs(bmle) > 1e-4)
  t0 <- as.integer((sum(abs(theta) > 1e-4) - n) / 2)
  LX <- chol(theta)%*%X
  test.temp <- (LX - LX%*%bmle)%*%diag(1/sqrt(omgmle))
  test.res <- apply(test.temp,2,tcrossprod)
  S <- matrix(rowSums(test.res),n,n)
  tracetrm <- sum(diag(S))
  if(!is.null(block_idx)){
    thetatrm <- -p*sum(sapply(block_idx, function(x){sum(log(eigen(as.matrix(theta[x,x]))$values))}))  
    fn <- 1 * (s0 + t0) 
  }else{
    thetatrm <- -p*sum(log(eigen(theta)$values))
    fn <- log(max(n,p)) * (s0 + t0) 
  }
  negloglikelihood <- n*sum(log(omgmle)) + thetatrm + tracetrm
  # fn <- log(max(n,p)) * (s0 + t0) 
  # fn <- 1 * (s0 + t0) 
  BIC <- negloglikelihood + fn
  return(list(BIC = BIC,
              negloglikelihood = negloglikelihood / 2,
              s0 = s0,
              penalty=fn,
              thetatrm=thetatrm
  ))
}

networkDAG_sol_path <- function(
    X,
    block_size,
    zeropos_list,
    block_idx=NULL,
    lambda_len=10,
    lambda2=100,
    lambda1_max_div = 200,
    lambda1_max_div2 = 200,
    maxIter=100
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  lambda.path <- get_lam_path(p, X, rho.est = rep(1,p), lambda_len, 
                              div = lambda1_max_div,
                              div2 = lambda1_max_div2)
  if(lambda_len == 1){
    lambda.path[1] = 0
  }
  
  
  # lambdaff.path <- lseq(5000, 100, 4)
  # if(lambda_len == 1){
  #   lambda.path[1] = 0
  # }
  
  saveRDS(lambda.path, file = "lambda_path.rds")
  # saveRDS(lambdaff.path, file = "lambdaff.path.rds")
  BICscores_main <- minrowcor_main <- rep(0, length(lambda.path))
  BICscores_1iter_main <- minrowcor_1iter_main <- rep(0, length(lambda.path))
  BICscores_baseline <- minrowcor_baseline <- rep(0, length(lambda.path))
  # BICscores_ff <- minrowcor_ff <- rep(0, length(lambda.path))
  for(k in 1:length(lambda.path)){
    set.seed(1)
    cat(paste0("===============================================================","\n"))
    cat(paste0('[INFO] Lambda k = ', k, "\n"))
    res <- run_bcd(
      X = X, 
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx,
      lambda1 = lambda.path[k],
      lambda2 = lambda2,
      maxIter = maxIter,
      tol = 1e-7
    )
    # check if max degree s < n
    if(any(apply(res$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
      cat("[INFO] lambda is too small, algorithm stops at k = ", k, "\n")
      break
    }
    
    # resff <- flipflop(
    #   X = X,
    #   block_size = block_size,
    #   zeropos_list = zeropos_list,
    #   block_idx = block_idx,
    #   lambda1 = lambdaff.path[k],
    #   lambda2 = lambda2
    # )
    
    # if(any(apply(resff$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
    #   cat("[INFO] lambdaff is too small, algorithm stops at k = ", k, "\n")
    #   break
    # }
    
    # saveRDS((resff), file = paste0("ff_lam_", k, ".rds"))
    saveRDS(res, file = paste0("main_lam_", k, ".rds"))
    # compute min correlation
    cor_est <- cor(t(chol(res$thetahat)%*%(X - X%*%res$bhat)))
    minrowcor_main[k] <- sum(abs(cor_est[upper.tri(cor_est)]))
    cor_est_1iter <- cor(t(chol(res$thetahat_1iter)%*%(X - X%*%res$bhat_1iter)))
    minrowcor_1iter_main[k] <- sum(abs(cor_est_1iter[upper.tri(cor_est_1iter)]))
    cor_est_baseline <- cor(t(X - X%*%res$bhat_1iter))
    minrowcor_baseline[k] <- sum(abs(cor_est_baseline[upper.tri(cor_est_baseline)]))
    # cor_est_ff <- cor(t(chol(resff$thetahat)%*%(X - X%*%resff$bhat)))
    # minrowcor_ff[k] <- sum(abs(cor_est_ff[upper.tri(cor_est_ff)]))
    #compute MLE
    mle_result <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat, 
      Lhat = chol(res$thetahat)
    )
    mle_result_1iter <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat_1iter,
      Lhat = chol(res$thetahat_1iter)
    )
    mle_result_baseline <- dag_mle_estimation(
      X = X, 
      Bhat = res$bhat_1iter,
      Lhat = diag(n)
    )
    # mle_result_ff <- dag_mle_estimation(
    #   X = X,
    #   Bhat = resff$bhat,
    #   Lhat = chol(resff$thetahat)
    # )
    saveRDS(mle_result, file = paste0("mainMLE_lam_", k, ".rds"))
    saveRDS(mle_result_1iter, file = paste0("mainMLE_1iter_lam_", k, ".rds"))
    saveRDS(mle_result_baseline, file = paste0("baselineMLE_lam_", k, ".rds"))
    # saveRDS(mle_result_ff, file = paste0("ffMLE_lam_", k, ".rds"))
    # compute BIC using MLE estimates
    BIC_result <- BIC_dag(
      X = X,
      block_idx = block_idx,
      bmle = mle_result$Bmle,
      omgmle = mle_result$omgmlesq,
      theta = res$thetahat
    )
    BICscores_main[k] <- BIC_result$BIC
    BIC_1iter_result <- BIC_dag(
      X = X,
      block_idx = block_idx,
      bmle = mle_result_1iter$Bmle,
      omgmle = mle_result_1iter$omgmlesq,
      theta = res$thetahat_1iter
    )
    BICscores_1iter_main[k] <- BIC_1iter_result$BIC
    BIC_baseline_result <- BIC_dag(
      X = X,
      block_idx = block_idx,
      bmle = mle_result_baseline$Bmle, 
      omgmle = mle_result_baseline$omgmlesq,
      theta = diag(n)
    )
    BICscores_baseline[k] <- BIC_baseline_result$BIC
    # BIC_ff_result <- BIC_dag(
    #   X = X,
    #   block_idx = block_idx,
    #   # bmle = mle_result_ff$Bmle,
    #   bmle = resff$bhat,
    #   # omgmle = mle_result_ff$omgmlesq,
    #   omgmle = resff$omega2hat,
    #   theta = resff$thetahat
    # )
    # BICscores_ff[k] <- BIC_ff_result$BIC
    saveRDS(BIC_baseline_result, paste0('BIC_baseline_result_', k, '.rds'))
    saveRDS(BIC_result, paste0('BIC_main_result_', k, '.rds'))
    saveRDS(BIC_1iter_result,  paste0('BIC_1iter_result_', k, '.rds'))
    # saveRDS(BIC_ff_result,  paste0('BIC_ff_result_', k, '.rds'))
  }
  saveRDS(BICscores_main, "BICscores_main.rds")
  saveRDS(minrowcor_main, "minrowcor_main.rds")
  saveRDS(BICscores_1iter_main, "BICscores_1iter_main.rds")
  saveRDS(minrowcor_1iter_main, "minrowcor_1iter_main.rds")
  saveRDS(BICscores_baseline, "BICscores_baseline.rds")
  saveRDS(minrowcor_baseline, "minrowcor_baseline.rds")
  # saveRDS(BICscores_ff, "BICscores_ff.rds")
  # saveRDS(minrowcor_ff, "minrowcor_ff.rds")
}

sim_newalgo_ordered <- function(
    args, 
    estimands,
    start_sim, 
    end_sim,
    lamLen=10,
    lambda2=1,
    lambda1_max_div=50,
    lambda1_max_div2=1000
){
  for(sim in start_sim:end_sim){
    dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
    setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
    con <- file("test.log")
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
    X_ <- sim_X(
      vers = sim, 
      n = args$n,
      p = estimands$realp,
      omg.sq = estimands$omg.sq,
      sig = estimands$sig, 
      b = estimands$b
    )
    saveRDS(X_, file = "X.rds")
    X <- X_$X
    networkDAG_sol_path(
      X = X, 
      block_size=args$block_size, 
      zeropos_list = estimands$zeropos_list,
      lambda_len = lamLen,
      lambda1_max_div = lambda1_max_div,
      lambda1_max_div2 = lambda1_max_div2,
      maxIter = 15,
      lambda2 = lambda2
    )
    sink() 
    sink(type="message")
    cat(paste0("Experiment ", sim, " completed. \n"))
    setwd("~/Documents/research/dag_network")
  }
}


get_shd_ordered <- function(
    kmainbic = NULL, 
    kmaincor = NULL, 
    kbenchbic = NULL, 
    kbenchcor = NULL,
    k1iterbic = NULL,
    k1itercor = NULL,
    kffbic = NULL,
    kffcor = NULL,
    pc_flag = F,
    ges_flag = F,
    sim,
    simID,
    bstar_adj, 
    s0, 
    b0,
    support_size,
    theta0,
    thresh = 0.1
){
  n <- dim(theta0)[1]
  p <- dim(b0)[1]
  allShdS <- vector(mode = "list", length = length(thresh))
  # 1iter -------------------------------------------------------------------
  # BIC
  main1iter_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmainbic, '.rds'))
  main1iter_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', k1itercor, '.rds'))
  main_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmainbic, '.rds'))
  main_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmaincor, '.rds'))
  baseline_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchbic, '.rds'))
  baseline_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchcor, '.rds'))
  
  for(j in 1:length(thresh)){
    output_ordered <- vector(mode = "list")
    bhattemp <- EbayesThresh::threshld(main1iter_best_bic$bhat_1iter, thresh[j])
    bhat_adj <- 1*(abs(main1iter_best_bic$bhat_1iter) > thresh[j])
    shdXmain1iter <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXmain1iter['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    shdXmain1iter['theta_l2err'] = norm(theta0 - main1iter_best_bic$thetahat_1iter, type = '2')^2 / support_size
    output_ordered$shdXmain1iter=shdXmain1iter
    
    #minCor
    bhattemp <- EbayesThresh::threshld(main1iter_best_cor$bhat_1iter, thresh[j])
    bhat_adj <- 1*(abs(main1iter_best_cor$bhat_1iter) > thresh[j])
    shdXmain1iterCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXmain1iterCor['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    shdXmain1iterCor['theta_l2err'] = norm(theta0 - main1iter_best_cor$thetahat_1iter, type = '2')^2 / support_size
    output_ordered$shdXmain1iterCor=shdXmain1iterCor  
    
    # main --------------------------------------------------------------------
    # BIC
    bhattemp <- EbayesThresh::threshld(main_best_bic$bhat, thresh[j])
    bhat_adj <- 1*(abs(main_best_bic$bhat) > thresh[j])
    shdXmain <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXmain['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    shdXmain['theta_l2err'] = norm(theta0 - main_best_bic$thetahat, type = '2')^2 / support_size
    output_ordered$shdXmain=shdXmain
    # minCor
    bhattemp <- EbayesThresh::threshld(main_best_cor$bhat, thresh[j])
    bhat_adj <- 1*(abs(main_best_cor$bhat) > thresh[j])
    shdXmainCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXmainCor['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    shdXmainCor['theta_l2err'] = norm(theta0 - main_best_cor$thetahat, type = '2')^2 / support_size
    output_ordered$shdXmainCor = shdXmainCor
    # baseline ----------------------------------------------------------------
    # BIC
    bhattemp <- EbayesThresh::threshld(baseline_best_bic$bhat_1iter, thresh[j])
    # baseline_best_bic <- readRDS(file = paste0(simID, '--', sim, '/baseline_lam_', kbenchbic, '.rds'))
    # bhat_adj <- 1*(abs(baseline_best_bic$bhat) > thresh[j])
    
    bhat_adj <- 1*(abs(baseline_best_bic$bhat_1iter) > thresh[j])
    shdXbaseline <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXbaseline['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    # shdXbaseline['beta_l2err'] = norm(b0 - baseline_best_bic$bhat, type = '2')^2 / s0
    shdXbaseline['theta_l2err'] = -100
    output_ordered$shdXbaseline=shdXbaseline
    # minCor
    bhattemp <- EbayesThresh::threshld(baseline_best_cor$bhat_1iter, thresh[j])
    bhat_adj <- 1*(abs(baseline_best_cor$bhat_1iter) > thresh[j])
    # baseline_best_cor <- readRDS(file = paste0(simID, '--', sim, '/baseline_lam_', kbenchcor, '.rds'))
    # bhat_adj <- 1*(abs(baseline_best_cor$bhat) > thresh[j])
    shdXbaselineCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
    shdXbaselineCor['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
    # shdXbaselineCor['beta_l2err'] = norm(b0 - baseline_best_cor$bhat, type = '2')^2 / s0
    shdXbaselineCor['theta_l2err'] = -100
    output_ordered$shdXbaselineCor = shdXbaselineCor
    # flipflop ----------------------------------------------------------------
    # BIC
    if(!is.null(kffbic)){
      ff_best_bic <- readRDS(file = paste0(simID, '--', sim, '/ff_lam_', kffbic, '.rds'))
      bhattemp <- EbayesThresh::threshld(ff_best_bic$bhat, thresh[j])
      bhat_adj <- 1*(abs(ff_best_bic$bhat) > thresh[j])
      shdXff <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
      shdXff['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
      shdXff['theta_l2err'] = norm(theta0 - ff_best_bic$thetahat, type = '2')^2 / support_size
      output_ordered$shdXff=shdXff
      # mincor
      ff_best_cor <- readRDS(file = paste0(simID, '--', sim, '/ff_lam_', kffcor, '.rds'))
      bhattemp <- EbayesThresh::threshld(ff_best_cor$bhat, thresh[j])
      bhat_adj <- 1*(abs(ff_best_cor$bhat) > thresh[j])
      shdXffCor <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
      shdXffCor['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
      shdXffCor['theta_l2err'] = norm(theta0 - ff_best_cor$thetahat, type = '2')^2 / support_size
      output_ordered$shdXffCor = shdXffCor  
    }
    
    if(pc_flag){
      pc_res <- readRDS(file = paste0(simID, '--', sim, '/pc_ordered.rds'))
      pc_adj_b = matrix(0,p,p)
      for(i in 2:p){
        if(length(pc_res[[i]]$coefficients) > 1){
          pc_adj_b[as.integer(names(pc_res[[i]]$coefficients)[-1]), 
                   i] = pc_res[[i]]$coefficients[-1]
        }
      }
      bhattemp <- EbayesThresh::threshld(pc_adj_b, thresh[j])
      bhat_adj <- 1*(abs(pc_adj_b) > thresh[j])
      shdXpc <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
      shdXpc['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
      shdXpc['theta_l2err'] = -100
      output_ordered$shdXpc=shdXpc
    }
    
    
    if(ges_flag){
      ges_res <- readRDS(file = paste0(simID, '--', sim, '/ges_ordered.rds'))
      ges_adj_b = matrix(0, p, p)
      for(i in 2:p){
        if(length(ges_res[[i]]$coefficients) > 1){
          ges_adj_b[as.integer(names(ges_res[[i]]$coefficients)[-1]), 
                    i] = ges_res[[i]]$coefficients[-1]
        }
      }
      bhattemp <- EbayesThresh::threshld(ges_adj_b, thresh[j])
      bhat_adj <- 1*(abs(ges_adj_b) > thresh[j])
      shdXges <- compute_SHD_dag(adj1 = bhat_adj, adj_true = bstar_adj, s0) %>% unlist()
      shdXges['beta_l2err'] = norm(b0 - bhattemp, type = '2')^2 / s0
      shdXges['theta_l2err'] = -100
      output_ordered$shdXges=shdXges
    }
    
    
    allShdS[[j]] = output_ordered
  }
  
  if(length(allShdS) == 1){
    return(allShdS[[1]])
  }
  return(allShdS)  
}


process_output_ordered <- function(
    simID='001', 
    estimands,
    args,
    start=1,
    num_sim=args$num_sim,
    ff_flag = F,
    pc_flag = F,
    ges_flag = F, 
    threshholds=0.1
){
  setwd(paste0('output/', simID))
  # args <- readRDS('args.rds')
  # estimands <- readRDS('estimands.rds')
  bstar_adj <- 1*(abs(estimands$b) > 0)
  for(sim in start:num_sim){
    # BCD results
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    minrowcor_main <- readRDS(paste0(simID, '--', sim, '/minrowcor_main.rds'))
    # BIC 1iter results
    BICscores_main1iter <- readRDS(paste0(simID, '--', sim, '/BICscores_1iter_main.rds'))  
    minrowcor_main1iter <- readRDS(paste0(simID, '--', sim, '/minrowcor_1iter_main.rds'))  
    # baseline 
    BICscores_baseline <- readRDS(paste0(simID, '--', sim, '/BICscores_baseline.rds')) 
    minrowcor_baseline <- readRDS(paste0(simID, '--', sim, '/minrowcor_baseline.rds'))
    # flipflop
    if(ff_flag){
      BICscores_ff <- readRDS(paste0(simID, '--', sim, '/BICscores_ff.rds'))
      minrowcor_ff <- readRDS(paste0(simID, '--', sim, '/minrowcor_ff.rds'))
      bestk_bic_ff <- which.min(BICscores_ff)
      bestk_cor_ff <- which.min(minrowcor_ff)
    }else{
      bestk_bic_ff <- NULL
      bestk_cor_ff <- NULL
    }
    
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_1iter <- which.min(BICscores_main1iter)
    bestk_cor_1iter <- which.min(minrowcor_main1iter)
    bestk_bic_baseline <- which.min(BICscores_baseline)
    bestk_cor_baseline <- which.min(minrowcor_baseline)
    
    SHD_stats <- get_shd_ordered(
      kmainbic = bestk_bic_main, 
      kmaincor = bestk_cor_main, 
      kbenchbic = bestk_bic_baseline, 
      kbenchcor = bestk_cor_baseline,
      k1iterbic = bestk_bic_1iter,
      k1itercor = bestk_cor_1iter,
      kffbic = bestk_bic_ff,
      kffcor = bestk_cor_ff,
      pc_flag = pc_flag,
      ges_flag = ges_flag,
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj, 
      support_size = sum(estimands$theta != 0),
      s0 = estimands$s0, 
      theta0 = estimands$theta,
      b0 = estimands$b,
      thresh = threshholds
    )
    
    testX <- sim_X(
      vers = sim+100, 
      n = args$n,
      p = estimands$realp,
      omg.sq = estimands$omg.sq,
      sig = estimands$sig, 
      b = estimands$b
    )    
    
    testllres <- testll(
      testX$X,
      simID, 
      sim,
      kmainbic = bestk_bic_main, 
      kmaincor = bestk_cor_main, 
      kbenchbic = bestk_bic_baseline, 
      kbenchcor = bestk_cor_baseline,
      k1iterbic = bestk_bic_1iter,
      k1itercor = bestk_cor_1iter,
      kffbic = bestk_bic_ff,
      kffcor = bestk_cor_ff,
      pc_flag = pc_flag,
      ges_flag = ges_flag
    )
    
    SHD_stats$shdXmain1iter['testll'] = testllres$bic_1iter
    SHD_stats$shdXmain1iterCor['testll'] = testllres$cor_1iter
    SHD_stats$shdXmain['testll'] = testllres$bic_main
    SHD_stats$shdXmainCor['testll'] = testllres$cor_main
    SHD_stats$shdXbaseline['testll'] = testllres$bic_baseline
    SHD_stats$shdXbaselineCor['testll'] = testllres$cor_baseline
    
    
    if(ff_flag){
      SHD_stats$shdXff['testll'] = testllres$bic_ff
      SHD_stats$shdXffCor['testll'] = testllres$cor_ff
    }
    
    if(pc_flag){
      SHD_stats$shdXpc['testll'] = testllres$pc
    }
    
    if(ges_flag){
      SHD_stats$shdXges['testll'] = testllres$ges
    }
    
    saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDstats_JMLR.rds"))
  }
  print(paste0("Sim", sim, "done", "SHDstats_JMLR.rds saved."))
  setwd("~/Documents/research/dag_network")
}


testll_unordered <- function(Xp){
  n = dim(Xp)[1]
  ges = GES_sol(Xp, originalX = Xp, thetahat = diag(n), decor = F, test = T)
  cat('[INFO] Running PC... \n')
  pc = pc_sol(Xp, originalX = Xp, thetahat = diag(n), decor = F, test = T)
  cat('[INFO] Running sparsebn... \n')
  sbn = sparsebn_sol(Xp, originalX = Xp, thetahat = diag(n), decor = F, test = T)
  # decorrelation -----------------------------------------------------------
  Xdecor <- get_Xdecor(Xp) # select lowest BIC on lambda path and decor
  gesdecor = GES_sol(
    X = Xdecor$X_decor,
    originalX = Xp,
    thetahat = Xdecor$thehat,
    decor = T,
    test = T
  )
  pcdecor = pc_sol(
    X = Xdecor$X_decor,
    originalX = Xp,
    thetahat = Xdecor$thehat,
    decor = T,
    test = T
  )
  sbndecor = sparsebn_sol(
    X = Xdecor$X_decor,
    originalX = Xp,
    thetahat = Xdecor$thehat,
    decor = T,
    test = T
  )
  return(list(
    gesll = -ges$negloglikelihood,
    pcll = -pc$negloglikelihood,
    sbnll = -sbn$negloglikelihood,
    gesdecorll = -gesdecor$negloglikelihood,
    pcdecorll = -pcdecor$negloglikelihood,
    sbndecorll = -sbndecor$negloglikelihood
  ))
}

testll <- function(
    testX,
    simID,
    sim,
    kmainbic = bestk_bic_main, 
    kmaincor = bestk_cor_main, 
    kbenchbic = bestk_bic_baseline, 
    kbenchcor = bestk_cor_baseline,
    k1iterbic = bestk_bic_1iter,
    k1itercor = bestk_cor_1iter,
    kffbic = NULL,
    kffcor = NULL,
    pc_flag = F,
    ges_flag = F
){
  output_testll <- vector(mode = "list")
  
  # 1iter
  mle1iter_best_bic <- readRDS(file = paste0(simID, '--', sim, '/mainMLE_1iter_lam_', k1iterbic, '.rds'))
  res1iter_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', k1iterbic, '.rds'))
  mle1iter_best_cor <- readRDS(file = paste0(simID, '--', sim, '/mainMLE_1iter_lam_', k1itercor, '.rds'))
  res1iter_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', k1itercor, '.rds'))
  
  testll <- BIC_dag(
    X = testX,
    bmle = mle1iter_best_bic$Bmle,
    omgmle = mle1iter_best_bic$omgmlesq,
    theta = res1iter_best_bic$thetahat_1iter
  )
  output_testll['bic_1iter'] = -testll$negloglikelihood
  testll <- BIC_dag(
    X = testX,
    bmle = mle1iter_best_cor$Bmle,
    omgmle = mle1iter_best_cor$omgmlesq,
    theta = res1iter_best_cor$thetahat_1iter
  )
  output_testll['cor_1iter'] = -testll$negloglikelihood
  
  # main
  mlemain_best_bic <- readRDS(file = paste0(simID, '--', sim, '/mainMLE_lam_', kmainbic, '.rds'))
  resmain_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmainbic, '.rds'))
  mlemain_best_cor <- readRDS(file = paste0(simID, '--', sim, '/mainMLE_lam_', kmaincor, '.rds'))
  resmain_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kmaincor, '.rds'))
  testll <- BIC_dag(
    X = testX,
    bmle = mlemain_best_bic$Bmle,
    omgmle = mlemain_best_bic$omgmlesq,
    theta = resmain_best_bic$thetahat
  )
  output_testll['bic_main'] = -testll$negloglikelihood
  testll <- BIC_dag(
    X = testX,
    bmle = mlemain_best_cor$Bmle,
    omgmle = mlemain_best_cor$omgmlesq,
    theta = resmain_best_cor$thetahat
  )
  output_testll['cor_main'] = -testll$negloglikelihood
  # baseline
  mlebaseline_best_bic <- readRDS(file = paste0(simID, '--', sim, '/baselineMLE_lam_', kbenchbic, '.rds'))
  # resbaseline_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchbic, '.rds'))
  mlebaseline_best_cor <- readRDS(file = paste0(simID, '--', sim, '/baselineMLE_lam_', kbenchcor, '.rds'))
  # resbaseline_best_cor <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', kbenchcor, '.rds'))
  testll <- BIC_dag(
    X = testX,
    bmle = mlebaseline_best_bic$Bmle,
    omgmle = mlebaseline_best_bic$omgmlesq,
    theta = diag(dim(testX)[1])
  )
  output_testll['bic_baseline'] = -testll$negloglikelihood
  testll <- BIC_dag(
    X = testX,
    bmle = mlebaseline_best_cor$Bmle,
    omgmle = mlebaseline_best_cor$omgmlesq,
    theta = diag(dim(testX)[1])
  )
  output_testll['cor_baseline'] = -testll$negloglikelihood
  
  if(!is.null((kffbic))){
    # KGLasso
    mleff_best_bic <- readRDS(file = paste0(simID, '--', sim, '/ffMLE_lam_', kffbic, '.rds'))
    resff_best_bic <- readRDS(file = paste0(simID, '--', sim, '/ff_lam_', kffbic, '.rds'))
    mleff_best_cor <- readRDS(file = paste0(simID, '--', sim, '/ffMLE_lam_', kffcor, '.rds'))
    resff_best_cor <- readRDS(file = paste0(simID, '--', sim, '/ff_lam_', kffcor, '.rds'))
    testll <- BIC_dag(
      X = testX,
      bmle = mleff_best_bic$Bmle,
      omgmle = mleff_best_bic$omgmlesq,
      theta = resff_best_bic$thetahat
    )
    output_testll['bic_ff'] = -testll$negloglikelihood
    testll <- BIC_dag(
      X = testX,
      bmle = mleff_best_cor$Bmle,
      omgmle = mleff_best_cor$omgmlesq,
      theta = resff_best_cor$thetahat
    )
    output_testll['cor_ff'] = -testll$negloglikelihood
  }
  
  if(pc_flag){
    # PC given order
    mle_pc <- readRDS(file = paste0(simID, '--', sim, '/pc_MLE_result.rds'))
    testll <- BIC_dag(
      X = testX,
      bmle = mle_pc$Bmle,
      omgmle = mle_pc$omgmlesq,
      theta = diag(dim(testX)[1])
    )
    output_testll['pc'] = -testll$negloglikelihood
  }
  
  if(ges_flag){
    mle_ges <- readRDS(file = paste0(simID, '--', sim, '/ges_MLE_result.rds'))
    testll <- BIC_dag(
      X = testX,
      bmle = mle_ges$Bmle,
      omgmle = mle_ges$omgmlesq,
      theta = diag(dim(testX)[1])
    )
    output_testll['ges'] = -testll$negloglikelihood
  }
  return(output_testll)
}

get_all_shd_ordered <- function(
    simID, 
    estimands,
    start=1,
    ff_flag = F,
    pc_flag = F,
    ges_flag = F,
    num_sim=10,
    nthrs = 50,
    stats_file_name='SHDstats.rds'
){
  #' 
  #' 
  setwd(paste0("output/",simID))
  thrs <- seq(0, 0.5,length.out = nthrs)
  # allShdS <- vector(mode = "list", length = length(thrs))
  bstar_adj <- 1*(abs(estimands$b) > 0)
  for(sim in start:num_sim){
    cat(paste0("Processing sim ", sim, " ===========\n"))
    allshd <- readRDS(paste0(simID, "--", sim, '/', stats_file_name))
    # baseline
    BICscores_baseline <- readRDS(paste0(simID, "--", sim, "/BICscores_baseline.rds"))
    minrowcor_baseline <- readRDS(paste0(simID, "--", sim, "/minrowcor_baseline.rds"))
    # BCD 1iter
    minrowcor_1iter <- readRDS(paste0(simID, "--", sim, "/minrowcor_1iter_main.rds"))
    BICscores_1iter <- readRDS(paste0(simID, "--", sim, "/BICscores_1iter_main.rds"))
    # BCD main
    BICscores_main <- readRDS(paste0(simID, "--", sim, "/BICscores_main.rds"))
    minrowcor_main <- readRDS(paste0(simID, "--", sim, "/minrowcor_main.rds"))
    
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_baseline <- which.min(BICscores_baseline)
    bestk_cor_baseline <- which.min(minrowcor_baseline)
    bestk_bic_1iter <- which.min(BICscores_1iter)
    bestk_cor_1iter <- which.min(minrowcor_1iter)
    
    # ff
    if (ff_flag){
      BICscores_ff <- readRDS(paste0(simID, "--", sim, "/BICscores_ff.rds"))
      minrowcor_ff <- readRDS(paste0(simID, "--", sim, "/minrowcor_ff.rds"))
      bestk_bic_ff = which.min(BICscores_ff)
      bestk_cor_ff = which.min(minrowcor_ff)  
    }else{
      bestk_bic_ff = NULL
      bestk_cor_ff = NULL
    }
    statsLen = length(allshd$shdXbaseline)
    # for(j in 1:length(thrs)){
    # cat("------------------------------\n")
    # cat(paste0("Compute threshold = ", thrs[j], "\n"))
    allShdS <- get_shd_ordered(
      kmainbic = bestk_bic_main,
      kmaincor = bestk_cor_main,
      kbenchbic = bestk_bic_baseline,
      kbenchcor = bestk_cor_baseline,
      k1iterbic = bestk_bic_1iter,
      k1itercor = bestk_cor_1iter,
      kffbic = bestk_bic_ff,
      kffcor = bestk_cor_ff,
      pc_flag = pc_flag,
      ges_flag = ges_flag,
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj,
      s0 = estimands$s0,
      support_size = sum(estimands$theta!=0),
      theta0 = estimands$theta,
      b0 = estimands$b,
      thresh = thrs
    )
    # }
    # allShdS <- readRDS(file = paste0(simID, '--', sim, "/allShdS.rds"))
    saveRDS(allShdS, file = paste0(simID, '--', sim, "/allShdS_jmlr.rds"))
    # match baseline edges E with BCD
    bestbenchBIC <- which.min(sapply(
      allShdS, 
      function(x) abs(x$shdXbaseline['pnum'] - allshd$shdXmain['pnum'])
    ))
    bestbenchCor <- which.min(sapply(
      allShdS, 
      function(x) abs(x$shdXbaselineCor['pnum'] - allshd$shdXmainCor['pnum'])
    ))
    cat(paste0("[INFP] thr chosen to be ", bestbenchBIC, " with value ", thrs[bestbenchBIC], "\n"))
    allshd$shdXbaseline[1:(statsLen-1)] <- allShdS[[bestbenchBIC]]$shdXbaseline
    allshd$shdXbaselineCor[1:(statsLen-1)] <- allShdS[[bestbenchCor]]$shdXbaselineCor
    
    # match ff edges E with BCD
    if(ff_flag){
      bestffBIC <- which.min(sapply(
        allShdS,
        function(x) abs(x$shdXff['pnum'] - allshd$shdXmain['pnum'])
      ))
      bestffCor <- which.min(sapply(
        allShdS,
        function(x) abs(x$shdXffCor['pnum'] - allshd$shdXmainCor['pnum'])
      ))  
      allshd$shdXff[1:(statsLen-1)] <- allShdS[[bestffBIC]]$shdXff
      allshd$shdXffCor[1:(statsLen-1)] <- allShdS[[bestffCor]]$shdXffCor
    }
    
    if(pc_flag){
      best_pc <- which.min(sapply(
        allShdS,
        function(x) abs(x$shdXpc['pnum'] - allshd$shdXmain['pnum'])
      ))
      allshd$shdXpc[1:(statsLen-1)] <- allShdS[[best_pc]]$shdXpc
    }
    
    if(ges_flag){
      best_ges <- which.min(sapply(
        allShdS,
        function(x) abs(x$shdXges['pnum'] - allshd$shdXmain['pnum'])
      ))
      allshd$shdXges[1:(statsLen-1)] <- allShdS[[best_ges]]$shdXges
    }
    
    saveRDS(allshd, paste0(simID, "--", sim, "/SHDclose_jmlr.rds"))
    cat("[INFO] Sim", sim, "is done. \n")
  }
  setwd("~/Documents/research/dag_network")
}

# get_all_shd_unordered <- function(
    #   simID, 
#   estimands,
#   num_sim
# ){
#   return()
# }


get_average_shd_ordered <- function(
    simID, nsim, 
    ff_flag=F,
    pc_flag=F,
    ges_flag=F
){
  SHDres <- readRDS(paste0("output/", simID, '/', simID, '--1/', "SHDclose_jmlr.rds"))
  num_statistic <- length(names(SHDres[[1]]))
  total <- data.frame(
    row.names = names(SHDres[[1]]), 
    shdXmain=rep(0, num_statistic),
    shdXbaseline=rep(0, num_statistic),
    shdXmain1iter=rep(0, num_statistic),
    shdXmainCor=rep(0, num_statistic),
    shdXbaselineCor=rep(0, num_statistic),
    shdXmain1iterCor=rep(0, num_statistic)
  )
  if(ff_flag){
    total$shdXff = rep(0, num_statistic)
    total$shdXffCor = rep(0, num_statistic)
  }
  
  if(pc_flag){
    total$shdXpc=rep(0, num_statistic)
  }
  if(ges_flag){
    total$shdXges=rep(0, num_statistic)
  }
  
  allshd_list = list()
  for (sim in 1:nsim){
    allshd <- as.data.frame(
      readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDclose_jmlr.rds"))
    )
    allshd <- allshd[, names(total)]
    allshd_list[[sim]] <- allshd
    total <- total + allshd
  }
  total <- round(total / nsim, 7)
  total$B0 = c(estimands$s0, rep(0, num_statistic-1))
  shd_se <- data.frame(apply(sapply(allshd_list, as.matrix, simplify="array"), c(1,2), sd)) / sqrt(nsim)
  saveRDS(shd_se, paste0("output/", simID, "/shd_se_jmlr.rds"))
  saveRDS(total, paste0("output/", simID, "/shd_average_jmlr.rds"))
}


get_average_shd_unordered <- function(simID, start, nsim){
  SHDres <- readRDS(paste0("output/", simID, '/', simID, '--1/', "SHDstats.rds"))
  num_statistic <- length(names(SHDres[[1]]))
  total <- data.frame(
    row.names = names(SHDres[[1]]), 
    shdXmain = rep(0, num_statistic),
    shd_pc = rep(0, num_statistic),
    shd_pc_decor=rep(0, num_statistic),
    shd_ges=rep(0, num_statistic),
    shd_ges_decor=rep(0, num_statistic),
    shd_sbn=rep(0, num_statistic),
    shd_sbn_decor=rep(0, num_statistic)
  )
  for (sim in start:nsim){
    allshd <- as.data.frame(
      readRDS(paste0("output/", simID,  '/', simID, "--", sim, "/SHDstats.rds"))
    )
    allshd <- allshd[, names(total)]
    total <- total + allshd
  }
  total <-  total / (nsim - start + 1)
  total$B0 = c(estimands$s0, rep(0,num_statistic-1))
  # saveRDS(total, paste0("output/", simID, "/shd_average.rds"))
  shd_se <- data.frame(apply(sapply(allshd_list, as.matrix, simplify="array"), c(1,2), sd)) / sqrt(nsim)
  saveRDS(shd_se, paste0("output/", simID, "/shd_se_jmlr.rds"))
  saveRDS(total, paste0("output/", simID, "/shd_average_jmlr.rds"))
}





GES_sol <- function(
    X, 
    originalX,
    thetahat,
    block_idx=NULL,
    decor=F,
    test=F
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  set.seed(100)
  fgsX <- tetradrunner(
    algoId = 'fges',
    df = X,
    scoreId = 'sem-bic',
    dataType = 'continuous',
    faithfulnessAssumed=TRUE,
    maxDegree=-1,
    verbose=TRUE
  )
  adjmat_fgesCPDAG_X <- get_adjmat_from_fges(
    fgsX$edges,
    p = p, 
    varnames = dimnames(X)[[2]]
  )
  pdag = getGraph(adjmat_fgesCPDAG_X)
  dag = pdag2dag(pdag)
  # Adjacency Matrix G:
  # G[i,j] = 1/2 if edge mark of edge i-j at j is head(kid)/tail(parent).
  dag_adj = showAmat(dag$graph)
  mle_result = get_mle_gespc(X = X, dag_adj = dag_adj)
  BIC_result_GES <- BIC_dag(
    X = originalX,
    block_idx = block_idx,
    bmle = mle_result$Bmle,
    omgmle = mle_result$omgmlesq,
    theta = thetahat
  )
  
  if(test){
    if(decor){
      saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG_decor_test.rds')  
      saveRDS(mle_result, file = "fGES_mle_result_decor_test.rds")
      saveRDS(BIC_result_GES, 'fGES_BIC_result_decor_test.rds')
    }else{
      saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG_test.rds')  
      saveRDS(mle_result, file = "fGES_mle_result_test.rds")
      saveRDS(BIC_result_GES, 'fGES_BIC_result_test.rds')
    }
    return(BIC_result_GES)
  }
  
  cat("If test is TRUE, you shouldn't see this!!! \n")
  
  if(decor){
    saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG_decor.rds')  
    saveRDS(mle_result, file = "fGES_mle_result_decor.rds")
    saveRDS(BIC_result_GES, 'fGES_BIC_result_decor.rds')
  }else{
    saveRDS(adjmat_fgesCPDAG_X, file = 'adjmat_fges_CPDAG.rds')  
    saveRDS(mle_result, file = "fGES_mle_result.rds")
    saveRDS(BIC_result_GES, 'fGES_BIC_result.rds')
  }
}


pc_sol_ordered <- function(
    X,
    thetahat,
    originalX,
    block_idx=NULL,
    decor=F,
    test=F
){
  
}



pc_sol <- function(
    X,
    thetahat,
    originalX,
    block_idx=NULL,
    decor=F,
    test=F
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  suffstat <- list(C = cor(X), n = n)
  res_pc <- pcalg::pc(
    suffStat = suffstat, 
    indepTest = gaussCItest, 
    alpha = 0.05,
    m.max = 5,
    labels = dimnames(X)[[2]]
  )
  adjmat_pc_CPDAG <- as(res_pc, "amat")
  dag_pc = pdag2dag(res_pc@graph)
  dag_adj_pc = showAmat(dag_pc$graph)
  mle_result = get_mle_gespc(X = X, dag_adj = dag_adj_pc)
  BIC_result_PC <- BIC_dag(
    X = originalX,
    block_idx = block_idx,
    bmle = mle_result$Bmle,
    omgmle = mle_result$omgmlesq,
    theta = thetahat
  )
  
  if(test){
    return(BIC_result_PC)
  }
  
  if(decor){
    saveRDS(adjmat_pc_CPDAG, file = 'adjmat_pc_CPDAG_decor.rds')
    saveRDS(mle_result, file = "pc_mle_result_decor.rds")
    saveRDS(BIC_result_PC, 'pc_BIC_result_decor.rds')
  }
  else{
    saveRDS(adjmat_pc_CPDAG, file = 'adjmat_pc_CPDAG.rds')  
    saveRDS(mle_result, file = "pc_mle_result.rds")
    saveRDS(BIC_result_PC, 'pc_BIC_result.rds')
  }
}

sparsebn_sol <- function(
    X,
    thetahat,
    originalX,
    block_idx=NULL,
    decor=F,
    test=F,
    target_nedges = 170
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  sbX <- sparsebnData(X, type = "continuous")
  sb_path <- estimate.dag(sbX, verbose = F)
  closestDAG = select(sb_path, edges = target_nedges)
  # sol_idx <- select.parameter(sb_path, sbX)
  # sol_base <- sb_path[[sol_idx]] # this is a DAG
  adjmat_sbbase <- get.adjacency.matrix(closestDAG) %>% as.matrix()
  #convert to CPDAG
  adjmat_sb_cpdag <- bnstruct::dag.to.cpdag(dag.adj.matrix = adjmat_sbbase)
  # get graph of DAG and compute MLE
  dag = getGraph(adjmat_sbbase)
  dag_adj = showAmat(dag)
  mle_result = get_mle_gespc(X = X, dag_adj = dag_adj)
  BIC_result_sbn <- BIC_dag(
    X = originalX,
    block_idx = block_idx,
    bmle = mle_result$Bmle,
    omgmle = mle_result$omgmlesq,
    theta = thetahat
  )
  
  if(test){
    return(BIC_result_sbn)
  }
  
  if(decor){
    saveRDS(adjmat_sb_cpdag, file = 'adjmat_sparsebn_CPDAG_decor.rds')
    saveRDS(mle_result, file = "sbn_mle_result_decor.rds")
    saveRDS(BIC_result_sbn, 'sbn_BIC_result_decor.rds')
  }
  else{
    saveRDS(adjmat_sb_cpdag, file = 'adjmat_sparsebn_CPDAG.rds')  
    saveRDS(mle_result, file = "sbn_mle_result.rds")
    saveRDS(BIC_result_sbn, 'sbn_BIC_result.rds')
  }
}

get_Xdecor <- function(Xp, best_bic_predefined = NULL){
  bic_score <- readRDS('BICscores_main.rds')
  bic_score_1iter <- readRDS('BICscores_1iter_main.rds')
  bic_score_baseline <- readRDS('BICscores_baseline.rds')
  
  best_bic <- which.min(bic_score)
  best_bic_1iter <- which.min(bic_score_1iter)
  best_bic_baseline <- which.min(bic_score_baseline)
  
  cat('[INFO] lambda index for the best BIC: ',
      best_bic, "|", paste0(best_bic_1iter, " (1iter)"),"\n")
  if(!is.null(best_bic_predefined)){
    best_res <- readRDS(paste0('main_lam_', best_bic_predefined, '.rds'))
    best_res_1iter <- readRDS(paste0('main_lam_', best_bic_predefined, '.rds'))  
  }else{
    best_res <- readRDS(paste0('main_lam_', best_bic, '.rds'))
    best_res_1iter <- readRDS(paste0('main_lam_', best_bic_1iter, '.rds'))  
  }
  # theta_l2err = norm(theta0 - main1iter_best_bic$thetahat_1iter, type = '2')^2 / support_size
  X_decor <- chol(best_res$thetahat) %*% Xp
  X_decor_1iter <- chol(best_res_1iter$thetahat) %*% Xp
  dimnames(X_decor) <- dimnames(Xp)
  dimnames(X_decor_1iter) <- dimnames(Xp)
  return(list(
    X_decor=X_decor,
    X_decor_1iter=X_decor_1iter,
    thehat = best_res$thetahat,
    thehat1iter = best_res_1iter$thetahat,
    best_bic = best_bic,
    best_bic_1iter=best_bic_1iter,
    best_bic_baseline=best_bic_baseline))
}

sim_newalgo_unordered <- function(
    args, 
    estimands, 
    start_sim=1, 
    end_sim=args$num_sim, 
    lambda1_max_div=50,
    lambda1_max_div2 = 200,
    lambda2 = 1,
    lamLen=15
){
  for(sim in start_sim:end_sim){
    dir.create(path = paste0("output/",args$setting, "/", args$setting, "--", sim))
    setwd(paste0("output/",args$setting, "/", args$setting, "--", sim))
    con <- file("test.log")
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
    X_ <- sim_X(vers = sim+20, 
                n = args$n,
                p = estimands$realp,
                omg.sq = estimands$omg.sq,
                sig = estimands$sig, 
                b = estimands$b)
    X <- X_$X
    XXp <- Permute_X(X, seed = sim)
    Xp <- XXp$Xperm
    saveRDS(X_, file = "X.rds")
    saveRDS(XXp, file = "Xp.rds")
    # Other methods -----------------------------------------------------------
    cat('[INFO] Running GES... \n')
    GES_sol(Xp, originalX = Xp, thetahat = diag(args$n), decor = F)
    cat('[INFO] Running PC... \n')
    pc_sol(Xp, originalX = Xp, thetahat = diag(args$n), decor = F)
    cat('[INFO] Running sparsebn... \n')
    sparsebn_sol(Xp, originalX = Xp, thetahat = diag(args$n), decor = F)
    # NetworkDAG --------------------------------------------------------------
    networkDAG_sol_path(
      X = Xp, 
      block_size=args$block_size, 
      zeropos_list = estimands$zeropos_list,
      lambda_len = lamLen, 
      lambda1_max_div = lambda1_max_div,
      lambda1_max_div2 = lambda1_max_div2,
      lambda2 = lambda2,
      maxIter = 15   
    ) 
    # decorrelation -----------------------------------------------------------
    Xdecor <- get_Xdecor(Xp) # select lowest BIC on lambda path and decor
    GES_sol(
      X = Xdecor$X_decor,
      originalX = Xp,
      thetahat = Xdecor$thehat,
      decor = T
    )
    pc_sol(
      X = Xdecor$X_decor,
      originalX = Xp,
      thetahat = Xdecor$thehat,
      decor = T
    )
    sparsebn_sol(
      X = Xdecor$X_decor,
      originalX = Xp,
      thetahat = Xdecor$thehat,
      decor = T
    )
    sink() 
    sink(type="message")
    setwd("~/Documents/research/dag_network")  
  }
}


get_shd_unordered <- function(
    bestk_bic_main, 
    sim,
    simID,
    bstar_adj, 
    theta_true,
    s0, 
    thresh = 0.1
){
  output_unordered <- vector(mode = "list")
  main_best_bic <- readRDS(file = paste0(simID, '--', sim, '/main_lam_', bestk_bic_main, '.rds'))
  bhat_adj_cpdag <- bnstruct::dag.to.cpdag((abs(main_best_bic$bhat) > thresh)+0)
  shdXmain <- unlist(compute_SHD_detail(bhat_adj_cpdag, bstar_adj, s0))
  theta_l2err = norm(theta_true - main_best_bic$thetahat, type = '2')^2 / sum(theta_true!=0)
  saveRDS(theta_l2err, 'theta_l2err.rds')
  cat('SHD for BCD is done. \n')
  # get PC shd --------------------------------------------------------------
  adjmat_pc_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_pc_CPDAG.rds'))
  adjmat_fges_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_pc_CPDAG_decor.rds'))
  shd_pc <- unlist(compute_SHD_detail(adjmat_pc_CPDAG, bstar_adj, s0))
  shd_pc_decor <- unlist(compute_SHD_detail(adjmat_fges_CPDAG_decor, bstar_adj, s0))
  cat('SHD for PC is done. \n')
  # get GES shd -------------------------------------------------------------
  adjmat_fges_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_fges_CPDAG.rds'))
  adjmat_fges_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_fges_CPDAG_decor.rds'))
  shd_ges <- unlist(compute_SHD_detail(adjmat_fges_CPDAG, bstar_adj, s0))
  shd_ges_decor <- unlist(compute_SHD_detail(adjmat_fges_CPDAG_decor, bstar_adj, s0))
  cat('SHD for GES is done. \n')
  # get sparsebn shd --------------------------------------------------------
  adjmat_sparsebn_CPDAG <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_sparsebn_CPDAG.rds'))
  adjmat_sparsebn_CPDAG_decor <- readRDS(
    file = paste0(simID, '--', sim, '/adjmat_sparsebn_CPDAG_decor.rds'))
  shd_sbn <- unlist(compute_SHD_detail(adjmat_sparsebn_CPDAG, bstar_adj, s0))
  shd_sbn_decor <- unlist(compute_SHD_detail(adjmat_sparsebn_CPDAG_decor, bstar_adj, s0))
  cat('SHD for sparsebn is done. \n')
  output_unordered$shdXmain=shdXmain
  output_unordered$shd_pc=shd_pc
  output_unordered$shd_pc_decor=shd_pc_decor
  output_unordered$shd_ges=shd_ges
  output_unordered$shd_ges_decor=shd_ges_decor
  output_unordered$shd_sbn=shd_sbn
  output_unordered$shd_sbn_decor=shd_sbn_decor
  return(output_unordered)
}


process_output_unordered <- function(simID = simID,start=1, nsim=10, thr = 0.1){
  setwd(paste0('output/', simID))
  estimands <- readRDS('estimands.rds')
  bstar_adj_cpdag <- bnstruct::dag.to.cpdag(1*(estimands$b != 0)) 
  p = dim(bstar_adj_cpdag)[1]
  dimnames(bstar_adj_cpdag) <- list(as.character(1:p), as.character(1:p))
  for(sim in start:nsim){
    cat(paste0('[INFO] Processing sim ', sim, '\n'))
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    bestk_bic_main <- which.min(BICscores_main)
    SHD_stats <- get_shd_unordered(
      bestk_bic_main = bestk_bic_main, 
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj_cpdag, 
      s0 = estimands$s0, 
      theta_true = estimands$theta,
      thresh = thr
    )
    # saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDstats.rds"))
    
    testX <- sim_X(
      vers = sim+100, 
      n = args$n,
      p = estimands$realp,
      omg.sq = estimands$omg.sq,
      sig = estimands$sig, 
      b = estimands$b
    ) 
    X = testX$X
    XXp <- Permute_X(X, seed = sim+1)
    Xp <- XXp$Xperm
    saveRDS(testX, file = paste0(simID, '--', sim, "/testX.rds"))
    saveRDS(XXp, file = paste0(simID, '--', sim, "/testXp.rds"))
    
    setwd(dir = paste0(simID, '--', sim))
    testllres <- testll_unordered(Xp)
    setwd('../')
    
    # SHD_stats <- readRDS('SHDstats.rds')    
    SHD_stats$shd_pc['testll'] = testllres$pcll
    SHD_stats$shd_pc_decor['testll'] = testllres$pcdecorll
    SHD_stats$shd_ges['testll'] = testllres$gesll
    SHD_stats$shd_ges_decor['testll'] = testllres$gesdecorll
    SHD_stats$shd_sbn['testll'] = testllres$sbnll
    SHD_stats$shd_sbn_decor['testll'] = testllres$sbndecorll
    SHD_stats$shdXmain['testll'] = NA
    saveRDS(SHD_stats, file = paste0(simID, '--', sim, "/SHDstats.rds"))
    
  }
  setwd("~/Documents/research/dag_network")
}

sample_sc_data <- function(
    full_log_vals,
    full_idx,
    size=rep(20, 7),
    seed=1
){
  set.seed(seed)  
  block_idx <- lapply(sample, full_idx, size=size)
  
  block_idx <- lapply(block_idx, sort)
  subsetID <- sort(unlist(block_idx))
  subsetXp <- full_log_vals[subsetID, ]
  carry = 0
  for(i in 1:length(block_idx)){
    block_idx[[i]][names(block_idx[[i]])] <- 1:length(block_idx[[i]]) + carry 
    carry <-  block_idx[[i]][length(block_idx[[i]])]
  }
  # subsetXp %>% dim()
  temp <- substr(dimnames(subsetXp)[[1]], start = 0, stop = 3) 
  cellnames <- sub(pattern = "_$",replacement =  "", x = temp)
  return(list(
    subsetXp=subsetXp,
    cellnames=cellnames,
    block_idx=block_idx
  ))
}

get_cell_block_idx <- function(cellnames){
  #' Return a list where each element is a vector of row indices belonging to the same block/cluster
  #' 
  unique_cell <- unique(cellnames)
  cellname_idx <- numeric(length = length(unique_cell))
  for(i in 1:length(unique_cell)){
    cellname_idx[i] <- which(cellnames == unique_cell[i]) %>% head(1)
  }
  cellname_idx <- c(cellname_idx, length(cellnames)+1)
  block_idx <- list()
  for(i in 1:length(unique_cell)){
    block_idx[[i]] <- seq(cellname_idx[i], cellname_idx[i+1]-1)
  }
  return(block_idx)
}



radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}



plot_cpdag <- function(cpdag, rescale = F, file){
  source("~/Documents/research/dag_network/R/newalgo/igraphplot2.R")
  # cpdag = fgesdag
  environment(plot.igraph2) <- asNamespace('igraph')
  environment(igraph.Arrows2) <- asNamespace('igraph')
  idx <- which(cpdag !=0, arr.ind=T)
  edgelist <- cbind(rownames(cpdag)[idx[,"row"]], colnames(cpdag)[idx[,"col"]])
  dir <- rep(T, dim(edgelist)[1])
  for(i in 1:length(dir)){
    if (cpdag[edgelist[i,2], edgelist[i,1]] == 1){
      dir[i] = F
    }
  }
  cat("Number of undirected edges is: ", sum(!dir) / 2)
  df <- data.frame(cbind(edgelist, dir))
  df$dir <- dir
  g <- graph_from_data_frame(df)  
  la <- layout.circle(g)
  lab.locs <- radian.rescale(x=1:length(V(g)), direction=-1, start=0)
  if(rescale){
    plot(g, 
         layout = la,
         vertex.label.degree = lab.locs,
         # layout = layout.fruchterman.reingold,
         edge.arrow.size=0.6, vertex.size=0.1, 
         vertex.label.dist=0.5,
         vertex.label.cex = 1.5)
  }else{
    setEPS()
    postscript(paste0(file, ".eps"))
    # par(mar = c(2,2,2,2))
    plot.igraph2(g, layout=la, 
                 edge.arrow.size= c(0.4,0.4)[as.factor(edge_attr(g, 'dir'))],
                 # edge.color = "black",
                 edge.color = c('grey', "grey")[as.factor(edge_attr(g, 'dir'))],
                 vertex.label.cex=0.7, 
                 vertex.label.degree = lab.locs,
                 vertex.label.dist=1.5,
                 vertex.label.family="Helvetica",
                 vertex.label.font=2,
                 vertex.shape="circle", 
                 vertex.size=2, 
                 vertex.label.color="black", 
                 edge.width=1, 
                 rescale = F,
                 asp = 0,
                 ylim=c(-1,1),xlim=c(-1,1)
    )
    dev.off()
    
    # plot(g, 
    #      rescale = rescale,
    #      asp = 0,
    #      ylim=c(-25,25),xlim=c(-30,30),
    #      # layout = layout.fruchterman.reingold,
    #      edge.arrow.size=0.6, vertex.size=0.1, 
    #      vertex.label.dist=0.2,
    #      vertex.label.cex = 1.5)
  }
}


reorder_data <- function(
    sub_grp,
    targetgene
){
  if(!assertthat::are_equal(dim(targetgene)[2], 51)){
    stop('Dimension of target gene is incorrect!!!')
  }
  # select the cells that appear in the clusters. Unnecessary if no cells were dropped.
  targetgene_small <- targetgene[(rownames(targetgene) %in% names(sub_grp)), ]
  group_idx = unique(sub_grp)
  block_rownames = vector(mode = 'list', length = length(group_idx))
  for(i in 1:length(block_rownames)){
    block_rownames[[i]] = names(sub_grp[which(sub_grp == group_idx[i])])
  }
  block_idx = vector(mode = 'list', length = length(group_idx))
  newdf <- targetgene_small[block_rownames[[1]],,drop=F]
  carry = 0
  block_idx[[1]] <- 1:length(block_rownames[[1]]) + carry 
  carry <-  length(block_idx[[1]])
  for(i in 2:length(block_idx)){
    newdf <- rbind(newdf, targetgene_small[block_rownames[[i]], ,drop=F])
    block_idx[[i]] <- 1:length(block_rownames[[i]]) + carry 
    carry <-  carry + length(block_rownames[[i]])
  }
  return(list(
    df = newdf,
    block_idx = block_idx
  ))
}


two_step_cluster <- function(
    othergenes,
    targetgene,
    sc_idx_full, 
    method = 'complete',
    corr_thr = c(0.820, 0.840, 0.770, 0.800, 0.830, 0.815, 0.780)
){
  if(!assertthat::are_equal(dim(targetgene)[2], 51)){
    stop('Dimension of target gene is incorrect!!!')
  }
  for(i in 1:length(sc_idx_full)){
    cell.cor <- othergenes[, sc_idx_full[[i]]] %>% cor(use="pairwise.complete.obs")
    cell.dist <- as.dist(1 - abs(cell.cor))
    cell.tree <- hclust(cell.dist, method=method)
    plot(cell.tree, 
         cex=0.3,
         cex.axis=1.5,
         ann=F,
         main= NA,
         hang = 0.03,
         cex.lab = 1.5)
    sub_grp <- cutree(cell.tree, h = 1-corr_thr[i])
    # saveRDS(table(sub_grp), file = paste0('sub_grp_', i, '.rds'))
    cat("[INFO] Cell type", i, ". Size of the largest cluster is: ", max(table(sub_grp)), "\n")
    targetgene_subset <- targetgene[(rownames(targetgene) %in% names(sub_grp)), ]  
    # reorder the rows here so that theta is block-diagonal
    df_tmp <- reorder_data(sub_grp = sub_grp, targetgene = targetgene_subset)  
    if(i == 1){
      resdf <- df_tmp$df
      reslist = df_tmp$block_idx
    }else{
      resdf <- rbind(resdf, df_tmp$df)  
      newlist <- lapply(df_tmp$block_idx, function(x) x + length(sc_idx_full[[i-1]]))
      reslist <- c(reslist, newlist)
    }
  }
  return(list(
    df = resdf,
    block_idx = reslist
  ))
}

get_mle_gespc <- function(X, dag_adj){
  n = dim(X)[1]
  p = dim(X)[2]
  bhat = matrix(0, p, p)
  omg_new_sq <- rep(0,p)
  names(omg_new_sq) = dimnames(dag_adj)[[1]]
  dimnames(bhat) = dimnames(dag_adj)
  if(!all(dimnames(bhat)[[1]] == dimnames(X)[[2]])){
    stop('Node names in X and dag_adj dont match!! \n')
  }
  for(node in colnames(X)){
    y = X[, node]
    par = names(which(dag_adj[node, ] == 2))
    if(length(par) != 0){
      x = X[, par]
      b = lm(y~x-1)
      bhat[par, node] = coef(b)  
      omg_new_sq[node] = sum(b$residuals^2) / n
    }else{
      omg_new_sq[node] = sum(y^2) / n
    }
  }
  return(list(
    Bmle = bhat,
    omgmlesq = omg_new_sq
  ))
}


get_roc_data_ordered <- function(simID, thrlb=0, thrub=0.8, thrlen=20){
  #' generate ROC curves
  #' 
  setwd(paste0('output/', simID))
  args <- readRDS('args.rds')
  estimands <- readRDS('estimands.rds')
  num_sim <- args$num_sim
  bstar_adj <- 1*(abs(estimands$b) > 0)
  
  threshes <- seq(thrlb, thrub, length.out = thrlen)
  allShdS <- vector(mode = "list", length = length(threshes))
  
  for(sim in 1:args$num_sim){
    cat("[INFO]: Processing sim ", sim, "\n")
    BICscores_main <- readRDS(paste0(simID, '--', sim, '/BICscores_main.rds'))  
    minrowcor_main <- readRDS(paste0(simID, '--', sim, '/minrowcor_main.rds'))
    BICscores_main1iter <- readRDS(paste0(simID, '--', sim, '/BICscores_1iter_main.rds'))
    minrowcor_main1iter <- readRDS(paste0(simID, '--', sim, '/minrowcor_1iter_main.rds'))
    BICscores_baseline <- readRDS(paste0(simID, '--', sim, '/BICscores_baseline.rds')) 
    minrowcor_baseline <- readRDS(paste0(simID, '--', sim, '/minrowcor_baseline.rds'))
    
    bestk_bic_main <- which.min(BICscores_main)
    bestk_cor_main <- which.min(minrowcor_main)
    bestk_bic_1iter <- which.min(BICscores_main1iter)
    bestk_cor_1iter <- which.min(minrowcor_main1iter)
    bestk_bic_baseline <- which.min(BICscores_baseline)
    bestk_cor_baseline <- which.min(minrowcor_baseline)
    
    BICscores_ff <- readRDS(paste0(simID, '--', sim, '/BICscores_ff.rds'))
    minrowcor_ff <- readRDS(paste0(simID, '--', sim, '/minrowcor_ff.rds'))
    bestk_bic_ff <- which.min(BICscores_ff)
    bestk_cor_ff <- which.min(minrowcor_ff)
    
    
    allShdS <- get_shd_ordered(
      kmainbic = bestk_bic_main, 
      kmaincor = bestk_cor_main, 
      kbenchbic = bestk_bic_baseline, 
      kbenchcor = bestk_cor_baseline,
      k1iterbic = bestk_bic_1iter,
      k1itercor = bestk_cor_1iter,
      sim = sim,
      simID = simID,
      bstar_adj = bstar_adj, 
      s0 = estimands$s0, 
      support_size = sum(estimands$theta!=0),
      theta0 = estimands$theta,
      b0 = estimands$b,
      kffbic = bestk_bic_ff,
      kffcor = bestk_cor_ff,
      pc_flag = T,
      ges_flag = T,
      thresh = threshes
    )
    saveRDS(allShdS, file = paste0(simID, '--', sim, "/allShdS_jmlr.rds"))
  }
  setwd("~/Documents/research/dag_network")
}



get_roc_plots <- function(simID, ymin=50, xmax=101, thrLen=20){
  #' generate ROC plot for BCD and baseline
  setwd(paste0("output/",simID))
  args <- readRDS('args.rds')
  estimands <- readRDS('estimands.rds')
  rocmain <- rocbench <- rocff <-  array(data = 0, dim = c(2, thrLen))
  for(sim in 1:args$num_sim){
    allShdS <- readRDS(paste0(simID, '--', sim, '/allShdS.rds'))  
    rocmain__ <- sapply(allShdS, function(x) c(FP =(x$shdXmain['pnum'] - x$shdXmain['TP']),
                                               TP = x$shdXmain['TP'] ))
    rocbench__ <- sapply(allShdS, function(x) c(FP =(x$shdXbaseline['pnum'] - x$shdXbaseline['TP']),
                                                TP = x$shdXbaseline['TP']))
    rocff__  <- sapply(allShdS, function(x) c(FP =(x$shdXff['pnum'] - x$shdXff['TP']),
                                              TP = x$shdXff['TP']))
    
    rocmain <- rocmain + rocmain__
    rocbench <- rocbench + rocbench__
    rocff <- rocff + rocff__
  }
  rocmain <- rocmain / thrLen
  rocbench <- rocbench / thrLen
  rocff <- rocff / thrLen
  
  
  # ymax <- 120
  # xmax <- min(estimands$s0,max(rocmain[1,]))
  xmin <- 0
  ymax <- min(estimands$s0, max(rocmain[2,]))
  # ymax = 83
  
  setEPS()
  postscript(paste0(args$setting, ".eps"))
  mar.default <- c(4,4,2,2)
  par(mar = mar.default + c(1, 1, 0, 0)) 
  plot(t(rocmain),
       xlim=c(xmin,xmax),
       ylim = c(ymin, ymax),
       type = "l",
       lty = 4,
       lwd = 4,
       col="red",
       # xlab = "False Positive",
       xlab = "",
       # ylab= "True Positive",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       # cex.axis = 2,
       # cex.lab = 3
  )
  points(t(rocbench),
         col = "blue",
         lty = 5,
         lwd = 4,
         type="l")
  points(t(rocff),
         col = "green",
         lty = 5,
         lwd = 4,
         type="l")
  axis(
    side = 1, 
    at = seq(xmin, xmax, length.out = 5),
    labels = round(seq(xmin, xmax, length.out = 5),1),
    cex.axis = 2.3)
  axis(
    side = 2,
    at = seq(ymin, ymax,length.out = 5), 
    labels = round(seq(ymin, ymax,length.out = 5),1),
    cex.axis = 2.3)
  legend("bottomright",
         legend = c("BCD", "Baseline", "KGLasso"),
         col = c("red", "blue",
                 "green"), 
         lty = c(5,5,5), 
         lwd = c(4,4,4),
         cex = 2)
  # legend("bottomright", legend = c("BCD"),
  #        col = c("red"), lty = c(4))
  dev.off()
  setwd("~/Documents/research/dag_network")
}

get_roc_plots_noff <- function(
    simID, ymin=50, 
    xmax=101,
    thrLen=20,
    start=1,
    nsim=10, 
    ff_flag=F
){
  #' generate ROC plot for BCD and baseline
  setwd(paste0("output/",simID))
  args <- readRDS('args.rds')
  estimands <- readRDS('estimands.rds')
  rocmain <- rocbench <- rocff <- array(data = 0, dim = c(2, thrLen))
  rocpc <- rocges <- array(data = 0, dim = c(2, thrLen))
  
  
  for(sim in start:nsim){
    # allShdS <- readRDS(paste0(simID, '--', sim, '/allShdS.rds'))
    allShdS <- readRDS(paste0(simID, '--', sim, '/allShdS_jmlr.rds'))  
    rocmain__ <- sapply(allShdS, function(x) c(FP =(x$shdXmain['pnum'] - x$shdXmain['TP']),
                                               TP = x$shdXmain['TP'] ))
    rocbench__ <- sapply(allShdS, function(x) c(FP =(x$shdXbaseline['pnum'] - x$shdXbaseline['TP']),
                                                TP = x$shdXbaseline['TP']))
    
    roc_pc_ <- sapply(
      allShdS, 
      function(x) c(FP =(x$shdXpc['pnum'] - x$shdXpc['TP']), 
                    TP = x$shdXpc['TP'])
    )
    
    roc_ges_ <- sapply(
      allShdS, 
      function(x) c(FP =(x$shdXges['pnum'] - x$shdXges['TP']), 
                    TP = x$shdXges['TP'])
    )
    rocmain <- rocmain + rocmain__
    rocbench <- rocbench + rocbench__
    rocpc <- rocpc + roc_pc_
    rocges <- rocges + roc_ges_
    
    if(ff_flag){
      rocff__ <- sapply(allShdS, function(x) c(FP =(x$shdXff['pnum'] - x$shdXff['TP']),
                                               TP = x$shdXff['TP']))
      rocff <- rocff + rocff__  
    }    
  }
  rocmain <- rocmain / (nsim - start + 1)
  rocbench <- rocbench /  (nsim - start + 1)
  rocpc <- rocpc /  (nsim - start + 1)
  rocges <- rocges /  (nsim - start + 1)
  if(ff_flag){
    rocff <- rocff / (nsim - start + 1)
  }
  
  # ymax <- 120
  # xmax <- min(estimands$s0,max(rocmain[1,]))
  xmin <- 0
  ymax <- min(estimands$s0, max(rocmain[2,]))
  # ymax = 83
  
  setEPS()
  postscript(paste0(args$setting, "_jmlr.eps"))
  mar.default <- c(4,4,2,2)
  par(mar = mar.default + c(1, 1, 0, 0)) 
  plot(t(rocmain),
       xlim=c(xmin,xmax),
       ylim = c(ymin, ymax),
       type = "l",
       lty = 4,
       lwd = 4,
       col="red",
       # xlab = "False Positive",
       xlab = "",
       # ylab= "True Positive",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       # cex.axis = 2,
       # cex.lab = 3
  )
  points(t(rocbench),
         col = "blue",
         lty = 5,
         lwd = 4,
         type="l")
  points(t(rocpc),
         col = "peachpuff4",
         lty = 5,
         lwd = 4,
         type="l")
  points(t(rocges),
         col = "purple",
         lty = 5,
         lwd = 4,
         type="l")
  if(ff_flag){
    points(t(rocff),
           col = "green",
           lty = 5,
           lwd = 4,
           type="l")
  }
  axis(
    side = 1, 
    at = seq(xmin, xmax, length.out = 5),
    labels = round(seq(xmin, xmax, length.out = 5),1),
    cex.axis = 2.3)
  axis(
    side = 2,
    at = seq(ymin, ymax,length.out = 5), 
    labels = round(seq(ymin, ymax,length.out = 5),1),
    cex.axis = 2.3)
  legend("bottomright",
         legend = c("BCD", "Baseline", 'KGLasso', "PC", "GES"),
         col = c("red", "blue", "green", "yellow", 'purple'), 
         lty = c(5,5,5,5,5), 
         lwd = c(4,4,4,4,4),
         cex = 2)
  # legend("bottomright", legend = c("BCD"),
  #        col = c("red"), lty = c(4))
  dev.off()
  setwd("~/Documents/research/dag_network")
}


get_shd_ji_diff <- function(simIDarr){
  setwd("~/Documents/research/dag_network")
  num_sims = length(simIDarr)
  shdgesdiff <- shdpcdiff <- shdsbdiff <- rep(0,num_sims * 10)
  JIgesdiff <- JIpcdiff <- JIsbdiff <- rep(0,num_sims * 10)
  
  skip = 0
  for(id in simIDarr){
    # 5 + 5
    cat('[INFO]: processing sim', id, '\n')
    setwd(paste0("output/",id))
    args <- readRDS('args.rds')
    simLen = as.numeric(args$num_sim)
    for (i in 1:simLen){
      # 10
      allshd <- readRDS(paste0(id, '--', i, '/SHDstats.rds'))  
      shdgesdiff[as.integer(i+skip)] <- allshd$shd_ges['myshd'] - allshd$shd_ges_decor['myshd']
      shdpcdiff[as.integer(i+skip)] <- allshd$shd_pc['myshd'] - allshd$shd_pc_decor['myshd']
      shdsbdiff[as.integer(i+skip)] <- allshd$shd_sbn['myshd'] - allshd$shd_sbn_decor['myshd']
      JIgesdiff[as.integer(i+skip)] <- allshd$shd_ges_decor['JI'] - allshd$shd_ges['JI']
      JIpcdiff[as.integer(i+skip)] <- allshd$shd_pc_decor['JI'] - allshd$shd_pc['JI']
      JIsbdiff[as.integer(i+skip)] <- allshd$shd_sbn_decor['JI'] - allshd$shd_sbn['JI']
    }
    setwd("~/Documents/research/dag_network")
    skip <- skip + simLen
  }
  
  return(list(
    shdgesdiff = shdgesdiff,
    shdpcdiff = shdpcdiff,
    shdsbdiff = shdsbdiff,
    JIgesdiff = JIgesdiff,
    JIpcdiff = JIpcdiff,
    JIsbdiff = JIsbdiff
  ))
}

plot_boxplot_shd_ja_bnlearn <- function(shd_data, simID, nsim=10){
  setwd('~/Documents/research/dag_network/output/simID')
  # thetas = c("Equalcorr", "Toeplitz", "Star", "ExpDecay", "AR")
  mydata <- data.frame(
    shddiff = c(
      shd_data$shdgesdiff, 
      shd_data$shdpcdiff,
      shd_data$shdsbdiff
    ),
    jadiff = c(
      shd_data$JIgesdiff,
      shd_data$JIpcdiff,
      shd_data$JIsbdiff
    ),
    id = rep(
      rep(
        rep(thetas, each=nsim), # 5 settings
        2), # (n < p) and (n ? p)
      nmethod), # 3 methods
    label = rep(rep(c("n<p", "n>p"), each = 50),nmethod),
    method = rep(c("GES", "PC", "SBN"), each = 2 * 5 * nsim)
  )
  
  for(theta in thetas){
    plot_shd_util(df = subset(mydata, (id == theta) & (shddiff < 1000)), theta_type = theta)
    plot_ja_util(df = subset(mydata, (id == theta) & (shddiff < 1000)), theta_type = theta)
  }
  setwd("~/Documents/research/dag_network")
}


plot_boxplot_shd_ja <- function(shd_data, nsim=10, nmethod=3){
  setwd('~/Documents/research/dag_network/output/')
  thetas = c("Equalcorr", "Toeplitz", "Star", "ExpDecay", "AR")
  mydata <- data.frame(
    shddiff = c(
      shd_data$shdgesdiff, 
      shd_data$shdpcdiff,
      shd_data$shdsbdiff
    ),
    jadiff = c(
      shd_data$JIgesdiff,
      shd_data$JIpcdiff,
      shd_data$JIsbdiff
    ),
    id = rep(
      rep(
        rep(thetas, each=nsim), # 5 settings
        2), # (n < p) and (n ? p)
      nmethod), # 3 methods
    label = rep(rep(c("n<p", "n>p"), each = 50),nmethod),
    method = rep(c("GES", "PC", "SBN"), each = 2 * 5 * nsim)
  )
  
  for(theta in thetas){
    plot_shd_util(df = subset(mydata, (id == theta) & (shddiff < 1000)), theta_type = theta)
    plot_ja_util(df = subset(mydata, (id == theta) & (shddiff < 1000)), theta_type = theta)
  }
  setwd("~/Documents/research/dag_network")
}

plot_shd_util <- function(
    df, 
    theta_type="Equal corr",
    labels=c("n<p" = "(100, 200)", "n>p" = "(300, 100)")
){
  setEPS()
  postscript(paste0(theta_type, "-shd.eps"))
  par(mar = c(2,2,2,2))
  plt <- ggplot(df, aes(x = label, y = shddiff)) + 
    geom_boxplot(aes(fill = method), width=1.5) + 
    xlab("Sample size (n, p)") + 
    ylab("Decrease in SHD") + 
    # facet_grid(.~label) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
    theme(
      strip.text.x = element_text(size = 25),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(colour = "black"),
      # axis.title=element_text(size=30,face="bold"),
      axis.title=element_blank(),
      axis.text = element_text(size=30),
      legend.text=element_text(size=30),
      # axis.text.x=element_blank(),
      legend.title = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "white")
    ) + 
    scale_x_discrete(labels=labels)
  print(plt)
  dev.off()
}

plot_ja_util <- function(df, theta_type="Equal corr", labels=c("n<p" = "(100, 200)", "n>p" = "(300, 100)")){
  scaleFUN <- function(x) sprintf("%.2f", x)
  setEPS()
  postscript(paste0(theta_type, "-ja.eps"))
  par(mar = c(1,1,1,1))
  plt <- ggplot(df, aes(x = label, y = jadiff)) + 
    geom_boxplot(aes(fill = method), width=1.5) + 
    xlab("Sample size (n, p)") + 
    ylab("Increase in Jaccard index") + 
    # facet_grid(.~label) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) + 
    theme(
      strip.text.x = element_text(size = 25),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(colour = "black"),
      # axis.title=element_text(size=30,face="bold"),
      axis.title=element_blank(),
      axis.text = element_text(size=30),
      legend.text=element_text(size=30),
      # axis.text.x=element_blank(),
      legend.title = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "white")
    ) + 
    scale_x_discrete(labels=labels) + 
    scale_y_continuous(labels = scaleFUN)
  print(plt)
  dev.off()
}

plot_testll <- function(df, theta_type="Equal corr"){
  # cols = gg_color_hue(3)
  setEPS()
  postscript(paste0(theta_type, "-3testll.eps"))
  par(mar = c(2,2,2,2))
  plt <- df %>% ggplot(aes(x=method, y=testll, fill = method)) + 
    geom_boxplot() + 
    xlab("Methods") + 
    ylab("Test data log likelihood") + 
    # scale_fill_manual(values = cols[1:2]) +
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) +
    theme(
      strip.text.x = element_text(size = 25),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(colour = "black"),
      # axis.title=element_text(size=30,face="bold"),
      axis.title=element_blank(),
      axis.text = element_text(size=30),
      legend.text=element_text(size=30),
      axis.text.x=element_blank(),
      legend.title = element_text(size = 30),
      strip.background = element_rect(color = "black", fill = "white")
    ) 
  print(plt)
  dev.off()
}

networkDAG_sol_path2 <- function(
    X,
    block_size,
    zeropos_list,
    block_idx=NULL,
    lambda_len=5,
    lambda2=100,
    lambdaff_max = 1000,
    lambdaff_min = 300,
    maxIter=50
){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  lambdaff.path <- lseq(lambdaff_max, lambdaff_min,  lambda_len)
  
  
  saveRDS(lambdaff.path, file = "lambdaff.path.rds")
  BICscores_ff <- minrowcor_ff <- rep(0, length(lambdaff.path))
  for(k in 1:length(lambdaff.path)){
    set.seed(1)
    cat(paste0("===============================================================","\n"))
    cat(paste0('[INFO] Lambda k = ', k, "\n"))
    
    resff <- flipflop(
      X = X,
      block_size = block_size,
      zeropos_list = zeropos_list,
      block_idx = block_idx,
      lambda1 = lambdaff.path[k],
      lambda2 = lambda2, 
      maxIter = maxIter
    )
    
    # if(any(apply(resff$bhat, 2, function(x) sum(abs(x) > 1e-4)) >= n)){
    #   cat("[INFO] lambdaff is too small, algorithm stops at k = ", k, "\n")
    # }
    saveRDS((resff), file = paste0("ff_lam_", k, ".rds"))
    # compute min correlation
    cor_est_ff <- cor(t(chol(resff$thetahat)%*%(X - X%*%resff$bhat)))
    minrowcor_ff[k] <- sum(abs(cor_est_ff[upper.tri(cor_est_ff)]))
    #compute MLE
    mle_result_ff <- dag_mle_estimation(
      X = X,
      Bhat = resff$bhat,
      Lhat = chol(resff$thetahat)
    )
    saveRDS(mle_result_ff, file = paste0("ffMLE_lam_", k, ".rds"))
    # compute BIC using MLE estimates
    BIC_ff_result <- BIC_dag(
      X = X,
      block_idx = block_idx,
      bmle = resff$bhat,
      omgmle = resff$omega2hat,
      theta = resff$thetahat
    )
    BICscores_ff[k] <- BIC_ff_result$BIC
    saveRDS(BIC_ff_result,  paste0('BIC_ff_result_', k, '.rds'))
  }
  saveRDS(BICscores_ff, "BICscores_ff.rds")
  saveRDS(minrowcor_ff, "minrowcor_ff.rds")
}