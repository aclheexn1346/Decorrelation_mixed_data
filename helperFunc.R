get_adjmat_from_fges <- function(edgelist, p, varnames){
  # return adjmatrix of cpdag
  myadjmat <- matrix(0, p, p)
  dimnames(myadjmat) <- list(varnames, varnames)
  for (i in 1:length(edgelist)) {
    par_name <- word(edgelist[i], 1)
    chil_name <- word(edgelist[i], -1)
    par_ind <- which(varnames == par_name)
    chil_ind <- which(varnames == chil_name)
    myadjmat[par_ind, chil_ind] <- 1
    if (grepl("<->", edgelist[i])  || grepl("---", edgelist[i]) ) {
      myadjmat[chil_ind, par_ind] <- 1
    }
  }
  return(myadjmat)
}

gen.B <- function(p, b.mag = .9, s0= 2*p, seed = 480, lower.thresh = .6){
  # Generate some random beta matrix with ordering such that it
  # is an upper diagonal matrix
  if(p <= 5 ) stop("p is too small!")
  set.seed(seed*2)
  invisible(capture.output(bb <- randomDag(seed = seed, numNodes = p,numEdges = s0)))
  b <-  get_adjmat_from_fges(bb$edges, length(bb$nodes), bb$nodes)
  b[b!=0] = runif(length(bb$edges), lower.thresh, b.mag)*(2*rbinom(length(bb$edges),1,0.5)-1)
  
  realp <- pp <- p
  dimnames(b) <- list(as.character(1:realp), as.character(1:realp))
  return(list(b=b, s0=s0, realp = realp, pp = pp))
}

sim_X <- function(vers, p, n, omg.sq, sig, b){
  #' simulate X from network DAG given its parameters
  #'
  #' \code{sim_X} returns X and E matrices generated from the network DAG
  #'
  #' @param
  #'
  set.seed(vers)
  eps_mat <- matrix(0, n, p)
  eps_mat[,1] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[1]*sig)
  eps_mat[,2] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[2]*sig)
  
  X <- matrix(0, n, p)
  
  X[,1] <- eps_mat[,1]
  X[,2] <- X[,1]*b[1,2] + eps_mat[,2]
  
  for(i in 3:p) {
    eps_mat[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[i]*sig) # Sigma is covariance between the n samples
    X[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat[,i]
    #if (i %% 50 == 0)
    #  cat("Getting ", i, "th column of X. \n" )
  }
  X <- apply(X, c(1, 2), binarize_func)
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat))
}

region_func = function(biv_value, c1, c2){
  if(biv_value[1] > c1 & biv_value[2] > c2){
    return(c(1,1))
  }
  if(biv_value[1] > c1 & biv_value[2] < c2){
    return(c(1,0))
  }
  if(biv_value[1] < c1 & biv_value[2] > c2){
    return(c(0,1))
  }
  if(biv_value[1] < c1 & biv_value[2] < c2){
    return(c(0,0))
  }
}

sim_X_LUM <- function(vers, p, n, omg.sq, sig, b){
  #' simulate X from network DAG given its parameters
  #'
  #' \code{sim_X} returns X and E matrices generated from the network DAG
  #'
  #' @param
  #' Create two epsilons I think and compare and whichever one is the one
  #' That gets 0 or 1 based on argmax
  #' LUM - Latent Utility Model
  set.seed(vers)
  eps_mat_1 <- matrix(0, n, p)
  eps_mat_1[,1] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[1]*sig)
  eps_mat_1[,2] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[2]*sig)
  
  
  X <- matrix(0, n, p)
  Z_mat <- matrix(0, n, p)
  
  trunc_vals = matrix(0, n, p)
  
  X[,1] <- as.numeric(eps_mat_1[,1] > 0)
  Z_mat[,1] <- eps_mat_1[,1]
  trunc_vals[,1] = rep(0,n)
  X[,2] <- as.numeric(X[,1]*b[1,2] + eps_mat_1[,2] > 0)
  Z_mat[,2] = X[,1]*b[1,2] + eps_mat_1[,2]
  trunc_vals[,2] = -X[,1]*b[1,2]
  
  
  for(i in 3:p) {
    eps_mat_1[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[i]*sig) # Sigma is covariance between the n samples
    X[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat_1[,i]
    trunc_vals[,i] = -rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*"))
    Z_mat[,i] = X[,i]
    X[,i] <- as.numeric(X[,i] > 0)
    
    #if (i %% 50 == 0)
    #  cat("Getting ", i, "th column of X. \n" )
  }
  # X <- apply(X, c(1, 2), binarize_func) # Just does it based on 0
  
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat_1, trunc_vals = trunc_vals, Z = Z_mat))
}

beta_est_loop4 = function(data, init_beta, init_epsilon ,block_sizes, loops, trueB){
  update_beta = init_beta
  eps_draw = init_epsilon
  rmse = numeric(loops)
  rmse_trueB = numeric(loops)
  btw_betas = numeric(loops)
  hold_uncordata = list()
  n = nrow(data)
  p = ncol(data)
  update_betas = list()

  for(i in 1:loops){
    # obtain -X\beta bounds
    trunc_vals = obtain_trunc_vals(data, update_beta)
    print(summary(c(update_beta)))
    # Estimate Sigma after running through one iteration
    if(i == 2){
      estimated_Sigma = Sig_Estimate_DAG(X = data, trunc_vals = trunc_vals, block_sizes = block_sizes)
      if(is.positive.definite(estimated_Sigma) == F){
        estimated_Sigma = nearPD(estimated_Sigma, corr = T, doSym = T)
      }
      # Checking eigen values for non-invertible blocks and adding regularizer of 0.05 if close to non-inv.
      cluster_number = length(block_sizes)
      for(j in 1:(cluster_number)){
        check_cov = estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]
        if(any(eigen(check_cov)$values < 0.001)){
          # Regularize covariance matrix if we get a non invertible matrix
          estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] = (estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] + 0.1*diag(block_sizes[j]))/(1 + 0.1)
        }
      }
      # Obtain cholesky decomposition
      L_hat = chol(solve(estimated_Sigma))
    }
    else if(i == 1){
      estimated_Sigma = diag(nrow = nrow(data))
      L_hat = chol(estimated_Sigma)
    }
    if(i %% 5 == 0){
      print(i)
    }
    
    # Cholesky decomposition of Theta = inverse(Sigma) such that L_hat %*% Sigma %*% t(L_hat) = I_n
    eps_draw = withTimeout({
      # Draw epsilon according to truncation values and estimated Sigma
      epsilon_draw(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = i)
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    # Obtain Z with updated epsilon
    Z = obtain_hidden_Z(X = data, beta = update_beta, epsilon = eps_draw)
    
    # Decorrelate data
    uncor_data = L_hat %*% Z
    if(i != 1){
      hold_uncordata[[i]] = uncor_data
    }
    
    LX = L_hat %*% data
    

    start = Sys.time()
    # obtain new beta and iterate between beta and drawing epsilon
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = update_beta)
    end = Sys.time()
    print(paste0("new beta step: ", end - start))
    update_betas[[i]] = new_beta_vals[[1]]
    btw_betas[i] = sqrt(mean(c(update_beta - new_beta_vals[[1]])^2))
    update_beta = new_beta_vals[[1]]
    rmse[i] = new_beta_vals[[2]]
    diff = c(update_beta - trueB)
    rmse_trueB[i] = sqrt(mean((diff[which(diff != 0)])^2))
  }
  return(list("beta" = update_beta, 
              "Sigma" = estimated_Sigma, 
              "Error" = rmse, 
              "Error_trueB" = rmse_trueB,
              "between" = btw_betas,
              "hold_data" = hold_uncordata,
              "est_betas" = update_betas))
}

beta_est_loop_double_cov_est = function(data, init_beta, init_epsilon ,block_sizes, loops, trueB){
  update_beta = init_beta
  eps_draw = init_epsilon
  rmse = numeric(loops)
  rmse_trueB = numeric(loops)
  btw_betas = numeric(loops)
  hold_uncordata = list()
  n = nrow(data)
  p = ncol(data)
  update_betas = list()
  # spec_trunc = list()
  # spec_eps = numeric(loops)
  # spec_beta = numeric(loops)
  # spec_Z = numeric(loops)
  # spec_uncor = numeric(loops)
  # spec_trunc_vals = numeric(loops)
  for(i in 1:loops){
    # if(n > p & i == 1){
    #   print("here")
    #   update_beta = probit_beta(data = data, init_beta = update_beta)[[1]]
    # }
    
    trunc_vals = obtain_trunc_vals(data, update_beta)
    print(summary(c(update_beta)))
    if(i == 2|i == 20){
      estimated_Sigma = Sig_Estimate_DAG(X = data, trunc_vals = trunc_vals, block_sizes = block_sizes)
      if(is.positive.definite(estimated_Sigma) == F){
        estimated_Sigma = nearPD(estimated_Sigma, corr = T, doSym = T)
      }
      # Checking eigen values for non-invertible blocks and adding regularizer of 0.05 if close to non-inv.
      cluster_number = length(block_sizes)
      for(j in 1:(cluster_number)){
        check_cov = estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]
        if(any(eigen(check_cov)$values < 0.001)){
          estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] = (estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] + 0.1*diag(block_sizes[j]))/(1 + 0.1)
        }
        # (sum(block_sizes[-(j:cluster_number)]))
        # sum(block_sizes[1:j])
      }
      L_hat = chol(solve(estimated_Sigma))
      if (i == 2){
        start_est_Sigma = estimated_Sigma
      }
    }
    else if(i == 1){
      estimated_Sigma = diag(nrow = nrow(data))
      L_hat = chol(estimated_Sigma)
    }
    if(i %% 5 == 0){
      print(i)
    }
    
    # Cholesky decomposition of Theta = inverse(Sigma) such that L_hat %*% Sigma %*% t(L_hat) = I_n
    eps_draw = withTimeout({
      # need to work on changing the epsilon draw based on the variable size blocks
      epsilon_draw(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = i)
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    Z = obtain_hidden_Z(X = data, beta = update_beta, epsilon = eps_draw)
    
    
    uncor_data = L_hat %*% Z
    if(i != 1){
      hold_uncordata[[i]] = uncor_data
    }
    
    LX = L_hat %*% data
    
    #true_uncor_data = L_true %*% Z
    #true_LX = L_true %*% data
    start = Sys.time()
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = update_beta)
    end = Sys.time()
    print(paste0("new beta step: ", end - start))
    update_betas[[i]] = new_beta_vals[[1]]
    btw_betas[i] = sqrt(mean(c(update_beta - new_beta_vals[[1]])^2))
    update_beta = new_beta_vals[[1]]
    # spec_eps[i] = eps_draw[38,43]
    # spec_Z[i] = Z[38,43]
    # spec_uncor[i] = uncor_data[38,43]
    # spec_beta[i] = update_beta[38,43]
    # spec_trunc_vals[i] = trunc_vals[38,43]
    # spec2_beta[i] = update_beta[2,3]
    rmse[i] = new_beta_vals[[2]]
    diff = c(update_beta - trueB)
    rmse_trueB[i] = sqrt(mean((diff[which(diff != 0)])^2))
    #summary(c(update_beta))
  }
  return(list("beta" = update_beta, 
              "Sigma" = estimated_Sigma, 
              "Error" = rmse, 
              "Error_trueB" = rmse_trueB,
              "between" = btw_betas,
              "hold_data" = hold_uncordata,
              "est_betas" = update_betas,
              "start_estimated_Sigma" = start_est_Sigma))
}

beta_est_loop_input_covariance = function(data, adj_mat, cov_est, init_epsilon ,block_sizes){
  eps_draw = init_epsilon
  update_beta = adj_mat
  n = nrow(data)
  p = ncol(data)
  estimated_Sigma = cov_est
  L_hat = chol(solve(estimated_Sigma))
  for(i in 1:2){
    trunc_vals = obtain_trunc_vals(data, update_beta)
    eps_draw = withTimeout({
      epsilon_draw(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = 2)
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    Z = obtain_hidden_Z(X = data, beta = update_beta, epsilon = eps_draw)
    
    
    uncor_data = L_hat %*% Z
    LX = L_hat %*% data
    
    #true_uncor_data = L_true %*% Z
    #true_LX = L_true %*% data
    start = Sys.time()
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = adj_mat)
    end = Sys.time()
    print(paste0("new beta step: ", end - start))
    #summary(c(update_beta))
    update_beta = new_beta_vals[[1]]
  }
  
  return(list("beta" = update_beta, 
              "Sigma" = estimated_Sigma))
}
Sig_Estimate_DAG_test = function(X, beta, block_sizes){
  data = X
  trunc_vals = obtain_trunc_vals(X, beta)
  c_sim = trunc_vals # [1,] for first row, etc.
  t_data = t(data)
  cluster_number = length(block_sizes)
  estimated_sig = diag(nrow(data))
  avg_corr = matrix(data = 0, nrow = max(block_sizes), ncol = max(block_sizes)) # dont know if I need or not
  # for(i in 1:cluster_size){
  #   sig[((i-1)*block_size + 1):(i*block_size),((i-1)*block_size + 1):(i*block_size)] = struc_matrix
  # }
  # result = sum_log_lik(data = t_data, rho_vec = rho, c1 = c1_sim, c2 = c2_sim)
  
  
  #num_pairs = t(combn(rows, 2))
  
  rho = seq(0,0.95, length = 10)
  j = 1
  num_pairs = t(combn(1:block_sizes[j], 2))
  for(i in 1:nrow(num_pairs)){
    pair_data = t_data[,num_pairs[i,]]
    c1_run = c_sim[num_pairs[i,1],]
    c2_run = c_sim[num_pairs[i,2],]
    result = sum_log_lik(data = pair_data, rho_vec = rho, c1 = c1_run, c2 = c2_run)
    if(length(rho[which(result == max(result))]) == 0){
      estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]] = mean(avg_corr)
      estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1]] = mean(avg_corr)
    }
    else{
      estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]] = rho[which(result == max(result))]
      estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1]] = rho[which(result == max(result))]
    }
  }
  if(length(which(is.na(estimated_sig[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]))) == 0){
    avg_corr[1:block_sizes[j], 1:block_sizes[j]] = avg_corr[1:block_sizes[j], 1:block_sizes[j]] + estimated_sig[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]/cluster_number
  }
  
  for(k in 1:cluster_number){
    if(length(which(is.na(estimated_sig[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])]))) > 0){
      estimated_sig[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] = avg_corr
    }
  }
  
  if(is.positive.definite(estimated_sig) == F){
    estimated_Sigma = nearPD(estimated_sig, corr = T, doSym = T)
  }
  else{
    estimated_Sigma = estimated_sig
  }
  # Checking eigen values for non-invertible blocks and adding regularizer of 0.05 if close to non-inv.
  cluster_number = length(block_sizes)
  for(j in 1:(cluster_number)){
    check_cov = estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]
    if(any(eigen(check_cov)$values < 0.001)){
      estimated_Sigma[1:block_sizes[j],1:block_sizes[j]] = estimated_Sigma[1:block_sizes[j],1:block_sizes[j]] + 0.1*diag(block_sizes[j])/(1 + 0.1)
    }
    # (sum(block_sizes[-(j:cluster_number)]))
    # sum(block_sizes[1:j])
  }
  return(estimated_Sigma)
}

new_beta = function(uncor_Z, LX, init_beta){
  beta_update = matrix(data = 0, nrow = nrow(init_beta), ncol = ncol(init_beta))
  max_parents = 0
  n = nrow(LX)
  p = ncol(LX)
  for(i in 1:ncol(init_beta)){
    #print(i)
    if(length(which(init_beta[,i] != 0)) > 0){
      parents = which(init_beta[,i] != 0)
      if(length(parents) > 1){
        model = glmnet(x = data.matrix(LX[,parents]), y = uncor_Z[,i], alpha = 0, lambda = 20/n+50/p)
        beta_estimate = coef(model)
        beta_update[parents,i] = beta_estimate[2:nrow(beta_estimate),1]
        
      }
      else{
        model = lm(uncor_Z[,i]~LX[,parents])
        beta_estimate = summary(model)$coefficients
        beta_update[parents,i] = beta_estimate[2:nrow(beta_estimate),1]
      }
    }
  }
  diff = (beta_update - init_beta)
  rmse_error = sqrt(mean(diff[which(diff != 0)]^2))
  #print(max_parents)
  return(list(beta_update, rmse_error))
}

diff_sig_rho = function(n, estimated_Sig, rho){
  sig_2 = block_diag_sep(n = n, block_size = 2, struc_matrix = matrix(c(1,rho,rho,1), nrow = 2, ncol = 2, byrow = T))
  difference = c((sig_2 - estimated_Sig)^2)
  non_zero = difference[which(difference != 0)]
  return(sqrt(mean(non_zero)))
}

diff_sig = function(n, estimated_Sig, cov_struc){
  sig_2 = block_diag_sep(n = n, block_size = 5, struc_matrix = cov_struc)
  difference = c((sig_2 - estimated_Sig)^2)
  non_zero = difference[which(difference != 0)]
  return(sqrt(mean(non_zero)))
}

sing_val_func = function(uncor_Z, LX, init_beta){
  sing_vals = list()
  for(i in 1:ncol(init_beta)){
    #print(i)
    if(length(which(init_beta[,i] != 0)) > 0){
      parents = which(init_beta[,i] != 0)
      if(any(is.na(uncor_Z[,i]))){
        uncor_Z_clean = uncor_Z[-which(is.na(uncor_Z[,i])),i]
        LX_clean = LX[-which(is.na(uncor_Z[,i])),parents]
        model = lm(uncor_Z_clean[,i]~LX_clean[,parents])
      }
      else{
        model = lm(uncor_Z[,i]~LX[,parents])
      }
      sing_vals[[i]] = min(svd(model.matrix(model))$d)
    }
  }
  return(sing_vals)
}

beta_true = function(data, init_beta, init_epsilon, true_sigma, block_size, loops){
  update_beta = init_beta
  eps_draw = init_epsilon
  abs_errors = numeric(loops)
  n = nrow(data)
  p = ncol(data)
  trunc_vals = obtain_trunc_vals(data, update_beta, update_epsilon)
  #summary(c(true_trunc_vals))
  
  # try fixing epsilon
  # Cholesky decomposition of Theta = inverse(Sigma) such that L_hat %*% Sigma %*% t(L_hat) = I_n
  eps_draw = epsilon_draw(Sigma_hat = true_sigma, data = data, trunc_vals = trunc_vals, block_size = block_size, prev_iter = eps_draw, iter_num = 2) # if iter_num = 1, then won't use previous epsilon
  L_true = chol(solve(true_sigma))
  
  #update_epsilon = expectation_epsilon(Sigma = estimated_Sigma, B = update_beta, X = data, block_size = 2, init_epsilon = update_epsilon)
  Z = obtain_hidden_Z(X = data, beta = update_beta, epsilon = eps_draw)
  #summary(c(update_epsilon))
  #summary(c(Z))
  
  #uncor_data = L_hat %*% Z
  #LX = L_hat %*% data
  
  true_uncor_data = L_true %*% Z
  true_LX = L_true %*% data
  new_beta_vals = new_beta(uncor_Z = true_uncor_data, LX = true_LX, init_beta = update_beta)
  update_beta = new_beta_vals[[1]]
  #summary(c(update_beta))
  return(list("beta" = update_beta))
}

epsilon_draw = function(Sigma_hat, data, trunc_vals, block_sizes, prev_iter, iter_num){
  # have a n by p truncation values.  Now just check each element in data
  a = obtain_bounds(data,trunc_vals)$a
  b = obtain_bounds(data,trunc_vals)$b
  eps_draws = matrix(0, nrow = nrow(data), ncol = ncol(data))
  cluster_number = length(block_sizes)
  # Clean up bounds
  # for(k in 1:ncol(data)){
  #   out_bounds = which(prev_iter[,k] < a[,k] | prev_iter[,k] > b[,k])
  #   if(length(out_bounds) > 0){
  #     for(l in 1:length(out_bounds)){
  #       bounds = c(a[out_bounds[l],k],b[out_bounds[l],k])
  #       prev_iter[out_bounds[l], k] = bounds[which.min(abs(bounds - prev_iter[out_bounds[l],k]))]
  #     }
  #   }
  # }
  # (sum(block_sizes[-(j:cluster_number)]))
  # sum(block_sizes[1:j])
  for(i in 1:cluster_number){
    
    init_index = (sum(block_sizes[-(i:cluster_number)]))+1
    end_index = sum(block_sizes[1:i])
    data_cluster = data[init_index:end_index,]
    a_bounds = a[init_index:end_index,]
    b_bounds = b[init_index:end_index,]
    Sigma_cluster = Sigma_hat[init_index:end_index,init_index:end_index]
    for(j in 1:ncol(data)){
      if(iter_num == 1){
        #print(dim(data_cluster))
        if(block_sizes[i] == 1){
          draws = tmvtnorm::rtmvnorm(n = 15, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = mean(draws)
        }
        else{
          draws = tmvtnorm::rtmvnorm(n = 15, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = colMeans(draws)
        }
        if(any(is.nan(eps_draws[init_index:end_index,j]))){
          eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
        }
      }
      else{
        if(block_sizes[i] == 1){
          draws = tmvtnorm::rtmvnorm(n = 15, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = mean(draws)
        }
        else{
          draws = tmvtnorm::rtmvnorm(n = 15, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = colMeans(draws)
        }
        if(any(is.nan(eps_draws[init_index:end_index,j]))){
          eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
        }
      }
      # if(any(is.nan(eps_draws[,j]))){
      #   eps_draws[init_index:end_index,j] = tmvtnorm::rtmvnorm(n = 10, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j])
      # }
    }
  }
  return(eps_draws)
}

epsilon_draw_single = function(Sigma_hat, data, trunc_vals, block_sizes, prev_iter, iter_num){
  # have a n by p truncation values.  Now just check each element in data
  a = obtain_bounds(data,trunc_vals)$a
  b = obtain_bounds(data,trunc_vals)$b
  eps_draws = matrix(0, nrow = nrow(data), ncol = ncol(data))
  cluster_number = length(block_sizes)
  # Clean up bounds
  # for(k in 1:ncol(data)){
  #   out_bounds = which(prev_iter[,k] < a[,k] | prev_iter[,k] > b[,k])
  #   if(length(out_bounds) > 0){
  #     for(l in 1:length(out_bounds)){
  #       bounds = c(a[out_bounds[l],k],b[out_bounds[l],k])
  #       prev_iter[out_bounds[l], k] = bounds[which.min(abs(bounds - prev_iter[out_bounds[l],k]))]
  #     }
  #   }
  # }
  # (sum(block_sizes[-(j:cluster_number)]))
  # sum(block_sizes[1:j])
  for(i in 1:cluster_number){
    
    init_index = (sum(block_sizes[-(i:cluster_number)]))+1
    end_index = sum(block_sizes[1:i])
    data_cluster = data[init_index:end_index,]
    a_bounds = a[init_index:end_index,]
    b_bounds = b[init_index:end_index,]
    Sigma_cluster = Sigma_hat[init_index:end_index,init_index:end_index]
    for(j in 1:ncol(data)){
      if(iter_num == 1){
        #print(dim(data_cluster))
        if(block_sizes[i] == 1){
          draws = tmvtnorm::rtmvnorm(n = 1, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = mean(draws)
        }
        else{
          draws = tmvtnorm::rtmvnorm(n = 1, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = colMeans(draws)
        }
        if(any(is.nan(eps_draws[init_index:end_index,j]))){
          eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
        }
      }
      else{
        if(block_sizes[i] == 1){
          draws = tmvtnorm::rtmvnorm(n = 1, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = mean(draws)
        }
        else{
          draws = tmvtnorm::rtmvnorm(n = 1, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 10)
          eps_draws[init_index:end_index,j] = colMeans(draws)
        }
        if(any(is.nan(eps_draws[init_index:end_index,j]))){
          eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
        }
      }
      # if(any(is.nan(eps_draws[,j]))){
      #   eps_draws[init_index:end_index,j] = tmvtnorm::rtmvnorm(n = 10, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j])
      # }
    }
  }
  return(eps_draws)
}


XB_vals <- function(X, beta){
  #' Obtain hidden variable Z from Beta, X, and epsilon
  #' Epsilon should be generated from the expectation
  
  
  XB = matrix(0, n, p)
  
  XB[,1] <- rep(0,n)
  XB[,2] = X[,1]*beta[1,2]
  
  
  for(i in 3:p) {
    XB[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, beta[1:i-1,i], "*"))
  }
  dimnames(XB) <- list(NULL, as.character(1:p))
  return(XB)
}


t_normal_likelihood = function(testdata, Sigma_estimate, est_beta){
  total_ll = 0
  trunc_vals = obtain_trunc_vals(testdata, est_beta)
  a_bounds = obtain_bounds(testdata, trunc_vals)$a
  b_bounds = obtain_bounds(testdata, trunc_vals)$b
  for(i in 1:ncol(testdata)){
    total_ll = total_ll + log(pmvnorm(lower = a_bounds[,i], upper = b_bounds[,i], sigma = Sigma_estimate))
  }
  return(total_ll)
}

truncated_normal_likelihood = function(testdata, Sigma_estimate, est_beta){
  total_ll = 0
  trunc_vals = obtain_trunc_vals(testdata, est_beta)
  a_bounds = obtain_bounds(testdata, trunc_vals)$a
  b_bounds = obtain_bounds(testdata, trunc_vals)$b
  num_pairs = t(combn(1:nrow(testdata), 2))
  est_Sig = Sigma_estimate
  t_testdata = t(testdata)
  for(i in 1:nrow(num_pairs)){
    pair_data = t_testdata[,num_pairs[i,]]
    c1_run = trunc_vals[num_pairs[i,1],]
    c2_run = trunc_vals[num_pairs[i,2],]
    result = sum_log_lik_4likelihood(data = pair_data, rho_vec = est_Sig[num_pairs[i,1],num_pairs[i,2]], c1 = c1_run, c2 = c2_run)
    total_ll = total_ll + result
  }
  
  return(total_ll)
}

sum_log_lik_4likelihood = function(data, rho_vec, c1, c2){
    # if(i %% 20 == 0){
    #   print(i)
    # }
    llk = 0
    rho = rho_vec
    #print(rho)
    for(j in 1:nrow(data)){
      #print(j)
      #print(llk)
      if(all(data[j,] == c(1,1))){
        lower = c(c1[j], c2[j])
        upper = c(Inf, Inf)
        exactSolutionCoeff <- 1/(1-pnorm(lower[2]))
        
        integrand <- function(y) pnorm((lower[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        
        llk = llk + log((1 - exactSolutionInteg*exactSolutionCoeff)*(1-pnorm(lower[2], mean = 0, sd = 1)))
      }
      else if(all(data[j,] == c(1,0))){
        lower = c(c1[j], -Inf)
        upper = c(Inf, c2[j])
        
        exactSolutionCoeff <- 1/(pnorm(upper[2]))
        
        integrand <- function(y) pnorm((lower[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        if(is.nan(exactSolutionInteg) == T){
          print("here")
        }
        
        llk = llk + log((1-exactSolutionInteg*exactSolutionCoeff)*pnorm(upper[2], mean = 0, sd = 1)) # can get NaN here
      }
      else if(all(data[j,] == c(0,1))){
        lower = c(-Inf, c2[j])
        upper = c(c1[j], Inf)
        
        exactSolutionCoeff <- 1/(1-pnorm(lower[2]))
        
        integrand <- function(y) pnorm((upper[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        
        llk = llk + log((exactSolutionInteg*exactSolutionCoeff)*(1-pnorm(lower[2], mean = 0, sd = 1)))
      }
      else if(all(data[j,] == c(0,0))){
        lower = c(-Inf, -Inf)
        upper = c(c1[j], c2[j])
        exactSolutionCoeff <- 1/(pnorm(upper[2]))
        
        integrand <- function(y) pnorm((upper[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        llk = llk + log(exactSolutionInteg*exactSolutionCoeff*pnorm(upper[2], mean = 0, sd = 1))
      }
      #llk2 = llk2 + log(TruncatedNormal::pmvnorm(mu = c(0,0), sigma = sigma_rho, lb = lower, ub = upper, type = "mc")[1])
    }
  return(llk)
}

obtain_trunc_vals <- function(X, beta){
  #' Obtain truncation values from Beta, X, and epsilon
  #' Epsilon should be generated from the expectation
  n = nrow(X)
  p = ncol(X)
  trunc_vals = matrix(0, n, p)
  trunc_vals[,1] = rep(0,n)
  trunc_vals[,2] = -X[,1]*beta[1,2]
  for(i in 3:p) {
    trunc_vals[,i] = -rowSums(sweep(X[,1:i-1], MARGIN = 2, beta[1:i-1,i], "*"))
  }
  return(trunc_vals)
}


obtain_hidden_Z <- function(X, beta, epsilon){
  #' Obtain hidden variable Z from Beta, X, and epsilon
  #' Epsilon should be generated from the expectation
  n = nrow(X)
  p = ncol(X)
  Z_mat <- matrix(0, n, p)
  
  trunc_vals = matrix(0, n, p)
  
  Z_mat[,1] <- epsilon[,1]
  trunc_vals[,1] = rep(0,n)
  Z_mat[,2] = X[,1]*beta[1,2] + epsilon[,2]
  trunc_vals[,2] = -X[,1]*beta[1,2]
  
  
  for(i in 3:p) {
    Z_mat[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, beta[1:i-1,i], "*")) + epsilon[,i]
    trunc_vals[,i] = -rowSums(sweep(X[,1:i-1], MARGIN = 2, beta[1:i-1,i], "*"))
  }
  dimnames(Z_mat) <- list(NULL, as.character(1:p))
  return(Z_mat)
}


obtain_bounds = function(data, trunc_vals){
  bounds_a = matrix(data = 0, nrow = nrow(data), ncol = ncol(data))
  bounds_b = matrix(data = 0, nrow = nrow(data), ncol = ncol(data))
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      if(data[i,j] == 0){
        bounds_a[i,j] = -Inf
        bounds_b[i,j] = trunc_vals[i,j]
      }
      else{
        bounds_a[i,j] = trunc_vals[i,j]
        bounds_b[i,j] = Inf
      }
    }
  }
  return(list("a" = bounds_a, "b" = bounds_b))
}

probit_beta = function(data, init_beta){
  beta_update = matrix(data = 0, nrow = nrow(init_beta), ncol = ncol(init_beta))
  colnames(beta_update) = paste("X", 1:nrow(beta_update), sep = "")
  rownames(beta_update) = paste("X", 1:nrow(beta_update), sep = "")
  data = data.frame(data)
  for(i in 1:ncol(init_beta)){
    #print(i)
    if(length(which(init_beta[,i] != 0)) > 0){
      parents = which(init_beta[,i] != 0)
      X_vars = paste("X", parents, sep="")
      y_vars = paste("X", i, sep = "")
      #mod = glmnet(x = data.matrix(data[,parents]), y = data[,i], family = binomial(link = "probit"), alpha = 1, lambda = 0.1)
      
      form = as.formula(paste(y_vars, paste(X_vars, collapse= "+"), sep = "~"))
      model = glm(form, family = binomial(link = "probit"), 
                 data = data)
      # if(any(svd(model.matrix(model))$d < 0.01)){
      #   print("Singular values close to 0 for Probit")
      # }
      beta_estimate = summary(model)$coefficients
      #beta_estimate = coef(mod)
      beta_update[names(beta_estimate[,1])[-1],i] = beta_estimate[2:nrow(beta_estimate),1]
    }
  }
  rmse_error = sqrt(mean((beta_update - init_beta)^2))
  return(list(beta_update, rmse_error))
}


block_diag_sep_var = function(n, min_block_size, max_block_size, cov_struc){
  block_sizes = c()
  sig = matrix(0, nrow = n, ncol = n)
  rows_left = n
  i = 0
  while(rows_left > min_block_size){
    i = i + 1
    sampled_block_size = sample(min_block_size:max_block_size, size = 1)
    while(rows_left - sampled_block_size < min_block_size){
      sampled_block_size = sample(min_block_size:max_block_size, size = 1)
      if(rows_left - sampled_block_size == 0){
        rows_left = rows_left - sampled_block_size
        block_sizes = c(block_sizes, sampled_block_size)
        cov_struc_chosen = sample(1:length(cov_struc), 1)
        sig[(sum(block_sizes[-i])+1):sum(block_sizes),(sum(block_sizes[-i])+1):sum(block_sizes)] = cov_struc[[cov_struc_chosen]](n = sampled_block_size)
        break
      }
    }
    if(rows_left == 0){
      return(list(sig, block_sizes))
    }
    cov_struc_chosen = sample(1:length(cov_struc), 1)
    rows_left = rows_left - sampled_block_size
    block_sizes = c(block_sizes, sampled_block_size)
    sig[(sum(block_sizes[-i])+1):sum(block_sizes),(sum(block_sizes[-i])+1):sum(block_sizes)] = cov_struc[[cov_struc_chosen]](n = sampled_block_size)
  }
  # block_sizes = c(block_sizes, rows_left)
  # print("here")
  # sig[(sum(block_sizes[-i])+1):sum(block_sizes),(sum(block_sizes[-i])+1):sum(block_sizes)] = cov_struc(n = rows_left)
  
  # for(i in 1:cluster_size){
  #   sig[((i-1)*block_size + 1):(i*block_size),((i-1)*block_size + 1):(i*block_size)] = struc_matrix
  # }
  return(list(sig, block_sizes))
}


Sig_Estimate_DAG = function(X, trunc_vals, block_sizes){
  data = X
  c_sim = trunc_vals # [1,] for first row, etc.
  t_data = t(data)
  cluster_number = length(block_sizes)
  estimated_sig = diag(nrow(data))
  avg_corr = matrix(data = 0, nrow = max(block_sizes), ncol = max(block_sizes)) # dont know if I need or not
  # for(i in 1:cluster_size){
  #   sig[((i-1)*block_size + 1):(i*block_size),((i-1)*block_size + 1):(i*block_size)] = struc_matrix
  # }
  # result = sum_log_lik(data = t_data, rho_vec = rho, c1 = c1_sim, c2 = c2_sim)
  
  
  #num_pairs = t(combn(rows, 2))

  rho = seq(0,0.95, length = 10)
  for(j in 1:cluster_number){
    #print(j)
    if(block_sizes[j] == 1){
      next
    }
    num_pairs = t(combn(1:block_sizes[j], 2))
    for(i in 1:nrow(num_pairs)){
      #print(paste0("i=",i))
      #print(paste0("j=", j))
      #print(c((sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]))
      pair_data = t_data[,c((sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2])]
      c1_run = c_sim[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],]
      c2_run = c_sim[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2],]
      result = sum_log_lik(data = pair_data, rho_vec = rho, c1 = c1_run, c2 = c2_run)
      if(length(rho[which(result == max(result))]) == 0){
        estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]] = mean(avg_corr)
        estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1]] = mean(avg_corr)
      }
      else{
        estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]] = rho[which(result == max(result))]
        estimated_sig[(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1]] = rho[which(result == max(result))]
      }
    }
    if(length(which(is.na(estimated_sig[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]))) == 0){
      avg_corr[1:block_sizes[j], 1:block_sizes[j]] = avg_corr[1:block_sizes[j], 1:block_sizes[j]] + estimated_sig[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])]/cluster_number
    }
  }
  for(k in 1:cluster_number){
    if(length(which(is.na(estimated_sig[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])]))) > 0){
      estimated_sig[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] = avg_corr
    }
  }
  return(estimated_sig)
}

sum_log_lik = function(data, rho_vec, c1, c2){
  llk_vect = numeric(length(rho_vec))
  
  #llk_vect_old = numeric(length(rho_vec))
  for(i in 1:length(rho_vec)){
    # if(i %% 20 == 0){
    #   print(i)
    # }
    llk = 0
    #llk2 = 0
    sigma_rho = matrix(c(1,rho_vec[i], rho_vec[i], 1), nrow = 2, ncol = 2, byrow = TRUE)
    rho = rho_vec[i]
    #print(rho)
    for(j in 1:nrow(data)){
      #print(j)
      if(all(data[j,] == c(1,1))){
        lower = c(c1[j], c2[j])
        upper = c(Inf, Inf)
        exactSolutionCoeff <- 1/(1-pnorm(lower[2]))
        
        integrand <- function(y) pnorm((lower[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        
        llk = llk + log((1 - exactSolutionInteg*exactSolutionCoeff)*(1-pnorm(lower[2], mean = 0, sd = 1)))
      }
      else if(all(data[j,] == c(1,0))){
        lower = c(c1[j], -Inf)
        upper = c(Inf, c2[j])
        
        exactSolutionCoeff <- 1/(pnorm(upper[2]))
        
        integrand <- function(y) pnorm((lower[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        if(is.nan(exactSolutionInteg) == T){
          print("here")
        }
        
        llk = llk + log((1-exactSolutionInteg*exactSolutionCoeff)*pnorm(upper[2], mean = 0, sd = 1)) # can get NaN here
      }
      else if(all(data[j,] == c(0,1))){
        lower = c(-Inf, c2[j])
        upper = c(c1[j], Inf)
        
        exactSolutionCoeff <- 1/(1-pnorm(lower[2]))
        
        integrand <- function(y) pnorm((upper[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        
        llk = llk + log((exactSolutionInteg*exactSolutionCoeff)*(1-pnorm(lower[2], mean = 0, sd = 1)))
      }
      else if(all(data[j,] == c(0,0))){
        lower = c(-Inf, -Inf)
        upper = c(c1[j], c2[j])
        exactSolutionCoeff <- 1/(pnorm(upper[2]))
        
        integrand <- function(y) pnorm((upper[1] - rho * y) / sqrt(1 - rho * rho)) * dnorm(y)
        upper_bound = upper[2]
        lower_bound = lower[2]
        exactSolutionInteg <- integrate(integrand, lower_bound, upper_bound)$value
        llk = llk + log(exactSolutionInteg*exactSolutionCoeff*pnorm(upper[2], mean = 0, sd = 1))
      }
      #llk2 = llk2 + log(TruncatedNormal::pmvnorm(mu = c(0,0), sigma = sigma_rho, lb = lower, ub = upper, type = "mc")[1])
    }
    llk_vect[i] = llk
    #llk_vect_old[i] = llk2
  }
  return(llk_vect)
}


init_epsilon_creation = function(n, p, sig, omg.sq){
  eps_mat <- matrix(0, n, p)
  for(i in 1:p){
    eps_mat[,i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[1]*sig)
  }
  return(eps_mat)
}

block_diag_sep = function(n, block_size, struc_matrix){
  cluster_size = ceiling(n/block_size)
  sig = matrix(0, nrow = n, ncol = n)
  for(i in 1:cluster_size){
    sig[((i-1)*block_size + 1):(i*block_size),((i-1)*block_size + 1):(i*block_size)] = struc_matrix
  }
  return(sig)
}

beta_estimate_mb = function(data, mb){
  df_data = data.frame(data)
  beta = matrix(data = 0, nrow = ncol(data), ncol = ncol(data))
  colnames(beta) = paste("X", 1:nrow(beta), sep = "")
  rownames(beta) = paste("X", 1:nrow(beta), sep = "")
  for(i in 1:length(mb)){
    X_vars = paste("X", mb[[i]], sep = "")
    y_vars = paste("X", i, sep = "")
    if(length(mb[[i]]) == 0){
      next
    }
    #mod = glmnet(x = data.matrix(data[,-i]), y = data[,i], family = binomial(link = "probit"), alpha = 1, lambda = lbd)
    form = as.formula(paste(y_vars, paste(X_vars, collapse= "+"), sep = "~"))
    model = glm(form, family = binomial(link = "probit"), 
                data = df_data)
    #coef_mod = coef(mod)
    #mb[[i]] = rownames(coef_mod)[which(coef_mod != 0)][-1]
    sum_model = summary(model)$coefficients
    beta[names(sum_model[,1])[-1],i] = sum_model[2:nrow(sum_model),1]
  }
  return(beta)
}


init_beta_0_mb = function(data, mb){
  df_data = data.frame(data)
  beta = matrix(data = 0, nrow = ncol(data), ncol = ncol(data))
  colnames(beta) = paste("X", 1:nrow(beta), sep = "")
  rownames(beta) = paste("X", 1:nrow(beta), sep = "")
  for(i in 1:length(mb)){
    if(length(mb[[i]]) == 0){
      next
    }
    beta[as.numeric(mb[[i]]),i] = 0.0001
  }
  return(beta)
}

init_beta_0_lasso = function(data, mb){
  df_data = data.frame(data)
  beta = matrix(data = 0, nrow = ncol(data), ncol = ncol(data))
  colnames(beta) = paste("X", 1:nrow(beta), sep = "")
  rownames(beta) = paste("X", 1:nrow(beta), sep = "")
  for(i in 1:length(mb)){
    if(length(mb[[i]]) == 0){
      next
    }
    beta[as.numeric(mb[[i]]),i] = 0.0001
  }
  return(beta)
}

toeplitz_struc = function(n, par_toep_min = 0.1, par_toep_max = 0.25){
  mat = matrix(0, nrow = n, ncol = n)
  par_toep = runif(1, par_toep_min, par_toep_max)
  for(i in 1:n){
    for(j in 1:n){
      mat[i,j] = par_toep^(abs(i - j)/5)
    }
  }
  return(mat)
}

const_struc = function(n, par_const_min = 0.4, par_const_max = 0.7){
  mat = matrix(0, nrow = n, ncol = n)
  par_const = runif(1, par_const_min, par_const_max)
  for(i in 1:n){
    for(j in 1:n){
      if(i == j){
        mat[i,j] = 1
      }
      else{
        mat[i,j] = par_const
      }
    }
  }
  return(mat)
}

neighborhood_selection = function(data, lbd){
  mb = list()
  for(i in 1:ncol(data)){
    mod = glmnet(x = data.matrix(data[,-i]), y = data[,i], family = binomial(link = "probit"), alpha = 1, lambda = lbd)
    coef_mod = coef(mod)
    mb[[i]] = rownames(coef_mod)[which(coef_mod != 0)][-1]
  }
  return(mb)
}

neighborhood_selection_parents_lasso = function(data, lbd){
  mb = list()
  for(i in 1:ncol(data)){
    mod = glmnet(x = data.matrix(data[,-i]), y = data[,i], family = binomial(link = "probit"), alpha = 1, lambda = lbd)
    coef_mod = coef(mod)
    mb[[i]] = rownames(coef_mod)[which(coef_mod != 0)][-1]
  }
  return(mb)
}

get_moral_edge = function(trueB, i){
  children = which(trueB[,i] != 0)
  parents = which(trueB[i,] != 0)
  spouse = c()
  for(j in 1:length(parents)){
    spouse = c(spouse, which(trueB[,parents[j]] != 0))
  }
  spouse = unique(spouse)
  spouse = spouse[-which(spouse == i)]
  edges = unique(c(children, parents, spouse))
  return(edges)
}

get_neighbors = function(trueB, i){
  children = which(trueB[,i] != 0)
  parents = which(trueB[i,] != 0)
  return(unique(c(children, parents)))
}

get_parents = function(trueB, i){
  parents = which(trueB[,i] != 0)
  return(parents)
}

get_children= function(trueB, i){
  children = which(trueB[i,] != 0)
  return(children)
}

get_mb = function(trueB, i){
  trueB[which(trueB != 0)] = 1
  ## Need to write own markov blanket code
  parents = which(trueB[,i] != 0)
  children = which(trueB[i,] != 0)
  spouses = c()
  if(length(children) != 0){
    for(j in 1:length(children)){
      spouses = c(spouses, which(trueB[,children[j]] != 0))
    }
  }
  markov_blanket = unique(c(parents, children, spouses))
  markov_blanket = markov_blanket[which(markov_blanket != i)]
  return(markov_blanket)
}

get_F1_mb = function(mb, trueB){
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for(i in 1:length(mb)){
    #print(i)
    children = as.numeric(mb[[i]])
    true = get_mb(trueB, i)
    if(length(true) == 0){
      next
    }
    true_pos = true_pos + length(which(true %in% children))
    #print(paste("% of true positive in node", i , length(which(true %in% children))/length(true)))
    false_pos = false_pos + length(children) - length(which(true %in% children))
    false_neg = false_neg + length(true) - length(which(true %in% children))
  }
  precision = true_pos/(true_pos + false_pos)
  recall = true_pos/(true_pos + false_neg)
  f1 = 2*precision*recall/(precision+recall)
  return(list("precision" = precision, "recall" = recall, "f1_score" = f1))
}


get_F1 = function(mb, trueB){
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for(i in 1:length(mb)){
    children = as.numeric(mb[[i]])
    true = get_neighbors(trueB, i)
    if(length(true) == 0){
      next
    }
    true_pos = true_pos + length(which(true %in% children))
    #print(paste("% of true positive in node", i , length(which(true %in% children))/length(true)))
    false_pos = false_pos + length(children) - length(which(true %in% children))
    false_neg = false_neg + length(true) - length(which(true %in% children))
  }
  precision = true_pos/(true_pos + false_pos)
  recall = true_pos/(true_pos + false_neg)
  f1 = 2*precision*recall/(precision+recall)
  return(list("precision" = precision, "recall" = recall, "f1_score" = f1))
}

get_F1_dag = function(mb, trueB){
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for(i in 1:length(mb)){
    children = as.numeric(mb[i,])
    true = get_children(trueB, i)
    if(length(true) == 0){
      next
    }
    true_pos = true_pos + length(which(true %in% children))
    #print(paste("% of true positive in node", i , length(which(true %in% children))/length(true)))
    false_pos = false_pos + length(children) - length(which(true %in% children))
    false_neg = false_neg + length(true) - length(which(true %in% children))
  }
  precision = true_pos/(true_pos + false_pos)
  recall = true_pos/(true_pos + false_neg)
  f1 = 2*precision*recall/(precision+recall)
  return(list("precision" = precision, "recall" = recall, "f1_score" = f1))
}

# the one using on 11/11
get_F1_dag_parents = function(mb, trueB){
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for(i in 1:length(mb)){
    parents = as.numeric(mb[[i]]$edges)
    true = as.numeric(trueB[[i]]$edges)
    if(length(true) == 0){
      next
    }
    true_pos = true_pos + length(which(true %in% parents))
    #print(paste("% of true positive in node", i , length(which(true %in% parents))/length(true)))
    false_pos = false_pos + length(parents) - length(which(true %in% parents))
    false_neg = false_neg + length(true) - length(which(true %in% parents))
  }
  precision = true_pos/(true_pos + false_pos)
  recall = true_pos/(true_pos + false_neg)
  f1 = 2*precision*recall/(precision+recall)
  return(list("precision" = precision, "recall" = recall, "f1_score" = f1))
}

get_F1_dag_parents_vote = function(vote_mat, trueB){
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for(i in 1:nrow(vote_mat)){
    parents = which(vote_mat[,i] != 0)
    true = get_parents(trueB, i)
    if(length(true) == 0){
      next
    }
    true_pos = true_pos + length(which(true %in% parents))
    #print(paste("% of true positive in node", i , length(which(true %in% parents))/length(true)))
    false_pos = false_pos + length(parents) - length(which(true %in% parents))
    false_neg = false_neg + length(true) - length(which(true %in% parents))
  }
  precision = true_pos/(true_pos + false_pos)
  recall = true_pos/(true_pos + false_neg)
  f1 = 2*precision*recall/(precision+recall)
  return(list("precision" = precision, "recall" = recall, "f1_score" = f1))
}


neighborhood_selection_mxm = function(data, test){
  mb = list()
  if(ncol(data) <= 100){
    max_k_num = 3
  }
  else{
    max_k_num = floor(nrow(data)/10)
  }
  for(i in 1:ncol(data)){
    mod = MXM::MMPC(target = data[,i], dataset = data[,-i], max_k = max_k_num, test = test, backward = T, threshold = 0.05)
    selected = mod@selectedVars
    selected[which(selected >= i)] = selected[which(selected >= i)]+1
    mb[[i]] = selected
  }
  return(mb)
}

neighborhood_selection_mmhc_mmpc = function(data){
  mb = list()
  skel = mmhc.skel(dataset = data, max_k = 7, test = "gSquare")$G
  for(i in 1:ncol(data)){
    mb[[i]] = which(skel[i,] != 0)
  }
  return(mb)
}

neighborhood_selection_mmhc_def = function(data){
  mb = list()
  skel = mmhc.skel(dataset = data, max_k = 7)$G
  for(i in 1:ncol(data)){
    mb[[i]] = which(skel[i,] != 0)
  }
  return(mb)
}

neighborhood_selection_pc = function(data){
  mb = list()
  skel = pcalg::pc(suffStat = list(C = cor(data), n = nrow(data)), indepTest = gaussCItest, alpha = 0.05, labels = colnames(data))@graph@edgeL
  for(i in 1:ncol(data)){
    mb[[i]] = skel[[i]]$edges
  }
  return(mb)
}

neighborhood_selection_mmhc = function(data){
  mb = list()
  skel = MXM::mmhc.skel(dataset = data, max_k = 7, test = "gSquare")$G
  for(i in 1:ncol(data)){
    mb[[i]] = which(skel[i,] != 0)
  }
  return(mb)
}

neighborhood_selection_mmhc_dag = function(data){
  mb = list()
  skel = pchc::mmhc(x = data, method = "cat", score = "bde")$dag$nodes
  for(i in 1:ncol(data)){
    mb[[i]] = skel[[i]]$mb
  }
  return(mb)
}

neighborhood_selection_mmhc_get_dag = function(data){ # to get parents
  parents = list()
  skel = pchc::mmhc(x = data, method = "cat", max_k = 7, score = "bde", alpha = 0.05)$dag$nodes
  for(i in 1:ncol(data)){
    parents[[i]] = skel[[i]]$parents
  }
  return(parents)
}

init_data_dag = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::mmhc(x = data, method = "cat", max_k = 7, score = "bde", alpha = 0.05)$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  return(parents)
}

init_data_dag_pchc = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::pchc(x = data, method = "cat", score = "bde")$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  return(parents)
}

init_data_dag_iamb = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  score = new("GaussL0penObsScore", data)
  skel = pcalg::ges(score)$essgraph$.in.edges
  for(i in 1:ncol(data)){
    #parent_list = sub('.', '', skel[[i]]$parents)
    parents[skel[[i]],i] = 1
  }
  return(parents)
}

init_data_dag_lasso = function(data){
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  n = nrow(data)
  p = ncol(data)
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    mod = glmnet(x = data.matrix(data[,-i]), y = data[,i], family = binomial(link = "probit"), alpha = 1, lambda = 0.2-0.0002*n+0.0000475*p)
    coef_mod = coef(mod)
    parents[rownames(coef_mod)[which(coef_mod != 0)][-1],i] = coef_mod[which(coef_mod != 0)][-1]
  }
  return(parents)
}

init_data_dag_pc = function(data){ # to get parents
  n <- nrow (data)
  V <- colnames(data) # labels aka node names
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = cor(data), n = n),
               indepTest = gaussCItest, ## indep.test: partial correlations
               alpha=0.05, labels = V)@graph@edgeL
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    parents[i, pc.fit[[i]]$edges] = 1
  }
  # parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  # colnames(parents) = rownames(parents) = 1:ncol(data)
  # for(i in 1:ncol(data)){
  #   parents[i,pc.D[[i]]$edges] = 1
  # }
  return(parents)
}

neighborhood_selection_dag_lasso = function(data){
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  n = nrow(data)
  p = ncol(data)
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    mod = glmnet(x = data.matrix(data[,-i]), y = data[,i], family = c("gaussian"), alpha = 1, lambda = 0.2-0.0002*n+0.0000475*p)
    coef_mod = coef(mod)
    parents[rownames(coef_mod)[which(coef_mod != 0)][-1],i] = 1
  }
  return(parents)
}



neighborhood_selection_pc_get_dag_decor = function(data){ # to get parents
  n <- nrow (data)
  V <- colnames(data) # labels aka node names
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = cor(data), n = n),
               indepTest = gaussCItest, ## indep.test: partial correlations
               alpha=0.05, labels = V)@graph@edgeL
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    parents[i, pc.fit[[i]]$edges] = 1
  }
  return(parents)
}


###### functions for copula PC ######

gaussCItestLocal <- function (x, y, S, suffStat) 
{
  # Using Fisher's z-transformation of the partial correlation, test for zero partial correlation. 
  #
  # Args:
  #   see function gaussCItest() in Package 'pcalg';
  #   'suffStat' contains an addtional element 'ESS.Mat', that is, the matrix of effective sample size. 
  #
  # Returns:
  #   the p-value of the current test.
  
  subMat <- suffStat$ESS.Mat[c(x,y,S),c(x,y,S)]
  ne <- mean(subMat[upper.tri(subMat)], na.rm = T)
  z <- zStat(x, y, S, C = suffStat$C, n = ne)
  2 * pnorm(abs(z), lower.tail = FALSE)
}

inferCopulaModel <- function (Y, n0 = dim(Y)[2] + 1, S0 = diag(dim(Y)[2])/n0, nsamp = 100, 
                              odens = max(1, round(nsamp/1000)), impute = any(is.na(Y)), 
                              plugin.threshold = 20, plugin.marginal = (apply(Y, 2, function(x) {
                                length(unique(x))
                              }) > plugin.threshold), verb = TRUE) 
{
  # Main function to perform inference for Gaussian copula models.
  #
  # Args: 
  #   Y, a n by p data matrix (n-observations, p-variables).
  #
  # Return:
  #   C.psamp, a collection of samples of the underlying correlation matrix.
  #
  # For details about the arguements and outputs, refer to function 'sbgcop.mcmc' in R package 'sbgcop',
  #   https://cran.r-project.org/web/packages/sbgcop/index.html.
  #
  # Author: Ruifei Cui
  
  require(sbgcop)
  ok_S0 <- all(eigen(S0)$val > 0) & dim(S0)[1] == dim(Y)[2] & 
    dim(S0)[2] == dim(Y)[2]
  ok_n0 <- (n0 >= 0)
  if (!ok_S0) {
    stop("Error: S0 must be a positive definite p x p matrix \n")
  }
  if (!ok_n0) {
    stop("Error: n0 must be positive \n")
  }
  
  vnames <- colnames(Y)
  Y <- as.matrix(Y)
  colnames(Y) <- vnames
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "average", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  S <- cov(Z)
  Y.pmean <- Y
  if (impute) {
    Y.pmean <- matrix(0, nrow = n, ncol = p)
  }
  LPC <- NULL
  C.psamp <- array(dim = c(p, p, floor(nsamp/odens)))
  Y.imp <- NULL
  if (impute) {
    Y.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  }
  dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 
                            1:floor(nsamp/odens))
  for (ns in 1:nsamp) {
    for (j in sample(1:p)) {
      Sjc <- S[j, -j] %*% solve(S[-j, -j]+0.01*diag(dim(S[-j,-j])[1]))
      sdj <- sqrt(S[j, j] - S[j, -j] %*% solve(S[-j, 
                                                 -j]+0.01*diag(dim(S[-j,-j])[1])) %*% S[-j, j])
      muj <- Z[, -j] %*% t(Sjc)
      if (!plugin.marginal[j]) {
        for (r in 1:Rlevels[j]) {
          ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
          lb <- suppressWarnings(max(Z[R[, j] == r - 
                                         1, j], na.rm = TRUE))
          ub <- suppressWarnings(min(Z[R[, j] == r + 
                                         1, j], na.rm = TRUE))
          Z[ir, j] <- qnorm(runif(length(ir), pnorm(lb, 
                                                    muj[ir], sdj), pnorm(ub, muj[ir], sdj)), 
                            muj[ir], sdj)
        }
      }
      ir <- (1:n)[is.na(R[, j])]
      Z[ir, j] <- rnorm(length(ir), muj[ir], sdj)
    }
    # relocate the mean to zero
    # added by Ruifei Cui
    Z = t( (t(Z)-apply(Z,2,mean)))
    
    S <- solve(rwish(solve(S0 * n0 + t(Z) %*% Z), n0 + 
                       n))
    if (ns%%odens == 0) {
      C <- S/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      lpc <- ldmvnorm(Z %*% diag(1/sqrt(diag(S))), 
                      C)
      LPC <- c(LPC, lpc)
      C.psamp[, , ns/odens] <- C
      if (impute) {
        Y.imp.s <- Y
        for (j in 1:p) {
          Y.imp.s[is.na(Y[, j]), j] <- quantile(Y[, 
                                                  j], pnorm(Z[is.na(Y[, j]), j], 0, sqrt(S[j, 
                                                                                           j])), na.rm = TRUE, type = 1)
        }
        Y.imp[, , ns/odens] <- Y.imp.s
        Y.pmean <- ((ns/odens - 1)/(ns/odens)) * Y.pmean + 
          (1/(ns/odens)) * Y.imp.s
      }
    }
    if (verb == TRUE & (ns%%(odens * 10)) == 0) {
      cat(round(100 * ns/nsamp), "percent done ", 
          date(), "\n")
    }
  }
  G.ps <- list(C.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp, 
               LPC = LPC)
  class(G.ps) <- "psgc"
  return(G.ps)
}



init_data_dag_copPC = function(data){ # to get parents
  n <- nrow (data)
  p = ncol(data)
  
  cop.obj <- inferCopulaModel(data, verb = T)
  # correlation matrix samples
  C_samples <- cop.obj$C.psamp[,,1:100]
  # average correlation matrix
  corr.cop <- apply(C_samples, c(1,2), mean)
  # local effective sample size
  less.cop <- ((1-corr.cop^2)^2)/apply(C_samples,c(1,2), var)
  
  n <- nrow (data)
  var.names <- colnames(data) # labels aka node names
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = corr.cop, n = n), 
                       indepTest = gaussCItest, labels = var.names, alpha = .05, conservative = T)@graph@edgeL
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    parents[i, pc.fit[[i]]$edges] = 1
  }
  # parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  # colnames(parents) = rownames(parents) = 1:ncol(data)
  # for(i in 1:ncol(data)){
  #   parents[i,pc.D[[i]]$edges] = 1
  # }
  return(parents)
}

neighborhood_selection_copPC_get_dag_decor = function(data){ # to get parents
  n <- nrow (data)
  p = ncol(data)
  var.names = colnames(data)
  cop.obj <- inferCopulaModel(data, verb = T)
  # correlation matrix samples
  C_samples <- cop.obj$C.psamp[,,1:100]
  # average correlation matrix
  corr.cop <- apply(C_samples, c(1,2), mean)
  # local effective sample size
  less.cop <- ((1-corr.cop^2)^2)/apply(C_samples,c(1,2), var)
  
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = corr.cop, n = n, ESS.Mat = less.cop), 
               indepTest = gaussCItestLocal, labels = var.names, alpha = .05, conservative = T)@graph@edgeL
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    parents[i, pc.fit[[i]]$edges] = 1
  }
  return(parents)
}


neighborhood_selection_tabu_get_dag_decor = function(data){ # to get parents
  n <- nrow (data)
  V <- colnames(data) # labels aka node names
  ## estimate CPDAG
  data = as.data.frame(data)
  
  tabu.fit <- tabu(data)$nodes
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = paste0(1:ncol(data))
  for(i in colnames(data)){
    parents[i, tabu.fit[[i]]$children] = 1
  }
  return(parents)
}


init_tabu_get_dag = function(data){ # to get parents
  n <- nrow (data)
  V <- colnames(data) # labels aka node names
  ## estimate CPDAG
  data = as.data.frame(data)
  data[] <- lapply(data, function(x) {
    if (is.numeric(x)) {
      as.factor(x)
    } else {
      x
    }
  })
  tabu.fit <- tabu(data, score = "bic")$nodes
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = paste0(1:ncol(data))
  for(i in colnames(data)){
    parents[i, tabu.fit[[i]]$children] = 1
  }
  return(parents)
}


###
init_iamb_get_dag = function(data){ # to get parents
  n <- nrow (data)
  V <- colnames(data) # labels aka node names
  ## estimate CPDAG
  data = as.data.frame(data)
  data[] <- lapply(data, function(x) {
    if (is.numeric(x)) {
      as.factor(x)
    } else {
      x
    }
  })
  iamb.fit <- bnlearn::iamb(x = data, test = "mi")$nodes
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = paste0(1:ncol(data))
  for(i in colnames(data)){
    if(length(iamb.fit[[i]]$children) > 0){
      parents[i, iamb.fit[[i]]$children] = 1
    }
    else{
      parents[i, iamb.fit[[i]]$nbr] = 1
    }
  }
  return(parents)
}

neighborhood_selection_iamb_get_dag_decor = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  data = as.data.frame(data)
  colnames(parents) = rownames(parents) = paste0(1:ncol(data))
  #skel = pchc::mmhc(x = data, max_k = 7, alpha = 0.05)$dag$nodes
  iamb.fit = bnlearn::iamb(x = data, test = "cor")$nodes
  for(i in 1:ncol(data)){
    #parent_list = sub('.', '', skel[[i]]$parents)
    if(length(iamb.fit[[i]]$children) > 0){
      parents[i, iamb.fit[[i]]$children] = 1
    }
    else{
      parents[i, iamb.fit[[i]]$nbr] = 1
    }
  }
  return(parents)
}


neighborhood_selection_mmhc_get_dag_decor = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::mmhc(x = data, max_k = 7, alpha = 0.05)$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  return(parents)
}

neighborhood_selection_pchc_get_dag_decor = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::pchc(x = data)$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  return(parents)
}

neighborhood_selection_ges_dag = function(data){
  mb = list()
  score <- new("GaussL0penObsScore", data)
  skel = pcalg::ges(score)$essgraph$.in.edges
  for(i in 1:ncol(data)){
    mb[[i]] = skel[[i]]
  }
  return(mb)
}