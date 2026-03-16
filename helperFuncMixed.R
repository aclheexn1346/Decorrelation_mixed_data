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

make_discrete_dense <- function(x, threshold = 4) {
  if (is.data.frame(x)) {
    # Process data frames
    x[] <- lapply(x, function(col) {
      nunique <- length(unique(col))
      if (nunique <= threshold) {
        col <- factor(col, levels = sort(unique(col)), ordered = TRUE)
        col <- as.integer(col) - 1
      }
      return(col)
    })
    return(as.data.frame(x))
    
  } else if (is.matrix(x)) {
    # Process matrices
    m <- apply(x, 2, function(col) {
      nunique <- length(unique(col))
      if (nunique <= threshold) {
        col <- factor(col, levels = sort(unique(col)), ordered = TRUE)
        col <- as.integer(col) - 1
      }
      return(col)
    })
    # ensure result is matrix, not character
    return(matrix(as.numeric(m), nrow = nrow(x), ncol = ncol(x),
                  dimnames = dimnames(x)))
  } else {
    stop("Input must be a data.frame or matrix")
  }
}


transform_vector <- function(vec) {
  # Check if the vector contains only 1s and 2s
  if (all(vec %in% c(1, 2))) {
    vec[vec == 1] <- 0
    vec[vec == 2] <- 1
  } else if (all(vec %in% c(0, 2))) {
    # If the vector contains only 0s and 2s, transform 2 -> 1
    vec[vec == 2] <- 1
  }
  return(vec)
}

sim_X_LUM_mixed = function(vers, p, n, omg.sq, sig, b, cutoffs, disc_prob = 0.5){
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
  
  #trunc_vals = matrix(0, n, p)
  cont_columns = rbinom(n = p, prob = disc_prob, size = 1)
  Z_mat[,1] <- eps_mat_1[,1]
  if(cont_columns[1] == 1){
    X[,1] <- as.numeric(cut(eps_mat_1[,1], cutoffs)) - 1 # cutoff
    if(length(levels(factor(X[,1]))) == (length(cutoffs) - 2)){
      X[,1] = transform_vector(X[,1])
    }
  }
  else{
    X[,1] = Z_mat[,1]
  }

  Z_mat[,2] = X[,1]*b[1,2] + eps_mat_1[,2]
  if(cont_columns[2] == 1){
    X[,2] <- as.numeric(cut(X[,1]*b[1,2] + eps_mat_1[,2], cutoffs)) - 1
    if(length(levels(factor(X[,2]))) == (length(cutoffs) - 2)){
      X[,2] = transform_vector(X[,2])
    }
  }
  else{
    X[,2] = Z_mat[,2]
  }
  
  
  for(i in 3:p) {
    eps_mat_1[, i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[i]*sig) # Sigma is covariance between the n samples
    Z_mat[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, b[1:i-1,i], "*")) + eps_mat_1[,i]
    if(cont_columns[i] == 1){
      X[,i] <- as.numeric(cut(Z_mat[,i], cutoffs)) - 1
      if(length(levels(factor(X[,i]))) == (length(cutoffs) - 2)){
        X[,i] = transform_vector(X[,i])
      }
    }
    else{
      X[,i] = Z_mat[,i]
    }

  }
  
  dimnames(X) <- list(NULL, as.character(1:p))
  return(list(X = X, eps_mat = eps_mat_1, Z = Z_mat, discCols = cont_columns))
}

tau_likelihood_func_mixed = function(delta, data, beta, discCols, variable){ # I think can use this for likelihood
  log_lik = 0
  tau = c(-Inf, cumsum(delta), Inf)
  if(discCols[variable] == 1){
    for(i in 1:nrow(data)){
      #print(paste0(i, " : ",log_lik))
      z_mean = as.numeric(data[i,] %*% beta[,variable])
      upper_tau = tau[data[i,variable]+2]
      lower_tau = tau[data[i,variable]+1]
      upper_lik = pnorm(upper_tau - z_mean)
      lower_lik = pnorm(lower_tau - z_mean)
      log_lik = log_lik + log(upper_lik - lower_lik)
    }
  }
  else{
    for(i in 1:nrow(data)){
      z_mean = as.numeric(data[i,] %*% beta[,variable])
      log_lik = log_lik + log(dnorm(data[i, variable], mean = z_mean, sd = 1))
    }
  }
  return(log_lik)
}

optimize_log_likelihood_mixed <- function(delta, data, beta, discCols, variable) {
  return(-tau_likelihood_func_mixed(delta, data, beta, discCols, variable))  # optim minimizes, so negate log-likelihood
}

make_discrete_dense <- function(x, threshold = 10) {
  if (is.data.frame(x)) {
    # Process data frames
    x[] <- lapply(x, function(col) {
      nunique <- length(unique(col))
      if (nunique <= threshold) {
        col <- factor(col, levels = sort(unique(col)), ordered = TRUE)
        col <- as.integer(col) - 1
      }
      return(col)
    })
    return(as.data.frame(x))
    
  } else if (is.matrix(x)) {
    # Process matrices
    m <- apply(x, 2, function(col) {
      nunique <- length(unique(col))
      if (nunique <= threshold) {
        col <- factor(col, levels = sort(unique(col)), ordered = TRUE)
        col <- as.integer(col) - 1
      }
      return(col)
    })
    # ensure result is matrix, not character
    return(matrix(as.numeric(m), nrow = nrow(x), ncol = ncol(x),
                  dimnames = dimnames(x)))
  } else {
    stop("Input must be a data.frame or matrix")
  }
}


EM_beta_given_tau_mixed = function(data, beta, tau, discCols, reps){
  beta_list = list(beta)
  n = nrow(data)
  p = ncol(data)
  diffs_list = list()
  for(rep in 1:reps){
    mu = data %*% beta
    cont_data = matrix(0, nrow = nrow(data), ncol = ncol(data))
    # Z - truncated normal expectation
    for(i in 1:n){
      for(j in 1:p){
        if(discCols[j] == 1){
          boundaries = c(-Inf, tau[[j]], Inf)
          cont_data[i,j] <- truncnorm::etruncnorm(a = boundaries[data[i,j]+1], b = boundaries[data[i,j]+2], mean = mu[i,j])
        }
        else{
          cont_data[i,j] = data[i,j]
        }
      }
    }
    # B - linear regress Maximization
    for(k in 1:p){
      #print(k)
      parents = which(beta[,k] != 0)
      if(length(parents) > 0){
        if(length(parents) > 1){
          model = glmnet(x = data.matrix(data[,parents]), y = cont_data[,k], alpha = 0, lambda = n*p / (10*(n + p)) + n*p / (sqrt(n * p)*10))
          beta_estimate = coef(model)
          beta[parents,k] = beta_estimate[2:nrow(beta_estimate),1]
          
        }
        else{
          #model = lm(cont_data[,k]~data[,parents])
          #beta_estimate = summary(model)$coefficients
          #beta[parents,k] = beta_estimate[2:nrow(beta_estimate),1]
          lambda = n*p / (10*(n + p)) + n*p / (sqrt(n * p)*10)
          y = cont_data[,k] - mean(cont_data[,k])
          x = data[,parents] - mean(data[,parents])
          #print(paste0("Condition number: ", kappa(t(x) %*% x)))
          numerator = sum(x*y)
          denom = sum(x^2) + lambda
          coeff = numerator/denom
          beta[parents,k] = coeff
        }
      }
    }
    # From here Feb 1: 5pm
    beta_list[[rep+1]] = beta
    #print(summary(beta[beta != 0]))
    diff = beta - beta_list[[rep]]
    #print(paste0("beta RMSE from previous iter: ",sqrt(mean(diff[diff != 0]^2))))
    diffs_list[[rep]] = sqrt(mean(diff[diff != 0]^2))
  }
  return(list(beta, diffs_list))
}

data_generation_mixed = function(n, p, cov_struc, reps, vers, thresholds){
  hold_data = list()
  hold_parents = list()
  hold_beta = list()
  hold_sigma = list()
  hold_block_size = list()
  for(rep in 1:reps){
    vers = vers + 1
    omg.sq = rep(1, p)
    block_obj = block_diag_sep_var(n = n, min_block_size = 5, max_block_size = 10, cov_struc = cov_struc)
    sig_2 = block_obj[[1]]
    block_sizes = block_obj[[2]]
    data_prop = rep(1,p)
    while(any(data_prop == 1)){
      vers = vers + 1
      B = gen.B(p = p, seed = vers)
      sim_result_nplarge = sim_X_LUM_mixed(vers, p, n, omg.sq, sig = sig_2, b = B$b, cutoffs = c(-Inf,thresholds,Inf))
      data = sim_result_nplarge$X
      for(j in 1:ncol(data)){
        data_prop[j] = length(unique(data[,j]))
      }
    }
    
    beta = B$b
    data = sim_result_nplarge$X
    init_beta = beta
    init_beta[init_beta != 0] = 1e-4
    hold_data[[rep]] = data
    hold_parents[[rep]] = init_beta
    hold_beta[[rep]] = B$b
    hold_sigma[[rep]] = sig_2
    hold_block_size[[rep]] = block_sizes
  }
  return(list("Data" = hold_data, "parents_mat" = hold_parents, "true_beta" = hold_beta, "true_sigma" = hold_sigma, "block_sizes" = hold_block_size, "discCols" = sim_result_nplarge$discCols))
}

sample_uniform <- function() {
  if (runif(1) < 0.5) {
    runif(1, -0.9, -0.6)
  } else {
    runif(1, 0.6, 0.9)
  }
}


genB_from_realDAG = function(B_mat){
  B_mat = apply(B_mat, c(1,2), function(x) {
    if (x == 1) sample_uniform() else 0
  })
  rownames(B_mat) = colnames(B_mat) = 1:nrow(B_mat)
  return(B_mat)
}

data_generation_mixed_realDAG = function(n, cov_struc, reps, vers, thresholds, Bmatrix){
  hold_data = list()
  hold_parents = list()
  hold_beta = list()
  hold_sigma = list()
  hold_block_size = list()
  p = ncol(Bmatrix)
  for(rep in 1:reps){
    vers = vers + 1
    omg.sq = rep(1, p)
    block_obj = block_diag_sep_var(n = n, min_block_size = 5, max_block_size = 15, cov_struc = cov_struc)
    sig_2 = block_obj[[1]]
    block_sizes = block_obj[[2]]
    data_prop = rep(1,p)
    while(any(data_prop == 1)){
      vers = vers + 1
      B = Bmatrix
      sim_result_nplarge = sim_X_LUM_mixed(vers, p, n, omg.sq, sig = sig_2, b = B, cutoffs = c(-Inf,thresholds,Inf))
      data = sim_result_nplarge$X
      for(j in 1:ncol(data)){
        data_prop[j] = length(unique(data[,j]))
      }
    }
    
    beta = B
    data = sim_result_nplarge$X
    init_beta = beta
    init_beta[init_beta != 0] = 1e-4
    hold_data[[rep]] = data
    hold_parents[[rep]] = init_beta
    hold_beta[[rep]] = B
    hold_sigma[[rep]] = sig_2
    hold_block_size[[rep]] = block_sizes
  }
  return(list("Data" = hold_data, "parents_mat" = hold_parents, "true_beta" = hold_beta, "true_sigma" = hold_sigma, "block_sizes" = hold_block_size, "discCols" = sim_result_nplarge$discCols))
}


EM_tau_beta_mixed = function(data, parent_set, discColumns, reps){
  beta_estimates = list()
  tau_estimates = list()
  n = nrow(data)
  p = ncol(data)
  beta_diffs = list()
  init_beta = parent_set
  #### Initial tau estimate with quantiles ####
  init_tau = vector("list", ncol(data))
  likelihood = numeric(reps)
  for(k in 1:p){
    if(discColumns[k] == 1){
      quants = table(data[,k])/length(data[,k])
      quants = quants[-length(quants)]
      init_tau[[k]] = qnorm(cumsum(quants))
    }
  }
  tau_estimates[[1]] = init_tau
  beta_estimates[[1]] = init_beta
  store_likelihood = matrix(data = 0, nrow = reps+1, ncol = p)
  ##############################################
  for(i in 1:reps){
    #print(i)
    # update beta given tau
    update_beta = EM_beta_given_tau_mixed(data, init_beta, init_tau, discCols = discColumns, reps = 5)
    beta_estimates[[i+1]] = update_beta[[1]]
    beta_diffs[[i]] = update_beta[[2]]
    # update tau given beta
    update_tau = vector("list", ncol(data))
    
    for(j in 1:p){
      #print(j)
      if (i == 1){
        store_likelihood[i,j] = tau_likelihood_func_mixed(diff(c(0, init_tau[[j]])), data = data, beta = update_beta[[1]], discCols = discColumns, variable = j)
      }

      # Generate data using the current vers value
      if(discColumns[j] == 1){
        if(length(unique(data[,j])) != 3){
          data[,j] = match(data[,j], sort(unique(data[,j]))) - 1
        }
        result <- optim(
          par = diff(c(0,init_tau[[j]])),
          fn = optimize_log_likelihood_mixed,
          data = data,
          beta = update_beta[[1]],
          discCols = discColumns,
          variable = j,
          method = "L-BFGS-B",
          lower = c(-Inf, rep(1e-5, length(diff(c(0,init_tau[[j]])))-1))  # Ensure increments are positive
        )
        optimized_delta <- result$par
        update_tau[[j]] = cumsum(optimized_delta)
      }

      # Extract optimized deltas and construct tau
      if(discColumns[j] == 1){
        store_likelihood[i+1,j] = tau_likelihood_func_mixed(diff(c(0, update_tau[[j]])), data = data, beta = update_beta[[1]], discCols = discColumns, variable = j)
      }
      else{
        store_likelihood[i+1,j] = tau_likelihood_func_mixed(diff(c(0, init_tau[[j]])), data = data, beta = update_beta[[1]], discCols = discColumns, variable = j)
      }
      likelihood[i] = likelihood[i] + store_likelihood[i,j]
    }
    tau_estimates[[i+1]] = update_tau
    init_tau = update_tau
    init_beta = update_beta[[1]]
  }
  return(list(update_tau, update_beta[[1]], beta_estimates, tau_estimates,store_likelihood, likelihood, beta_diffs))
}

get_bounds <- function(observed, thresholds) {
  n <- length(thresholds)
  if (observed == 0) {
    # Case for observed value 0: lower bound is -Inf, upper bound is first threshold
    lower_bound <- -Inf
    upper_bound <- thresholds[1]
  } else if (observed > n-1) {
    # Case for observed value exceeding number of thresholds: last threshold to Inf
    lower_bound <- thresholds[n]
    upper_bound <- Inf
  } else {
    # General case: bounds are between consecutive thresholds
    lower_bound <- thresholds[observed]
    upper_bound <- thresholds[observed + 1]
  }
  return(c(lower_bound, upper_bound))
}

Sig_Estimate_DAG_mixed = function(X, trunc_vals, block_sizes, discCols){
  data = X
  t_data = t(data)
  cluster_number = length(block_sizes)
  estimated_sig = diag(nrow(data))
  
  avg_corr = matrix(data = 0, nrow = max(block_sizes), ncol = max(block_sizes)) # dont know if I need or not
  for(j in 1:cluster_number){
    #print(j)
    if(block_sizes[j] == 1){
      next
    }
    num_pairs = t(combn(1:block_sizes[j], 2))
    indices = (sum(block_sizes[-(j:cluster_number)]))
    for(i in 1:nrow(num_pairs)){
      col1 = indices+num_pairs[i,1]
      col2 = indices+num_pairs[i,2]
      cols = c(col1,col2)
      pair_data = t_data[,cols]
      c1_run = mapply(get_bounds, t_data[,col1], trunc_vals, SIMPLIFY = FALSE)
      c2_run = mapply(get_bounds, t_data[,col2], trunc_vals, SIMPLIFY = FALSE)

      log_lik_wrapper <- function(rho) {
        sum_log_lik_mixed_singrho(data = pair_data, rho = rho, c1 = c1_run, c2 = c2_run, discCols = discCols)
      }
      opt_result <- tryCatch({
        optimize(log_lik_wrapper, interval = c(0, 0.99), maximum = TRUE)
      }, error = function(e) {
        cat("Error during rho optimization:", conditionMessage(e), "\n")
        return(list(maximum = NA, objective = NA))
      })
      
      max_rho <- opt_result$maximum
      
      estimated_sig[indices+num_pairs[i,1],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,2]] = max_rho
      estimated_sig[indices+num_pairs[i,2],(sum(block_sizes[-(j:cluster_number)]))+num_pairs[i,1]] = max_rho
      
    }
    
  }
  return(estimated_sig)
}

beta_est_loop = function(data, init_beta, init_epsilon ,block_sizes, loops, trueB){
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
    trunc_vals = obtain_trunc_vals(data, update_beta)
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
          estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] = (estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] + 0.1*diag(block_sizes[j]))/(1 + 0.1)
        }
      }
      L_hat = chol(solve(estimated_Sigma))
    }
    else if(i == 1){
      estimated_Sigma = diag(nrow = nrow(data))
      L_hat = chol(estimated_Sigma)
    }
    eps_draw = withTimeout({
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
    start = Sys.time()
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = update_beta)
    end = Sys.time()
    # print(paste0("new beta step: ", end - start))
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


beta_est_loop_mixed = function(data, init_beta, init_epsilon ,block_sizes, loops, trueB, discColumns){
  # Run EM_tau_beta
  
  # est_parents = init_data_dag_mixed(data = data)
  output_thresholds <- EM_tau_beta_mixed(
    data = data, 
    parent_set = init_beta,
    discColumns = discColumns, # fix discCols
    reps = 5
  )
  hold_data = data
  hold_likelihoods = output_thresholds[[6]]
  est_sigma = Sig_Estimate_DAG_mixed(X = data, trunc_vals = output_thresholds[[1]], block_sizes = block_sizes, discCols = discColumns)
  #est_sigma = diag(nrow(data))
  # Start from here to use sigma estimate and obtained beta estimates to recover cont
  if(is.positive.definite(est_sigma) == F){
    estimated_Sigma = nearPD(est_sigma, corr = T, doSym = T)
  }
  else{
    estimated_Sigma = est_sigma
  }
  # Checking eigen values for non-invertible blocks and adding regularizer of 0.05 if close to non-inv.
  block_sizes = block_sizes
  trunc_vals = output_thresholds[[1]]
  update_beta = output_thresholds[[2]]
  
  cluster_number = length(block_sizes)
  for(k in 1:(cluster_number)){
    check_cov = estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])]
    if(any(eigen(check_cov)$values < 0.001)){
      estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] = (estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] + 0.1*diag(block_sizes[k]))/(1 + 0.1)
    }
  }
  L_hat = chol(solve(estimated_Sigma))
  #####
  # data and do de-correlation method at least on n = 100, p = 100, 300
  eps_draw = matrix(0, nrow = nrow(data), ncol = ncol(data))
  btw_betas = numeric(20)
  rmse = numeric(20)
  rmse_betas = numeric(20)
  hold_uncor_data = list()
  # From here Feb 3rd 2pm
  
  for(j in 1:20){
    
    #print(j)
    eps_draw = withTimeout({
      epsilon_draw_mixed(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = j, discCols = discColumns)
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    Z = obtain_hidden_Z_mixed(X = data, beta = update_beta, epsilon = eps_draw, discCols = discColumns)
    
    
    uncor_data = L_hat %*% Z
    if(j != 1){
      hold_uncor_data[[j]] = uncor_data
    }
    
    LX = L_hat %*% data
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = update_beta)
    btw_betas[j] = sqrt(mean(c(update_beta - new_beta_vals[[1]])^2))
    update_beta = new_beta_vals[[1]]
    rmse[j] = new_beta_vals[[2]]
    #diff = c(update_beta - data_gen$true_beta[[1]])
    #rmse_betas[j] = sqrt(mean((diff[which(diff != 0)])^2))
  }
  return(list(hold_data, hold_uncor_data, hold_likelihoods, estimated_Sigma))
}

beta_est_loop_input_covariance_mixed = function(data, adj_mat, cov_est, init_epsilon ,block_sizes, discColumns){
  eps_draw = init_epsilon
  update_beta = adj_mat
  n = nrow(data)
  p = ncol(data)
  estimated_Sigma = cov_est
  L_hat = chol(solve(estimated_Sigma))
  output_thresholds <- EM_tau_beta_mixed(
    data = data, 
    parent_set = adj_mat,
    discColumns = discColumns,
    reps = 5
  )
  for(i in 1:2){
    trunc_vals = output_thresholds[[1]]
    eps_draw = withTimeout({
      epsilon_draw_mixed(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = 2, discCols = discColumns)
      
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    Z = obtain_hidden_Z_mixed(X = data, beta = update_beta, epsilon = eps_draw, discCols = discColumns)
    
    
    uncor_data = L_hat %*% Z
    LX = L_hat %*% data
    
    #true_uncor_data = L_true %*% Z
    #true_LX = L_true %*% data
    start = Sys.time()
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = adj_mat)
    end = Sys.time()
    #print(paste0("new beta step: ", end - start))
    #summary(c(update_beta))
    update_beta = new_beta_vals[[1]]
  }
  
  return(list("beta" = update_beta, 
              "Sigma" = estimated_Sigma,
              "thresholds" = output_thresholds[[1]]))
}


sum_log_lik_mixed_singrho = function(data, rho, c1, c2, discCols){
  llk = 0
  sigma_rho = matrix(c(1,rho, rho, 1), nrow = 2, ncol = 2, byrow = TRUE)
  for(j in 1:nrow(data)){
    if(discCols[j] == 1){
      lower_bounds = c(c1[[j]][1], c2[[j]][1])
      upper_bounds = c(c1[[j]][2], c2[[j]][2])
      mean = c(0,0)
      cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
      cdf <- pmvnorm(lower = lower_bounds, upper = upper_bounds, mean = mean, sigma = cov_matrix)
      llk = llk + log(cdf)
    }
    else{
      mean = c(0,0)
      cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
      pdf = log(dmvnorm(matrix(data[j,], nrow = 1, ncol = 2), mean = c(0,0), sigma = cov_matrix))
      llk = llk + pdf
    }
  }
  return(llk)
}

compute_shd_cpdag <- function(cpdag_est, cpdag_true) {
  
  # adjacency matrices
  A_est  <- as(cpdag_est, "matrix")
  A_true <- as(cpdag_true, "matrix")
  
  # skeletons
  Skel_est  <- (A_est != 0) | (t(A_est) != 0)
  Skel_true <- (A_true != 0) | (t(A_true) != 0)
  
  # 1. Skeleton differences (additions + deletions)
  shd_skel <- sum(Skel_est != Skel_true) / 2
  
  # 2. Orientation differences (only where skeleton exists in both)
  common_edges <- Skel_est & Skel_true
  
  orient_diff <- 0
  p <- nrow(A_est)
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j && common_edges[i, j]) {
        # check direction mismatch
        if (A_est[i, j] != A_true[i, j] ||
            A_est[j, i] != A_true[j, i]) {
          orient_diff <- orient_diff + 1
        }
      }
    }
  }
  
  shd_skel + orient_diff
}


obtain_hidden_Z_mixed <- function(X, beta, epsilon, discCols){
  #' Obtain hidden variable Z from Beta, X, and epsilon
  #' Epsilon should be generated from the expectation
  n = nrow(X)
  p = ncol(X)
  Z_mat <- matrix(0, n, p)
  

  if(discCols[1] == 1){
    Z_mat[,1] <- epsilon[,1]
  }
  else{
    Z_mat[,1] = X[,1]
  }
  
  if(discCols[2] == 1){
    Z_mat[,2] = X[,1]*beta[1,2] + epsilon[,2]
  }
  else{
    Z_mat[,2] = X[,2]
  }
  
  
  for(i in 3:p) {
    if(discCols[i] == 1){
      Z_mat[,i] <- rowSums(sweep(X[,1:i-1], MARGIN = 2, beta[1:i-1,i], "*")) + epsilon[,i]
    }
    else{
      Z_mat[,i] = X[,i]
    }
  }
  dimnames(Z_mat) <- list(NULL, as.character(1:p))
  return(Z_mat)
}


init_epsilon_creation = function(n, p, sig, omg.sq){
  eps_mat <- matrix(0, n, p)
  for(i in 1:p){
    eps_mat[,i] <- mvrnorm(1, mu = rep(0, n), Sigma = omg.sq[1]*sig)
  }
  return(eps_mat)
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
  return(list(sig, block_sizes))
}

init_data_dag_mixed = function(data){ # do it with mxm package and get hc
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::mmhc(x = data, max_k = 7, alpha = 0.05)$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  #bnlearn_graph <- empty.graph(colnames(skeleton$G))
  #amat(bnlearn_graph) <- skeleton$G
  #df = data.frame(data)
  #colnames(df) = colnames(skeleton$G)
  #pchc::mmhc(data, skel = skeleton)
  #parents[1:ncol(data),1:ncol(data)] = skeleton$G
  return(parents)
}

init_data_dag_mixed_MXM = function(data){
  # Suppose your data is in df
  df <- data  
  
  # Set threshold for "few unique levels"
  # For example, 10 or fewer unique values
  threshold <- 10  
  
  df <- as.data.frame(lapply(df, function(col) {
    nunique <- length(unique(col))
    
    if ((is.factor(col) || is.character(col)) && nunique <= threshold) {
      # convert to ordered factor
      factor(col, levels = sort(unique(col)), ordered = TRUE)
    } else {
      # keep numeric or high-cardinality factor as-is
      col
    }
  }))
  
  output = MXM::mmhc.skel(df, test = NULL)$G
  return(output)
}


gaussCItestLocal <- function (x, y, S, suffStat) 
{
  # Using Fisher's z-transformation of the partial correlation, test for zero partial correlation. 
  #
  # Args:
  #   see function gaussCItest() in Package 'pcalg'
  #   suffStat contains an addtional element 'ESS.Mat', that is, the matrix of effective sample size 
  #
  # Returns:
  #   the p-value of the current test
  
  subMat <- suffStat$ESS.Mat[c(x,y,S),c(x,y,S)]
  ne <- mean(subMat[upper.tri(subMat)], na.rm = T)
  z <- zStat(x, y, S, C = suffStat$C, n = ne)
  2 * pnorm(abs(z), lower.tail = FALSE)
}

init_data_dag_mixed_pc = function(data){ # do it with mxm package and get hc
  het.obj <- hetcor(data)
  ## heterogeneous correlations
  het.corr <- het.obj$correlations
  ## standard deviations
  het.sd <- het.obj$std.errors
  ## effective sample size
  het.ess <- ((1-het.corr^2)/het.sd)^2
  
  #### 3. Call the PC Algorithm for Causal Discovery on Mixed data ####
  
  ##
  pc.fit <- pc(suffStat = list(C = het.corr, n = nrow(data), ESS.Mat = het.ess), 
                  indepTest = gaussCItestLocal, alpha = 0.01, p = ncol(data))
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  for(i in 1:ncol(data)){
    parents[i, pc.fit@graph@edgeL[[i]]$edges] = 1
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


### change this
beta_est_loop_multi = function(data, init_beta, init_epsilon ,block_sizes, loops, trueB){
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
    trunc_vals = obtain_trunc_vals(data, update_beta)
    #print(summary(c(update_beta)))
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
          estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] = (estimated_Sigma[((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j]),((sum(block_sizes[-(j:cluster_number)]))+1):sum(block_sizes[1:j])] + 0.1*diag(block_sizes[j]))/(1 + 0.1)
        }
        # (sum(block_sizes[-(j:cluster_number)]))
        # sum(block_sizes[1:j])
      }
      L_hat = chol(solve(estimated_Sigma))
    }
    else if(i == 1){
      estimated_Sigma = diag(nrow = nrow(data))
      L_hat = chol(estimated_Sigma)
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
    #print(paste0("new beta step: ", end - start))
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
              "est_betas" = update_betas))
}

get_eps_bounds_mixed = function(data, trunc_vals, discCols){
  a = matrix(0, nrow = nrow(data), ncol = ncol(data))
  b = matrix(0, nrow = nrow(data), ncol = ncol(data))
  for(j in 1:length(trunc_vals)){
    if(discCols[j] == 1){
      thresholds = c(-Inf, trunc_vals[[j]], Inf)
      if(length(unique(data[,j])) != 3){
        data[,j] = match(data[,j], sort(unique(data[,j]))) - 1
      }
      for(i in 1:nrow(data)){
        a[i,j] = thresholds[data[i,j]+1]
        b[i,j] = thresholds[data[i,j]+2]
      }
    }
  }
  return(list("a" = a, "b" = b))
}

epsilon_draw_mixed = function(Sigma_hat, data, trunc_vals, block_sizes, prev_iter, iter_num, discCols){
  # have a n by p truncation values.  Now just check each element in data
  a = get_eps_bounds_mixed(data,trunc_vals,discCols)$a
  b = get_eps_bounds_mixed(data,trunc_vals, discCols)$b
  eps_draws = matrix(0, nrow = nrow(data), ncol = ncol(data))
  cluster_number = length(block_sizes)
  for(i in 1:cluster_number){
    #print(paste0("i is:", i))
    
    init_index = (sum(block_sizes[-(i:cluster_number)]))+1
    end_index = sum(block_sizes[1:i])
    data_cluster = data[init_index:end_index,]
    a_bounds = a[init_index:end_index,]
    b_bounds = b[init_index:end_index,]
    Sigma_cluster = Sigma_hat[init_index:end_index,init_index:end_index]
    for(j in 1:ncol(data)){
      if(discCols[j] == 1){
        #print(paste0("j is:", j))
        if(iter_num == 1){
          #print(dim(data_cluster))
          if(block_sizes[i] == 1){
            draws = tmvtnorm::rtmvnorm(n = 20, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 20)
            eps_draws[init_index:end_index,j] = mean(draws)
          }
          else{
            draws = tmvtnorm::rtmvnorm(n = 20, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 20)
            eps_draws[init_index:end_index,j] = colMeans(draws)
          }
          if(any(is.nan(eps_draws[init_index:end_index,j]))){
            eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
          }
        }
        else{
          if(block_sizes[i] == 1){
            draws = tmvtnorm::rtmvnorm(n = 20, mean = 0, sigma = Sigma_cluster, lower = a_bounds[j], upper = b_bounds[j],algorithm = "gibbs", burn.in.samples = 20)
            eps_draws[init_index:end_index,j] = mean(draws)
          }
          else{
            draws = tmvtnorm::rtmvnorm(n = 20, mean = rep(0,nrow(data_cluster)), sigma = Sigma_cluster, lower = a_bounds[,j], upper = b_bounds[,j],algorithm = "gibbs", burn.in.samples = 20)
            eps_draws[init_index:end_index,j] = colMeans(draws)
          }
          if(any(is.nan(eps_draws[init_index:end_index,j]))){
            eps_draws[init_index:end_index,j] = prev_iter[init_index:end_index,j]
          }
        }
      }
    }
  }
  return(eps_draws)
}

obtain_hidden_Z <- function(X, beta, epsilon, discCols){
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
####
######

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
        #print(paste0("Condition number: ", kappa(t(LX[,parents]) %*% LX[,parents])))
        
        model = glmnet(x = data.matrix(LX[,parents]), y = uncor_Z[,i], alpha = 0, lambda = n*p / (10*(n + p)) + n*p / (sqrt(n * p)*10))
        beta_estimate = coef(model)
        beta_update[parents,i] = beta_estimate[2:nrow(beta_estimate),1]
        
      }
      else{
        #model = lm(uncor_Z[,i]~LX[,parents])
        #beta_estimate = summary(model)$coefficients
        lambda = n*p / (10*(n + p)) + n*p / (sqrt(n * p)*10)
        y = uncor_Z[,i] - mean(uncor_Z[,i])
        x = LX[,parents] - mean(LX[,parents])
        #print(paste0("Condition number: ", kappa(t(x) %*% x)))
        numerator = sum(x*y)
        denom = sum(x^2) + lambda
        beta = numerator/denom
        beta_update[parents,i] = beta
      }
    }
  }
  diff = (beta_update - init_beta)
  rmse_error = sqrt(mean(diff[which(diff != 0)]^2))
  #print(max_parents)
  return(list(beta_update, rmse_error))
}

t_normal_likelihood_mixed = function(testdata, Sigma_estimate, est_beta, trunc_vals, discColumns){
  total_ll = 0
  Xbeta = testdata %*% est_beta
  trunc_vals = get_eps_bounds_mixed(testdata, trunc_vals, discColumns)
  a_bounds = trunc_vals$a
  b_bounds = trunc_vals$b
  for(i in 1:ncol(testdata)){
    if(discColumns[i] == 1){
      total_ll = total_ll + log(pmvnorm(lower = a_bounds[,i], upper = b_bounds[,i], mean = Xbeta[,i], sigma = Sigma_estimate))
    }
    else{
      total_ll = total_ll + log(dmvnorm(matrix(testdata[,i], nrow = 1, ncol = length(testdata[,i])), mean = Xbeta[,i], sigma = Sigma_estimate))
    }
  }
  return(total_ll)
}


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

neighborhood_selection_mmhc_get_dag_decor = function(data){ # to get parents
  parents = matrix(0,nrow = ncol(data), ncol = ncol(data))
  colnames(parents) = rownames(parents) = 1:ncol(data)
  skel = pchc::mmhc(x = data, max_k = 7, alpha = 0.05)$dag$nodes
  for(i in 1:ncol(data)){
    parents[skel[[i]]$parents,i] = 1
  }
  return(parents)
}


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

