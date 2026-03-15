source("libraries.R")
source("helperFuncMixed.R")

mixed_data_decor = function(data, block_sizes){
  # Identify discrete columns if fewer that 5 levels in a row
  discColumns = numeric(ncol(data))
  for(k in 1:ncol(data)){
    if(length(unique(data[,k])) < 5){
      discColumns[k] = 1
    }
  }
  
  # Fill gaps if there are missing levels in discrete variables
  # by reducing base level to 0 with increments of 1
  # i.e. 0,0,0,2,2,2 -> 0,0,0,1,1,1
  dense_data = make_discrete_dense(data)
  
  # Initial estimate of parents
  est_parents = init_data_dag_mixed(data = dense_data)
  
  # Pre-estimation of thresholds T and beta
  output_thresholds <- EM_tau_beta_mixed(
    data = dense_data, 
    parent_set = est_parents,
    discColumns = discColumns,
    reps = 5
  )
  
  # Pre-estimation of row-covariance Sigma
  est_sigma = Sig_Estimate_DAG_mixed(X = dense_data, trunc_vals = output_thresholds[[1]], block_sizes = block_sizes, discCols = discColumns)

  data = dense_data
  trunc_vals = output_thresholds[[1]]
  update_beta = output_thresholds[[2]]
  
  # Check if PD, if not move to nearest PD matrix
  if(is.positive.definite(est_sigma) == F){
    estimated_Sigma = nearPD(est_sigma, corr = T, doSym = T)
  }
  else{
    estimated_Sigma = est_sigma
  }
  
  # Checking eigen values for non-invertible blocks and adding regularizer if close to non-inv.
  cluster_number = length(block_sizes)
  for(k in 1:(cluster_number)){
    check_cov = estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])]
    if(any(eigen(check_cov)$values < 0.001)){
      estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] = (estimated_Sigma[((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k]),((sum(block_sizes[-(k:cluster_number)]))+1):sum(block_sizes[1:k])] + 0.1*diag(block_sizes[k]))/(1 + 0.1)
    }
  }
  # Obtain Cholesky de-correlation factor of precision matrix
  L_hat = chol(solve(estimated_Sigma))
  #####
  
  eps_draw = matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  # List containing de-correlated datasets for consensus approach
  list_uncor_data = list()
  
  # De-correlation procedure iterated 20 times
  for(j in 1:21){
    print(j)
    
    # Draws of epsilon subject to estimated thresholds
    eps_draw = withTimeout({
      epsilon_draw_mixed(Sigma_hat = estimated_Sigma, data = data, trunc_vals = trunc_vals, block_sizes = block_sizes, prev_iter = eps_draw, iter_num = j, discCols = discColumns)
    }, timeout=600, onTimeout="silent")
    if(is.null(eps_draw)){
      print("Epsilon draw long run-time")
      return(NULL)
    }
    
    # Recovery of latent Z for discret variables
    Z = obtain_hidden_Z_mixed(X = data, beta = update_beta, epsilon = eps_draw, discCols = discColumns)
    
    # De-correlate data with estimated covariance
    uncor_data = L_hat %*% Z
    if(j != 1){
      list_uncor_data[[j-1]] = uncor_data
    }
    
    LX = L_hat %*% data
    # Obtain update of beta
    new_beta_vals = new_beta(uncor_Z = uncor_data, LX = LX, init_beta = update_beta)
    update_beta = new_beta_vals[[1]]
  }
  
  # Dataset for average approach
  avg_uncor_data <- Reduce("+", list_uncor_data) / length(list_uncor_data)
  return(list(list_uncor_data, avg_uncor_data))
}

# load("test_data.RData")
# test_decor = mixed_data_decor(test_data, block_sizes = rep(5,20))

