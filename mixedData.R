source("libraries.R")
source("multicategory/codes/helperFuncMixed.R")


mixed_data_gen = function(n, B_mat, cov_struc, reps, vers, thresholds){
  B_mat = genB_from_realDAG(B_mat)
  repeat {
    tryCatch({
      # Generate data using the current vers value
      data_gen <- data_generation_realDAG(
        n = n, 
        cov_struc = cov_struc, 
        reps = 1, 
        vers = vers, 
        thresholds = thresholds,
        Bmatrix = B_mat
      )
      # Run EM_tau_beta
      est_parents = init_data_dag_mixed(data = data_gen$Data[[1]])
      output_thresholds <- EM_tau_beta_mixed(
        data = data_gen$Data[[1]], 
        parent_set = est_parents,
        discColumns = data_gen$discCols,
        reps = 5
      )
      
      # If EM_tau_beta runs without errors, break out of the loop
      if (!is.null(output_thresholds)) break
      
    }, error = function(e) {
      # Log the error and resample vers
      cat("Error encountered: ", conditionMessage(e), "\n")
      vers <<- sample(1000:10000, size = 1)
      cat("Resampled vers: ", vers, "\n")
    })
  }
  discColumns = data_gen$discCols
  hold_data = data_gen$Data[[1]]
  hold_likelihoods = output_thresholds[[6]]
  est_sigma = Sig_Estimate_DAG_mixed(X = data_gen$Data[[1]], trunc_vals = output_thresholds[[1]], block_sizes = data_gen$block_sizes[[1]], discCols = discColumns)
  # Start from here to use sigma estimate and obtained beta estimates to recover cont
  if(is.positive.definite(est_sigma) == F){
    estimated_Sigma = nearPD(est_sigma, corr = T, doSym = T)
  }
  else{
    estimated_Sigma = est_sigma
  }
  # Checking eigen values for non-invertible blocks and adding regularizer of 0.05 if close to non-inv.
  block_sizes = data_gen$block_sizes[[1]]
  data = data_gen$Data[[1]]
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
    
    print(j)
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
    diff = c(update_beta - data_gen$true_beta[[1]])
    rmse_betas[j] = sqrt(mean((diff[which(diff != 0)])^2))
  }
  return(list(hold_data, hold_uncor_data, hold_likelihoods, rmse, true_beta = data_gen$true_beta[[1]]))
}

MMHC_pre_post_f1 = function(estdag_truebeta, estdag_data, estdag_uncor_data){
  # Convert true DAG to CPDAG
  init_f1_vec = numeric(8)
  avg_f1_vec = numeric(8)
  consensus_f1_vec = numeric(8)
  for(j in 1:8){
    print(j)
    true_beta = estdag_truebeta[[j]]
    true_beta[which(true_beta != 0)] = 1
    colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
    g1 <- as(t(true_beta), "graphNEL")
    true_cpdag <- dag2cpdag(g1)
    ############################
    # Getting initial CPDAG estimate from data
    ###################################
    disc_data = estdag_data[[j]]
    disc_data = make_discrete_dense(disc_data)
    init_dag = init_data_dag_mixed(disc_data)
    init_g1 <- as(t(init_dag), "graphNEL")
    init_cpdag <- dag2cpdag(init_g1)
    init_f1 = get_F1_dag_parents(init_cpdag@edgeL, true_cpdag@edgeL)
    init_f1_vec[j] = init_f1$f1_score
    # # Getting post CPDAG from init which uses last uncor_z
    # #############################
    # # Getting post CPDAG from Consensus estimate
    consensus_data = estdag_uncor_data[[j]]
    consensus_mat = matrix(0, nrow = ncol(disc_data), ncol = ncol(disc_data))
    sampled_datasets = sample(2:20, size = 10)
    for(i in sampled_datasets){
      #print(i)
      parents = neighborhood_selection_mmhc_get_dag_decor(consensus_data[[i]])
      parents = dag.to.cpdag(parents)
      consensus_mat = consensus_mat + parents
    }
    consensus_mat[which(consensus_mat <= 4)] = 0
    consensus_mat[which(consensus_mat >= 5)] = 1
    consensus_g1 <- as(t(consensus_mat), "graphNEL")
    consensus_cpdag <- dag2cpdag(consensus_g1)
    consensus_f1 = get_F1_dag_parents(consensus_cpdag@edgeL, true_cpdag@edgeL)
    consensus_f1_vec[j] = consensus_f1$f1
    # #############################################
    
    avg_mat = matrix(data = 0, nrow = nrow(disc_data), ncol = ncol(disc_data))
    for(k in 2:20){
      avg_mat = avg_mat + estdag_uncor_data[[j]][[k]]
    }
    avg_mat = avg_mat/19
    avg_dag = neighborhood_selection_mmhc_get_dag_decor(data = avg_mat)
    avg_beta_g1 <- as(t(avg_dag), "graphNEL")
    avg_beta_cpdag <- dag2cpdag(avg_beta_g1)
    avg_beta_f1 = get_F1_dag_parents(avg_beta_cpdag@edgeL, true_cpdag@edgeL)
    #true_beta_f1_list[[j]] = true_beta_f1
    avg_f1_vec[j] = avg_beta_f1$f1_score
    
  }
  return(list("consensus_f1" = consensus_f1_vec, "avg_f1" = avg_f1_vec, "init_f1" = init_f1_vec))
}

# struc_diff_1_mixed = mixed_data_gen(n = 100, B_mat = B_mat, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0, thresholds = c(-1,1))
output = mixed_data_gen(n = 100, B_mat = B_mat, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0, thresholds = c(-1,0,1))



f1_estdag_data_mixed = list()
f1_estdag_uncor_data_mixed = list()
f1_estdag_truebeta_mixed = list()
for(i in 1:8){
  print(paste0("We are on iteration: ", i))
  output = mixed_data_gen(n = 100, B_mat = B_mat, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0, thresholds = c(-1,0,1))
  f1_estdag_data_mixed[[i]] = output[[1]]
  f1_estdag_uncor_data_mixed[[i]] = output[[2]]
  f1_estdag_truebeta_mixed[[i]] = output$true_beta
  f1_data_iter = list(output[[1]], output[[2]], output$true_beta)
  save(f1_data_iter, file = paste0("100_mixed_data_",i,".RData"))
}

f1_data = list(f1_estdag_data_mixed, f1_estdag_uncor_data_mixed, f1_estdag_truebeta_mixed)
save(f1_data, file = paste0("100_mixed_data.RData"))





