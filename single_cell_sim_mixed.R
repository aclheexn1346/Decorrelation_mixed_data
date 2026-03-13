source("helperFuncMixed.R")
# load("single_cell_multicategorical_5_2025/othergenes.RData")
# load("single_cell_multicategorical_5_2025/block_sizes_alldata.RData")
# load("single_cell_multicategorical_5_2025/tot_sample_alldata.RData")
# load("single_cell_multicategorical_5_2025/allData_df_mixed.RData")
source("libraries.R")
source("helperFunc.R")
source("single_cell_funcs.R")
library(dplyr)
library(mclust)
targetgene <- readRDS("single_cell_data/sig_genes_log_val.rds")
sig_genes_log_val <- readRDS("single_cell_data/sig_genes_log_val_full.rds")
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene %>% dim()
set.seed(12)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 2000) %>% as.matrix() %>% t()
# use 2000 genes to cluster the 1018 cells


# default clustering ------------------------------------------------------
sc_block_idx_full <- readRDS("single_cell_data/single_cell_block_idx_full.rds")
for(i in 1:length(sc_block_idx_full)){
  names(sc_block_idx_full[[i]]) <- colnames(targetgene)[sc_block_idx_full[[i]]]
}
# Xp <- readRDS(file = "data/single_cell_data/sig_genes_log_val.rds")
Xp <- t(targetgene)
Xp %>% dim()
# randomly sample 20 cells from each cell type and merge the indices

full_log_vals = Xp

fit_gaussian_vs_mixture <- function(df, max_components = 3) {
  new_df = df
  discretes = 0
  conts = 0
  for(col in 1:ncol(df)){
    x = df[,col]
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    n <- length(x)
    
    if (n < 10) return(NA)  # not enough data
    
    # --- Single Gaussian log-likelihood
    mu <- mean(x)
    sigma2 <- var(x)
    ll_gauss <- sum(dnorm(x, mean = mu, sd = sqrt(sigma2), log = TRUE))
    bic_gauss <- ll_gauss - 0.5 * 2 * log(n)   # 2 params: mean, variance
    
    # --- Gaussian Mixture log-likelihood (up to max_components)
    gmm <- tryCatch(Mclust(x, G = 1:max_components, verbose = FALSE), 
                    error = function(e) NULL)
    
    if (is.null(gmm)) return(NA)
    
    bic_best <- max(gmm$BIC, na.rm = TRUE)   # already penalized
    k_best <- gmm$G
    if (bic_best > bic_gauss){
      discretized = FALSE
      for (cluster in 1:k_best){
        # check each cluster to see if it is normal
        # If any not normal, then discretize
        # If all normal, keep continuous
        cluster_data = x[gmm$classification == cluster]
        if(length(unique(cluster_data)) < 3){
          next
        }
        if(shapiro.test(cluster_data)$p.value < 0.05){
          new_df[,col] = gmm$classification - 1
          discretes = discretes + 1
          discretized = TRUE
          break
        }
      }
      if(discretized == FALSE){
        new_df[,col] = x
        conts = conts + 1
      }
      
    }
    else{
      new_df[,col] = x
      conts = conts + 1
    }
  }
  print(paste0("Number of discrete variables: ", discretes))
  print(paste0("Number of continuous variables: ", conts))
  
  return(new_df)
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


mixed_full_log_vals = fit_gaussian_vs_mixture(full_log_vals)

pearson_distance = function(X) {
  corr_mat <- cor(t(X), method = "pearson", use = "pairwise.complete.obs")
  dist_mat <- 1 - corr_mat
  dist_mat[dist_mat > 1] <- 1
  diag(dist_mat) <- 0
  as.dist(dist_mat)
}

othergenes_clean <- othergenes[!grepl("^HFF", rownames(othergenes)), ]
d = pearson_distance(scale(othergenes_clean))

# d <- dist(scale(othergenes), method = "euclidean")
hc1 <- hclust(d, method = "complete")
plot(hc1, cex = 0.6, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4, xlab = "", main = "Cell-type Clustering Dendrogram")


sub_grp <- cutree(hc1, k=100)
sub_grp %>% table()

clusters <- sub_grp[!grepl("^HFF", names(sub_grp))]
uniq_clusters = unique(clusters)
#cluster_size = 15 to 30 sampling from full data
tot_sample = c()
block_sizes = c()
for(i in uniq_clusters){
  cluster_sub_grp = sub_grp[which(sub_grp == i)]
  cluster_size = length(cluster_sub_grp)
  # cluster_sample = sample(cluster_sub_grp, cluster_size)
  block_sizes = c(block_sizes, cluster_size)
  tot_sample = c(tot_sample, cluster_sub_grp)
}

obtain_data_from_sample = function(tot_sample, orig_data_values){
  extracted_cell_names = names(tot_sample)
  extracted_samples = orig_data_values[extracted_cell_names,]
  df = data.frame(extracted_samples)
  df$cell_name = extracted_cell_names
  df$cluster_num = tot_sample
  df = df[c(52, 53, 1:51)]
  return(df)
}
df_mixed = obtain_data_from_sample(tot_sample, mixed_full_log_vals)
save(df_mixed, file = "sample_mixed.RData")



# estimating covariance from background genes
background_gene_cov_estimation = function(bgr_gene_df, n = 400, tot_sample, block_sizes){
  remove_columns = c() # Need to have a decent number of factors
  for(j in 1:ncol(bgr_gene_df)){
    if(length(levels(as.factor(bgr_gene_df[names(tot_sample),j]))) < 10){
      remove_columns = c(remove_columns, j)
    }
  }
  sample_from_bgr = setdiff(1:ncol(bgr_gene_df), remove_columns)
  sampled = sample(sample_from_bgr, size = n)
  extracted_cell_names = names(tot_sample)
  bgr_extracted_samples = bgr_gene_df[extracted_cell_names,sampled]
  bgr_data = fit_gaussian_vs_mixture(bgr_extracted_samples)
  init_epsilon = matrix(0, nrow = nrow(bgr_data), ncol = ncol(bgr_data))
  remove_column_indices = c() # remove any columns with less than two factors since every data is the same
  for(i in 1:ncol(bgr_data)){
    if(length(levels(as.factor(bgr_data[,i]))) < 2){
      remove_column_indices = c(remove_column_indices, i)
    }
  }
  if(length(remove_column_indices) > 0){
    bgr_data = bgr_data[,-remove_column_indices]
  }
  bgr_data = as.matrix(bgr_data)
  rownames(bgr_data) = 1:dim(bgr_data)[1]
  colnames(bgr_data) = 1:dim(bgr_data)[2]
  bgr_data = make_discrete_dense(bgr_data)
  est_parents = init_data_dag_mixed(data = bgr_data)
  est_parents[which(est_parents != 0)] = 0.0001
  
  init_epsilon = matrix(0, nrow = nrow(bgr_data), ncol = ncol(bgr_data))
  
  discCols = numeric(ncol(bgr_data))
  for(k in 1:ncol(bgr_data)){
    if(length(unique(bgr_data[,k])) < 4){
      discCols[k] = 1
    }
  }
  # Getting estimated covariance and the uncorrelated data for the concensus cpdag estimate
  beta_est_ident_lr <- beta_est_loop_mixed(data = as.matrix(bgr_data), init_beta = est_parents, init_epsilon = init_epsilon, block_sizes = block_sizes, loops = 21, trueB = init_beta, discColumns = discCols)
  est_Sigma = beta_est_ident_lr[[4]]
  
  return(est_Sigma)
}

bgr_est_Sigma = background_gene_cov_estimation(bgr_gene_df = othergenes_clean, n = 100, tot_sample = tot_sample, block_sizes = block_sizes)

 # Data preparation
data = df_mixed[,c(-1,-2)]
gene_names = names(data)
data = as.matrix(data)
rownames(data) = 1:dim(data)[1]
colnames(data) = 1:dim(data)[2]
init_epsilon = matrix(0, nrow = nrow(data), ncol = ncol(data))


k <- 10
total <- sum(block_sizes)
target <- total / k

cum_sum <- cumsum(block_sizes)

# Greedily cut when cumulative size exceeds multiples of the target
cut_points <- sapply(1:(k-1), function(i) which.min(abs(cum_sum - i * target)))

group_id <- cut(seq_along(block_sizes),
                breaks = c(0, cut_points, length(block_sizes)),
                labels = FALSE)

# Check how balanced they are
group_sizes = c(tapply(block_sizes, group_id, sum))


# Cross validation likelihood by cluster

CV_likelihood_all = function(data, block_sizes, bgr_Estimated_Sigma, group_sizes, group_id){
  data = make_discrete_dense(data)
  
  likelihood_folds_init_2 = rep(0,10)
  likelihood_folds_concensus_2 = rep(0, 10)
  likelihood_folds_concensus_nocov_2 = rep(0, 10)
  #likelihood_folds_complete_nocov = rep(0, length(block_sizes))
  #likelihood_folds_complete_cov = rep(0, length(block_sizes))
  
  for(cluster_number in 10:length(group_sizes)){
    print(cluster_number)
    init_index = (sum(group_sizes[0:(cluster_number-1)]))+1
    end_index = sum(group_sizes[1:cluster_number])
    test_data = data[init_index:end_index,]
    train_data = data[-(init_index:end_index),]
    test_data = make_discrete_dense(test_data, threshold = 4)
    train_data = make_discrete_dense(train_data, threshold = 4)
    remove_column_indices = c() # remove any columns with less than two factors since every data is the same
    for(i in 1:ncol(train_data)){
      if(length(levels(as.factor(train_data[,i]))) < 2){
        remove_column_indices = c(remove_column_indices, i)
      }
    }
    if(length(remove_column_indices) > 0){
      train_data = train_data[,-remove_column_indices]
      test_data = test_data[,-remove_column_indices]
    }
    train_data = as.matrix(train_data)
    rownames(train_data) = 1:dim(train_data)[1]
    colnames(train_data) = 1:dim(train_data)[2]
    
    init_dag = init_data_dag_mixed(data = train_data)
    init_beta = init_dag
    init_beta[which(init_beta != 0)] = 0.0001
    
    init_epsilon = matrix(0, nrow = nrow(train_data), ncol = ncol(train_data))
    
    discCols = numeric(ncol(data))
    for(k in 1:ncol(data)){
      if(length(unique(data[,k])) < 5){
        discCols[k] = 1
      }
    }
    # Getting estimated covariance and the uncorrelated data for the concensus cpdag estimate
    beta_est_ident_lr <- beta_est_loop_mixed(data = as.matrix(train_data), init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_sizes[-which(group_id == cluster_number)], loops = 21, trueB = init_beta, discColumns = discCols[-remove_column_indices])
    
    # getting the final cpdag from the initial data
    #init_dag = init_data_dag(train_data)
    init_g1 <- as(t(init_dag), "graphNEL")
    init_cpdag <- dag2cpdag(init_g1)
    init_adj_mat = (as(init_cpdag, "matrix"))
    
    # Getting the final cpdag from the concensus data
    concensus_data = beta_est_ident_lr[[2]]
    concensus_mat = matrix(0, nrow = ncol(train_data), ncol = ncol(train_data))
    sampled_datasets = sample(2:20, size = 10)
    for(i in sampled_datasets){
      #print(i)
      parents = neighborhood_selection_mmhc_get_dag_decor(concensus_data[[i]])
      parents = dag2cpdag(parents)
      concensus_mat = concensus_mat + parents
    }
    concensus_mat[which(concensus_mat <= 4)] = 0
    concensus_mat[which(concensus_mat >= 5)] = 1
    concensus_g1 <- as(t(concensus_mat), "graphNEL")
    concensus_cpdag = dag2cpdag(concensus_g1)
    concensus_adj_mat = (as(concensus_cpdag, "matrix"))
    
    # Getting the final beta using the final cpdag structure obtained
    # First is using the final cpdag from the concensus and using the estimated covariance to make one last beta estimate
    concensus_adj_mat[concensus_adj_mat != 0] = 0.0001
    conc_final = beta_est_loop_input_covariance_mixed(data = train_data, adj_mat = concensus_adj_mat, cov_est = nearPD(beta_est_ident_lr[[4]]), init_epsilon = init_epsilon, block_sizes = block_sizes[-which(group_id == cluster_number)], discColumns = discCols[-remove_column_indices])
    conc_final_beta = conc_final[[1]]
    conc_thresholds = conc_final[[3]]
    
    # Using the initial cpdag and assumption of identity matrix but obtaining the final beta using same method as above
    init_final = beta_est_loop_input_covariance_mixed(data = train_data, adj_mat = init_beta, cov_est = diag(nrow(train_data)), init_epsilon = init_epsilon, block_sizes = block_sizes[-which(group_id == cluster_number)], discColumns = discCols[-remove_column_indices])
    init_final_beta = init_final[[1]]
    init_thresholds = init_final[[3]]

    # tau_likelihood_func_mixed with delta = diff(c(0, init_tau[[j]]))
    conc_ll2 = t_normal_likelihood_mixed(testdata = test_data, Sigma_estimate = bgr_Estimated_Sigma[init_index:end_index, init_index:end_index], est_beta = conc_final_beta, trunc_vals = conc_thresholds, discColumns = discCols[-remove_column_indices])
    conc_nocov_ll2 = t_normal_likelihood_mixed(testdata = test_data, Sigma_estimate = diag(nrow(test_data)), est_beta = conc_final_beta, trunc_vals = conc_thresholds, discColumns = discCols[-remove_column_indices])
    init_ll2 = t_normal_likelihood_mixed(testdata = test_data, Sigma_estimate = diag(nrow(test_data)), est_beta = init_final_beta, trunc_vals = init_thresholds, discColumns = discCols[-remove_column_indices])

    
    likelihood_folds_concensus_2[cluster_number] = conc_ll2[1]/(nrow(test_data)*ncol(test_data))
    likelihood_folds_concensus_nocov_2[cluster_number] = conc_nocov_ll2[1]/(nrow(test_data)*ncol(test_data))
    likelihood_folds_init_2[cluster_number] = init_ll2[1]/(nrow(test_data)*ncol(test_data))
    ll_list = list(init_ll2[1]/(nrow(test_data)*ncol(test_data)), conc_nocov_ll2[1]/(nrow(test_data)*ncol(test_data)), conc_ll2[1]/(nrow(test_data)*ncol(test_data)))
    save(ll_list, file = paste0("loglikelihood_normalized_fold_alldata_", cluster_number,".RData"))
  }
  return(list(
    "init_mvt" = likelihood_folds_init_2,
    "consensus_mvt" = likelihood_folds_concensus_2,
    "consensus_nocov_mvt" = likelihood_folds_concensus_nocov_2,
    "Final_consensus_CPDAG" = concensus_adj_mat))
}

# stuck on cluster 4 idk why. look into that one
CV_res_mixed = CV_likelihood_all(data = data, block_sizes = block_sizes, bgr_Estimated_Sigma = bgr_est_Sigma,group_sizes = group_sizes, group_id = group_id)


single_cell_mixed_train_all = function(data, block_sizes, bgr_Estimated_Sigma){
  data = make_discrete_dense(data)
  train_data = data
  train_data = make_discrete_dense(train_data, threshold = 4)
  remove_column_indices = c() # remove any columns with less than two factors since every data is the same
  for(i in 1:ncol(train_data)){
    if(length(levels(as.factor(train_data[,i]))) < 2){
      remove_column_indices = c(remove_column_indices, i)
    }
  }
  if(length(remove_column_indices) > 0){
    train_data = train_data[,-remove_column_indices]
  }
  train_data = as.matrix(train_data)
  rownames(train_data) = 1:dim(train_data)[1]
  colnames(train_data) = 1:dim(train_data)[2]
  init_dag = init_data_dag_mixed(data = train_data)
  init_beta = init_dag
  init_beta[which(init_beta != 0)] = 0.0001
  
  init_epsilon = matrix(0, nrow = nrow(train_data), ncol = ncol(train_data))
  
  discCols = numeric(ncol(data))
  for(k in 1:ncol(data)){
    if(length(unique(data[,k])) < 5){
      discCols[k] = 1
    }
  }
    # Getting estimated covariance and the uncorrelated data for the concensus cpdag estimate
  beta_est_ident_lr <- beta_est_loop_mixed(data = as.matrix(train_data), init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_sizes, loops = 21, trueB = init_beta, discColumns = discCols)
    
    # getting the final cpdag from the initial data
  #init_dag = init_data_dag(train_data)
  init_g1 <- as(t(init_dag), "graphNEL")
  init_cpdag <- dag2cpdag(init_g1)
  init_adj_mat = (as(init_cpdag, "matrix"))

  # Getting the final cpdag from the concensus data
  concensus_data = beta_est_ident_lr[[2]]
  concensus_mat = matrix(0, nrow = ncol(train_data), ncol = ncol(train_data))
  sampled_datasets = sample(2:20, size = 10)
  for(i in sampled_datasets){
    #print(i)
    parents = neighborhood_selection_mmhc_get_dag_decor(concensus_data[[i]])
    parents = dag2cpdag(parents)
    concensus_mat = concensus_mat + parents
  }
  concensus_mat[which(concensus_mat <= 4)] = 0
  concensus_mat[which(concensus_mat >= 5)] = 1
  concensus_g1 <- as(t(concensus_mat), "graphNEL")
  concensus_cpdag = dag2cpdag(concensus_g1)
  concensus_adj_mat = (as(concensus_cpdag, "matrix"))
  
  # Getting the final beta using the final cpdag structure obtained
  # First is using the final cpdag from the concensus and using the estimated covariance to make one last beta estimate
  concensus_adj_mat[concensus_adj_mat != 0] = 0.0001
  conc_final = beta_est_loop_input_covariance_mixed(data = train_data, adj_mat = concensus_adj_mat, cov_est = nearPD(beta_est_ident_lr[[4]]), init_epsilon = init_epsilon, block_sizes = block_sizes, discColumns = discCols)
  conc_final_beta = conc_final[[1]]
  conc_thresholds = conc_final[[3]]
  
    # Using the initial cpdag and assumption of identity matrix but obtaining the final beta using same method as above
  init_final = beta_est_loop_input_covariance_mixed(data = train_data, adj_mat = init_beta, cov_est = diag(nrow(train_data)), init_epsilon = init_epsilon, block_sizes = block_sizes, discColumns = discCols)
  init_final_beta = init_final[[1]]
  init_thresholds = init_final[[3]]
  
  return(list("Final_consensus_CPDAG" = concensus_adj_mat))
}


final_consensus_mat_mixed = single_cell_mixed_train_all(data = data, block_sizes = block_sizes, bgr_Estimated_Sigma = bgr_est_Sigma)
final_consensus_mat_mixed$Final_consensus_CPDAG[final_consensus_mat_mixed$Final_consensus_CPDAG != 0] = 1
rownames(final_consensus_mat_mixed$Final_consensus_CPDAG) = colnames(final_consensus_mat_mixed$Final_consensus_CPDAG) = gene_names
final_cpdag = final_consensus_mat_mixed$Final_consensus_CPDAG
save(final_cpdag, file = "alldata_cpdag.RData")


final_cpdag_edgecount = final_cpdag
final_cpdag_edgecount[final_cpdag_edgecount == 1] = 0
final_cpdag_edgecount = final_cpdag_edgecount + edge_count
save(final_cpdag_edgecount, file = "single_cell_multicategorical_5_2025/single_cell_mixed_conf_edges/alledges.RData")
# adjMat should be a square matrix with rownames and colnames

# Get indices of edges present in the final_cpdag
edge_indices <- which(final_cpdag == 1, arr.ind = TRUE)

# Create a data frame with edges and their counts
edge_counts_df <- data.frame(
  from = rownames(final_cpdag)[edge_indices[, 1]],
  to   = colnames(final_cpdag)[edge_indices[, 2]],
  count = final_cpdag_edgecount[edge_indices]
)

# Preview the result
print(edge_counts_df)

# Make a new column that defines an edge key in both directions
edge_counts_df$pair <- paste(edge_counts_df$from, edge_counts_df$to, sep = "_")
reverse_pairs <- paste(edge_counts_df$to, edge_counts_df$from, sep = "_")

# Check if the reverse pair exists
edge_counts_df$type <- ifelse(reverse_pairs %in% edge_counts_df$pair, "undirected", "directed")

# Remove duplicate undirected edges (keep only one direction)
# edge_counts_df_unique <- edge_counts_df[!(edge_counts_df$type == "undirected" & edge_counts_df$from > edge_counts_df$to), ]

# Split into directed and undirected
directed <- subset(edge_counts_df, type == "directed")
undirected <- subset(edge_counts_df, type == "undirected")

# View result
print("Directed edges:")
print(directed)

print("Undirected edges:")
print(undirected)

write.csv(undirected, file = "single_cell_multicategorical_5_2025/single_cell_mixed_conf_edges/undirected_final_cpdag_conf.csv")

get_edge_lists <- function(adjMat) {
  adjMat <- as.matrix(adjMat)
  
  # All edges
  edges <- which(adjMat == 1, arr.ind = TRUE)
  edge_list <- data.frame(
    Parent = rownames(adjMat)[edges[, 1]],
    Child  = colnames(adjMat)[edges[, 2]],
    stringsAsFactors = FALSE
  )
  
  # Check which edges are reciprocated
  is_undirected <- mapply(function(p, c) adjMat[c, p] == 1, 
                          edge_list$Parent, edge_list$Child)
  
  # Undirected edges (keep each pair only once)
  undirected <- edge_list[is_undirected, ]
  if (nrow(undirected) > 0) {
    # Order pairs so A-B and B-A collapse to one
    undirected <- undirected[!duplicated(t(apply(undirected, 1, sort))), ]
  }
  
  # Directed edges (not reciprocated)
  directed <- edge_list[!is_undirected, ]
  
  return(list(
    Directed = directed,
    Undirected = undirected
  ))
}

edgeList_df = get_edge_lists(final_cpdag)
save(edgeList_df, file = "single_cell_multicategorical_5_2025/alldata_cpdag_edgelist.RData")

directed = edgeList_df$Directed
undirected = edgeList_df$Undirected

write.csv(weighted_edge_list, file = "single_cell_multicategorical_5_2025/confidence_edges.csv")

#######  Importing data for bootstrap ########

source("libraries.R")
source("helperFunc.R")
source("single_cell_funcs.R")
library(dplyr)
targetgene <- readRDS("single_cell_data/sig_genes_log_val.rds")
sig_genes_log_val <- readRDS("single_cell_data/sig_genes_log_val_full.rds")
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene %>% dim()
set.seed(12)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 2000) %>% as.matrix() %>% t()
# use 2000 genes to cluster the 1018 cells


# default clustering ------------------------------------------------------
sc_block_idx_full <- readRDS("single_cell_data/single_cell_block_idx_full.rds")
for(i in 1:length(sc_block_idx_full)){
  names(sc_block_idx_full[[i]]) <- colnames(targetgene)[sc_block_idx_full[[i]]]
}
# Xp <- readRDS(file = "data/single_cell_data/sig_genes_log_val.rds")
Xp <- t(targetgene)
Xp %>% dim()
# randomly sample 20 cells from each cell type and merge the indices

full_log_vals = Xp

# targetgene will be the 51 genes and othergenes are the background genes
# # clustering --------------------------------------------------------------

B = 10
frac <- 0.5
n = dim(full_log_vals)[1]
for(i in 1:B){
  idx <- sample(1:n, size = floor(frac*n), replace = FALSE)
  othergenes_samp = othergenes[idx,]
  targetgene_samp = targetgene[idx,]
  # clustering
  d <- dist(scale(othergenes), method = "euclidean")
  hc1 <- hclust(d, method = "complete")
  # getting sample
  sub_grp <- cutree(hc1, k=100)
  clusters <- sub_grp[!grepl("^HFF", names(sub_grp))]
  uniq_clusters = unique(clusters)
  
  # get total sample and put into tot_sample and get blocksizes
  tot_sample = c()
  block_sizes = c()
  for(i in uniq_clusters){
    cluster_sub_grp = sub_grp[which(sub_grp == i)]
    cluster_size = length(cluster_sub_grp)
    # cluster_sample = sample(cluster_sub_grp, cluster_size)
    block_sizes = c(block_sizes, cluster_size)
    tot_sample = c(tot_sample, cluster_sub_grp)
  }
  
  
  
}



####### Plotting graph of estimated GRN #########
g = graph_from_data_frame(all_edges)
deg <- degree(g, mode = "all")  # or mode="in"/"out" if you prefer

ord <- order(deg, decreasing = TRUE)
ord[1] = 11
ord[26] = 17
# 11 <-> 17
# 34 <-> 8
ord[17] = 8
ord[23] = 34
# 6 <-> 14
ord[21] = 14
ord[29] = 6
# 28 <-> 37
g = permute(g, ord)

library(igraph)
E(g)$type <- all_edges$type
lay <- layout_in_circle(g) 
E(g)$color <- "black" 
E(g)$width <- 1 
E(g)$arrow.size <- ifelse(E(g)$type == "undirected", 0, 0.2)
V(g)$size <- 3
V(g)$color <- "black" 
V(g)$frame.color <- NA 
n <- vcount(g) 
angle <- 360 * (seq_len(n) - 1) / n 
V(g)$label.degree <- angle * pi / 180 
V(g)$label.dist <- 1.75 
V(g)$label.cex <- .75
arrow_sizes <- ifelse(E(g)$type == "undirected", 0, 0.2)
plot(g, 
     layout = lay, 
     vertex.label = V(g)$name, 
     vertex.label.color = "black",
     edge.color = ifelse(E(g)$type == "directed", "black", NA),
     edge.width = E(g)$width,
     edge.arrow.size = 1,
     asp = 1, margin = -.1, 
     vertex.label.degree=lab.locs)

undirected_edges <- E(g)[E(g)$type == "undirected"]
for(e in undirected_edges) {
  from <- ends(g, e)[1]
  to   <- ends(g, e)[2]
  segments(
    lay[which(V(g)$name == from), 1], lay[which(V(g)$name == from), 2],
    lay[which(V(g)$name == to), 1],   lay[which(V(g)$name == to), 2],
    col = "DodgerBlue", lwd = 1
  )
}
title("Estimated Single-cell GRN")
