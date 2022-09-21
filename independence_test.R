
#' Independence test between L and R amino-acids given C and (phi,psi)
#'
#' For a given central amino-acid C and a sample of (phi, psi) points on the Ramachandran map, the function tests the
#' indepedence of left and right neighbors identities (L and R) for the fixed C and a set of fixed values of (phi, psi)
#' computed after discretizing the Ramachandran space.
#'
#' @param central_name The name of the central residue C.
#' @param N_bin The number of bins of the discretized space (ignored if N_rep is not NULL).
#' @param N_n For discretization method I, the minimum number of points that a space division must contain to be kept. For methods II and III, 
#' the exact number of points that the space division must contain.
#' @param ball If FALSE, the grid division is a square (discretization method I). If TRUE, the grid division is a torus-ball (discretization methods II or III).
#' @param N_rep If NULL, the grid is homogeneous and deterministically chosen (discretization methods I and II). Else, 
#' the grid is random and its representative points are chosen by sub-sampling N_rep points from the sample of all points having C as
#' central residue. 
#' @param dist_tol If ball is TRUE (for discretization methods II and II), the maximum radius that a subdivision (ball) with N_n points
#' must have to be accepted.
#' @param aa_list A list of the allowed values of L and R.
#' @param aa_list_left If not NULL, a list of the allowed values of L (it ignores aa_list for L when not NULL).
#' @param aa_list_right If not NULL, a list of the allowed values of R (it ignores aa_list for R when not NULL).
#' @param set_path The path to the directory where the data (of the central residue C) is contained. 
#' @param plot Logical: whether to plot the sample and the performed discretization on the torus.
#' @param sim_pvalue The parameter simulate.p.value of R's chisq.test function.
#' @param only_trans Logical: whether to filter CIS conformations and keep only TRANS ones. 
#'
#' @return A dataframe with one row per (phi, psi) value of the discretization (i.e. one row per performed test) containing as columns:
#' \item{phi_grid}{The value of phi for the given division representative point.}
#' \item{psi_grid}{The value of psi for the given division representative point.}
#' \item{p_value}{The p-value for the division-specific independence test (L and R independent given C and phi_grid, psi_grid).}
#' \item{D}{The test statistic realization for the division-specific independence test (L and R independent given C and phi_grid, psi_grid).}
#' \item{mat_dim}{The dimension of the contingence matrix for the division-specific independence test (L and R independent given C and phi_grid, psi_grid).}
#' \item{mat_size}{The number of points included in the contigence matrix for the division-specific independence test (L and R independent given C and phi_grid, psi_grid).}
#' \item{p_value_bonf}{The Bonferroni-corrected p_value.}
#' \item{p_value_sidak}{The Sidak-corrected p_value.}
#' 
#' @export


independence_test <- function(central_name, N_bin, N_n = NULL, ball = FALSE, N_rep = NULL, dist_tol = 1e-2,
                              aa_list, aa_list_left = NULL, aa_list_right = NULL,
                              set_path = data_path_all, plot = TRUE, sim_pvalue = FALSE, only_trans = TRUE){
  
  
  
  if(is.null(aa_list_left)){aa_list_left <- aa_list}
  if(is.null(aa_list_right)){aa_list_right <- aa_list}
  
  # Import central amino-acid data
  data_central <- get(load(file.path(set_path, paste(c(central_name,'angles.RData'),collapse='_'))))
  
  # Remove tripeptides with PRO and GLY as left or right neighbors
  data_central <- data_central[-which(data_central$Res1 == 'PRO' | data_central$Res3 == 'PRO' | data_central$Res1 == 'GLY' | data_central$Res3 == 'GLY'), ]
  
  # Remove cis conformations
  if(only_trans == TRUE){data_central <- data_central[which(data_central$CIS.0..TRANS.1. == 1), ]}
  
  # Needed variables for coding
  data_central$torus_dist <- NA
  data_central$out <- as.factor(rep(0, nrow(data_central)))
  levels(data_central$out) <- 0:1
  
  # Choose homogeneous or random grid
  if(is.null(N_rep) == TRUE){  # Homogeneous grid
    step <- 2*pi/N_bin 
    phi_grid <- seq(-pi, pi, length = N_bin + 1)
    psi_grid <- seq(-pi, pi, length = N_bin + 1)
    pgrid <- as.data.frame(expand.grid(phi_grid[1:N_bin], psi_grid[1:N_bin]))}else{ # Random grid
      ind <- sample(1:nrow(data_central),size = N_rep)
      pgrid <- data_central[ind, c('Phi_res_2','Psi_res_2')]
    }
  colnames(pgrid) <- c('phi_grid','psi_grid')
  
  pgrid$p_value <- NA;  pgrid$D <- NA
  pgrid$mat_dim <- NA;  pgrid$mat_size <- NA
  n_empty <- 0
  
  # For each point on the grid:
  #1. Define the associated division (set of points with equivalent angle values in terms of discretization),
  #2. Perform the chi-square independence test for the set of points in the division.
  
  for(i in 1:nrow(pgrid)){
    
    # Choose how to build the divisions (squares defined by the grid or balls defined by the torus distance)
    
    if(ball == FALSE){ # Square grid
      data_bin <- data_central[which(data_central$Phi_res_2 >= pgrid[i,1] & data_central$Phi_res_2 < pgrid[i,1] + step &
                                       data_central$Psi_res_2 >= pgrid[i,2] & data_central$Psi_res_2 < pgrid[i,2] + step), ]
      
      # Dataset for test performing (for contingency matrix building)
      
      if(nrow(data_bin) < N_n){ # Remove bins with less than N_n points. This option can be unset by letting N_n = Inf.
        data_central$out[which(data_central$Phi_res_2 >= pgrid[i,1]
                               & data_central$Phi_res_2 < pgrid[i,1] + step & data_central$Psi_res_2 >= pgrid[i,2]
                               & data_central$Psi_res_2 < pgrid[i,2]+step)] <- 1
        next} 
      
    }else{ 
      
      # Bins are balls defined by torus distance with center at each grid point
      # Every ball is defined as including exactly N_n points and a maximum radius
      
      point <- c(pgrid[i,1], pgrid[i,2]) # The grid point (the center of the ball)
      
      # Computing torus distance between all the points and the ball center
      dist_vector <- sqrt(pmin(abs(point[1] - data_central[ ,c('Phi_res_2')]), 2*pi - abs(point[1] - data_central[ ,c('Phi_res_2')]))^2 + 
                            pmin(abs(point[2] - data_central[ ,c('Psi_res_2')]), 2*pi - abs(point[2] - data_central[ ,c('Psi_res_2')]))^2)
      
      # Take the N_n points closer to the grid point and
      # remove balls with radius > tolerance
      
      if(max(sort(dist_vector)[1:N_n]) > dist_tol){next}else{
        
        data_bin <- data_central[which(dist_vector %in% sort(dist_vector)[1:N_n]), ] # Dataset for test performing (for contingency matrix building)
        data_central$torus_dist[which(dist_vector %in% sort(dist_vector)[1:N_n])] <- i # Categorical variable for color plotting
      }}
    
    # Contingency matrix building
    # Counting the number of tripeptides in each division
    
    left_right_list <- as.data.frame(expand.grid(aa_list_left, aa_list_right)) # L-R list for the fixed central amino-acid
    colnames(left_right_list) <- c('left','right')
    left_right_list$left <- as.character(left_right_list$left); left_right_list$right <- as.character(left_right_list$right)
    left_right_list$N_trip <- NA # Number of points with the corresponding L and R
    
    # Filling N_trip variable in left_right_list
    for(k in 1:nrow(left_right_list)){
      left_right_list$N_trip[k] <- nrow(data_bin[which(data_bin$Res1 == left_right_list$left[k] & data_bin$Res3 == left_right_list$right[k]), ])
    }
    
    # Transforming left_right_list into a (contingency) matrix
    cont_mat <- matrix(data = left_right_list$N_trip,nrow = length(aa_list_left), ncol = length(aa_list_right), byrow = FALSE)
    colnames(cont_mat) <- aa_list_right; rownames(cont_mat) <- aa_list_left
    
    #Removing empty rows/columns
    if(sum(cont_mat) == 0){
      
      cat('\n No points in this grid division\n')
      n_empty <- n_empty + 1
      next}else{
        
        while(nrow(cont_mat)*ncol(cont_mat) != 0){
          
          nxt<-FALSE
          sum_col <- apply(cont_mat, FUN = sum, MARGIN = 2)
          
          if(any(sum_col == 0) == TRUE){
            
            if(length(which(sum_col == 0)) >= ncol(cont_mat) - 1){nxt <- TRUE; break}
            
            cont_mat <- cont_mat[-which(sum_col == 0), ]
            cont_mat <- cont_mat[ ,-which(sum_col == 0)]}
          
          sum_row <- apply(cont_mat, FUN = sum, MARGIN = 1)
          if(any(sum_row == 0) == TRUE){
            
            if(length(which(sum_row == 0)) >= nrow(cont_mat) - 1){nxt <- TRUE; break}
            cont_mat <- cont_mat[-which(sum_row == 0), ]
            cont_mat <- cont_mat[ ,-which(sum_row == 0)]
            
          }
          if(any(sum_row == 0) == FALSE & any(sum_col == 0) == FALSE){break}
        }
        
        if(nrow(cont_mat)*ncol(cont_mat) == 0 | nxt == TRUE){next}else{
          
          test <- chisq.test(cont_mat, simulate.p.value = sim_pvalue) # Performs chi-square independence test
          # Test is performed using the asymptotic distribution of the Pearson's chi-square statistic
          
          pgrid$p_value[i] <- test$p.value  #test p-value
          pgrid$D[i] <- test$statistic #test statistic value
          pgrid$mat_dim[i] <- min(dim(cont_mat)) #Contingency matrix dimension
          pgrid$mat_size[i] <- sum(cont_mat) #Contingency matrix number of points
        }}
  }
  
  # Multiple-test corrections
  N_test <- nrow(pgrid)*length(aa_list) # All the test (for each central amino-acid) must be performed using the same discretization parameters
  pgrid$p_value_bonf <- p.adjust(pgrid$p_value, method = 'bonferroni') # Bonferroni adjusted p-values
  pgrid$p_value_holm <- p.adjust(pgrid$p_value, method = 'holm') # Holm adjusted p-values
  
  # Plotting S^1\times S^1 space with square grid discretization
  plot_1 <- ggplot2::ggplot(data = data_central, aes(x = Phi_res_2, y = Psi_res_2, col = out, group = out))+
      ggplot2::geom_point(size = 0.005)+
      ggplot2::xlab('Phi')+
      ggplot2::ylab('Psi')+
      ggplot2::ggtitle(central_name)+
      ggplot2::theme(legend.position = 'none')
  
  # Plotting S^1\times S^1 space with torus ball discretization
  plot_2 <- ggplot2::ggplot(data = data_central, aes(x = Phi_res_2, y = Psi_res_2, col = as.factor(torus_dist)))+
    ggplot2::geom_point(size = 0.005)+
    ggplot2::xlab('Phi')+
    ggplot2::ylab('Psi')+
    ggplot2::xlim(-pi,pi)+
    ggplot2::ylim(-pi,pi)+
    #ggplot2::ggtitle(central_name)+
    ggplot2::theme(legend.position = 'none')
  
  if(plot & !ball){
    print(plot_1+
            ggplot2::geom_vline(size = 0.2, colour = 'darkblue', linetype = 'dashed', xintercept = phi_grid)+
            ggplot2::geom_hline(size = 0.2, colour = 'darkblue', linetype = 'dashed', yintercept = psi_grid))}
  
  if(plot & ball & is.null(N_rep)){
    print(plot_2+
            ggplot2::geom_vline(size = 0.2, colour = 'darkblue', linetype = 'dashed', xintercept = phi_grid)+
            ggplot2::geom_hline(size = 0.2, colour = 'darkblue', linetype = 'dashed', yintercept = psi_grid))}
  
  if(plot & ball & !is.null(N_rep)){print(plot_2)}
  
  #If no number-of-points restrictions are imposed, the number of empty regions after discretizing:
  #cat('\n',n_empty,'empty divisions out of',N_bin^2,'\n')
  
  return(pgrid)
  
}
