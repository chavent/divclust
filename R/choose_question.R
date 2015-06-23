#' @title Choice of the best binary question 
#' @description  This function finds the best binary question do divide a cluster A  
#' into to subclusters such that the bipartition (A_l, A_l_c) has maximum between-clusters inertia.
#' @details This function works for both categorical, numerical and mixed data. 
#' This is the core function of the divclust algorithm.
#' We are seeking the binary question which gives the best 
#' bipartition. A binary question is defined with a cutting variable (quantitative or qualitative), and 
#' a cutting value. For quantitative variable, the cutting value is a real number. For qualitative, the
#' cutting value is
#' @param X the data matrix of dimension (nxp) where p is equal to the number of numerical variables
#' plus the number of categories. This matrix is used to construct the binary questions
#' @param Z the numerical data matrix of dimension (nxk) used to compute the inertia criterion 
#' (the matrix of the principal components for instance)
#' @param indices vector of indices for the cluster A to divide.
#' @param vec_quali vector containing the number of categories for each modalities (according to the
#' categories observed in A_l)
#' @param w weights vector
#' @param D diagonal distance matrix coefficients
#' @param vec_order vector containing TRUE if the categories of the variable are ordered
#' @return \item{inert}{the between-clusters inertia of the bipartition (A_l, A_l_c)}
#' @return \item{A_l}{the vector of indices of the cluster A_l}
#' @return \item{A_l}{the vector of indices of the cluster A_l_c}
#' @return \item{cut_ind}{the index of the cutting variables}
#' @return \item{cut_val}{a list with : \itemize{
#' \item $type: the type of the cutting variable (quantitative or qualitative)
#' \item $value: a real value if the variable is quantitative and a biparition of categories 
#'  if the variable is qualitative. 
#'  The first cluster of this bipartition of categories is given by $value$bipart whereas the second
#' cluster is given by $value$bipart_c}}
#' @keywords internal  
#' @export 

choose_question <- function (X, Z, indices, vec_quali = c(), w = rep(1. / nrow(Z), nrow(Z)), 
                             D = rep(1, ncol(Z)),vec_order) {
  l = length(indices)
  #if (l <= 1)
  #   stop(' choose_question : cannot divide singleton or empty sets')
  p <- dim(X)[2]
  inert_max <- 0
  second_inert <- 0
  c_max <- 0
  j_max <- -1
  A_l_max <- c()
  A_l_c_max <- c()
  nb_quanti <- p - sum(vec_quali) 
  nb_quali <- length(vec_quali)
  if (l > 1) {
    # first we list all the quantitative variables
    for (j in seq_len(nb_quanti)) {
      cut_vals <- compute_cut_values(X[indices, j])
      inerts <- sapply(cut_vals, function(l) { between_cluster_inert(Z, indices[X[indices, j] <= l], indices[X[indices, j] > l], w, D) }) 
      l_max <- which.max(inerts) 
      if ( length(l_max)== 1 && inerts[l_max] >=inert_max) { # beware, think to add a empty test
        denom <- inert_1D(X, j, indices, w) 
        nume <- between_cluster_inert_1D(X, j , indices[X[indices, j] <= cut_vals[l_max]], indices[X[indices, j] > cut_vals[l_max]], w)
        second_inert_candidat <- nume / denom
        if (inerts[l_max] == inert_max) { # beware, may be a tolerance test would be better
          if (second_inert_candidat >second_inert ) {
            c_max <- list(type = "quanti", value = cut_vals[l_max])
            j_max <- j
            second_inert <- second_inert_candidat
            inert_max <- inerts[l_max]
            A_l_max = indices[X[indices, j_max] <= c_max$value]
            A_l_c_max = indices[X[indices, j_max] > c_max$value]
          }
        } # equality case end
        else {
          second_inert <- second_inert_candidat
          j_max <- j 
          c_max <- list(type = "quanti", value = cut_vals[l_max])
          inert_max <- inerts[l_max]
          A_l_max = indices[X[indices, j_max] <= c_max$value]
          A_l_c_max = indices[X[indices, j_max] > c_max$value]
        } 
      }
    } 
    
    
    start_offsets = c(0, cumsum(vec_quali)) + nb_quanti 
    
    for(j in seq_len(nb_quali)) {
      nb_mod = vec_quali[[j]] 
      to_keep = seq_len(nb_mod)[apply(X[indices, start_offsets[[j]] + seq_len(nb_mod)],2,sum) != 0]
      if (vec_order[j]==TRUE)
      {
        biparts <- list()
        if (length(to_keep) >1)
          biparts = lapply(1:(length(to_keep)-1),function(x){to_keep[1:x]})
      } else 
        biparts = bipart(length(to_keep))
    
      nb_biparts = length(biparts)
      
      select_quali <- function(bp_index) {
        A_l <- indices[apply(as.matrix(X[indices, start_offsets[[j]] + to_keep[biparts[[bp_index]]]] == 1), 1, any)]
        A_l_c <- indices[ apply(as.matrix(X[indices, start_offsets[[j]] + to_keep[biparts[[bp_index]]]] == 0), 1, all)]
        return(between_cluster_inert(Z, A_l, A_l_c, w, D))}
      
      
      inerts <- lapply(seq_len(nb_biparts), select_quali)
      l_max <- which.max(inerts)
      if (length(l_max) >=1 && inerts[[l_max]] > inert_max) {
        inert_max <- inerts[[l_max]]
        value <- to_keep[biparts[[l_max]]]
        bipart_c <- setdiff(to_keep, value) 
        c_max <- list(type = "quali", value = list(bipart=value,bipart_c = bipart_c))
        j_max <- j 
        A_l_max <- indices[apply(as.matrix(X[indices, start_offsets[[j]] + c_max$value$bipart] == 1), 1, any)]
        A_l_c_max <- indices[apply(as.matrix(X[indices, start_offsets[[j]] + c_max$value$bipart] == 0), 1, all)]
      }
      
    }
  }
  return (list(inert = inert_max, A_l = A_l_max, A_l_c = A_l_c_max, cut_ind = j_max, cut_val = c_max))
}

