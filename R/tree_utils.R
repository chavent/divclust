#' @title list_leaves
#' @description 
#' This function finds terminal clusters, their monothetic description and the terminal partition. 
#' If the tree is complete, the terminal partition is the partition of singleton (one observation in each cluster) 
#' @param cluster the divclust object
#' @return leaves a list with the list of observations of a cluster ($class) and its monthetic description ($path) 
#' @return classes a vector of integers indicating the cluster to which each observation is allocated.
#' @keywords internal
#' @export
list_leaves <- function(cluster) { 
  l <- list()
  f <- fifo_create()
  #fifo_add(f, list(node = cluster$tree, path = ""))
  fifo_add(f, list(node = cluster$tree, path = NULL))
  while ( !fifo_is_empty(f)) {
    z <- fifo_remove(f)
    node <- z$node
    path <- z$path
    
    sentence <- make_path(node, path)
    path_l <- sentence$l
    path_r <- sentence$r
    #sentence <- make_sentences(cluster, node)
    #path_l <- paste(path, sentence$l, sep=', ')
    #path_r <- paste(path, sentence$r, sep=', ')
    if (is.null(node$l$l)) {
      l[[length(l) + 1]] <-  list(class = node$v$A_l, path = path_l) 
    } else {
      fifo_add(f, list(node = node$l, path = path_l))
    }
    if (is.null(node$r$l)) {
      l[[length(l) + 1]]  <- list(class = node$v$A_l_c, path = path_r)
    } else {
      fifo_add(f, list(node = node$r, path = path_r))
    }  
  }
  c <- rep(0, length(l)) 
  k <- 1 
  for (i in l) { 
    for (j in i$class) {
      c[j] <- k 
    } 
    k <- k + 1 
  }
  return(list(leaves = l, classes = c))
}

#' @title Builds two necessary conditions
#' @description Builds the monothetic description (path) of the two clusters of a bipartition
#' @param node current node
#' @param path of the current node 
#' @return the path of the two clusters of the biparition of the current node (left and rigth)
#' @keywords internal
#' @export

make_path <- function(node,path)
{
  if (node$v$cut_val$type=="quanti") {
    z <- Intervals(c(-Inf,node$v$cut_val$value))
    z_comp <- Intervals(c(node$v$cut_val$value,Inf))
  }
  if (node$v$cut_val$type=="quali") {
    z <- node$v$cut_val$value$bipart
    z_comp <- node$v$cut_val$value$bipart_c
  }
  
  if (is.null(path)) {
    j <- 1
    path_l <- list()
    path_r <- list()
  }
  else {
    path_l <- path
    path_r <- path
    vec_type <- unlist(lapply(path,function(x){node$v$cut_val$type==x$type}))
    vec_cut_ind <- unlist(lapply(path,function(x){node$v$cut_ind==x$cut_ind}))
    if (2 %in% (vec_type+vec_cut_ind))  {
      j <- which(vec_type+vec_cut_ind==2)
      if (node$v$cut_val$type=="quanti") {
        z <- interval_intersection(z,path[[j]]$val)
        z_comp <- interval_intersection(z_comp,path[[j]]$val)
      }
      if (node$v$cut_val$type=="quali") {
        z <- intersect(z,path[[j]]$val)
        z_comp <- intersect(z_comp,path[[j]]$val)
      }
    }
    else j <- length(vec_type)+1
  }
  
  path_l[[j]]<- list(type=node$v$cut_val$type,cut_ind=node$v$cut_ind,val = z)
  path_r[[j]] <- list(type=node$v$cut_val$type,cut_ind=node$v$cut_ind,val = z_comp)
  return(list(left = path_l, right = path_r))
}

#' @title Description of the leaves
#' @description Builds the monothetic description of the leaves
#' @param leaves a list of leaves
#' @param cluster the divclust object containing the tree 
#' @return The monothetic description of the leaves
#' @keywords internal
#' @export
#' 

make_description <- function(leaves,cluster)
{  
  mod_quali <- cluster$mod_quali
  cnames_quanti <- cluster$cnames_quanti
  cnames_quali <- cluster$cnames_quali
  description <- list()
  # suppress the two first character of a string
  cut_1 <- function(s) { Reduce(paste0,sapply(strsplit(s,""), identity)[-c(1,2)])} 
  
  for (j in 1:length(leaves))
  {
    description[[j]] <- character()
    for (k in 1:length(leaves[[j]]$path))
    {
      x <-leaves[[j]]$path[[k]]
      if (x$type =="quanti") {
        A=paste0("[",round(x$val@.Data[1],digits=2)," ; ",round(x$val@.Data[2],digits=2),"[")
        sentence <- paste0(cnames_quanti[x$cut_ind], " = ", A)
      }
      if (x$type =="quali") {
        select_mod <- mod_quali[[x$cut_ind]][x$val]
        select_str <- Reduce(function (x,y) {paste(x,y, sep=",")}, select_mod)
        symbol_in <- " = "
        sentence <- paste0(cnames_quali[[x$cut_ind]], symbol_in, "{", select_str, "}")  
      } 
      description[j] <- paste(description[[j]], sentence, sep=' , ')
    }    
    description[j] <- cut_1(description[[j]])
  }
  return(description)
} 




#' @title Builds two binary questions
#' @description Builds the two binary questions of a bipartition
#' @param cluster the divclust object containing the tree
#' @param node current node for which we wants to build the questions
#' @return a pair of the two complementaries binary questions (left and rigth)
#' @keywords internal
#' @export

make_sentences <- function(cluster, node) {
  
  mod_quali <- cluster$mod_quali
  cnames_quanti <- cluster$cnames_quanti
  cnames_quali <- cluster$cnames_quali
  index <-  node$v$cut_ind
  
  if (node$v$cut_val$type =="quanti") {
    #value <-  node$v$cut_val$value
    #cut_val <- value * cluster$sig_quanti[index] + cluster$mu_quanti[index]
    cut_val <-  node$v$cut_val$value
    sentence <- paste0(cnames_quanti[index], " < ", round(cut_val, digits=2))
    sentence_comp <- paste0(cnames_quanti[index], " >= ", round(cut_val, digits=2))
  }   
  
  else if (node$v$cut_val$type =="quali") {
    value <-  node$v$cut_val$value
    select_mod <- mod_quali[[index]][value$bipart]
    select_comp <- setdiff(mod_quali[[index]][value$bipart_c], select_mod)
    select_str <- Reduce(function (x,y) {paste(x,y, sep=",")}, select_mod)
    select_comp_str <- Reduce(function (x,y) {paste(x,y, sep=",")}, select_comp)
    #symbol_in <- " \u220A "
    symbol_in <- " = "
    sentence <- paste0(cnames_quali[[index]], symbol_in, "{", select_str, "}") 
    sentence_comp <- paste0(cnames_quali[[index]], symbol_in,"{", select_comp_str, "}") 
  }
  
  else { 
    stop("unknown cut_val$type : must be quali or quanti") 
  }
  return(list(left = sentence, right = sentence_comp))
}

#' @title Binary transformation of categorical data
#' @description Organizes categorical data into an indicator matrix and gives a description of the categorical variables 
#' @keywords internal
#' @param frame a data frame with n rows and p categorical colums (columns are of class factor)
#' @return \item{Y}{an indicator matrix with n rows and m columns where m is the total number of categories} 
#' @return \item{vec_quali}{a vector of size p containing the number of categories of each variable}
#' @return \item{mod_quali}{a list with the set of categories of each variable}
#' @export 
format_cat_data <- function(frame) {
  #X<-as.matrix(frame)
  X <- frame
  nb_vars <- ncol(X)
  cnames_quali <- colnames(frame)
  #mod_quali <- apply(X, 2, unique)
  mod_quali <- lapply(X,function(x){levels(x)})
  vec_quali <- sapply(mod_quali, length)
  nb_cols <- sum(vec_quali)  
  offsets <- c(0, cumsum(vec_quali))
  # builds the disjonctif table
  Y<-matrix(0, nrow(X), nb_cols)
  for(var in seq_len(nb_vars))
    for(mod in seq_len(length(mod_quali[[var]])))
      Y[ X[,var] == mod_quali[[var]][[mod]], offsets[var] + mod] <-1
  return(list(Y = Y, vec_quali = vec_quali, mod_quali = mod_quali, cnames_quali = cnames_quali))
}


#' @title Generalized Singular Value Decomposition
#' @description Performs the generalized singular value decomposition \eqn{A=UDV} with weights M on the columns and N on the rows.
#' @param A a n times p numerical matrix.
#' @param M a vector of size p with diag(M)  the metric on \eqn{R^p}.
#' @param N a vector of size n with diag(N)  the metric on \eqn{R^n}.
#' @return a 3-tuple \eqn{(U,V,d)} with the left singular vectors, right singular vectors and singular values.
#' @export
gsvd <- function(A, M, N) {
  if (!is.vector(N))
    stop('N must be a vector')
  if (!is.vector(M))
    stop('M must be a vector')
  z_t <- diag(sqrt(N)) %*% A %*% diag(sqrt(M))
  sig_t <- svd(z_t)
  u <- diag(1./sqrt(N)) %*% sig_t$u 
  v <- diag(1./sqrt(M)) %*% sig_t$v
  return (list(u = u, v = v, d = sig_t$d))
}

#' @title 1D inertia
#' @description Computes the 1D inertia of a column of a matrix data on a subset of rows of a data matrix.
#' @param Z matrix data
#' @param compo column number
#' @param indices vector representing the subset of rows
#' @param w weight vector
#' @param D diagonal distance matrix
#' @keywords internal
#' @export

inert_1D <- function (Z, compo, indices, w = rep(1. / nrow(Z), nrow(Z)), D = rep(1, ncol(Z))) {
  if (length(indices) <= 1)
    return (0)
  subZ <- Z[indices, compo]
  subw <- w[indices]
  g <- sum(subZ * subw)/sum(subw)
  subZ  <- subZ - g
  inert <- D[compo] * sum(subw * subZ^2)
  return(inert)
}


#' @title 1D between-clusters inertia
#' @description Computes 1D between-clusters inertia of a bipatition \eqn{(C_1, C_2)} of a subset of rows of a data matrix.
#' @param Z matrix data
#' @param compo columns number
#' @param indices indices vectors of the first class (\eqn{C_1})
#' @param indices_c indices vectors of the second class (\eqn{C_2})
#' @param w weights vectors
#' @param D diagonal distance matrix vector 
#' @keywords internal
#' @export
#' 
between_cluster_inert_1D <- function (Z, compo, indices, indices_c, w = rep(1. / nrow(Z), nrow(Z)), D = rep(1, ncol(Z))) 
{ 
  # checks there is not a empty class
  if (is.null(indices) | is.null(indices_c)) stop('between_cluster_inert : one of the two cluster is empty!')
  mu <- sum(w[indices])
  g <- sum(w[indices] * Z[indices, compo]) / mu
  
  mu_c <- sum(w[indices_c])
  g_c <- sum(w[indices_c] * Z[indices_c, compo]) / mu_c 
  
  inert <- D[compo] * (g - g_c)^2
  inert <- (mu_c * mu) / (mu_c + mu) * inert
  return(inert)
}

#' @title Between-clusters inertia
#' @description Performs the between-clusters inertia of a biparition \eqn{(C_1, C_2)} of a subset of rows of a data matrix.
#' @param Z a numerical data matrix.
#' @param indices indices vectors of the first class (\eqn{C_1})
#' @param indices_c indices vectors of the second class (\eqn{C_2})
#' @param w weights vector
#' @param D diagonal distance matrix vector 
#' @keywords internal
#' @export

between_cluster_inert <- function (Z, indices, indices_c, w = rep(1. / nrow(Z), nrow(Z)), D = rep(1, ncol(Z))) 
{ 
  # checks there is not a empty class
  if (is.null(indices) | is.null(indices_c)) stop('between_cluster_inert : one of the two cluster is empty!')
  if (length(indices) > 1) {
    mu <- sum(w[indices])
    g <- colSums(w[indices] * Z[indices,]) / mu
  }
  else {
    mu <-  w[indices]
    g <- Z[indices,]
  }
  if (length(indices_c) > 1) {
    mu_c <- sum(w[indices_c])
    g_c <- colSums(w[indices_c] * Z[indices_c,]) / mu_c 
  }
  else {
    mu_c <-  w[indices_c]
    g_c <- Z[indices_c,]
  }
  inert <- sum( D * (g - g_c)^2)
  inert <- (mu_c * mu) / (mu_c + mu) * inert
  return(inert)
}


#' @title Computes the cut values for a numerical variable
#' @description Function which computes the cut values of the binary questions defined 
#' with a numerical variables (a column of the data matrix).
#' @param col the column vector in the matrix
#' @return a vector of cut values
#' @keywords internal
#' @seealso \link{choose_question}
#' @export
compute_cut_values <- function(col) {
  cuts <- unique(sort(col))
  if (length(cuts) >= 2) {
    cuts <- 0.5*(cuts[-1]+cuts[1:(length(cuts)-1)])
  }
  return (cuts)
}

#' @export
#' @title computes_height
#' @param tree an internal tree
#' @keywords internal
computes_height <- function (tree) {
   if (!is.null(tree)) {
      return (1+max(computes_height(tree$l), computes_height(tree$r)))
   }
   else {
       return(0)
   }
}


