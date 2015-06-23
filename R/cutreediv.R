#' @title Cut the tree
#' @description 
#' This function cuts the tree into several cluters by specifying the desired number of clusters. 
#' @param tree the divclust object
#' @param K an integer with the desired number of clusters.
#' @return \item{clusters}{the list of observations in each  cluster}
#' @return \item{description}{the monothetic description of each cluster}
#' @return \item{which_cluster}{a vector of integers indicating the cluster of each observation} 
#' @return \item{B}{the proportion of inertia explained by the partition (between-cluster inertia/total inertia)}
#' @return \item{leaves}{an internal list of \strong{leaves}}
#' @export
#' @examples
#' data(protein) # pure quantitatives data
#' tree <- divclust(protein) # full clustering
#' p_5 <- cutreediv(tree,K=5) # partition in 5 clusters
#' p_5
#' 
#' data(dogs) # pure qualitative data
#' tree <- divclust(dogs) # full clustering
#' p_4 <- cutreediv(tree,K=4) # partition in 4 clusters

#' data(wine) # mixed data 
#' data <- wine[,1:29]
#' tree <- divclust(data) # full clustering
#' p_4 <- cutreediv(tree, 4) # 
#' 

cutreediv <- function(tree,K) { 
  cl <- match.call()
  cluster <- tree
  if (!inherits(cluster, "divclust")) 
    stop("use only with \"divclust\" objects")
  if (!(K %in% 1:(length(cluster$clusters)-1))) stop(paste("K must be an integer between 0 and",length(cluster$clusters)-1))
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
    #if (is.null(node$l$l)) {
    if (node$l$v$inert <= cluster$height[K]) {   
      l[[length(l) + 1]] <-  list(class = node$v$A_l, path = path_l) 
    } else {
      fifo_add(f, list(node = node$l, path = path_l))
    }
    if (node$r$v$inert <= cluster$height[K]) {
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
  part <- list()
  part$description <- make_description(l,cluster)
  names(part$description) <- paste("C", 1:K, sep = "")
  part$clusters <- lapply(l, function(x) {cluster$rnames[x$class]})
  names(part$clusters) <- paste("C", 1:K, sep = "")
  part$which_cluster <- c
  part$B <- sum(cluster$height[1:(K-1)])/cluster$T
  part$leaves <- l
  part$call <- cl
  class(part) <- "cutreediv"
  return(part)
}