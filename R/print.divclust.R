#' @export 
#' @title Prints a 'divclust' object
#' @description This is a method for the function print for objects of the class \code{divclust}.
#' @param x  object of class divclust
#' @param ... further arguments passed from other methods. They are ignored in this function.
 
print.divclust <- function(x, ...){
  if (!inherits(x, "divclust")) 
    stop("use only with \"divclust\" objects")
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "") 
  cat("\n")
  
  res <- matrix("",7,2)
  colnames(res) <-c("name","description")
  res[1,] <- c("$tree", "an internal tree")
  res[2,] <- c("$clusters", "list of observations in each cluster")
  res[3,] <- c("$description", "monothetic description of each cluster")
  res[4,] <- c("$which_cluster", "cluster memberships")
  res[5,] <- c("$height", "the clustering height")
  res[6,] <- c("$B", "the proportion of explained inertia")
  res[7,] <- c("$vec_quali", "number of categories of qualitative variables")
  print(res)
}