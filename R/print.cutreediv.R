#' @export 
#' @title Prints a 'cutdiv' object
#' @description This is a method for the function print for objects of the class \code{cutdiv}.
#' @param x  object of class cutdiv
#' @param ... further arguments passed from other methods. They are ignored in this function.
 
print.cutreediv <- function(x, ...){
  if (!inherits(x, "cutreediv")) 
    stop("use only with \"cutdiv\" objects")
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "") 
  cat("\n")
  
  res <- matrix("",4,2)
  colnames(res) <-c("name","description")
  res[1,] <- c("$clusters", "list of observations in each cluster")
  res[2,] <- c("$description", "monothetic description of each cluster")
  res[3,] <- c("$which_cluster", "cluster memberships")
  res[4,] <- c("$B", "the proportion of explained inertia")
  print(res)
}