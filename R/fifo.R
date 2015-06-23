fifo_create <- function() {
   e <- new.env()
   e$l <- list()
   return(e)
}

fifo_add <- function(f, a) {
    f$l <- append(f$l, list(a))
}

fifo_remove <- function(f) {
  if (length(f$l) >= 1) {
     a <- f$l[[1]]
     f$l[[1]] <- NULL
     return(a)
  }
  else {
    stop("fifo_remove : no element to remove")
  }
}

fifo_is_empty <- function(f) {
   return (length(f$l) == 0)
}
