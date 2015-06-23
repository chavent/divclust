#' tests if two lists of vectors represent or not the same partition
#' @param l1 first partition
#' @param l2 first partition
#' @return a logical 
#' @keywords internal
# complexity : O(min(|l1|,|l2|)

partitionsequal <-function (l1,l2) {
   n1 = length(l1)
   n2 = length(l2)
   if (n1 == 0 & n2 == 0) {
        return(TRUE)
   }
   else if(n1 ==0 & n2 >= 1) {
       return(FALSE)
   }
   else if(n1 >=1 & n2 == 0) {
        return(FALSE)
   }
   else {
        pos = Position(function(l) {setequal(l1[[1]],l)},l2)
       if (!is.na(pos))
          {
           l1[[1]] <- NULL
           l2[[pos]] <-NULL
           return(partitionsequal(l1,l2))
          }
       else  {
        return(FALSE)  
      } 
  }
}   
