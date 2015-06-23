#' @title Splits in quantitative and qualitative datasets
#' @description Splits a dataframe in two data frames. The first one contains the numerical columns
#' and the second one contains the categorical columns.
#' @param base the data
#' @return \item{data_quanti}{the numerical data frame}
#' @return \item{data_quali}{the categorical data frame}
#' @examples
#' data(wine)
#' split_mix(wine)$data_quanti
#' split_mix(wine)$data_quali
#' @export

split_mix <- function (base) 
{
  type <- NULL
  base <- data.frame(base, check.names = T)
  j_quant <- unlist(lapply(base,is.numeric))
  
  if (sum(j_quant)!=0)
    data_quanti <- base[, j_quant,drop=FALSE] 
  else data_quanti=NULL
  if (sum(!j_quant)!=0)
    data_quali <- base[, !j_quant,drop=FALSE]
  else data_quali=NULL
  return(list(data_quanti=data_quanti,data_quali=data_quali))
}