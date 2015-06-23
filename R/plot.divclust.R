#' @title Dendrogram with binary questions
#' @method plot divclust
#' @description Plot the  dendrogram produced by divclust algorithm with the binary questions
#' associated with each split. The number of binary questions drawn on the dendrogram can be changed
#' and if the text of the question is to long, is is printed on the console.
#' @param x  object of class divclust
#' @param nqbin an integer between 0 and K-1 (K is the number of leaves in the tree). Indicates
#' the number of binary questions drawn on plot. This parameter is used to plot the binary questions
#' only of the top levels of the dendrogram. The default value is 3.
#' @param label If TRUE, the labels of the observations are drawn.
#' @param ... further arguments passed from other methods
#' @export
#' @seealso  \link{divclust}
#' @examples
#' data(protein) # pure quantitatives data
#' c_tot <- divclust(protein) # full clustering
#' plot(c_tot)
#' plot(c_tot,nqbin=4) # the text of the 4th question is printed in the console
#' c_4 <- divclust(protein, K=4) # stops the clustering to 4 clusters  
#' plot(c_4)
#'
#' data(dogs) # pure qualitative data
#' c_tot <- divclust(dogs) # full clustering
#' plot(c_tot,nqbin=4) # the text of the 4th question is printed in the console
#'
#' data(wine) # mixed data 
#' data <- wine[,1:29]
#' c_tot <- divclust(data) # full clustering
#' plot(c_tot)

plot.divclust <- function(x, nqbin=3,label=TRUE,...) {
  if (!inherits(x, "divclust")) 
    stop("use only with \"divclust\" objects")
  plot.new() 
  if ((x$kmax-1) < nqbin) nqbin <- x$kmax-1
  if (!(nqbin %in% 1:(length(x$clusters)-1))) stop(paste("nqbin must be an integer between 0 and",length(x$clusters)-1))
  usr <- par("usr")
  x_min <- usr[1]
  x_max <- usr[2]
  y_min <- usr[3]
  y_max <- usr[4]
  delta_x <- diff(usr[1:2])
  delta_y <- diff(usr[3:4])
  delta_inert_max <- x$tree$v$inert
  height_unit = 0.7 * delta_y / delta_inert_max
  deep <- computes_height(x$tree) 
  width <- 0.5 * delta_x / (1 - 0.25^(deep+1))
  draw(x, x_min + 0.25 * delta_x, y_min + 0.9 * delta_y , width, height_unit, nqbin,label,delta_x/5,...)
}


#' @export
#' @title Draw the tree
#' @description Draw a tree representing the hierarical structure 
#' @param cluster object returned by divclust function
#' @param x_pos initial abcissa for drawing
#' @param y_pos initial ordinate for drawing
#' @param width elementary horizontal length for the tree
#' @param height elementary vertically length for the tree
#' @param threshold cut-off value to decide if the text is written in the tree
#' @param ... further arguments passed from other methods
#' @keywords internal


draw <- function(cluster, x_pos, y_pos, width, height, nqbin, label,threshold,...) {
  
  rnames <- cluster$rnames
  bool_name <- (cluster$kmax == length(cluster$cluster)) 
  width_init <- width
  height_init <- height 
  chh <- par()$cxy[ 1 ]  ##  character height

  rec_draw <- function(tree, x_pos, y_pos, width, height,...) {
    if (!is.null(tree)) {
      segments(x_pos, y_pos, x_pos + width , y_pos)
      if (((nqbin!=0) & (tree$v$inert %in% cluster$height[1:nqbin]))){
        sentence <- make_sentences(cluster, tree)
        if (width < threshold) {
          k <- which(cluster$height[1:nqbin]==tree$v$inert)
          text(x_pos, y_pos + chh, paste0("Q",k,"=True"),adj=c(0.5,0),...)
          text(x_pos + width, y_pos + chh, paste0("Q",k,"=False"),adj=c(0.5,0),...)
          #print(paste(paste0("Q",k),":",sentence$left, "or",sentence$right ))
          message(paste(paste0("Q",k),":",sentence$left, "or",sentence$right ))
        } else {
          text(x_pos, y_pos + chh, sentence$left,,adj=c(0.5,0),...)
          text(x_pos + width, y_pos + chh, sentence$right,,adj=c(0.5,0),...)
        }
      }
    } 
    
    if (!is.null(tree$l$l)) {
      h_l = tree$v$inert -  tree$l$v$inert
      rec_draw(tree$l, x_pos - width/4., y_pos - height * h_l, width / 2., height,...)  
    }
    else {
      h_l <- tree$v$inert
      if (label==TRUE) {
        if ((bool_name)) {
          sentence = rnames[tree$v$A_l[1]]
        } 
        else {
          sentence = paste0("C",
                           as.character(cluster$which_cluster[tree$v$A_l[1]]))
        }
        text(x_pos-0.03*width_init, y_pos - height * h_l-0.02, sentence, srt = 270, pos=4,...)
      }
    } 
    segments(x_pos, y_pos, x_pos, y_pos - height * h_l)
    
    if (!is.null(tree$r$l)) {
      h_r = tree$v$inert -  tree$r$v$inert
      rec_draw(tree$r, x_pos + width - width/4, y_pos - height * h_r, width / 2., height,...)
    }
    else {
      h_r <- tree$v$inert
      if (label==TRUE) {
        if (bool_name) {
          sentence = rnames[tree$v$A_l_c[1]]
        }
        else {
          sentence = paste("C",
                           as.character(cluster$which_cluster[tree$v$A_l_c[1]]),sep="")
        }
        text(x_pos+width-0.03*width_init, y_pos - height * h_r-0.02, sentence, srt = 270, pos = 4,...)
      }
    }
    segments(x_pos + width, y_pos, x_pos + width , y_pos -  height * h_r)
  }
  rec_draw(cluster$tree, x_pos, y_pos, width, height,...) 
}


