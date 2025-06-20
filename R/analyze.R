#' Function to plot heatmap of results
#' 
#' @param SEEK.out Output of run_SEEK.
#' @param vocabulary Vector of vocabulary in the order corresponding to the input corpus matrix.
#' @param filter.value Vocabulary words with no off-diagonal value greater than or equal to filter.value will be omitted (default is 0)
#' @returns A pheatmap object.
#' @examples
#' plot_results(SEEK.out, vocabulary, filter.value = 0.5)
#' @export
plot_results <- function(SEEK.out, vocabulary, filter.value = 0, font.size = 8){
  O <- SEEK.out$O
  rownames(O) <- vocabulary
  colnames(O) <- vocabulary
  
  tmp <- O
  tmp[tmp < filter.value] <- 0
  diag(tmp) <- 0
  keep.idx <- rowSums(tmp) > 0
  heatmap.m <- O[keep.idx,keep.idx]
  
  pm <- pheatmap::pheatmap(heatmap.m,
                           color = colorRampPalette(c("white", "navy"))(100),
                           breaks = seq(0,1,by=0.01),
                           border_color = NA,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           fontsize = font.size,
                           legend = FALSE)
  return(pm)
}


#' Function to extract diagonal values from ensembled joint hallmark projection matrix
#' 
#' @param SEEK.out Output of run_SEEK.
#' @param vocabulary Vector of vocabulary in the order corresponding to the input corpus matrix.
#' @returns A dataframe containing the vocabulary and diagonal values.
#' @examples
#' extract_diagonals(SEEK.out, vocabulary)
#' @export
extract_diagonals <- function(SEEK.out,vocabulary){
  O <- SEEK.out$O
  df <- data.frame(word = vocabulary, value = diag(O))
  return(df)
}


#' Function to assess the stability of a given model with respect to the ensemble matrix
#' 
#' @param SEEK.out Output of run_SEEK.
#' @param vocabulary Vector of vocabulary in the order corresponding to the input corpus matrix.
#' @param K K of interest.
#' @param n.top.words Number of top words per topic to analyze (default is 10).
#' @returns A pheatmap object and the stability measures of the block of top words from each topic.
#' @examples
#' assess_stability(SEEK.out, vocabulary, 6, n.top.words = 10)
#' @export
assess_stability <- function(SEEK.out,vocabulary,K,n.top.words = 10,font.size=8){
  set.seed(1)
  Ks <- SEEK.out$Ks
  B <- SEEK.out$loadings.matrices[[which(Ks == K)]]
  O <- SEEK.out$O
  rownames(O) <- vocabulary
  colnames(O) <- vocabulary
  
  top.words <- sapply(1:ncol(B),extract.top.words,B,n.top.words)
  rownames(top.words) <- vocabulary
  TWT <- top.words %*% t(top.words)
  TWT0 <- TWT
  diag(TWT0) <- 0
  order <- hclust(dist(TWT0))$order
  keep.idx <- intersect(order,which(diag(TWT) > 0))
  
  heatmap.m <- O[keep.idx,keep.idx]
  pm <- pheatmap::pheatmap(heatmap.m,
                           color = colorRampPalette(c("white", "navy"))(50),
                           treeheight_row = 0,
                           treeheight_col = 0,
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           fs = font.size,
                           legend = FALSE)
  
  # average weight within each block
  block.stability <- lapply(1:K,get.block.stability,heatmap.m,top.words)
  
  return(list(plot = pm,
              block.stability = block.stability))
  
}


# helper function to get n top words for each topic from the given loadings matrix
extract.top.words <- function(i,B,n){
  v <- B[,i]
  new <- rep(0,length(v))
  new[order(v, decreasing = TRUE)[1:n]] <- 1
  return(new)
}

# helper function to calculate the stability of a given topic
get.block.stability <- function(i,m,top.words){
  block.words <- row.names(top.words[top.words[,i]==1,])
  idx <- which(rownames(m) %in% block.words)
  block <- m[idx,idx]
  return(list(block = block.words,
              stability = mean(block)))
}







