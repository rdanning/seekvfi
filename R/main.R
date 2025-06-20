#' Main function to run SEEK-VEC
#' 
#' @param D A vocabulary x documents count matrix.
#' @param Ks List of candidate values of K.
#' @param method Topic modeling method to use - one of c("ts" (default),"nmf","ctm","lda").
#' @param threshold A proportion between 0 and 1 or a number of top words per topic (default is 0.1).
#' @returns A SEEK.out object containing the input Ks; the topic, loadings, hallmark, and joint hallmark projection matrices for each value of K; the normalized eigenvector corresponding to the consensus matrix, and the ensembled joint hallmark projection matrix O.
#' @examples
#' run_SEEK(D, 3:12, method = "ts", threshold = 0.2)
#' run_SEEK(D, 3:12, method = "ts", threshold = 10)
#' @export
run_SEEK <- function(D, Ks, method = "ts", threshold = 0.1){
  
  # check for valid inputs
  check.inputs(D,Ks,method,threshold)
  
  # run topic models and convert output to joint hallmark projection matrices
  topic.matrices <- sapply(unique(Ks), get.topic.matrix, D, method)
  loadings.matrices <- lapply(topic.matrices, prop.table, 1)
  hallmark.matrices <- lapply(loadings.matrices, get.hallmark.matrix, threshold)
  hp.matrices <- lapply(hallmark.matrices, XXT)
  
  # run zero-reduce procedure to make ensembling more efficient
  zero.reduce.out <- zero.reduce(hp.matrices)
  reduced.hp.matrices <- zero.reduce.out$reduced.hps
  nonzero.idx <- zero.reduce.out$keep.idx
  
  # ensemble zero-reduced matrices
  SE.out <- SE.eigenscore(reduced.hp.matrices)
  O.reduced <- SE.out$O
  u <- SE.out$u
  
  # re-insert all-zero rows and columns
  O <- merge.reduced.O(O.reduced,nonzero.idx,nrow(D))
  
  # return output
  return(list(Ks = Ks,
              topic.matrices = topic.matrices,
              loadings.matrices = loadings.matrices,
              hallmark.matrices = hallmark.matrices,
              hp.matrices = hp.matrices,
              u = u,
              O = O))
}

# function to check for valid inputs
check.inputs <- function(D, Ks, method, threshold){
  if(sum(is.na(D)) > 0){
    stop("Values in counts matrix must be non-missing")
  }
  if(min(D) < 0){
    stop("Values in counts matrix must be non-negative")
  }
  if(!(all(Ks > 0)) | !(all(Ks == floor(Ks)))){
    stop("Candidate values of K must be positive integers")
  }
  if(!(tolower(method) %in% c("ts","nmf","lda","ctm"))){
    stop("Topic modeling method must be one of c('ts','nmf','lda','ctm')")
  }
  if(threshold < 0){
    stop("Threshold must be a loadings proportion between 0 and 1 or a number of top words 1 or greater")
  }
  if(threshold > 1 & (floor(threshold) != threshold)){
    stop("Number of top words must be an integer")
  }
}

# function to run topic modeling with the given method
get.topic.matrix <- function(K, D, method){
  print(paste0("Fitting a topic model with ",K," topics using ",toupper(method)))
  p <- nrow(D)
  if(method == "ts"){
    return(TopicScore::topic_score(K, prop.table(D,2))$A_hat)
  }
  if(method == "nmf"){
    return(fastTopics::fit_topic_model(t(D),k=K)$F)
  }
  if(method == "ctm"){
    ctm.fit <- topicmodels::CTM(t(D),k=K)
    ctm.topics <- tidytext::tidy(ctm.fit, matrix = "beta")
    ctm.topics$term <- rep(1:p,each=K)
    A.ctm <- tidyr::pivot_wider(ctm.topics, names_from = topic, values_from = beta)
    A.ctm <- A.ctm[,-which(names(A.ctm) == "term")]
    return(as.matrix(A.ctm))
  }
  if(method == "lda"){
    lda.fit <- topicmodels::LDA(t(D),k=K)
    lda.topics <- tidytext::tidy(lda.fit, matrix = "beta")
    lda.topics$term <- rep(1:p,each=K)
    A.lda <- tidyr::pivot_wider(lda.topics, names_from = topic, values_from = beta)
    A.lda <- A.lda[,-which(names(A.lda) == "term")]
    return(as.matrix(A.lda))
  }
}

# function to binarize loadings matrix
get.hallmark.matrix <- function(X, threshold){
  if(threshold < 1){
    vals <- as.numeric(X)
    threshold <- vals[order(vals,decreasing = TRUE)][threshold * ncol(X)]
    H <- X
    H[H < threshold] <- 0
    H[H >= threshold] <- 1
  } else{
    H <- sapply(1:ncol(X),extract.top.words,X,threshold)
  }
  return(H)
}

# function to select top n words from each column of a matrix
extract.top.words <- function(i,X,n){
  v <- X[,i]
  new <- rep(0,length(v))
  new[order(v, decreasing = TRUE)[1:n]] <- 1
  return(new)
}

# helper function to self transpose multiply
XXT <- function(X){
  XXT <- X %*% t(X)
  return(XXT)
}

# function to remove rows and columns of words that are uninformative across all candidate models
reduce.O <- function(i,O.list,keep.idx){
  O <- O.list[[i]]
  return(O[keep.idx,keep.idx])
}

# function to return list of zero-reduced hp matrices and the corresponding indices
zero.reduce <- function(hp.list){
  rs.df <- Reduce('+', hp.list)
  all.rs <- rowSums(rs.df)
  keep.idx <- which(all.rs != 0)
  reduced.hps <- lapply(1:length(hp.list), reduce.O, hp.list, keep.idx)
  return(list(reduced.hps = reduced.hps,
              keep.idx = keep.idx))
}

# helper function to re-insert zeroes
pad.cols <- function(i,O.reduced,v){
  if(v[i] == 0){
    return(rep(0,nrow(O.reduced)))
  } else{
    return(as.numeric(O.reduced[,v[i]]))
  }
}

# function to re-insert all-zero rows and columns
merge.reduced.O <- function(O.reduced,keep.idx,n){
  tf <- 1:n %in% keep.idx
  df <- data.frame(idx = 1:n,
                   tf = tf)
  v <- ave(df$idx, df$tf, FUN = seq_along)
  v[!(df$tf)] <- 0
  O.tmp <- sapply(1:n,pad.cols,O.reduced,v)
  O <- sapply(1:n,pad.cols,t(O.tmp),v)
  return(O)
}
