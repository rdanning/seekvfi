#' Main function to run SEEK-VFI
#'
#' @param D A genes x cells expression count matrix with gene names as the rownames.
#' @param Ks List of candidate values of K.
#' @param im iter.max parameter for kmeans (default 100)
#' @param ns start parameter for kmeans (default 1)
#' @param seed Random seed (default 1)
#' @returns A data frame containing the SEEK-VFI scores of each gene
#' @examples
#' run_seekvfi(D, 3:12)
#' @export
run_seekvfi <- function(counts, Ks, ns = 1, im = 100, seed = 1){

  set.seed(seed)

  # check for valid inputs
  check.inputs(counts,Ks)

  # column-normalize
  D <- prop.table(counts,2)

  # run topic models and convert output to joint hallmark projection matrices
  SVD.out <- run_svd(D,max(Ks))
  topic.matrices <- sapply(unique(Ks),
                           get.topic.matrix,
                           D,SVD.out,ns,im,seed)
  loadings.matrices <- lapply(topic.matrices, prop.table, 1)
  sparsity.vectors <- lapply(loadings.matrices, get.svs)

  # ensemble sparsity vectors
  scores <- get.scores(sparsity.vectors)

  # return output
  return(data.frame(gene = rownames(counts),
                    SV.score = scores))
}

# function to check for valid inputs
check.inputs <- function(D, Ks){
  if(sum(is.na(D)) > 0){
    stop("Values in counts matrix must be non-missing")
  }
  if(min(D) < 0){
    stop("Values in counts matrix must be non-negative")
  }
  if(is.null(rownames(D))){
    stop("Matrix rownames must contain gene names")
  }
  if(!(all(Ks > 0)) | !(all(Ks == floor(Ks)))){
    stop("Candidate values of K must be positive integers")
  }
}

# function to run topic modeling with the given K
get.topic.matrix <- function(K, D, SVD.out, ns, im, seed){
  print(paste0("Fitting a topic model with ",K," topics"))
  return(run_TopicScore(K,D,SVD.out,ns,im,seed))
}

# helper function to extract sparsity
get.sparsity <- function(i,m){
  v <- m[i,]
  s <- sum(v*v)
  return(s)
}

# function to extract sparsity vectors
get.svs <- function(m){
  sv <- sapply(1:nrow(m),get.sparsity,m)
  return(sv)
}


