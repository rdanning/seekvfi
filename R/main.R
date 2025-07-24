#' Main function to run SEEK-VFI
#'
#' @param D A genes x cells expression count matrix with gene names as the rownames.
#' @param Ks List of candidate values of K.
#' @returns A data frame containing the SEEK-VFI scores of each gene
#' @examples
#' seekvfi(D, 3:12)
#' @export
run_seekvfi <- function(D, Ks){

  # check for valid inputs
  check.inputs(D,Ks)

  # run topic models and convert output to joint hallmark projection matrices
  SVD.out <- run_svd(D,max(Ks))
  topic.matrices <- sapply(unique(Ks), get.topic.matrix, D, SVD.out)
  loadings.matrices <- lapply(topic.matrices, prop.table, 1)
  sparsity.vectors <- lapply(loadings.matrices, get.svs)

  # ensemble sparsity vectors
  scores <- get.scores(sparsity.vectors)

  # return output
  return(data.frame(gene = rownames(D),
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
get.topic.matrix <- function(K, D, SVD.out){
  print(paste0("Fitting a topic model with ",K," topics"))
  return(run_TopicScore(K, prop.table(D,2),SVD.out))
}

# helper function to extract sparsity, rescaled from [1/K,1] to [0,1]
get.sparsity <- function(i,m){
  v <- m[i,]
  K <- ncol(m)
  s <- (sum(v*v)-(1/K))/(1-(1/K))
  return(s)
}

# function to extract sparsity vectors
get.svs <- function(m){
  sv <- sapply(1:nrow(m),get.sparsity,m)
  return(sv)
}


