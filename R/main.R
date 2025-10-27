#' Main function to run SEEK-VFI
#'
#' @param D A genes x cells expression count matrix with gene names as the rownames.
#' @param Ks List of candidate values of K.
#' @param parallel Whether to parallelize within Topic-SCORE code: one of c(FALSE, "cluster", "multisession"). Default value is FALSE.
#' @returns A data frame containing the SEEK-VFI scores of each gene
#' @examples
#' seekvfi(D, 3:12)
#' @export
run_seekvfi <- function(counts, Ks, parallel = FALSE, maxSize = 8 * 1024^3){

  # plan future if running in parallel
  if(parallel == "cluster"){
    future::plan(future::cluster)
    options(future.globals.maxSize = maxSize)
    mapper1 <- furrr::future_map
    mapper2 <- furrr::future_map2_dbl
  }
  if(parallel == "multicore"){
    future::plan(future::multicore)
    options(future.globals.maxSize = maxSize)
    mapper1 <- furrr::future_map
    mapper2 <- furrr::future_map2_dbl
  }
  if(parallel == "multisession"){
    future::plan(future::multisession)
    options(future.globals.maxSize = maxSize)
    mapper1 <- furrr::future_map
    mapper2 <- furrr::future_map2_dbl
  } else{
    mapper1 <- purrr::map
    mapper2 <- purrr::map2_dbl
  }

  # check for valid inputs
  check.inputs(counts,Ks)

  # column-normalize
  D <- prop.table(counts,2)

  # run topic models and convert output to joint hallmark projection matrices
  SVD.out <- run_svd(D,max(Ks))
  topic.matrices <- mapper1(unique(Ks), run_TopicScore, D, SVD.out, mapper2)
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
get.topic.matrix <- function(K, D, SVD.out, mapper2){
  print(paste0("Fitting a topic model with ",K," topics"))
  return(run_TopicScore(K, D, SVD.out, mapper2))
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


