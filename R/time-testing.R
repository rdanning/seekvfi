#' @export
time_seekvfi <- function(counts, Ks, parallel = FALSE, maxSize = 8 * 1024^3){
  
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
  
  start.time <- Sys.time()
  
  # column-normalize
  D <- prop.table(counts,2)
  
  # run topic models and convert output to joint hallmark projection matrices
  SVD.out <- run_svd(D,max(Ks))
  topic.matrices <- mapper1(unique(Ks), run_TopicScore, D, SVD.out, mapper2)
  loadings.matrices <- lapply(topic.matrices, prop.table, 1)
  sparsity.vectors <- lapply(loadings.matrices, get.svs)
  
  # ensemble sparsity vectors
  scores <- get.scores(sparsity.vectors)
  
  end.time <- Sys.time()
  
  # return output
  return(as.numeric(difftime(end.time,start.time,units="secs")))
}
