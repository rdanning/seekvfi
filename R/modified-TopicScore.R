#' Lightly-modified TopicScore code
#' Adapted from https://github.com/cran/TopicScore/blob/master/R/topic_score.R
#' Optimizes running multiple models by calculating the singular values once; parallelizes some of the vertex hunting steps
#' TopicScore paper: Ke, Z. T. & Wang, M. Using SVD for Topic Modeling. Journal of the American Statistical Association 119, 434–449. http://dx.doi.org/10.1080/01621459.2022.2123813 (Oct. 2022).


vertices_est <- function(R,K0,m,num_start,mapper2){
  K <- dim(R)[2] + 1

  obj <- kmeans(R,m,iter.max=100,nstart=num_start)
  theta <- as.matrix(obj$centers)
  theta_original <- theta

  inner <- theta%*%t(theta)
  distance <- diag(inner)%*%t(rep(1,m)) + rep(1,m)%*%t(diag(inner)) - 2*inner
  top2 <- which(distance==max(distance),arr.ind=TRUE)[1,]
  theta0 <- as.matrix(theta[top2,])
  theta <- as.matrix(theta[-top2,])

  if (K0 > 2){
    for (k0 in 3:K0){
      inner <- theta%*%t(theta)
      distance <- rep(1,k0-1)%*%t(diag(inner))-2*theta0%*%t(theta)
      ave_dist <- colMeans(distance)
      index <- which(ave_dist==max(ave_dist))[1]
      theta0 <- rbind(theta0, theta[index,])
      theta <- as.matrix(theta[-index,])
    }
    theta <- theta0
  }

  comb <- combn(1:K0, K)
  max_values <- rep(0, dim(comb)[2])

  is <- 1:dim(comb)[2]
  js <- 1:K0
  grid <- matrix(mapper2(rep(is,times = length(js)),
                   rep(js, each = length(is)),
                   simplex_dist_parallel,
                   theta,
                   comb),
                 nrow = length(is))
  max_values <- apply(grid,1,max)

  min_index <- which(max_values == min(max_values))

  return(list(V=theta[comb[,min_index[1]],], theta=theta_original))

}

simplex_dist_parallel <- function(i,j, theta0, comb){
  theta <- as.matrix(theta0[j,])
  V <- as.matrix(theta0[comb[,i],])

  VV <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))%*%V
  D <- VV%*%t(VV)
  d <- VV%*%(theta-V[dim(V)[1],])

  A <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))
  b0 <- c(rep(0,dim(V)[1]-1),-1)

  obj <- quadprog::solve.QP(D, d, A, b0)
  return(sum((theta-V[dim(V)[1],]) ^2)+ 2*obj$value)
}

run_svd <- function(D, max.K, Mquantile=0){
  print("Computing SVD")
  M <- rowMeans(D)
  M_trunk <- pmin(M,quantile(M,Mquantile,na.rm = TRUE),na.rm=TRUE)
  obj <- RSpectra::svds(sqrt(M_trunk^(-1))*D, max.K)
  Xi <- obj$u
  return(list(Xi=Xi,
              M_trunk=M_trunk))
}

run_TopicScore <- function(K, D, SVD.out, mapper2, Mquantile=0, num_start = 1){

  Xi <- SVD.out$Xi
  M_trunk <- SVD.out$M_trunk

  K0 <- ceiling(1.5*K)
  m <- 10*K

  p <- dim(D)[1]
  n <- dim(D)[2]

  Xi <- Xi[,1:K]

  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(Xi[,2:K],2,function(x) x/Xi[,1])

  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m,num_start,mapper2)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta

  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,0)
  Pi <- prop.table(Pi,1)

  #Step 4
  A_hat <- sqrt(M_trunk)*Xi[,1]*Pi

  #Step 5
  A_hat <- prop.table(A_hat,2)

  return(A_hat)
}






