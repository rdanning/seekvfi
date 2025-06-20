# get element of similarity matrix
get.G.ij <- function(i,cv1.do,cv2.do,sv.list,total.combos){
  pct <- 100*round(i/total.combos, 3)
  print(paste0("Calculating similarity matrix: ",pct,"% complete"))
  A <- sv.list[[cv1.do[i]]]
  B <- sv.list[[cv2.do[i]]]
  g <- cov(A,B)
  return(g)
}


# helper function to fill similarity matrix across the diagonal
fill.Gijs <- function(i,combos,Gij.df){
  row <- combos[i,]
  if(row$Var1 <= row$Var2){
    return(Gij.df[Gij.df$Var1 == row$Var1 & Gij.df$Var2 == row$Var2,]$Gij)
  } else{
    return(Gij.df[Gij.df$Var1 == row$Var2 & Gij.df$Var2 == row$Var1,]$Gij)
  }
}

# construct similarity matrix
get.similarity.matrix <- function(sv.list){
  n <- length(sv.list)
  combos <- expand.grid(1:n,1:n)
  cv1 <- combos$Var1
  cv2 <- combos$Var2
  do.vec <- as.vector(combos[1]<=combos[2])
  cv1.do <- cv1[do.vec]
  cv2.do <- cv2[do.vec]
  calc.Gijs <- sapply(1:sum(do.vec), get.G.ij, cv1.do, cv2.do, sv.list, sum(do.vec))
  Gij.df <- data.frame(Var1 = cv1.do,
                       Var2 = cv2.do,
                       Gij = calc.Gijs)
  Gijs <- sapply(1:nrow(combos),fill.Gijs,combos,Gij.df)
  G <- matrix(Gijs,nrow=n,ncol=n)
  return(G)
}

# helper function to scale each matrix by the eigenvector element
scale.sv <- function(i,sv.list,u){
  return(sv.list[[i]]*u[i])
}

# ensemble sparsity vectors with eigenscore
get.scores <- function(sv.list){
  G <- get.similarity.matrix(sv.list)
  u <- abs(eigen(G)$vectors[,1])
  u.norm <- u/sum(u)
  scaled.svs <- lapply(1:length(sv.list),scale.sv,sv.list,u.norm)
  scores <- Reduce('+',scaled.svs)
  return(scores)
}








