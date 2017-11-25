
## SCORE ##########################################################################

SCORE<-function(A,K,mode="collapse"){
  # Input:
  #   the igraph A, and number of communities K
  # Output:
  #   the cluster vector indicating the community of each node

  # pre-process data
  dat = igraph_to_matrix(A,K,mode)
  A = dat$A
  gc = dat$gc
  direct = dat$direct
  eta = dat$eta

  # Set parameter
  n = gorder(A)
  T_n = log(n) #threshold

  # Calculate R_star
  R_star<-matrix(0,nrow=n,ncol=K-1)
  for(i in 1:n){
    for(k in 1:(K-1)){
      r = eta[i,k+1]/eta[i,1]
      if (  r < -T_n )
      { R_star[i,k] = -T_n }
      if (  r > T_n )
      { R_star[i,k] = T_n }
      else {
        R_star[i,k]<- r
      }
    }
  }

  # K-means step
  k<-kmeans(R_star,centers=K)
  community = matrix(k$cluster,gorder(A),1)
  rownames(community) = V(A)$label
  return(list(community = community, direct = direct, gc = gc, R = R_star))

}

## mixedSCORE #####################################################################

mixedSCORE <- function(A, K, L_theo = FALSE, mode = "collapse"){
  # Input:
  #   the adjacency matrix A, and number of communities K
  #   L_theo indicates whether estimate L by the theoretical approach
  # Output:
  #   the estimated membership matrix Pi

  # pre-process data
  dat = igraph_to_matrix(A,K,mode)
  A = dat$A
  gc = dat$gc
  direct = dat$direct
  eta = dat$eta
  lambda = dat$lambda
  
  # Set parameter
  n = gorder(A)
  T_n = log(n) #threshold

  # Calculate R_star
  R_star<-matrix(0,nrow=n,ncol=K-1)
  for(i in 1:n){
    for(k in 1:(K-1)){
      r = eta[i,k+1]/eta[i,1]
      if (  r < -T_n )
      { R_star[i,k] = -T_n }
      if (  r > T_n )
      { R_star[i,k] = T_n }
      else {
        R_star[i,k]<- r
      }
    }
  }

  #  Estimate L by the theoretical approach
  if (L_theo == TRUE){
    l_star = K+1
    while(kmeans(R_star,centers=(l_star+1))$tot.within/kmeans(R_star,centers=l_star)$tot.within 
          >= 1/log(log(n))){
      l_star=l_star+1
    }
    V = vertice_hunting(R_star,l_star)$matrix
  }
  #  Estimate L by the practical approach
  if (L_theo == FALSE){
    s <- c()
    v_old= vertice_hunting(R_star,K)$matrix
    for (l in (K+1):(3*K)){
      VH = vertice_hunting(R_star,l)
      # obtain d_l and (v1_(l), ..., vk_(l))
      d_l = VH$d
      v_new = VH$matrix
      # calculate delta_l
      delta_l = delta(v_old, v_new, K)
      # update s
      s <- c(s, delta_l/(1+d_l))
      v_old = v_new
    }
    # this is for picking the larger index if there is a tie
    l_star = K+1+length(s)-which.min(s[length(s):1])
    V <- vertice_hunting(R_star,l_star)$matrix
  }

  # Calculate b1
  b_1 <- matrix(0,1,K)
  for (k in 1:K){
    if (K == 2){
      b_1[1,k] <- 1/sqrt(lambda[1] + lambda[2]* matrix(V[k,],1,K-1) %*% matrix(V[k,],K-1,1))
    }
    if (K > 2){
      b_1[1,k] <- 1/sqrt(lambda[1] + matrix(V[k,],1,K-1) %*% diag(lambda[2:K]) %*% matrix(V[k,],K-1,1))
    }
  }

  # Recover B and P
  B = matrix(0,K,K)
  B[,1]=b_1
  for (l in 2:K){
    for(k in 1:K){
      B[k,l] = V[k,l-1]*b_1[k]
    }
  }
  P = B %*% diag(lambda) %*% t(B)
  
  # Calculate Pi
  Pi_hat = matrix(0,n,K)
  for (i in 1:n){
    w_i = solve(rbind(t(V),1),c(R_star[i,],1))
    Pi_star = (w_i/b_1)*(w_i/b_1>0)
    Pi_hat[i,] = Pi_star/sum(Pi_star)
  }
  rownames(Pi_hat) = V(A)$label
  colnames(Pi_hat) = 1:K
  
  return(list(member = Pi_hat, P = P, direct = direct, gc = gc,
              B = B, R = R_star, vertice = V, l = l_star))
}

## topicSCORE #####################################################################

topicSCORE <- function(D,K){
  # Input:
  #   the text corpus matrix D, and number of topics K
  # Output:
  #   the estimated topic matrix A

  # Set parameter
  n=ncol(D)
  p=nrow(D)
  L=10*K
  T_n = 2*log(n)
  M_diag = as.vector((D %*% matrix(1,n,1))/n)

  # Calculate K singular values and left singular vector
  s = svd(diag(1/sqrt(M_diag)) %*% D)
  sigma = s$d[1:K]
  Xi = s$u[,1:K]

  # Calculate R_star
  R_star<-matrix(0,nrow=p,ncol=K-1)
  for(i in 1:n){
    for(k in 1:(K-1)){
      r = Xi[i,k+1]/Xi[i,1]
      if (  r < -T_n )
      { R_star[i,k] = -T_n }
      if (  r > T_n )
      { R_star[i,k] = T_n }
      else {
        R_star[i,k]<- r
      }
    }
  }

  # Calculate the vertex matrix V, and Pi
  V = vertice_hunting(R_star, L)
  Ver = t(V$matrix)
  V_hat = rbind(1,Ver)
  Pi = matrix(0,p,K)
  for (j in 1:p){
    Pi_star = Solve(V_hat,c(1,R_star[j,]))
    Pi_star = Pi_star*(Pi_star>0)
    Pi[j,] = Pi_star/sum(abs(Pi_star))
  }

  # Calculate A_star and A_hat
  A_star = diag(sqrt(M_diag)) %*% diag(Xi[,1]) %*% Pi
  A_hat = matrix(0,p,K)
  for (i in 1:K){
    A_hat[,i] = A_star[,i]/sum(abs(A_star[,i]))
  }
  return(list(topic = A_hat, weight = Pi, R = R_star, vertice = V))

}

## Ancillary functions ###########################################################

igraph_to_matrix <- function(A,K,mode){
  
  require(igraph)
  # convert to undirected network if necessary
  if(is.directed(A) == FALSE){
    direct = NULL
  }
  if(is.directed(A) == TRUE){
    A = as.undirected(A,mode)
    direct = mode
  }
  adj_matrix <- as_adjacency_matrix(A,sparse=FALSE)
  
  # extract giant component if necessary
  eigen = eigen(adj_matrix,symmetric=TRUE)
  lambda = eigen$values[1:K]
  eta = eigen$vectors[,1:K]
  
  n_zero = sum(eta[,1]==0) #number of zero component
  if(n_zero > 0){ # extract giant component
    c = clusters(A)
    A_giant <- induced.subgraph(A,which(c$membership == which.max(c$csize)))
    gc = V(A)$label[setdiff(1:gorder(A),which(c$membership == which.max(c$csize)))]
    adj_matrix <- as_adjacency_matrix(A_giant,sparse=FALSE)
    A = A_giant
    # calculate eigendecomposition
    eigen = eigen(adj_matrix,symmetric=TRUE)
    lambda = eigen$values[1:K]
    eta = eigen$vectors[,1:K]
  }
  if(n_zero == 0){
    gc = NULL
  }
  
  return(list(A = A, gc = gc, direct = direct,eta = eta, lambda =lambda))
}


vertice_hunting <- function(R_star,L){
  # Inputs:
  #   R_star is the matrix constructed from the eigenvectors of A by SCORE
  #   L is the tuning parameter indicating the number of clusters
  # Output:
  #   The estimated vertices (v1,...,vK)^T and d_L

  require(quadprog)
  require(utils)
  require(limSolve)
  require(lpSolve)

  # apply k-means method to R_star assuming L clusters, write the centers m
  cluster<-kmeans(R_star,centers=L)
  m<-cluster$centers
  K <- ncol(R_star)+1
  
  # if K = 2, the convex hull can always enclose all the points
  if(K==2){ 
    d_L = 0
    V = matrix(c(max(m),min(m)),2,1)
  }

  if(K>2){
    # Get K vertices out of the L cluster centers
    comb <-combn(L,K)
    combin = ncol(comb)
    
    d_max <- matrix(0,1,combin) # list of the max dist
    d_L<- 0                     # minimum of max dist
    
    for(i in 1:combin){
      c<-c()
      for (j in 1:K) {c<-rbind(c,m[comb[j,i],])}
      V<-matrix(c,K,K-1) # get vertex matrix

      for(k in 1:L){
        f.obj <- 1
        f.con <- rbind(t(V),matrix(1,1,K),diag(K))
        f.dir <- c(matrix("=",1,K),matrix(">=",1,K))
        f.rhs <- c(m[k,],1,matrix(0,1,K))
        result<-lp ("min", f.obj, f.con, f.dir, f.rhs)

        if (result$status==2){ #outside the hull
          A <- t(V)
          B <- matrix(m[k,],K-1,1)
          E <- matrix(1,1,K)
          G <- diag(K)
          H <- matrix(0,K,1)
          out<- lsei (A = A, B = B, E = E, F = 1, G = G, H = H,
                      Wx = NULL, Wa = NULL, type = 1, tol = sqrt(.Machine$double.eps),
                      tolrank = NULL, fulloutput = FALSE, verbose = TRUE)

          if (d_max[i] < sqrt(out$solutionNorm) ){
            d_max[i] <- sqrt(out$solutionNorm)
          }
        }
      }
      # find min of the max dist
      d_L<- min(d_max)
      i<-which.min(d_max)
      c<-c()
      for (j in 1:K) {c<-rbind(c,m[comb[j,i],])}
      V<-matrix(c,K,K-1)
    }
  }
  return(list(matrix = V, d = d_L))
}

delta <- function(v_old, v_new, K){
  # Input:
  #   v_old and v_new are the estimated vertices matrix
  #   K is the number of communities
  # Output:
  #   delta_L

  require(gtools)
  perm = permutations(K,K)
  
  d <- c()
  for (i in 1:nrow(perm)){
    A<-matrix(0,K,K)
    for(j in 1:K){ #generate permutation matrix
      A[j,perm[i,j]]<-1
    }
    d <- c(d, max(rowSums((A %*% v_new - v_old)^2)))
  }
  
  return(min(d))
}
