generate_sparse_ising = function(n,p,K,true_K,B,seed,filename){
  #This function generate W_k matrices, alpha parameters, theta_jj parameters.
  #Then, it generates B replications of (p x n) response matrices from the Ising similarity regression model
  set.seed(seed)
  Y = array(data = NA,dim = c(p,n,B))
  W = array(data = NA,dim = c(p,p,K))
  
  
  alpha = rep(0,K)
  
  #randomly generate true values for alpha uniformly from [-0.4,-0.3] U [0.3,0.4]
  alpha[1:true_K] = runif(true_K,min = -0.4,max = -0.3)
  for(k in 1:true_K){
    sign = sample(c(-1,1),size = 1)
    alpha[k] = alpha[k] * sign
  }
  
  
  #Generate K similarity matrices where w_{j1,j2,k} = exp(-d_{j1,j2,k}^2) with d_{j1,j2,k} i.i.d U(p^{-1/2}, p^{1/2})
  for(k in 1:K){
    d = exp( - (runif(p*(p-1)/2, min = p^(-0.5) ,max=p^(0.5)) ^2))
    W[,,k][upper.tri(W[,,k])] = d
    W[,,k][lower.tri(W[,,k])] =  t(W[,,k])[lower.tri(W[,,k])]
    diag(W[,,k]) = 0
  }
  
  
  
  
  #randomly generate theta_{jj} uniformly from [-1,-0.5] U [0.5,1] for j=1,...,p u
  theta_jj = runif(p,min = -1, max = 0)
  theta_jj[theta_jj > -0.5] = theta_jj[theta_jj > -0.5] + 1
  
  #Theta matrix as a linear combination of W_k
  Theta = array(data = 0, dim= c(p,p))
  for(k in 1:K){
    Theta = Theta + alpha[k] * W[,,k]
  }
  
  #Inverse temperature beta always = 1
  Beta = 1
  for(b in 1:B){
    Y[,,b] = t(IsingSampler(n, Theta, theta_jj, Beta, 1000/p,
                            responses = c(0L,1L), method = "MH"))
  }
  
  
  save(Y,W,alpha,theta_jj,Theta,Beta,n,p,K,true_K,B,seed,file = filename)
}
