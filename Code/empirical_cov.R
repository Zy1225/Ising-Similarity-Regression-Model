empirical_cov = function(Y_mat,W,theta_jj,alpha){
  p = dim(W)[1]; K = dim(W)[3]; n = dim(Y_mat)[2]
  if(is.na(K)){K=1}
  empirical_hessian = Matrix(0,nrow = p+K, ncol = p+K)
  for(i in 1:n){
    empirical_hessian = empirical_hessian - hessian_i(Y_mat,W,theta_jj,alpha,i)
  }
  
  empirical_cov_score = array(0, dim = c(p+K,p+K))
  for(i in 1:n){
    empirical_cov_score = empirical_cov_score + matrix(score_i(Y_mat,W,theta_jj,alpha,i), nrow = p+K) %*% matrix(score_i(Y_mat,W,theta_jj,alpha,i), nrow = 1)
  }
  term = solve(empirical_hessian)
  return(term %*% empirical_cov_score %*% term)
}

Ximat = function(W,alpha){
  p = dim(W)[1]; K = dim(W)[3]
  if(is.na(K)){K=1}
  TW = array(data=0, dim = c(p,p))
  if(K > 1){
    for (i in 1:K){
      TW = TW + alpha[i] * W[,,i]
    }
  }else{
    TW = TW + alpha[1] * W
  }
  return(TW)
}

Dl_i_Dtheta_j = function(Y_mat,W,theta_jj,alpha,i,j){
  term = exp(theta_jj[j] + Ximat(W,alpha)[j,] %*% Y_mat[,i] )
  return(c(Y_mat[j,i] - (term)/(1 + term)))
}

Dl_i_Dtheta = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]
  vec = array(NA, dim = p)
  for(j in 1:p){
    vec[j] = Dl_i_Dtheta_j(Y_mat,W,theta_jj,alpha,i,j)
  } 
  return(vec)
}

Dl_i_Dalpha_k = function(Y_mat,W,theta_jj,alpha,i,k){
  K = dim(W)[3]
  if(is.na(K)){K=1}
  if(K>1){
    return(c(t(W[,,k] %*% Y_mat[,i]) %*% Dl_i_Dtheta(Y_mat,W,theta_jj,alpha,i)))
  }else{
    return(c(t(W %*% Y_mat[,i]) %*% Dl_i_Dtheta(Y_mat,W,theta_jj,alpha,i)))
  }
}

Dl_i_Dalpha = function(Y_mat,W,theta_jj,alpha,i){
  K = dim(W)[3]
  if(is.na(K)){K=1}
  vec = array(NA, dim = K)
  for(k in 1:K){
    vec[k] = Dl_i_Dalpha_k(Y_mat,W,theta_jj,alpha,i,k)
  }
  return(vec)
}

score_i = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]; K = dim(W)[3]
  if(is.na(K)){K=1}
  vec = array(NA, dim = p+K)
  vec[1:p] = Dl_i_Dtheta(Y_mat,W,theta_jj,alpha,i)
  vec[(p+1):(p+K)] = Dl_i_Dalpha(Y_mat,W,theta_jj,alpha,i)
  return(vec)
}

#Hessian submatrices
D2l_i_Dtheta_Dtheta = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]
  out_mat = array(0, dim = c(p,p))
  for(j in 1:p){
    term = exp(theta_jj[j] + Ximat(W,alpha)[j,] %*% Y_mat[,i] )
    out_mat[j,j] = ( - term) / (( 1 + term)^2)
  }
  return(out_mat)
}

D2l_i_Dalpha_Dtheta = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]; K = dim(W)[3]
  if(is.na(K)){K=1}
  out_mat = array(0, dim = c(p,K))
  for(j in 1:p){
    term = exp(theta_jj[j] + Ximat(W,alpha)[j,] %*% Y_mat[,i] )
    if(K > 1){
      for(k in 1:K){
        out_mat[j,k] = ( - (c(W[j,,k] %*% Y_mat[,i])) * term) / (( 1 + term)^2)
      }
    }else{
      out_mat[j,1] = ( - (c(W[j,] %*% Y_mat[,i])) * term) / (( 1 + term)^2)
    }
  }
  return(out_mat)
}

D2l_i_Dtheta_Dalpha = function(Y_mat,W,theta_jj,alpha,i){
  return(t(D2l_i_Dalpha_Dtheta(Y_mat,W,theta_jj,alpha,i))) 
}

D2l_i_Dalpha_Dalpha = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]; K = dim(W)[3]
  if(is.na(K)){K=1}
  out_mat = array(0, dim = c(K,K))
  if(K > 1){
    for(k1 in 1:K){
      for(k2 in 1:K){
        out_mat[k1,k2] = 0
        for(j in 1:p){
          term = exp(theta_jj[j] + Ximat(W,alpha)[j,] %*% Y_mat[,i] )
          out_mat[k1,k2] = out_mat[k1,k2]  + (- (c(W[j,,k1] %*% Y_mat[,i])) * (c(W[j,,k2] %*% Y_mat[,i]))* term ) / ((1 + term)^2)
        }
      }
    }
  }else{
    out_mat[1,1] = 0
    for(j in 1:p){
      term = exp(theta_jj[j] + Ximat(W,alpha)[j,] %*% Y_mat[,i] )
      out_mat[1,1] = out_mat[1,1]  + (- (c(W[j,] %*% Y_mat[,i])) * (c(W[j,] %*% Y_mat[,i]))* term ) / ((1 + term)^2)
    }
  }
  return(out_mat)
}

hessian_i = function(Y_mat,W,theta_jj,alpha,i){
  p = dim(W)[1]; K = dim(W)[3]
  if(is.na(K)){K = 1}
  out_mat = array(NA,  dim = c(p+K,p+K))
  out_mat[1:p,1:p] = D2l_i_Dtheta_Dtheta(Y_mat,W,theta_jj,alpha,i)
  out_mat[1:p, ((p+1):(p+K))] = D2l_i_Dalpha_Dtheta(Y_mat,W,theta_jj,alpha,i)
  out_mat[((p+1):(p+K)), 1:p] = D2l_i_Dtheta_Dalpha(Y_mat,W,theta_jj,alpha,i)
  out_mat[((p+1):(p+K)),((p+1):(p+K))] = D2l_i_Dalpha_Dalpha(Y_mat,W,theta_jj,alpha,i)
  return(Matrix(out_mat,sparse = T))
}

