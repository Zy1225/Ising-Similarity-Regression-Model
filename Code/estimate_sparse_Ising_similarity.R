estimate_sparse_Ising_similarity = function(Y,W,fold_id, type = 'penalized', true_K = NULL){
  n = dim(Y)[2]; p = dim(Y)[1]; K = dim(W)[3]
  
  if(is.na(K)){K=1}
  intercept_matrix = diag(p) %x% rep(1,n)
  
  covariate_matrix = array(data=NA, dim=c(n*p,K))
  if(K==1){
    covariate_matrix[,1] = c(t(W%*% Y))
  }else{
    for(k in 1:K){
      covariate_matrix[,k] = c(t(W[,,k] %*% Y))
    }
  }
  
  
  glm_Y = c(t(Y))
  glm_X = cbind(intercept_matrix,covariate_matrix)
  
  #Fit unpenalized logistic regression to obtain the unpenalized estimate to construct the adaptive weights
  logistic_reg = glm(glm_Y~glm_X + 0,family = binomial,control = glm.control(maxit = 100))
  unpenalized_alpha = logistic_reg$coefficients[(p+1): (p+K)]
  ad.weight <- 1/abs(unpenalized_alpha)
  
  if(type == 'penalized'){
    #Fit penalized logistic regression using adaptive lasso penalty with block cross-validation
    cv.adalasso <- cv.glmnet(glm_X, glm_Y, family = "binomial", alpha=1, standardize= TRUE, intercept = FALSE, penalty.factor= c(rep(0,p),ad.weight),foldid = fold_id, nfolds = 10)
    est_parameters = predict(cv.adalasso,type="coefficients",s=cv.adalasso$lambda.min)[2:(p+K+1)]
    hat_theta_jj = est_parameters[1:p]
    hat_alpha = est_parameters[(p+1): (p+K)]
    best_lambda = cv.adalasso$lambda.min
  }
  
  if(type == 'unpenalized'){
    hat_theta_jj = logistic_reg$coefficients[1:p]
    hat_alpha = logistic_reg$coefficients[(p+1): (p+K)]
    best_lambda = NULL
  }
  
  if(type == 'oracle'){
    glm_X = cbind(intercept_matrix,covariate_matrix[,1:true_K])
    oracle_logistic_reg = glm(glm_Y~glm_X + 0,family = binomial,control = glm.control(maxit = 100))
    hat_theta_jj = oracle_logistic_reg$coefficients[1:p]
    hat_alpha = c(oracle_logistic_reg$coefficients[(p+1):(p+true_K)], rep(0,K-true_K))
    best_lambda = NULL
  }
  
  return(list(hat_theta_jj = hat_theta_jj,hat_alpha = hat_alpha,best_lambda = best_lambda))
}