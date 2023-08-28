rm(list = ls())
library(IsingSampler)
library(tidyr)
library(glmnet)


here::i_am("Simulation/simulation_script.R")
library(here)
source(here("Code","generate_sparse_ising.R"))
source(here("Code","estimate_sparse_Ising_similarity.R"))

#Generating Simulation Data
ns = c(50,100,200,400)
ps = c(10,25,50,100,200)
K = 20
true_K = 5
seed = 12
B = 1000

for(n in ns){
  for(p in ps){
    filename = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    generate_sparse_ising(n,p,K,true_K,B,seed,filename)
  }
}



#Penalized Pseudo-likelihood estimator
for(n in ns){
  for(p in ps){
    filename = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    file_result = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_result.csv"))
    load(filename)
    
    for(b in 1:B){
      tryCatch({
        est_fit = estimate_sparse_Ising_similarity(Y[,,b],W, fold_id = rep(rep(1:10, each = (n /10)), times = p))
        output = t(as.matrix(c(b, est_fit$hat_alpha, est_fit$hat_theta_jj)))
        if(b == 1){
          write.table(output,  file=file_result,  sep=",", row.names=F, col.names=F)
        }else{
          write.table(output,  file=file_result, append = T, sep=",", row.names=F, col.names=F)
        }
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}




#Oracle Estimator
for(n in ns){
  for(p in ps){
    filename = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    file_result = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_oracle_result.csv"))
    load(filename)
    
    for(b in 1:B){
      tryCatch({
        est_fit = estimate_sparse_Ising_similarity(Y[,,b],W, fold_id = rep(rep(1:10, each = (n /10)), times = p), type = "oracle", true_K = 5)
        output = t(as.matrix(c(b, est_fit$hat_alpha, est_fit$hat_theta_jj)))
        if(b == 1){
          write.table(output,  file=file_result,  sep=",", row.names=F, col.names=F)
        }else{
          write.table(output,  file=file_result, append = T, sep=",", row.names=F, col.names=F)
        }
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}



#Unpenalized pseudo-likelihood estimator
for(n in ns){
  for(p in ps){
    filename = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    file_result = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_unpenalized_result.csv"))
    load(filename)
    
    for(b in 1:B){
      tryCatch({
        est_fit = estimate_sparse_Ising_similarity(Y[,,b],W, fold_id = rep(rep(1:10, each = (n /10)), times = p), type = "unpenalized", true_K = NULL)
        output = t(as.matrix(c(b, est_fit$hat_alpha, est_fit$hat_theta_jj)))
        if(b == 1){
          write.table(output,  file=file_result,  sep=",", row.names=F, col.names=F)
        }else{
          write.table(output,  file=file_result, append = T, sep=",", row.names=F, col.names=F)
        }
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}



#Compute MSE_theta and MSE_alpha for three types of estimator
alpha_mse_vec = NULL
theta_jj_mse_vec = NULL
theta_jj_non_scaled_mse_vec = NULL

alpha_mse_vec_oracle = NULL
theta_jj_mse_vec_oracle = NULL
theta_jj_non_scaled_mse_vec_oracle = NULL

alpha_mse_vec_unpenalized = NULL
theta_jj_mse_vec_unpenalized = NULL
theta_jj_non_scaled_mse_vec_unpenalized = NULL
for(n in ns){
  for(p in ps){
    file_data = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    file_result = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_result.csv"))
    file_result_oracle = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_oracle_result.csv"))
    file_result_unpenalized = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_unpenalized_result.csv"))
    
    
    load(file_data)
    alpha = alpha
    theta_jj = theta_jj
    
    
    
    result = read.csv(file_result,header = F, sep = ",")
    result_alpha = result[,2:(K+1)]
    result_theta_jj = result[,(K+2):(K+1+p)]
    
    result_oracle = read.csv(file_result_oracle,header = F, sep = ",")
    result_alpha_oracle = result_oracle[,2:(K+1)]
    result_theta_jj_oracle = result_oracle[,(K+2):(K+1+p)]
    
    result_unpenalized = read.csv(file_result_unpenalized,header = F, sep = ",")
    result_alpha_unpenalized = result_unpenalized[,2:(K+1)]
    result_theta_jj_unpenalized = result_unpenalized[,(K+2):(K+1+p)]
    
    
    
    alpha_mse = (norm(as.matrix(result_alpha - matrix(rep(alpha, B), nrow = B, byrow = TRUE)), type = "F")^2)/(K*B)
    theta_jj_mse = (norm(as.matrix(result_theta_jj - matrix(rep(theta_jj,B),nrow = B, byrow = TRUE)), type = "F")^2)/(p*B)
    
    
    alpha_mse_oracle = (norm(as.matrix(result_alpha_oracle - matrix(rep(alpha, B), nrow = B, byrow = TRUE)), type = "F")^2)/(K*B)
    theta_jj_mse_oracle = (norm(as.matrix(result_theta_jj_oracle - matrix(rep(theta_jj,B),nrow = B, byrow = TRUE)), type = "F")^2)/(p*B)
    
    alpha_mse_unpenalized = (norm(as.matrix(result_alpha_unpenalized - matrix(rep(alpha, B), nrow = B, byrow = TRUE)), type = "F")^2)/(K*B)
    theta_jj_mse_unpenalized = (norm(as.matrix(result_theta_jj_unpenalized - matrix(rep(theta_jj,B),nrow = B, byrow = TRUE)), type = "F")^2)/(p*B)
    
    
    alpha_mse_vec = c(alpha_mse_vec,alpha_mse)
    theta_jj_mse_vec = c(theta_jj_mse_vec,theta_jj_mse)
    
    alpha_mse_vec_oracle = c(alpha_mse_vec_oracle,alpha_mse_oracle)
    theta_jj_mse_vec_oracle = c(theta_jj_mse_vec_oracle,theta_jj_mse_oracle)
    
    alpha_mse_vec_unpenalized = c(alpha_mse_vec_unpenalized,alpha_mse_unpenalized)
    theta_jj_mse_vec_unpenalized = c(theta_jj_mse_vec_unpenalized,theta_jj_mse_unpenalized)
  }
}
mse_results = expand.grid(c(10,25,50,100,200),c(50,100,200,400))
colnames(mse_results) = c("p","n")
mse_results$p = as.factor(mse_results$p)
mse_results$type = "Penalized"

mse_results$alpha_mse = alpha_mse_vec
mse_results$theta_jj_mse = theta_jj_mse_vec

mse_results_oracle = expand.grid(c(10,25,50,100,200),c(50,100,200,400))
colnames(mse_results_oracle) = c("p","n")
mse_results_oracle$p = as.factor(mse_results_oracle$p)
mse_results_oracle$type = "Oracle"

mse_results_oracle$alpha_mse = alpha_mse_vec_oracle
mse_results_oracle$theta_jj_mse = theta_jj_mse_vec_oracle

mse_results_unpenalized = expand.grid(c(10,25,50,100,200),c(50,100,200,400))
colnames(mse_results_unpenalized) = c("p","n")
mse_results_unpenalized$p = as.factor(mse_results_unpenalized$p)
mse_results_unpenalized$type = "Unpenalized"

mse_results_unpenalized$alpha_mse = alpha_mse_vec_unpenalized
mse_results_unpenalized$theta_jj_mse = theta_jj_mse_vec_unpenalized

final_mse_results = rbind(mse_results,mse_results_oracle,mse_results_unpenalized)

alpha_mse_wide = spread(final_mse_results[,c(1,2,3,4)], key = type, value = alpha_mse)
colnames(alpha_mse_wide) = c("p","n","Oracle_alpha_mse","Penalized_alpha_mse","Unpenalized_alpha_mse")

theta_jj_mse_wide = spread(final_mse_results[,c(1,2,3,5)], key = type, value = theta_jj_mse)
colnames(theta_jj_mse_wide) = c("p","n","Oracle_theta_jj_mse","Penalized_theta_jj_mse","Unpenalized_theta_jj_mse")

final_scaled_MSE_table = cbind(alpha_mse_wide , theta_jj_mse_wide[,3:5])
final_scaled_MSE_table = final_scaled_MSE_table[order(final_scaled_MSE_table[,1],final_scaled_MSE_table[,2]),]
final_scaled_MSE_table[,2:8] = round(final_scaled_MSE_table[,2:8],3)
rownames(final_scaled_MSE_table) = NULL

#Compute TPR and FPR of penalized pseudo-likelihood estimator
TPR = NULL
FPR = NULL

for(n in ns){
  for(p in ps){
    file_data = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,".Rdata"))
    file_result = here("Simulation",paste0("Sparse_Ising_n",n,"_","p",p,"_result.csv"))
    
    load(file_data)
    alpha = alpha
    theta_jj = theta_jj
    
    
    result = read.csv(file_result,header = F, sep = ",")
    result_alpha = result[,2:(K+1)]
    
    TPR = c(TPR,mean(result_alpha[,1:true_K] !=0))
    FPR = c(FPR,mean(result_alpha[,(true_K + 1):K]!=0))
    
  }
}
selection_results = expand.grid(ps,ns)
colnames(selection_results) = c("p","n")
selection_results$p = as.factor(selection_results$p)
selection_results$TPR = TPR
selection_results$FPR = FPR
selection_results = selection_results[order(selection_results[,1],selection_results[,2]),] 
selection_results[,3:4] = round(selection_results[,3:4],3)
rownames(selection_results) = NULL

TPR_mat = matrix(selection_results[,3], nrow = 5, ncol = 4, byrow = TRUE)
FPR_mat =  matrix(selection_results[,4], nrow = 5, ncol = 4, byrow =TRUE)
colnames(TPR_mat) = colnames(FPR_mat) = c("n=50","n=100","n=200","n=400")
rownames(TPR_mat) = rownames(FPR_mat) = c("p=10","p=25","p=50","p=100","p=200")



#Final Results
final_scaled_MSE_table
TPR_mat
FPR_mat

