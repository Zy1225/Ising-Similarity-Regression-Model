rm(list=ls())
library(glmnet)
library(network)
library(ggplot2)
library(sna)
library(ggpubr)
library(GGally)


here::i_am("Application/application.R")
library(here)
source(here("Code","estimate_sparse_Ising_similarity.R"))
source(here("Code","empirical_cov.R"))

#Read data
vote_data = read.csv(here("Application",'Roll_Call_Voting_Congress_117.csv'),header = TRUE,check.names = FALSE)
senator_data = read.csv(here("Application",'Congress_117_Senator_Data.csv'), header = TRUE)
senator_wiki_data = read.csv(here("Application",'Congress_117_Senator_wiki.csv'),header = FALSE)

#Extract senator's full name/party/state from US Senate roll call voting data
full_names = colnames(vote_data)[16:115]
full_names = sapply(strsplit(full_names," \\("),"[",1)

party = colnames(vote_data)[16:115]
party= sapply(strsplit(party," \\("),"[",2)
party = sapply(party,function(x)substring(x, 1, 1))
names(party) = NULL

state = colnames(vote_data)[16:115]
state = sapply(strsplit(state,"\\-"),"[",2)
state = sapply(state,function(x)substring(x, 1, 2))
names(state) = NULL

#Remove procedural votes with Yea/Nay proportion outside of [0.3,0.7]
proportion = vote_data$Yea / (vote_data$Yea + vote_data$Nay)
argue = (proportion >= 0.3 & proportion <= 0.7)
argue_vote_data = vote_data[argue,]
Y = argue_vote_data[,c(16:115)]

check = c()
for(bill in 1:nrow(Y)){
  check = c(check,unique(Y[bill,]))
}
unique(check)


#Convert Yea/Nay or Guilty/Not Guilty into binary 0/1
for(bill in 1:nrow(Y)){
  for(s in 1:ncol(Y)){
    if(Y[bill,s] == 'Yea'){
      Y[bill,s] = 1
    }
    if(Y[bill,s] == 'Nay'){
      Y[bill,s] = 0
    } 
    if(Y[bill,s] == 'Guilty'){
      Y[bill,s] = 1
    }
    if(Y[bill,s] == 'Not Guilty'){
      Y[bill,s] = 0
    }
  }
}

#Impute missing votes
temp_Y = Y
temp_Y[temp_Y == 'Not Voting'] = NA
temp_Y[temp_Y == 'Present'] = NA

set.seed(1)
for(bill in 1:nrow(Y)){
  for(s in 1:ncol(Y)){
    if(!Y[bill,s] %in% c(0,1)){
      if(party[s]!='I'){
        #Imputing for Democrat or Republican
        Y[bill,s] = ifelse(mean(as.numeric(temp_Y[bill, party == party[s]]), na.rm = T) > 0.5, 1,0)
      }else{
        #Imputing for Independent
        Y[bill,s] = sample(c(0,1),1,prob = c(0.5,0.5))
      }
    }
  }
}

#Constructing Similarity Matrices for State/Party/Class
W_state = matrix(NA,nrow = 100, ncol = 100)
W_party = matrix(NA,nrow = 100, ncol = 100)
W_class = matrix(NA,nrow = 100, ncol = 100)

for(j1 in 1:100){
  for(j2 in 1:j1){
    if(j1 == j2){
      W_state[j1,j2] = 0
      W_party[j1,j2] = 0
      W_class[j1,j2] = 0
    }else{
      if(senator_data$State[j1] == senator_data$State[j2]){
        W_state[j1,j2] = 1
      }else{
        W_state[j1,j2] = 0
      }
      
      if(senator_data$Party[j1] == senator_data$Party[j2]){
        W_party[j1,j2] = 1
      }else{
        W_party[j1,j2] = 0
      }
      
      if(senator_data$Class[j1] == senator_data$Class[j2]){
        W_class[j1,j2] = 1
      }else{
        W_class[j1,j2] = 0
      }
    }
  }
}

W_state[upper.tri(W_state)] =  t(W_state)[upper.tri(W_state)]
W_party[upper.tri(W_party)] =  t(W_party)[upper.tri(W_party)]
W_class[upper.tri(W_class)] =  t(W_class)[upper.tri(W_class)]

# Construct similarity matrices using data from wikipedia, i.e. for Age/Gender/Lawyer/Executive/Business/Farmer/Army/Teacher/Professor

#Manually reformat the names in the wikipedia to match up with US Senate website
senator_wiki_data$V2 = sapply(senator_wiki_data$V2,function(x)substring(x, 1, nchar(x)-1))

senator_wiki_data$V2[senator_wiki_data$V2 == "Chris Murphy"] = "Christopher Murphy"
senator_wiki_data$V2[senator_wiki_data$V2 == "Tom Carper"] = "Thomas Carper"
senator_wiki_data$V2[senator_wiki_data$V2 == "Chris Coons"] = "Christopher Coons"
senator_wiki_data$V2[senator_wiki_data$V2 == "Jim Risch"] = "James  Risch"
senator_wiki_data$V2[senator_wiki_data$V2 == "Dick Durbin"] = "Richard Durbin"
senator_wiki_data$V2[senator_wiki_data$V2 == "Chuck Grassley"] = "Charles Grassley"
senator_wiki_data$V2[senator_wiki_data$V2 == "Ed Markey"] = "Edward Markey"
senator_wiki_data$V2[senator_wiki_data$V2 == "Jacky Rosen"] = "Jacklyn Rosen"
senator_wiki_data$V2[senator_wiki_data$V2 == "Bob Menendez"] = "Robert Menendez"
senator_wiki_data$V2[senator_wiki_data$V2 == "Ben Ray LujÃ¡n"] = "Ben Lujan"
senator_wiki_data$V2[senator_wiki_data$V2 == "Chuck Schumer"] = "Charles Schumer"
senator_wiki_data$V2[senator_wiki_data$V2 == "Thom Tillis"] = "Thomas Tillis"
senator_wiki_data$V2[senator_wiki_data$V2 == "Sherrod Brown"] = "Sherrod  Brown"
senator_wiki_data$V2[senator_wiki_data$V2 == "Bob Casey Jr."] = "Bob Casey"
senator_wiki_data$V2[senator_wiki_data$V2 == "Pat Toomey"] = "Patrick Toomey"
senator_wiki_data$V2[senator_wiki_data$V2 == "Jack Reed"] = "John Reed"
senator_wiki_data$V2[senator_wiki_data$V2 == "Bernie Sanders"] = "Bernard Sanders"
senator_wiki_data$V2[senator_wiki_data$V2 == "Tim Kaine"] = "Timothy Kaine"
senator_wiki_data$V2[senator_wiki_data$V2 == "Joe Manchin"] = "Joseph Manchin"
senator_wiki_data$V2[senator_wiki_data$V2 == "Shelley Moore Capito"] = "Shelley Capito"

reordered_senator_wiki_data = senator_wiki_data[match(full_names,senator_wiki_data$V2),]
sum(reordered_senator_wiki_data$V2 == full_names)

#Age
age = reordered_senator_wiki_data$V4
age = sapply(strsplit(age,"age "),"[",2)
age = sapply(strsplit(age,"\\)"),"[",1)
age = as.numeric(age)
age = scale(age)[,1]

W_age = matrix(NA, nrow = 100, ncol =100)
for(j1 in 1:100){
  for(j2 in 1:j1){
    W_age[j1,j2] = exp(-(age[j1] - age[j2])^2)
  }
}
W_age[upper.tri(W_age)] =  t(W_age)[upper.tri(W_age)]
diag(W_age) = 0

#Gender
gender = reordered_senator_wiki_data$V10
W_gender = matrix(NA,nrow = 100, ncol = 100)
for(j1 in 1:100){
  for(j2 in 1:j1){
    if(j1 == j2){
      W_gender[j1,j2] = 0
    }else{
      if(gender[j1] == gender[j2]){
        W_gender[j1,j2] = 1
      }else{
        W_gender[j1,j2] = 0
      }
    }
  }
}
W_gender[upper.tri(W_gender)] =  t(W_gender)[upper.tri(W_gender)]
diag(W_gender) = 0

#Create a similarity matrix for each of the most common job - Lawyer/Executive/Business/Farmer/Army/Teacher/Professor
lawyer_jobs = c('Lawyer','United States Attorney','U.S. Attorney','Texas Solicitor General','Deputy Assistant United States Attorney, U.S. Department of Justice Office of Legislative Affairs','Assistant United States Attorney')
executive_jobs = c("Corporate executive","Nonprofit organization executive")
army_jobs = c('U.S. Navy officer','United States Navy Reserve officer','United States Army officer','U.S. Marine Corps officer',
              'U.S. Air Force officer/Judge Advocate','Specialist Fourth Class, U.S. Army','Petty officer third class, U.S. Navy',
              'Member, United States Army Reserve','Marine Corps Reserve Sergeant','Marine Corps Officer','Captain, U.S. Army Reserve',
              'Army Reserve officer','Army officer','Army National Guard Officer','Army National Guard officer','Air Force Reserve officer')
teacher_jobs = c('Teacher','Music teacher')
professor_jobs = c('Professor','College professor','College professor and lecturer','University president')


W_lawyer = matrix(NA,nrow = 100,ncol = 100)
W_executive = matrix(NA,nrow = 100,ncol = 100)
W_business = matrix(NA,nrow = 100,ncol = 100)
W_farmer = matrix(NA,nrow = 100,ncol = 100)
W_army = matrix(NA,nrow = 100,ncol = 100)
W_teacher = matrix(NA,nrow = 100,ncol = 100)
W_professor = matrix(NA,nrow = 100,ncol = 100)

for(j1 in 1:100){
  for(j2 in 1:100){
    job1 = unlist(strsplit(reordered_senator_wiki_data$V5[j1],";"))
    job1[length(job1)] = substring(job1[length(job1)],1,nchar(job1[length(job1)])-1)
    
    job2 = unlist(strsplit(reordered_senator_wiki_data$V5[j2],";"))
    job2[length(job2)] = substring(job2[length(job2)],1,nchar(job2[length(job2)])-1)
    
    if( (sum(job1 %in% lawyer_jobs) > 0) &  (sum(job2 %in% lawyer_jobs) > 0) ){
      W_lawyer[j1,j2] = 1
    }else{
      W_lawyer[j1,j2] = 0
    }
    
    if( (sum(job1 %in% executive_jobs) > 0) &  (sum(job2 %in% executive_jobs) > 0) ){
      W_executive[j1,j2] = 1
    }else{
      W_executive[j1,j2] = 0
    }
    
    if( (sum(job1 %in% "Businessman") > 0) &  (sum(job2 %in% "Businessman") > 0) ){
      W_business[j1,j2] = 1
    }else{
      W_business[j1,j2] = 0
    }
    
    if( (sum(job1 %in% "Farmer") > 0) &  (sum(job2 %in% "Farmer") > 0) ){
      W_farmer[j1,j2] = 1
    }else{
      W_farmer[j1,j2] = 0
    }
    
    
    if( (sum(job1 %in% army_jobs) > 0) &  (sum(job2 %in% army_jobs) > 0) ){
      W_army[j1,j2] = 1
    }else{
      W_army[j1,j2] = 0
    }
    
    if( (sum(job1 %in% teacher_jobs) > 0) &  (sum(job2 %in% teacher_jobs) > 0) ){
      W_teacher[j1,j2] = 1
    }else{
      W_teacher[j1,j2] = 0
    }
    
    if( (sum(job1 %in% professor_jobs) > 0) &  (sum(job2 %in% professor_jobs) > 0) ){
      W_professor[j1,j2] = 1
    }else{
      W_professor[j1,j2] = 0
    }
    
    
  }
}
diag(W_lawyer) = 0
diag(W_executive) = 0
diag(W_business) = 0
diag(W_farmer) = 0
diag(W_army) = 0
diag(W_teacher) = 0
diag(W_professor) = 0

# Construct similarity matrices using data from Twitter - Tweets/Followers
twitter_df = read.table(file = here("Application",'twitter_df.csv') , header = T, sep = ",")
sum(twitter_df$name == reordered_senator_wiki_data$V2)

#Number of tweets
tweets = twitter_df$tweets
tweets = scale(tweets)[,1]
W_tweets = matrix(NA, nrow = 100, ncol =100)
for(j1 in 1:100){
  for(j2 in 1:j1){
    W_tweets[j1,j2] = exp(-(tweets[j1] - tweets[j2])^2)
  }
}
W_tweets[upper.tri(W_tweets)] =  t(W_tweets)[upper.tri(W_tweets)]
diag(W_tweets) = 0


#Number of followers
followers = twitter_df$followers
followers = scale(followers)[,1]
W_followers = matrix(NA, nrow = 100, ncol =100)
for(j1 in 1:100){
  for(j2 in 1:j1){
    W_followers[j1,j2] = exp(-(followers[j1] - followers[j2])^2)
  }
}
W_followers[upper.tri(W_followers)] =  t(W_followers)[upper.tri(W_followers)]
diag(W_followers) = 0


#A symmetric adjacency matrix based on follower-followee relationship on twitter among senators, where W_{ij} = W_{ji} = 1 if i follows j OR j follows i
W_follow_adj = read.table(here("Application",'W_follow_adj.csv'),header  = T, sep = ",")
W_follow_adj = as.matrix(W_follow_adj)
colnames(W_follow_adj) = NULL

#### Fit Ising Similarity Regression Model ####
Y = t(as.matrix(sapply(Y, as.numeric)))
K=15
W = array(NA, dim = c(100,100,K))
W[,,1] = W_state
W[,,2] = W_party
W[,,3] = W_class
W[,,4] = W_age
W[,,5] = W_gender
W[,,6] = W_lawyer
W[,,7] = W_executive
W[,,8] = W_business
W[,,9] = W_farmer
W[,,10] = W_army
W[,,11] = W_teacher
W[,,12] = W_professor
W[,,13] = W_tweets
W[,,14] = W_followers
W[,,15] = W_follow_adj

n = dim(Y)[2] ; p =dim(Y)[1]

#Initialize fold ID for block cross-validation
fold_id = rep(rep(1:10, each = ceiling(n /10)), times = p)
extra = ceiling(n /10) * 10 - n
extra_indices = which(fold_id == 10) [(which(fold_id == 10) -1 ) %% ceiling(n /10) >= (ceiling(n/10) -extra )]
fold_id = fold_id[-extra_indices]

#Fit Ising similarity Regression Model
real_data_result = estimate_sparse_Ising_similarity(Y = Y, W =  W, fold_id = fold_id)
real_data_result$hat_alpha

#Compute 95% CI based on empirical sandwich covariance matrix as described in Section S.3 of the paper
post_sandwich_cov = empirical_cov(Y_mat = Y,W[,,which(real_data_result$hat_alpha!=0)],real_data_result$hat_theta_jj,real_data_result$hat_alpha[which(real_data_result$hat_alpha!=0)])
post_lower_ci = real_data_result$hat_alpha[which(real_data_result$hat_alpha!=0)] - qnorm(0.975) * sqrt(diag(post_sandwich_cov))[-c(1:100)]
post_upper_ci = real_data_result$hat_alpha[which(real_data_result$hat_alpha!=0)] + qnorm(0.975) * sqrt(diag(post_sandwich_cov))[-c(1:100)]


lower_vec = rep(0,K); upper_vec = rep(0,K)
lower_vec[which(real_data_result$hat_alpha!=0)] = post_lower_ci
upper_vec[which(real_data_result$hat_alpha!=0)] = post_upper_ci

data.frame(similarity_measure = c('State','Party','Class', 'Age','Gender','Lawyer','Executive','Business','Farmer','Army','Teacher','Professor','Tweets','Followers','Twitter Follower-Followee Relationship'),
           hatalpha = round(real_data_result$hat_alpha,3), 
           lower = round(lower_vec,3), upper = round(upper_vec,3))

#Histogram of hat_theta_jj's
hist_hat_theta_jj= ggplot(data = data.frame(hat_theta_jj = real_data_result$hat_theta_jj), aes(x = hat_theta_jj)) + geom_histogram(fill = 'grey', color='black') + theme_bw() + ylab('Frequency') + 
  xlab(bquote(hat(theta)[jj])) + ggtitle(bquote(hat(theta)[jj]~"of 100 Senators")) +  theme(plot.title = element_text(hjust = 0.5))
hist_hat_theta_jj

#Weighted Graphs for State/Party/Twitter follower-followee relationship/Theta
theta_mat = Ximat(W,real_data_result$hat_alpha)
num_sub_senator = 20
length(unique(c(which(theta_mat >= tail(sort(theta_mat), n=20)[1], arr.ind=TRUE))))
sub_senators = sort(unique(c(which(theta_mat >= tail(sort(theta_mat), n=20)[1], arr.ind=TRUE))))

set.seed(1)
sub_W_mat = W[sub_senators, sub_senators, c(15,1,2)]
sub_alpha = real_data_result$hat_alpha[c(15,1,2)]

party[party == 'R'] = 'Republican'
party[party == 'D'] = 'Democrat'
#Color for party: Democrat-Cyan, Republican-Purple
color_party = rep(NA,100)
color_party[party == "Democrat"] = "cyan"
color_party[party == "Republican"] = "purple"
color_party = color_party[sub_senators]
names(color_party) = party[sub_senators]


#Twitter follower-followee relationship 
twitter_net = as.network(sub_W_mat[,,1], directed =  FALSE, matrix.type = 'adjacency', print.adj = T, ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(twitter_net) = full_names[sub_senators]
twitter_net %v% "state" = state[sub_senators]
twitter_net %v% "party" = party[sub_senators]
twitter_net %v% "color" = color_party
x = gplot.layout.fruchtermanreingold(twitter_net, NULL)
twitter_net %v% "x" = x[, 1]
twitter_net %v% "y" = x[, 2]

#State
state_net = as.network(sub_W_mat[,,2], directed =  FALSE, matrix.type = 'adjacency', print.adj = T, ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(state_net) = full_names[sub_senators]
state_net %v% "state" = state[sub_senators]
state_net %v% "party" = party[sub_senators]
state_net %v% "color" = color_party
state_net %v% "x" = x[, 1]
state_net %v% "y" = x[, 2]

#Party
party_net = as.network(sub_W_mat[,,3], directed =  FALSE, matrix.type = 'adjacency', print.adj = T, ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(party_net) = full_names[sub_senators]
party_net %v% "state" = state[sub_senators]
party_net %v% "party" = party[sub_senators]
party_net %v% "color" = color_party
party_net %v% "x" = x[, 1]
party_net %v% "y" = x[, 2]


state_plot = ggnet2(state_net, mode = c("x","y"),label = "state", label.color = "color", alpha = 0, edge.size = sub_alpha[2]) + ggtitle('State') + theme(plot.title = element_text(hjust = 0.5))
party_plot = ggnet2(party_net, mode = c("x","y"),label = "state", label.color = "color", alpha = 0, edge.size = sub_alpha[3]) + ggtitle('Party') + theme(plot.title = element_text(hjust = 0.5))
twitter_plot = ggnet2(twitter_net,mode = c("x","y") ,label = "state", label.color = "color", alpha = 0, edge.size = sub_alpha[1]) + ggtitle('Twitter') + theme(plot.title = element_text(hjust = 0.5))

#HatTheta
theta_mat = Ximat(W,real_data_result$hat_alpha)
sub_theta_mat = theta_mat[sub_senators,sub_senators]
theta_net = as.network(sub_theta_mat, directed = FALSE, matrix.type = 'adjacency', ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(theta_net) = full_names[sub_senators]
theta_net %v% "state" = state[sub_senators]
theta_net %v% "party" = party[sub_senators]
theta_net %v% "color" = color_party
theta_net %v% "x" = x[, 1]
theta_net %v% "y" = x[, 2]

set.edge.attribute(theta_net, "edge_linetype", ifelse(theta_net %e% "weight" > 0, 1, 3))
set.edge.attribute(theta_net, "edge_alpha_median", ifelse(theta_net %e% "weight" > median(theta_mat), 1,0))
theta_net %e% 'weight' =  abs(theta_net %e% 'weight')

theta_df = as.data.frame(theta_net)
for(i in 1:nrow(theta_df)){
  if(theta_df$edge_alpha_median[i] == 0){
    theta_net[theta_df$.tail[i], theta_df$.head[i]] = 0
    theta_net[theta_df$.head[i], theta_df$.tail[i]] = 0
  }
}
theta_df = as.data.frame(theta_net)


####Checking if the edges from Theta after removing State and Party edges, is a subset of Twitter edges
theta_edge = as.edgelist(theta_net)
theta_edge_c = paste0(theta_edge[,1],',',theta_edge[,2])
twitter_edge = as.edgelist(twitter_net)
twitter_edge_c = paste0(twitter_edge[,1],',',twitter_edge[,2])
party_edge = as.edgelist(party_net)
party_edge_c = paste0(party_edge[,1],',',party_edge[,2])
state_edge = as.edgelist(state_net)
state_edge_c = paste0(state_edge[,1],',',state_edge[,2])


#theta_edge_c_removing_state_party =theta_edge_c[ (!theta_edge_c %in% state_edge_c )& (!theta_edge_c %in% party_edge_c) ]
theta_edge_c_removing_state_party =theta_edge_c[ !(theta_edge_c %in% union(state_edge_c, party_edge_c)) ]
mean(theta_edge_c_removing_state_party %in%  twitter_edge_c)

### Setting edges in theta_net that is between pairs of senators from different party AND different state to a different color 
edge_color_mat = matrix(0.5, nrow = 20, ncol =20)
for(i in 1:length(strsplit(theta_edge_c_removing_state_party, ','))){
  row_ind = as.numeric(strsplit(theta_edge_c_removing_state_party, ',')[[i]][1])
  col_ind = as.numeric(strsplit(theta_edge_c_removing_state_party, ',')[[i]][2])
  edge_color_mat[row_ind,col_ind] = edge_color_mat[col_ind,row_ind] = 1
}
diag(edge_color_mat) = 0

edge_color_net = as.network(edge_color_mat , directed = FALSE, matrix.type = 'adjacency', ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(edge_color_net) = full_names[sub_senators]
edge_color_net %v% "state" = state[sub_senators]
edge_color_net %v% "party" = party[sub_senators]
edge_color_net %v% "color" = color_party
edge_color_net %v% "x" = x[, 1]
edge_color_net %v% "y" = x[, 2]

#
final_theta_net  = as.network(sub_theta_mat, directed = FALSE, matrix.type = 'adjacency', ignore.eval = FALSE, names.eval = 'weight')
network.vertex.names(final_theta_net) = full_names[sub_senators]
final_theta_net %v% "state" = state[sub_senators]
final_theta_net %v% "party" = party[sub_senators]
final_theta_net %v% "color" = color_party
final_theta_net %v% "x" = x[, 1]
final_theta_net %v% "y" = x[, 2]

set.edge.attribute(final_theta_net, "edge_linetype", ifelse(final_theta_net %e% "weight" > 0, 1, 3))
set.edge.attribute(final_theta_net, "edge_alpha_median", ifelse(final_theta_net %e% "weight" > median(theta_mat), 1,0))
final_theta_net %e% 'weight' =  abs(final_theta_net %e% 'weight')
set.edge.attribute(final_theta_net, "edge_color", ifelse(edge_color_net %e% "weight" == 1, 'red','grey50'))

theta_plot= ggnet2(final_theta_net, mode = c("x","y"), color = 'party', color.legend = 'Party', label = "state", label.color = "color", alpha = 0, edge.size = 'weight',edge.lty = 'edge_linetype', edge.alpha = final_theta_net %e% 'edge_alpha_median', edge.color = 'edge_color') + ggtitle('Theta') + theme(plot.title = element_text(hjust = 0.5))+ 
  scale_color_manual(name = 'Party',values = c('cyan','purple')) 


ggarrange(state_plot,party_plot,twitter_plot,theta_plot, common.legend = T, legend = 'bottom')
