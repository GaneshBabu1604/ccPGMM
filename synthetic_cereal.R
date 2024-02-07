#Loading the scripts to implement PGMM(UUU), constrained PGMM and consensus clustering.
source('PGMM.R')
source('constrained_PGMM.R')
source('consensus_Clustering.R')

# Dataset with synthetic cereal
#Loading Threshold based labels for Nine hyperspectral images (three for each cereal type)
cluster_labels = read.csv('Threshold_Labels_Full_Data.csv')
cluster_labels = cluster_labels$X0

# Loading Original Cereal
original_cereal_data = read.csv('Data.csv')
Data = original_cereal_data
sim_data = 1
Q = 2
seed_list = c(100,200,300,400,500)
#Generating background data from factanal of known background pixels (based on threshold)
set.seed(seed_list[sim_data])
A = factanal(original_cereal_data[which(cluster_labels==1),2:102],Q,scores = 'regression')
means = apply(original_cereal_data[which(cluster_labels==1),2:102],2,mean)
new_means = rep(means,length(which(cluster_labels==1)))
new_means = matrix(new_means,101,length(which(cluster_labels==1)))
errors = t(mvrnorm(length(which(cluster_labels==1)),rep(0,101),diag(A$uniquenesses)))
scores = mvrnorm(length(which(cluster_labels==1)),rep(0,Q),diag(rep(1,Q)))
M = A$loadings%*% t(scores)
A_background = t(means + M + errors)
Data[which(cluster_labels==1),2:102] = A_background

#Generating wheat data from factanal of known wheat pixels (based on threshold)
A = factanal(original_cereal_data[which(cluster_labels==2),2:102],Q,scores = 'regression')
means = apply(original_cereal_data[which(cluster_labels==2),2:102],2,mean)
new_means = rep(means,length(which(cluster_labels==2)))
new_means = matrix(new_means,101,length(which(cluster_labels==2)))
errors = t(mvrnorm(length(which(cluster_labels==2)),rep(0,101),diag(A$uniquenesses)))
scores = mvrnorm((length(which(cluster_labels==2))),rep(0,Q),diag(rep(1,Q)))
M = A$loadings%*% t(scores)
A_wheat = t(means + M + errors)
Data[which(cluster_labels==2),2:102] = A_wheat

#Generating corn data from factanal of known corn pixels (based on threshold)
A = factanal(original_cereal_data[which(cluster_labels==3),2:102],Q,scores = 'regression')
means = apply(original_cereal_data[which(cluster_labels==3),2:102],2,mean)
new_means = rep(means,length(which(cluster_labels==3)))
new_means = matrix(new_means,101,length(which(cluster_labels==3)))
errors = t(mvrnorm(length(which(cluster_labels==3)),rep(0,101),diag(A$uniquenesses)))
scores = mvrnorm(length(which(cluster_labels==3)),rep(0,Q),diag(rep(1,Q)))
M = A$loadings%*% t(scores)
A_corn = t(means + M + errors)
Data[which(cluster_labels==3),2:102] = A_corn

#Generating rice data from factanal of known rice pixels (based on threshold)
A = factanal(original_cereal_data[which(cluster_labels==4),2:102],Q,scores = 'regression',lower = '0.0053')
means = apply(original_cereal_data[which(cluster_labels==4),2:102],2,mean)
new_means = rep(means,length(which(cluster_labels==4)))
new_means = matrix(new_means,101,length(which(cluster_labels==4)))
errors = t(mvrnorm(length(which(cluster_labels==4)),rep(0,101),diag(A$uniquenesses)))
scores = mvrnorm(length(which(cluster_labels==4)),rep(0,Q),diag(rep(1,Q)))
M = A$loadings%*% t(scores)
A_rice = t(means + M + errors)
Data[which(cluster_labels==4),2:102] = A_rice

original_cereal_data$X<-NULL
Data$X<-NULL

# Loading the larger(43.5%) selected constrained pixels based on the location in the HI image.
background_constraints_selected = read.csv('Full_data_background_constraints_random_larger.csv')
wheat_constraints_selected = read.csv('Full_data_wheat_constraints_random_larger.csv')
corn_constraints_selected = read.csv('Full_data_corn_constraints_random_larger.csv')
rice_constraints_selected = read.csv('Full_data_rice_constraints_random_larger.csv')
background_constraints_selected = background_constraints_selected$X0
wheat_constraints_selected = wheat_constraints_selected$X0
corn_constraints_selected = corn_constraints_selected$X0
rice_constraints_selected = rice_constraints_selected$X0
B_selected_constraints = list(background_constraints_selected,wheat_constraints_selected,corn_constraints_selected,rice_constraints_selected)

# Loading the smaller(24.7%) selected constrained pixels based on the location in the HI image.
background_constraints_selected = read.csv('Full_data_background_constraints_random_smaller.csv')
wheat_constraints_selected = read.csv('Full_data_wheat_constraints_random_smaller.csv')
corn_constraints_selected = read.csv('Full_data_corn_constraints_random_smaller.csv')
rice_constraints_selected = read.csv('Full_data_rice_constraints_random_smaller.csv')
background_constraints_selected = background_constraints_selected$X0
wheat_constraints_selected = wheat_constraints_selected$X0
corn_constraints_selected = corn_constraints_selected$X0
rice_constraints_selected = rice_constraints_selected$X0
B_selected_constraints = list(background_constraints_selected,wheat_constraints_selected,corn_constraints_selected,rice_constraints_selected)


# Parallel run
library(parallel)  
library(doParallel)
no_cores <- detectCores() - 1
no_cores
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
clusterEvalQ(cl, {
  library(LaplacesDemon)
  library(ggplot2)
  library(caTools)
  library(mclust)
  library(pgmm)
  library(fossil)
  library(MASS)
  library(IMIFA)
  library(MixSim)
})

# Fitting PGMM and constrained-PGMM on different settings of M & d
N = nrow(Data)
G = 4
Q = 2
e = 25000
# Generating randomly selected M subsets of features. Features per subset d
seed = seed_list[sim_data]
for(d in c(10,20)){
  set.seed(seed*1000)
  Y = list()
  for(M in c(10,25,50,100)){
    X = list()
    for (i in 1:M){#Change the Data object to just columns
      cols = sample(1:ncol(Data),d,replace = F)
      x = list(cols)
      X = append(X,list(x))
    }
    Y = append(Y,list(X))
  }
  print(paste('simulation',as.character(sim_data),'begins'))
  M_index = 1
  for(M in c(10,25,50,100)){
    for (constrained in c(0,1)){
      start_time = Sys.time()
      print(c(M,d,constrained))
      if(constrained == 0){
        #constrained_pgmmAECM with no constraints
        result <- foreach(m = Y[[M_index]]) %dopar% constrained_pgmmAECM_implementer(Data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = NULL)
      }
      else{
        #constrained_pgmmAECM with constraints
        result <- foreach(m = Y[[M_index]]) %dopar% constrained_pgmmAECM_implementer(Data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = B_selected_constraints)
      }
      print('Models computed')
      end_time = as.numeric(Sys.time() - start_time,units = "mins")
      print(end_time)
      parallel_list = list(no_models = M,constraint = constrained,parallel_time = end_time)
      n_result = 1
      for(s in result){
        location = paste('syntheticcereal_Data_',as.character(sim_data),'_M',as.character(M),'_d',as.character(d),'_cons',as.character(constrained),'_',as.character(n_result),'.RData',sep = '')
        save(s,file = location)
        n_result=n_result+1
      }
      location = paste('syntheticcereal_Data_',as.character(sim_data),'_M',as.character(M),'_d',as.character(d),'_cons',as.character(constrained),'_creation_time.RData',sep = '')
      save(parallel_list,file=location)
    }
    M_index = M_index + 1
  }    
}


stopCluster(cl)
# Consenus clustering with equal (e = 0) and entropy based weights (e = 1) and saving the consensus cluster solution for c-PGMM (c = 0) and cc-PGMM (c = 1)
for(d in c(10,20)){
  for(M in c(10,25,50,100)){
    for(c in c(0,1)){
      for(e in c(0,1)){
        print(c(c,e))
        loc = paste('syntheticcereal_Data_',as.character(sim_data),sep ='')
        sim_mat_hclust(loc,d,M,c,e,G,N)        
      }
    }    
  }
}

# Computing ARI
final_results = data.frame()
for(d in c(10,20)){
  for(M in c(10,25,50,100)){
    loc = paste('syntheticcereal_Data_',as.character(sim_data),sep ='')
    for(c in c(0,1)){
      for(e in c(0,1)){
        saved_loc = paste(loc,'_M',as.character(M),'_d',as.character(d),'_cons',as.character(c),'_entropy',as.character(e),'_consensus_cluster.RData',sep = '')
        load(saved_loc)
        final_results = rbind(final_results,cbind(sim_data,c,e,d,M,adjustedRandIndex(result[[1]],cluster_labels),result[[8]],result[[9]]))
      }
    }    
  }
}

# Fitting PGMM on p = 101 features
full_pgmm_kmeans = pgmmAECM(data.frame(Data),G = 4, Q = 1, epoch = 100000,zstart = 1)

#Fitting Mclust on p = 101 features
full_mclust = Mclust(Data,G = 4)

#Fitting DBScan
library(dbscan)
set.seed(100)
kNNdistplot(as.matrix(Data), minPts = 202)
abline(h=20, lt=2, col="blue")
full_dbscan <- dbscan(as.matrix(Data), minPts = 202, eps=20)