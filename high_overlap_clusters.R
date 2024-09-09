#Loading the scripts to implement PGMM(UUU), constrained PGMM and consensus clustering.
source('PGMM.R')
source('constrained_PGMM.R')
source('consensus_Clustering.R')

#Loading the threshold labels for three of the nine hyperspectral images (one for each cereal type)
cluster_labels = read.csv('Threshold_Labels_Subset_Data.csv')
cluster_labels = cluster_labels$X0

#Generating dataset with heavily overlapping clusters
sim_data = 1
seed_list = c(100,200,300,400,500)
set.seed(seed_list[sim_data])
#Number of factors
Q = 3
#Number of observations
N = length(cluster_labels)
#Number of variables
P = 101
#Number of clusters
G = 4
#Mean Matrix
mu = array(dim = c(G,P))
#Loadings Matrix
lambda = array(dim = c(P,Q,G))
#Uniqueness
psi = array(0,c(P,P,G))
temp_simulated_data = matrix(nrow = N,ncol = P)
correlation_span = 5
correlation_strength = 0.03
#Heavily overlapping data
#mu_g generator. Generated from a normal distribution.
mean_list = c(-1.25,0,1.25,2.5)
sd_list = c(1.5,1.5,1.5,1.5)

for(g in 1:G){
  mu[g,] <- rnorm(P,mean = mean_list[g],sd = sd_list[g])
}

#loading matrix (lambda_g) and uniqueness (psi_g) generation
for(g in 1:G){
  start_p = 1
  end_p = correlation_span
  while(start_p<=P){
    temp_mean = runif(1,0.3,0.9)
    lambda[start_p:end_p,,g] = matrix(rnorm(((end_p-start_p)+1)*Q,mean = temp_mean,sd = correlation_strength),nrow = ((end_p-start_p)+1),ncol = Q)
    if(end_p-start_p == 0){
      psi[start_p:end_p,start_p:end_p,g] = rgamma(((end_p-start_p)+1),runif(1,0,0.1),1)
    }
    else{
      psi[start_p:end_p,start_p:end_p,g] = diag(rgamma(((end_p-start_p)+1),runif(1,0,0.1),1))
    }
    start_p = start_p + correlation_span
    end_p = end_p + correlation_span
    if(end_p>P){
      end_p = P
    }
  }
}

#Generating the heavily overlapping data from four heavily overlapping Factor Analyzers(One for each cereal type and background)
for(g in 1:G){
  scores = mvrnorm(length(which(cluster_labels==g)),rep(0,Q),diag(Q))
  errors = t(mvrnorm(length(which(cluster_labels==g)),rep(0,P),psi[,,g]))
  lambda_scores = lambda[,,g]%*%t(scores)
  temp_simulated_data[which(cluster_labels==g),] =mu[g,]+t(lambda_scores+errors)
}

# Loading the selected constrained pixels based on the location in the HI image.
background_simulated_constraints_selected = read.csv('background_constraints_subset_larger.csv')
wheat_simulated_constraints_selected = read.csv('wheat_constraints_subset_larger.csv')
corn_simulated_constraints_selected = read.csv('corn_constraints_subset_larger.csv')
rice_simulated_constraints_selected = read.csv('rice_constraints_subset_larger.csv')
background_simulated_constraints_selected = background_simulated_constraints_selected$X0
wheat_simulated_constraints_selected = wheat_simulated_constraints_selected$X0
corn_simulated_constraints_selected = corn_simulated_constraints_selected$X0
rice_simulated_constraints_selected = rice_simulated_constraints_selected$X0
B_selected_constraints = list(background_simulated_constraints_selected,wheat_simulated_constraints_selected,corn_simulated_constraints_selected,rice_simulated_constraints_selected)

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
G = 4
Q = 1
e = 15000
# Generating randomly selected M subsets of features. Features per subset d
seed = seed_list[sim_data]
for(d in c(10,20)){
  set.seed(seed)
  Y = list()
  for(M in c(10,25,50,100)){
    X = list()
    for (i in 1:M){#Change the Data object to just columns
      cols = sample(1:ncol(temp_simulated_data),d,replace = F)
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
        #pgmmAECM is constrained_pgmmAECM with no constraints
        result <- foreach(m = Y[[M_index]]) %dopar% constrained_pgmmAECM_implementer(temp_simulated_data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = NULL)
      }
      else{
        #constrained_pgmmAECM with constraints
        result <- foreach(m = Y[[M_index]]) %dopar% constrained_pgmmAECM_implementer(temp_simulated_data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = B_selected_constraints)
      }
      print('Models computed')
      end_time = as.numeric(Sys.time() - start_time,units = "mins")
      print(end_time)
      parallel_list = list(no_models = M,constraint = constrained,parallel_time = end_time)
      n_result = 1
      for(s in result){
        location = paste('heavilyoverlappingcluster_Data_',as.character(sim_data),'_M',as.character(M),'_d',as.character(d),'_cons',as.character(constrained),'_',as.character(n_result),'.RData',sep = '')
        save(s,file = location)
        n_result=n_result+1
      }
      location = paste('heavilyoverlappingcluster_Data_',as.character(sim_data),'_M',as.character(M),'_d',as.character(d),'_cons',as.character(constrained),'_creation_time.RData',sep = '')
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
        loc = paste('heavilyoverlappingcluster_Data_',as.character(sim_data),sep ='')
        sim_mat_hclust(loc,d,M,c,e,G,N)        
      }
    }    
  }
}

# Computing ARI
final_results = data.frame()
for(d in c(10,20)){
  for(M in c(10,25,50,100)){
    loc = paste('heavilyoverlappingcluster_Data_',as.character(sim_data),sep ='')
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
full_pgmm_kmeans = pgmmAECM(data.frame(temp_simulated_data),G = 4, Q = 1, epoch = 100000,zstart = 1)

#Fitting Mclust on p = 101 features
full_mclust = Mclust(temp_simulated_data,G = 4)

#Fitting DBScan
library(dbscan)
set.seed(100)
kNNdistplot(as.matrix(temp_simulated_data), minPts = 202)
abline(h=17, lt=2, col="blue")
full_dbscan <- dbscan(as.matrix(temp_simulated_data), minPts = 202, eps=17)
adjustedRandIndex(full_dbscan$cluster,cluster_labels)
