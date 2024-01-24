source('PGMM.R')
source('constrained_PGMM.R')
source('consensus_Clustering.R')

# Dataset with well separated cluster
cluster_labels = read.csv('Threshold_Labels_Full_Data.csv')
cluster_labels = cluster_labels$X0

# Original Cereal
Data = read.csv('Data.csv')
Data$X<-NULL

# Loading the larger selected constrained pixels based on the location in the HI image.
background_constraints_selected = read.csv('Full_data_background_constraints_random_larger.csv')
wheat_constraints_selected = read.csv('Full_data_wheat_constraints_random_larger.csv')
corn_constraints_selected = read.csv('Full_data_corn_constraints_random_larger.csv')
rice_constraints_selected = read.csv('Full_data_rice_constraints_random_larger.csv')
background_constraints_selected = background_constraints_selected$X0
wheat_constraints_selected = wheat_constraints_selected$X0
corn_constraints_selected = corn_constraints_selected$X0
rice_constraints_selected = rice_constraints_selected$X0
B_selected_constraints = list(background_constraints_selected,wheat_constraints_selected,corn_constraints_selected,rice_constraints_selected)


# Loading the smaller selected constrained pixels based on the location in the HI image.
background_constraints_selected = read.csv('Full_data_background_constraints_random_smaller.csv')
wheat_constraints_selected = read.csv('Full_data_wheat_constraints_random_smaller.csv')
corn_constraints_selected = read.csv('Full_data_corn_constraints_random_smaller.csv')
rice_constraints_selected = read.csv('Full_data_rice_constraints_random_smaller.csv')
background_constraints_selected = background_constraints_selected$X0
wheat_constraints_selected = wheat_constraints_selected$X0
corn_constraints_selected = corn_constraints_selected$X0
rice_constraints_selected = rice_constraints_selected$X0
B_selected_constraints = list(background_constraints_selected,wheat_constraints_selected,corn_constraints_selected,rice_constraints_selected)


# Fitting PGMM and constrained-PGMM on different settings of r & d
constrained_pgmmAECM_implementer<-function(x,G=4,Q=2,epoch = 10000, zstart = 1,zlist = c(), seed = 123456, B = NULL){
  start_time = Sys.time()
  new <- try(constrained_pgmmAECM(data.frame(x),G=G,Q=Q,epoch = epoch, zstart = zstart, zlist = zlist,seed = seed, B=B),silent = T)
  if(is.list(new)==TRUE){
    prod_probability = new$probability_score*log(new$probability_score)
    prod_probability[is.na(prod_probability)] = 0
    weight <- -1*sum(colSums(prod_probability))
    weighted_B <- (1/weight)
    end_time = as.numeric(Sys.time() - start_time,units='mins')
    return(list(new$probability_score,new$converged,new$predicted_cluster,weighted_B,end_time,new$log_likelihood_list,new$initial_cluster)) 
  }
  else{
    return(new)
  }
}

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
library(ggplot2)

sim_data = 2
G = 4
Q = 2
e = 25000
# Generating randomly selected subsets r of features. Features per subset d
seed = 200
for(d in c(20)){
  set.seed(seed*1000)
  Y = list()
  for(r in c(25)){
    X = list()
    for (i in 1:r){#Change the Data object to just columns
      cols = sample(1:ncol(Data),d,replace = F)
      x = list(cols)
      X = append(X,list(x))
    }
    Y = append(Y,list(X))
  }
  print(paste('simulation',as.character(sim_data),'begins'))
  r_index = 1
  for(r in c(25)){
    for (constrained in c(1)){
      start_time = Sys.time()
      print(c(r,d,constrained))
      if(constrained == 0){
        #pgmmAECM is constrained_pgmmAECM with no constraints
        result <- foreach(m = Y[[r_index]]) %dopar% constrained_pgmmAECM_implementer(Data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = NULL)
      }
      else{
        #constrained_pgmmAECM with constraints
        result <- foreach(m = Y[[r_index]]) %dopar% constrained_pgmmAECM_implementer(Data[,m[[1]]],G=G,Q=Q,epoch = e, zstart = 1, B = B_selected_constraints)
      }
      print('Models computed')
      end_time = as.numeric(Sys.time() - start_time,units = "mins")
      print(end_time)
      parallel_list = list(no_models = r,constraint = constrained,parallel_time = end_time)
      n_result = 1
      for(s in result){
        location = paste('original_Data_',as.character(sim_data),'_r',as.character(r),'_d',as.character(d),'_cons',as.character(constrained),'_',as.character(n_result),'.RData',sep = '')
        save(s,file = location)
        n_result=n_result+1
      }
      location = paste('original_Data_',as.character(sim_data),'_r',as.character(r),'_d',as.character(d),'_cons',as.character(constrained),'_creation_time.RData',sep = '')
      save(parallel_list,file=location)
    }
    r_index = r_index + 1
  }    
}

stopCluster(cl)
# Consenus clustering with equal (e = 0) and entropy based weights (e = 1) and saving the consensus cluster solution for c-PGMM (c = 0) and cc-PGMM (c = 1)
r = 25
N = 28039
for(c in c(1)){
  for(e in c(0)){
    print(c(c,e))
    loc = paste('original_Data_',as.character(sim_data),sep ='')
    sim_mat_hclust(loc,d,r,c,e,G,N)        
  }
}

# Computing ARI
r = 25
loc = paste('original_Data_',as.character(sim_data),sep ='')
final_results = data.frame()
for(c in c(1)){
  for(e in c(0)){
    saved_loc = paste(loc,'_r',as.character(r),'_d',as.character(d),'_cons',as.character(c),'_entropy',as.character(e),'_consensus_cluster.RData',sep = '')
    load(saved_loc)
    final_results = rbind(final_results,cbind(sim_data,c,e,d,r,adjustedRandIndex(result[[1]],cluster_labels),result[[8]],result[[9]]))
  }
}

#Fitting Mclust on p = 101 features
library(mclust)
set.seed(100)
start_time = Sys.time()
out = Mclust(Data,G = 4,modelNames = c('VVE'))
end_time = as.numeric(Sys.time() - start_time,units = "mins")
full_mclust = list(model_output = out,time_taken = end_time)


# Fitting PGMM on p = 101 features
full_pgmm_kmeans = pgmmAECM(Data,G = 4, Q = 2, epoch = 100000,zstart = 1)

#Fitting DBScan
library(dbscan)
set.seed(100)
kNNdistplot(as.matrix(Data), minPts = 202)
abline(h=0.19, lt=2, col="blue")
full_dbscan <- dbscan(as.matrix(Data), minPts = 202, eps=0.19)
adjustedRandIndex(full_dbscan$cluster,cluster_labels)
