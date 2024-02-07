library(sys)
sim_mat_creator <- function(file_loc,d,n_models,c,e,N){
  #Input:
  #file_loc : Location of the saved results of constrained PGMM applied to each subset.
  #d : Number of variables per data subset.
  #n_models: Number of data subset.
  #c: 0 or 1. 1, if constraints are implemented.
  #e: 0 or 1. 1, if entropy weights should be used in consensus approach.
  #N: Number of observations.
  #Output:
  #Returns a list of consensus similarity matrix, the number of converged models, list of models which converged, number of data subsets,
  #number of variables per subset, c, e, time taken if models are implemented in series and the time taken if the models are implemented in parallel.
  start_time = Sys.time()
  models_series_time = 0
  parallel_build_time = 0
  series_build_time = 0
  fin_sim = matrix(0,nrow = N,ncol = N)
  converged = 0
  converged_list = c()
  counter = 0
  weights = 0
  temp = 1
  load_location = paste(file_loc,'_M',as.character(n_models),'_d',as.character(d),'_cons',as.character(c),'_creation_time.RData',sep = '')
  load(load_location)
  models_parallel_time = parallel_list$parallel_time
  #current_result = head(result,r)
  for (i in c(1:n_models)){
    print(i)
    load_location = paste(file_loc,'_M',as.character(n_models),'_d',as.character(d),'_cons',as.character(c),'_',as.character(i),'.RData',sep = '')
    load(load_location)
    if (length(s)==7){
      if(s[[2]] == 'converged'){
        converged = converged + 1
        converged_list = c(converged_list,i)
      }
      print('Computing Sim')
      temp_sim = s[[1]]%*%t(s[[1]])
      diag(temp_sim) = 1
      if(e==0){
        print('Adding Sim')
        fin_sim = fin_sim + temp_sim
        temp_sim = NULL
      }
      else{
        fin_sim = fin_sim + (temp_sim*s[[4]])
        temp_sim = NULL
      }
      weights = weights + s[[4]]
      models_series_time = models_series_time + s[[5]]
    }
  }
  if(e==0){
    fin_sim = fin_sim/n_models
  }
  else{
    fin_sim = fin_sim/weights
  }
  end_time = as.numeric(Sys.time() - start_time,units = 'mins')
  series_build_time = models_series_time + end_time
  parallel_build_time = models_parallel_time + end_time
  print('A')
  sim_mat_list = list(fin_sim,converged,converged_list,n_models,d,c,e,series_build_time,parallel_build_time)
  return(sim_mat_list)
}
#Performing consensus clustering
sim_mat_hclust <- function(file_loc,d,n_models,c,e,G,N){
  #Input:
  #file_loc : Location of the saved results of constrained PGMM applied to each subset.
  #d : Number of variables per data subset.
  #n_models: Number of data subset.
  #c: 0 or 1. 1, if constraints are implemented.
  #e: 0 or 1. 1, if entropy weights should be used in consensus approach.
  #G: Number of clusters.
  #N: Number of observations.
  #Output:
  #Saves the consensus cluster solution in the specified location(file_loc)
  #The saved list contains the consensus cluster solution, the number of converged models, array of models which converged, number of data subsets,
  #number of variables per subset, c, e, time taken if models are implemented in series and the time taken if the models are implemented in parallel.
  sim_mat_list = sim_mat_creator(file_loc,d,n_models,c,e,N)
  start_time = Sys.time()
  fin_dis_sim = 1 - sim_mat_list[[1]]
  consensus_hclust_data = hclust(as.dist(fin_dis_sim),method = "complete")
  groups = cutree(consensus_hclust_data, k=G)
  end_time = as.numeric(Sys.time() - start_time,units = 'mins')
  sim_mat_list[[8]] = sim_mat_list[[8]] + end_time
  sim_mat_list[[9]] = sim_mat_list[[9]] + end_time
  save_location = paste(file_loc,'_M',as.character(n_models),'_d',as.character(d),'_cons',as.character(c),'_entropy',as.character(e),'_consensus_cluster.RData',sep = '')
  result = list(groups,sim_mat_list[[2]],sim_mat_list[[3]],sim_mat_list[[4]],sim_mat_list[[5]],sim_mat_list[[6]],sim_mat_list[[7]],sim_mat_list[[8]],sim_mat_list[[9]])
  save(result,file = save_location)
}