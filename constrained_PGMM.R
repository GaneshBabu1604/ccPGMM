library(LaplacesDemon)
library(ggplot2)
library(caTools)
library(mclust)
library(pgmm)
library(fossil)
library(MASS)
library(sys)
library(MixSim)
library(IMIFA)


#Creating a Parsimonious Gaussian Mixture Model(UUU) Alternating Expectation - Conditional Maximization Algorithm:
constrained_pgmmAECM <- function(x,G = 1,Q = 2,epoch=100, convergence_value = 1e-6,zstart=1,zlist=c(),B=NULL,seed=123456){
  #Input:
  #x : Dataframe(N,P). Where N is the number of observations for creating the clusters and P is the features of each observation.
  #epoch : int. Number of iterations.
  #G: int. Number of clusters.
  #Q: int. Number of factors.
  #convergence_value: convergence score. log-likelihood ratio is used to assess convergence.
  #zstart: int. kmeans starting values (1) or user defined starting values(2)
  #zlist: vector. User defined starting values. Applicable only when zstart = 2
  #B: list of lists. Positive and negative constraints. Observations within a list are positively constrained.
  #   Observations in different list are negatively constrained.
  #seed: int. seed for kmeans clustering. Applicable only when zstart = 1
  
  #Output:
  #predicted_cluster: final cluster membership of observations.
  #initial_cluster: initial cluster membership of observations.
  #predicted_score: matrix of posterior probability of cluster membership.
  #log_likelihood_list: vector of log-likelihood.
  #converged: whether the model is converged or not.
  #bic: Bayesian information criterion.
  #time: time taken to fit the model.
  
  #computing the total number of training pixels
  start_time = Sys.time()
  N=nrow(x)
  #computing the number of features for training pixels
  P=ncol(x)
  #Creating a initial cluster with kmeans clustering
  if (zstart == 1){
    set.seed(seed)
    init_cluster = kmeans(x,G,nstart = 5)$cluster
  }
  else if(zstart==2){
    init_cluster = zlist
  }
  #Creating the initial data with initial clusters for maximum log-likelihood estimation
  initial_data=cbind(x,init_cluster)
  colnames(initial_data)=c(colnames(x),"init_cluster")
  #Computing the initial parameters
  tau = matrix(0,G,1)
  mu = matrix(0,G,P)
  lambda = array(0,c(P,Q,G))
  psi = array(0,c(P,P,G))
  sigma = array(0,c(P,P,G))
  beta = matrix(0,Q,P)
  theta = matrix(0,Q,Q)
  s = array(0,c(P,P,G))
  Log_Likelihood = c()
  predicted_cluster=c()
  z=matrix(0,N,G)
  z_hat=matrix(0,N,G)
  z_hat=unmap(initial_data[,"init_cluster"])
  covw_computed = covw(x,z_hat)
  #Initialization
  for (g in 1:G){
    #Computing the probability of belonging to each class
    tau[g,]=nrow(initial_data[initial_data[,'init_cluster']==g,])/N
    #Computing the mean of each class
    mu[g,]=covw_computed$mean[,g]
    #Computing the covariance matrix of each class
    s[,,g]=covw_computed$S[,,g]
    eigen_g=eigen(s[,,g])
    for (q in 1:Q){
      for (p in 1:P){
        if (eigen_g$vectors[p,q] < 0){
          lambda[p,q,g] = -exp((0.5*log(eigen_g$values[q])+log(-eigen_g$vectors[p,q])))
        }
        else{
          lambda[p,q,g] = exp((0.5*log(eigen_g$values[q])+log(eigen_g$vectors[p,q])))
        }
      }
    }
    diag(psi[,,g]) = abs(diag(s[,,g]-lambda[,,g]%*%t(lambda[,,g])))
    # Updating sigma_g based on latent factors
    sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
  }
  # AECM Algorithm
  for (e in 0:epoch){
    #Step 1
    #Conditional Maximization
    covw_computed = covw(x,z_hat)
    for (g in 1:G){
      #Computing the probability of belonging to each class
      tau[g,]=(sum(z_hat[,g]))/N
      #Computing the mean of each class
      mu[g,]=covw_computed$mean[,g]
    }
    #Creating a temporary variable
    z_old = z_hat
    #Expectation
    if(e!=0){
      #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
      for (g in 1:G){
        z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
      }
      z = constrained_z_implementer(z,B,G,N)
      #Computing Log-Sum-Exp
      z_max = apply(z,1,max)
      z_sum=log(apply(exp(z-z_max),1,sum))+z_max
      z_hat=exp(z-z_sum)
    }
    # Step 2:
    # Update s,beta,theta and Conditional Maximization
    covw_computed = cov_weighted(x,G,mu,z_hat)
    for (g in 1:G){
      s[,,g]=covw_computed[,,g]
      # Using woodbury identity for computing the Inverse of lambda%*%t(lambda) + psi
      psi_inv = solve(psi[,,g])
      beta = t(lambda[,,g]) %*% (psi_inv - (psi_inv %*% lambda[,,g] %*% solve(diag(Q)+t(lambda[,,g]) %*% psi_inv %*% lambda[,,g]) %*% t(lambda[,,g]) %*% psi_inv))
      theta = diag(Q) - (beta %*% lambda[,,g]) + (beta %*% s[,,g] %*% t(beta))
      # Estimating new lambda
      lambda[,,g] = s[,,g] %*% t(beta) %*% solve(theta)
      # Estimating new psi
      diag(psi[,,g]) = abs(diag(s[,,g]-lambda[,,g]%*%beta%*%s[,,g]))
      # Estimating new sigma based on latent factors
      sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
    }
    # Expectation
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    for (g in 1:G){
      z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log = TRUE)
    }
    temp_z = z
    z = constrained_z_implementer(z,B,G,N)
    #Computing Log-Sum-Exp
    temp_z_max = apply(temp_z,1,max)
    z_max = apply(z,1,max)
    temp_z_sum = log(apply(exp(temp_z-temp_z_max),1,sum))+temp_z_max
    z_sum=log(apply(exp(z-z_max),1,sum))+z_max
    temp_z_hat = exp(temp_z-temp_z_sum)
    z_hat=exp(z-z_sum)
    #Computing the log likelihood of the prediction
    Log_Likelihood =c(Log_Likelihood,sum(temp_z_sum))
    if(e>3){
      l_k_1 = (Log_Likelihood[length(Log_Likelihood)] - Log_Likelihood[length(Log_Likelihood)-1])/Log_Likelihood[length(Log_Likelihood)-1]
      if (abs(l_k_1) < convergence_value){
        converged = 'converged'
        print('Cluster Solution Converged')
        break
      }          
    }
    if (e == epoch){
      converged = 'not converged'
      print('Cluster Solution not Converged')
      break
    }
  }
  bic_model = 2*tail(Log_Likelihood,1) - PGMM_dfree(Q,P,G,'UUU') * log(N)
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  predicted_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  end_time = as.numeric(Sys.time() - start_time,units = "mins")
  return(list(predicted_cluster=predicted_cluster,initial_cluster = init_cluster,probability_score=z_hat,converged = converged,bic = bic_model,log_likelihood_list = Log_Likelihood,time = end_time))
}

cov_weighted = function(x,G=1,mu,z_hat){
  s = array(0,c(ncol(x),ncol(x),G))
  for(g in 1:G){
    x_temp = sweep(x,2,mu[g,],'-')
    x_temp_2 = sweep(x_temp,1,z_hat[,g],'*')
    s[,,g] = (t(x_temp_2) %*% as.matrix(x_temp))/sum(z_hat[,g])
  }
  return(s)
}

#Incorporating the constraints in to posterior probability computation
constrained_z_implementer<-function(z,B,G,N){
  #Input:
  #z: matrix of probability score based on estimated parameters (N x G)
  #B: list of arrays. Positive and negative constraints. Observations within a array are positively constrained.
  #   Observations in different array are negatively constrained.
  #G: int. Number of clusters.
  #N: int. Number of observations(pixels).
  if(is.null(B)!=TRUE){
    B_g_log_exp_sum=list()
    B_g_exp_sum = list()
    iter = 1
    for(b in B){
      temp=c()
      temp_2=c()
      for(g in 1:G){
        maximum = max(z[b,g])
        summed = log(sum(exp(z[b,g]-maximum)))+maximum
        temp = c(temp,summed)
        temp_2 = c(temp_2,exp(summed))
      }
      B_g_log_exp_sum[[iter]] = temp
      B_g_exp_sum[[iter]] = temp_2
      iter=iter+1
    }
    z_temp=matrix(0,N,G)
    g_list = 1:G
    b_list = 1:length(B)
    for(b in b_list){
      z_temp = constrained_z(z_temp,B,B_g_log_exp_sum,B_g_exp_sum,b_list,b,g_list,1)
      z[B[[b]],]=z_temp[B[[b]],]
    }
  }  
  return(z)
}

constrained_z<-function(z,B,B_g_log_sum,B_g_exp_sum,b_list,b_list_val,g_list,iter=1){
  temp = 0
  for(b in b_list[b_list_val]){
    b_list_2=b_list[b_list!=b]
    for(g in g_list){
      if(iter==1){
        temp = 0
      }
      g_list_2 = g_list[g_list!=g]
      if(length(b_list_2)>0){
        received_temp = constrained_z(z,B,B_g_log_sum,B_g_exp_sum,b_list_2,1,g_list_2,iter+1)
        temp = temp + ((B_g_exp_sum[[b]][g])*received_temp)
      }
      else{
        temp = temp + B_g_exp_sum[[b]][g]
      }
      if(iter==1){
        z[B[[b]],g] = log(temp)
      }
    }
  }
  if(iter==1){
    return(z)
  }
  else{
    return(temp)
  }
}

#Implementing constrained PGMM
constrained_pgmmAECM_implementer<-function(x,G=4,Q=2,epoch = 10000, zstart = 1,zlist = c(), seed = 123456, B = NULL){
  #Input:
  #x : Dataframe(N,P). Where N is the number of observations for creating the clusters and P is the features of each observation.
  #epoch : int. Number of iterations.
  #G: int. Number of clusters.
  #Q: int. Number of factors.
  #zstart: int. kmeans starting values (1) or user defined starting values(2)
  #zlist: vector. User defined starting values. Applicable only when zstart = 2
  #B: list of lists. Positive and negative constraints. Observations within a list are positively constrained.
  #   Observations in different list are negatively constrained.
  #seed: int. seed for kmeans clustering. Applicable only when zstart = 1
  
  #Output:
  #Returns a list of posterior probability(probability_score), predicted cluster (maximum a posteriori values), Entropy weights, time taken to fit, log likelihood and the initial_cluster.
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

#Creating constraints list
constrained_pixels <- function(n_cluster = 3, n_images = 4, images_size = list(),constrained_image = c(),constrained_image_x_y = list()){
  #Input:
  #n_cluster : int. Number of clusters G
  #n_images : int. Number of images.
  #image_size: list. Size of images(n_x,n_y). n_x - number of pixels along x axis, n_y - number of pixels along y axis. 
			   #Example : list(c(143,230),c(152,215),c(136,229),c(138,181))
  #constrained_image: array. Images where constraint information available. 1 if constraint information present, else 0. 
					  #Example : if n_images = 4, constrained_image = c(1,1,0,1) means images 1, 2, and 4 has known constraints.
  #constrained_image_x_y: list of lists. Each list (corresponds to each image) has multiple array specifying the blocks of pixels which belongs to a separate cluster.
					      #Example : list(list(c(1,1,35,150,230),c(1,110,143,1,50),c(2,90,96,55,70)),list(c(1,1,35,1,40),c(1,1,35,1,40),c(3,24,32,108,118)),list(),list())
						  #Here Image 1 has three blocks of pixels and Image 2 has three blocks of pixels with constraint information. Image 3 and 4 has no constraint information (empty list).
						  #Here c(1,1,35,150,230) refers to (g,x_start,x_end,y_start,y_end). The pixels within block (x_start,x_end,y_start,y_end) belong to cluster 1.
  #Output:
  #Returns a constraint list
  constraints_list = list()
  for(c in c(1:n_cluster)){
    constraints_list[[c]] = list()
  }
  pixels_n = 0
  for(i in c(1:n_images)){
    if(i!=1){
      pixels_n = pixels_n + (image_size[[i-1]][1]*images_size[[i-1]][2])
    }
    if(constrained_image[[i]]==1){
      for(j in c(1:n_cluster)){
        for(k in constrained_image_x_y[[i]]){
          if(j == k[1]){
            for(x in c(k[2]:k[3])){
              temp_val = ((x-1)*image_size[[i]][2])+pixels_n
              for(y in c(k[4]:k[5])){
                pixel_val = temp_val+y
                constraints_list[[j]] = append(constraints_list[[j]],pixel_val)
              }
            }
          }
        }
      }      
    }
  }
  for(c in c(1:n_cluster)){
    constraints_list[[c]] = unlist(constraints_list[[c]])
  }
  return(constraints_list)
}
#Testing the constrained_pixels function.
n_cluster = 3
n_images = 4
image_size = list(c(143,230),c(152,215),c(136,229),c(138,181))
constrained_image = c(1,1,0,0)
constrained_image_x_y= list(list(c(1,1,35,150,230),c(1,110,143,1,50),c(2,90,96,55,70),c(2,90,98,134,150),c(2,57,63,108,120)),list(c(1,1,35,1,40),c(1,1,35,1,40),c(3,24,32,108,118),c(3,20,28,135,145),c(3,124,131,128,143),c(3,56,66,144,150)),list(),list())
constraints_list = constrained_pixels(n_cluster,n_images,image_size,constrained_image,constrained_image_x_y)