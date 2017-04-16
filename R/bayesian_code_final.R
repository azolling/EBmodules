require(evd)

#' trans function
#' 
#' this function is the inverse of the logit. Useful for estimating p in the likelihood optimization
#' @param q parameter that we want to apply inverse logit on 
#' @keywords logit
#' @export
#' @examples
#' trans(rnorm(100,0,1))

trans <- function(q){return(exp(q)/(1+exp(q)))}


#' logit function
#' 
#' This function computes the logit of a quantity
#' @param p a vector with values between 0 and 1 which we want to take logit.
#' @keywords logit
#' @export
#' @examples
#' logit(runif(100,0,1))

logit <- function(p){return(log(p/(1-p)))}



#' fisher function
#' 
#' The fisher transformation of a vector of correlations
#' @param r is a vector of correlations
#' @export
#' @return a vector of normaly distributed data

fisher <- function(r){
	z <- 0.5*log((1+r)/(1-r))
	if(length(which(abs(r)==1))!=0){z[which(abs(r)==1)] <- 50}
	return(z)
	}


#' fisher_inv function
#' 
#' Inverse of the fisher transformation, 
#' @param z a normal vector 
#' @export
#' @return a vector of correlations.

fisher_inv <- function(z){
	return((exp(2*z)-1)/(exp(2*z)+1))
}

#' likelihood_logit function
#' 
#' This function computes the likelihood for the empirical bayesian model. It uses the logit transformation to estimate p (simpler) 
#' @param param is a vector of length L+2, where the L first vectors are sigma (the standard deviations for each of the different type1 studies), 
#' q=logit(p)=param[L+1], where p is the probability that theta is 0 (trans(q)=p) and 
#' tau=param[L+2], the sd of theta.
#' @param Z a vector of z values (see bayesian_modules)
#' @return the value of the likelihood
#' @export

likelihood_logit <- function(param,Z){
	L <- length(Z)
	b <- ((1/(param[L+2]^2))+sum(1/(param[1:L]^2)))
	w2 <- 1
	w1_part <- 0
	for(l in 1:L){
		w2 <- w2*dnorm(Z[[l]],mean=0,sd=sqrt(param[l]^2))	
		w1_part <- w1_part+(Z[[l]]/(param[l]^2))
	}
	w1 <- w2*(trans(param[L+1])/(sqrt(b*param[L+2]^2)))*exp((1/(2*b))*(w1_part^2))
	f <- -sum(log(w1+((1-trans(param[L+1]))*w2)),na.rm=TRUE)
	return(f)	
}



#' median_theta_function
#' 
#' Function to compute the posterior median of theta
#' @param mean the mean of theta 
#' @param sd the standard deviation of theta
#' @param p_z the posterior probability of being 0
#' @return the posterior median of theta
#' @export

median_theta_function <- function(mean,sd,p_z){
	med_theta <- rep(0,length(mean))
	mean_select <- mean[which(p_z>0.5)]
	p_z_select <- p_z[which(p_z>0.5)]
	bornes<-cbind(p_z_select*pnorm(0,mean_select,sd),p_z_select*pnorm(0,mean_select,sd)+(1-p_z_select))
	select1 <- which(bornes[,1]>0.5)
	med_theta[which(p_z>0.5)][select1] <- (sd*qnorm(1/(2*p_z_select[select1])))+mean_select[select1]
	select2 <- which(bornes[,2]<0.5)
	med_theta[which(p_z>0.5)][select2] <- (sd*qnorm(((2*p_z_select[select2])-1)/(2*p_z_select[select2])))+mean_select[select2]
	return(med_theta)
	}
	

#' bayesian modules function
#' 
#' This function computes the common correlation matrix from a set of correlation matrices
#' @param R is a list of correlation matrices of length L (the number of type 1 studies). The correlation matrices are of size pxp and are all of the same size. 
#' @param select is the number of items we want to select in order to estimate the parameters, default value is set to 10^5
#' @return  param a list of the estimated parameters
#' @return med_theta the posterior median of theta
#' @return mean_theta the posterior mean of theta
#' @return sd_theta the posterior standard deviation of theta
#' @return post_theta the posterior probability that theta is 0
#' @return gene_names the names of the genes included in the analysis
#' @export

bayesian_modules <- function(R,select=100000){
	gene_names <- rownames(R[[1]])
	## Fisher transformation of the correlations
	L <- length(R)
	R_vector <- lapply(R,function(x){x[upper.tri(x,diag=FALSE)]}) 
	# transform correlations matrices into vectors. we don't need the diag as it is 1 so not interesting and we know the matrix is symmetric
	Z <- lapply(R_vector,fisher)	
	Z_mat <- matrix(unlist(Z),nrow=L,byrow=TRUE)	
	G <- dim(Z_mat)[2]
	# estimation of the parameters, param=(sigma1,...sigmaL,q=logit(p),tau)
	# initialization of the param
	param <- c(rep(0.5,L),-1.2,0.5)
	#selection of some pairwise correlation of Z to estimate the parameters
	s <- sample(1:G,size=min(select,G))
	numb_it <- as.integer(100) # number of iterations for optim
	while(!is.finite(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(param,Z=lapply(Z,"[",s)))){
		s <- sample(1:G,size=min(select,G))
		}
	system.time(MLE <- optim(param,likelihood_logit,Z=lapply(Z,"[",s)))
	while(MLE$convergence!=0){
		if(MLE$convergence==1){
			numb_it <- as.integer(numb_it*10)
			if(!is.na(numb_it)){
			 while(!is.finite(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(param,Z=lapply(Z,"[",s)))){
			    s <- sample(1:G,size=min(select,G))
			 }
			  if(!is.finite(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(MLE$par,Z=lapply(Z,"[",s)))){
			system.time(MLE <- optim(param,likelihood_logit,Z=lapply(Z,"[",s),control=list(maxit=as.integer(numb_it))))
			  }
			  else {system.time(MLE <- optim(MLE$par,likelihood_logit,Z=lapply(Z,"[",s),control=list(maxit=as.integer(numb_it))))}
			} else {MLE$convergence <- 10}
		}
		if(MLE$convergence==10){
			select <- select*10	
			numb_it <- min(10^8,as.integer(numb_it*10))
			s <- sample(1:G,size=min(select,G))
			while(!is.finite(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(param,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(param,Z=lapply(Z,"[",s)))){
			s <- sample(1:G,size=min(select,G))
			}
			if(numb_it==10^8) {
			  cat("SANN, ")
			  numb_it <- 1000
			  if(!is.finite(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(MLE$par,Z=lapply(Z,"[",s)))){
			  MLE <- optim(param,likelihood_logit,Z=lapply(Z,"[",s),method="SANN",control=list(maxit=as.integer(numb_it)))}
			  else {MLE <- optim(MLE$par,likelihood_logit,Z=lapply(Z,"[",s),method="SANN",control=list(maxit=as.integer(numb_it)))} 
			}
      else{
        if(!is.finite(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.na(likelihood_logit(MLE$par,Z=lapply(Z,"[",s))) || is.nan(likelihood_logit(MLE$par,Z=lapply(Z,"[",s)))){
        system.time(MLE <- optim(param,likelihood_logit,Z=lapply(Z,"[",s),control=list(maxit=as.integer(numb_it))))}
        else system.time(MLE <- optim(MLE$par,likelihood_logit,Z=lapply(Z,"[",s),control=list(maxit=as.integer(numb_it))))
      }  
		}
	}
	p <- trans(MLE$par[L+1])
	sigma <- abs(MLE$par[1:L])
	tau <- abs(MLE$par[L+2])
	param <- c(sigma,p,tau)
	w2_log <- 0
	w1_part <- 0
	for(l in 1:L){
		w2_log <- w2_log+log(dnorm(Z[[l]],mean=0,sd=sqrt(sigma[l]^2)))
		w1_part <- w1_part+Z[[l]]/(sigma[l]^2)
	}
	w2 <- exp(w2_log)
	b <- min((1/(tau^2))+sum(1/(sigma^2)),10^300)
	w1 <- w2*(p/(sqrt(b*tau^2)))*exp((1/(2*b))*(w1_part^2))
	w2 <- (1-p)*w2
	p_z <- w1/(w1+w2)
	p_z[which(w2==0)] <- 1
	mean_theta <- w1_part/b
	mean_theta[which(mean_theta>10^300)] <- 10^300
	med_theta <- median_theta_function(mean_theta,1/sqrt(b),p_z)
	
	results <- list(param,med_theta,mean_theta,1/sqrt(b),p_z,gene_names)
	names(results) <- c("param","med_theta","mean_theta","sd_theta","post_theta","gene_names")
	return(results)		
}


#' lik_mixgpd_onexi function
#' 
#' This function computes the likelihood  of a mixture of two GPD distributions with a common xi parameter
#' @param sigma1, parameter of the first GPD distribution
#' @param sigma2, parameter of the second GPD distribution
#' @param xi common shape parameter for both distributions
#' @param p proportion of mixture of each distribution
#' @return value of the likelihood
#' @export

lik_mixgpd_onexi <- function(param,x){
	sigma1 <- sqrt(param[1]^2)
	sigma2 <- sqrt(param[2]^2)
	xi <- param[3]
	p <- trans(param[4])
	return(-sum(log(p*dgpd(x,scale=sigma1,shape=xi)+(1-p)*dgpd(x,scale=sigma2,shape=xi))))
}



#' scattered_extremes_onexi function
#' 
#' This function is used to discriminate between scattered and non scattered genes. It fits a mixture of GPD distributions with a common xi.
#' @details here we consider that xi1=xi2 for a better estimation of this parameter
#' @param T statistic whose tail/most extremes values should be the scattered genes for the most part
#' @param quantile_T is the threshold that we choose for the tail
#' @param hessian bolean to indicate whether we want to compute the standard deviation of the parameters
#' @return scattered gives the name of the scattered genes, the value of the parameters and if hessian is TRUE, their standard error
#' @export

scattered_extremes_onexi <- function(T,quantile_T=0.9,hessian=FALSE){
	u <- quantile(T,p=quantile_T)
	x <- T[which(T>u)]
	param <- c(0.5,1,-0.2,logit(0.8))
	while(is.na(lik_mixgpd_onexi(param=param,x=x-u))||is.nan(lik_mixgpd_onexi(param=param,x=x-u))||!is.finite(lik_mixgpd_onexi(param=param,x=x-u))){
	  cat("change param")
	  param<-c(runif(2,0.1,1),-0.2,logit(runif(1,0.1,0.9)))
	  }
	
	MLE <- optim(param,lik_mixgpd_onexi,x=x-u,control=list(maxit=10000),hessian=hessian)
	while(MLE$convergence!=0){
		if(MLE$convergence==1){
		#cat("try again \n")
		MLE <- optim(MLE$par,lik_mixgpd_onexi,x=x-u,control=list(maxit=10000),hessian=hessian)
		}
		if(MLE$convergence==10){
			#cat("SANN \n")
			MLE <- optim(MLE$par,lik_mixgpd_onexi,x=x-u,control=list(maxit=10000),method="SANN",hessian=hessian)			
		}
		}
		sigma1 <- abs(MLE$par[1])
		sigma2 <- abs(MLE$par[2])
		xi <- MLE$par[3]
		p <- trans(MLE$par[4])
		if(hessian) {se <- sqrt(diag(solve(MLE$hessian)))
			names(se) <- c("se_sigma1","se_sigma2","se_xi","se_q")}else se <- NA
			param <- c(sigma1,sigma2,xi,p)
		
		names(param) <- c("sigma1","sigma2","xi","p")
	post2 <- (1-p)*dgpd(x-u,scale=sigma2,shape=xi)/((1-p)*dgpd(x-u,scale=sigma2,shape=xi)+p*dgpd(x-u,scale=sigma1,shape=xi))
	post1 <- (p*dgpd(x-u,scale=sigma1,shape=xi))/((1-p)*dgpd(x-u,scale=sigma2,shape=xi)+p*dgpd(x-u,scale=sigma1,shape=xi))
	max_x <- which(x==max(x))
	if(post1[max_x]>post2[max_x]) post <- post1 else post <- post2
	results <- list(names(which(post>0.5)),param,se)
	names(results) <- c("scattered","param","se")
	return(results)
}

#' scattered_by_kmeans function
#' 
#' Function to find scattered genes by kmeans with k=2
#' @param T the statistic from which we want to discriminate between scattered and non scattered genes
#' @return scattered the names of the scattered genes
#' @export

scattered_by_kmeans <- function(T){
  results_kmeans <- kmeans(as.matrix(T),2)
  cluster_max <- results_kmeans$cluster[which(T==max(T))]
  scattered <- which(results_kmeans$cluster==cluster_max)
  return(scattered)
}


#' EB_module function
#' 
#' This function is the main function of the package. It computes the common gene modules for a collection of correlation matrices
#' @param R is a list of lenght L1 (number of studies) of correlation matrices of the same size
#' @param gap whether to use the gap statistic to determine the optimal number of clusters
#' @param k if gap is TRUE then k=NULL otherwise, k is an integer indicating the number of desired modules (scattered genes do not count as a module)
#' @param scattered: bolean whether the function should look for scattered genes before looking for modules.
#' @return common_r the common correlation matrix estimated through the empirical bayesian model
#' @return modules a list of modules with the scattered genes grouped inside a module (called -1)
#' @return bestk the optimal number of modules according to a modified version of the GAP statistic (if gap is set to TRUE)
#' @export

EB_module <- function(R,gap=FALSE,k=500,scattered=FALSE,method_scattered=c("GPD","kmeans"),Kmin=NULL,Kmax=NULL,Nperm=100,dicho_f=2){
	#start by computing the common correlation matrix for all studies
	n_gene <- dim(R[[1]])[1]
	names_gene <- rownames(R[[1]])
	theta_tilde <- bayesian_modules(R,select=100000)$med_theta
	R_tilde <- fisher_inv(theta_tilde)
	common_R <- theta_matrix <- matrix(nrow=n_gene,ncol=n_gene)
	rownames(common_R) <- colnames(common_R) <- rownames(theta_matrix) <- colnames(theta_matrix) <- names_gene
	common_R[upper.tri(common_R,diag=FALSE)] <- R_tilde
	diag(common_R) <- rep(1,n_gene)
	common_R[lower.tri(common_R,diag=FALSE)] <- t(common_R)[lower.tri(common_R,diag=FALSE)]
	
	#identify scattered genes
	theta_matrix[upper.tri(theta_matrix,diag=FALSE)] <- theta_tilde
	theta_matrix[lower.tri(theta_matrix,diag=FALSE)] <- t(theta_matrix)[lower.tri(theta_matrix,diag=FALSE)]
	#treat as NA genes that have 0 correlations
	
	extremes_stat <- -log(apply(abs(theta_matrix),1,var,na.rm=TRUE))
	names(extremes_stat) <- rownames(theta_matrix)
	if(scattered){
	  if(method_scattered=="GPD") {scattered_extr <- scattered_extremes_onexi(extremes_stat,quantile_T=0.9)$scattered}
	  else if(method_scattered=="kmeans") {scattered_extr <- scattered_by_kmeans(extremes_stat)}
	  else {stop("invalid method for scattered genes")}
	}
	else scattered_extr <- NULL
	
	#perform clustering to detect modules
	#remove scattered genes from the matrix to cluster
	if(length(scattered_extr)!=0){
		R_to_cluster <- common_R[setdiff(rownames(common_R),scattered_extr),setdiff(rownames(common_R),scattered_extr)]
		}
		else R_to_cluster <- common_R
		
	#cluster
	tree <- hclust(as.dist(1-R_to_cluster),method="ward.D")
	if(!gap){
		if((k==0)||(length(k)==0)){stop("ERROR the number of clusters k is null")}
		else {cluster <- cutree(tree,k)
		      bestk <- k}
	}
	if(gap){
	  bestk <- GAP_dicho(R=R_to_cluster,Nperm=Nperm,Kmin=Kmin,Kmax=Kmax,dicho_f=dicho_f)
	  cluster <- cutree(tree,bestk)
	}
	#prepare cluster to return (ie modules+scattered genes)
	cluster_scattered <- c(cluster,rep(-1,length(scattered_extr)))
	names(cluster_scattered) <-c(names(cluster),scattered_extr) 
	cluster_scattered_list <- split(names(cluster_scattered),cluster_scattered)
	result <- list(common_R,cluster_scattered_list,bestk)
	names(result) <- c("common_R","modules","bestk")
	return(result)
}




#' Wk function
#' 
#' function to compute the average sum of pairwise distances across clusters as defined in the paper from Tibshirani
#' @param R is a correlation matrix
#' @param cluster is the results of a clustering procedure, a vector of length number of genes indicating the cluster to which the corresponding gene belongs.
#' @return Wk the average sum of pairwise distances across clusters
#' @export

Wk <- function(R,cluster){
	p <- length(cluster)
	genes <- seq_len(p)
	return(0.5*sum(vapply(split(genes,cluster),function(l){
		R_select <- R[l,l]
		sum(as.dist(1-R_select)/nrow(R_select))
	},0)))	
}


#' hclust_fun function
#' 
#' function that performs hierarchical clustering and outputs the right format for our gap function
#' @param R a correlation matrix
#' @param k the number of clusters
#' @param met the method for hierarchical clustering, one of ward, average, complete or single
#' @return a vector where each gene is attributed a cluster
#' @export

hclust_fun <- function(R,k,met=c("ward","average","complete","single")){
	clust <- hclust(as.dist(1-R),method=met)
	cluster_id <- cutree(clust,k=k) 
	toreturn <- list(cluster_id)
	names(toreturn) <- "cluster"
	return(toreturn)
}

#' kmeans_fun function
#' 
#' function that performs kmeans clustering and outputs the right format for our gap function
#' @param R a correlation matrix
#' @param k the number of clusters
#' @param met default parameter set to NULL
#' @return a vector where each gene is attributed a cluster
#' @export


kmeans_fun <- function(R,k,met=NULL){
	toreturn <- list(kmeans((1-R),k)$cluster)
	names(toreturn) <- "cluster"
return(toreturn)
}



#' GAP_dicho function
#' 
#' This function is a modified version of the GAP statistic for very large datasets. 
#' @param R is a correlation matrix
#' @param Kmin gives the minimum number of clusters that we want to obtain
#' @param Kmax gives the maximum number of clusters that we want to obtain
#' @param Nperm is the number of resampling that we want to do
#' @param method can only be one of the proposed clustering algorithm. it must return a list with one of the elements named "cluster", which is a vector of length p (number of genes) and containing the id of the cluster for each gene
#' @param bestk_method the way to compute the best k as describe in the GAP package
#' @param dicho_f the factor for the dichotomy, A small value will lead to several rounds.
#' @return bestk the optimal number of clusters.
#' @export

GAP_dicho <- function(R,Kmin,Kmax,Nperm=100,method="ward.D",bestk_method="Tib",dicho_f=2){
	p <- nrow(R)
	if(Kmax>p) stop("error, Kmax too big!")
	if(Kmin<2) stop("error, Kmin is too small")
	genes <- rownames(R)
	true_tree <- hclust(as.dist(1-R),method=method)
	delta <- round((Kmax-Kmin)/20,digits=0)
	Kset <- seq(Kmin,Kmax,delta)
	while(delta!=0){
	  #cat("\n delta=",delta,"\n")
	  logW_resamp <- matrix(0,nrow=length(Kset),ncol=Nperm)
	  logW <- ElogW <- SE.sim <- numeric(length(Kset))
	  counting <- 1
		  for(k in Kset){
		    logW[(counting)] <- log(Wk(R,cutree(true_tree,k=k)))
		    counting <- counting+1
	    }
	    for(n in 1:Nperm){
	      #if(is.element(n,seq(0,Nperm,10))) {cat("n=",n,", ")}
	    	R_vec <- R[upper.tri(R,diag=FALSE)]
		    R_resamp_vec <- sample(R_vec,length(R_vec))
		    R_resamp <- matrix(nrow=p,ncol=p)
		    diag(R_resamp) <- rep(1,p)
		    R_resamp[upper.tri(R_resamp,diag=FALSE)] <- R_resamp_vec
		    R_resamp[lower.tri(R_resamp,diag=FALSE)] <- t(R_resamp)[lower.tri(R_resamp,diag=FALSE)]
		    tree_resamp <- hclust(as.dist(1-R_resamp),method=method)
		    counting <- 1
		      for(k in Kset){
			      if(k==1) cluster_resamp <- rep.int(1,p)
			      else cluster_resamp <- cutree(tree_resamp,k=k)
			      logW_resamp[(counting),n] <- log(Wk(R_resamp,cluster_resamp))
			      counting <- counting+1
		      }
		  }
	  ElogW <- rowMeans(logW_resamp)
	  SE.sim <- sqrt((1+(1/Nperm))*apply(logW_resamp,1,var))
	  gap <- ElogW-logW
	  if(bestk_method=="Tib"){
	    bestk <- Kmin
	    counting <- 1
	    for(k in Kset[-1]){
		    if(gap[counting]>=gap[counting+1]-SE.sim[counting+1]) {
			    bestk <- Kset[counting]
			    break;	
		    }
	      counting <- counting+1
	    }
	  }
	  else if (bestk_method=="first_max"){bestk <- K[which(diff(gap,differences=1)<=0)[1]]}
	  else if (bestk_method=="global_max"){bestk <- K[which.max(gap)]}
	  if(delta==1) {delta <- 0} else delta <- ceiling(delta/dicho_f) 
	  Kset <- sort(unique(c(bestk-seq(0,10,1)*delta,bestk+seq(0,10,1)*delta)))
	  Kset <- Kset[Kset>Kmin & Kset<Kmax]
}
	return(bestk)
}



#' RandIndex function
#' 
#' This function computes the rand index for two clusters clusX and clusY.
#' @param clusX and clusY lists of length "number of clusters" and trying to cluster the same number of genes
#' @return the value of the rand index
#' @export

RandIndex <- function(clusX,clusY){
	#need first to construct the matrix of agreements A
	if(sum(unlist(lapply(clusX,length)))!=sum(unlist(lapply(clusY,length))))stop("error both clusters do not cluster the same objects")
	p <- sum(unlist(lapply(clusX,length)))
	A <- matrix(nrow=length(clusX),ncol=length(clusY))
	gene <- seq.int(p)
	a <- n <-  0
	for(i in 1:length(clusX)){
		for(j in 1:length(clusY)){
			A[i,j] <- length(intersect(clusX[[i]],clusY[[j]]))
			a <- a+choose(A[i,j],2)	
			n <- n+A[i,j]#=p??
		}
	}
	sumrow <- apply(A,1,sum)
	b <- sum(choose(sumrow,2))
	sumcol <- apply(A,2,sum)
	c <- sum(choose(sumcol,2))
	ARI <- (a-(b*c/choose(n,2)))/((0.5*(b+c))-(b*c/choose(n,2)))
	return(ARI)
}

#' consensus_cluster function
#' 
#' This function performs consensus clustering.
#' @param R a correlation matrix
#' @param kmin the minimum number of clusters
#' @param kmax the maximum number of clusters
#' @param clusterAlg the clustering algorithm, one of the following kmeans, hclust_ward,hclust_average,hclust_complet,hclust_single
#' @param B the number of resampling steps
#' @param prop_resamp the proportion of the data that we want to resample
#' @return best k (number of clusters) chosen by this method.
#' @export


consensus_cluster <- function(R,kmin=2,kmax,clusterAlg="hclust_ward",B=10, prop_resamp=0.8){
	if(!isSymmetric(R)) stop("error, R is not a symmetric matrix")
	if(is.null(rownames(R))) rownames(R) <- colnames(R) <- 1:dim(R)[1]
	M <- I_mat <-  vector("list",kmax)
	A <- rep(0,kmax)
	for(k in 2:kmax){
		M_clust <- vector("list",B)
		I_mat <- vector("list",B)
		for(b in 1:B){
			gene_resamp <- sample(1:dim(R)[1],round(prop_resamp*dim(R)[1]))
			R_resamp <- R[gene_resamp,gene_resamp]
			rownames(R_resamp) <- colnames(R_resamp) <- rownames(R[gene_resamp,gene_resamp])
			if(clusterAlg=="hclust_ward") clust <- hclust_fun(R_resamp,k=k,met="ward")$cluster
			if(clusterAlg=="hclust_average") clust <- hclust_fun(R_resamp,k=k,met="average")$cluster
			if(clusterAlg=="hclust_complet") clust <- hclust_fun(R_resamp,k=k,met="complete")$cluster
			if(clusterAlg=="hclust_single") clust <- hclust_fun(R_resamp,k=k,met="single")$cluster
			if(clusterAlg=="kmeans") clust <- kmeans_fun(R_resamp,k=k)$cluster
			#compute the connectivity matrix of this cluster
			if(is.null(names(clust))) names(clust) <- rownames(R_resamp)
			clust_list <- split(names(clust),clust)
			M_clust[[b]] <- I_mat[[b]] <- matrix(0,nrow=p,ncol=p)
			rownames(M_clust[[b]]) <- colnames(M_clust[[b]]) <- rownames(I_mat[[b]]) <- colnames(I_mat[[b]]) <- rownames(R)
			I_mat[[b]][names(clust),names(clust)] <- matrix(1,nrow=length(clust),ncol=length(clust))
			for(i in 1:length(clust_list)){
				M_clust[[b]][clust_list[[i]],clust_list[[i]]] <- matrix(1,nrow=length(clust_list[[i]]),ncol=length(clust_list[[i]]))
			}		
	}
	M[[k]] <- Reduce("+",M_clust)/Reduce("+",I_mat)
	M_order <-  sort(M[[k]][upper.tri(M[[k]],diag=FALSE)])
	for(i in 2:length(M_order)){
		CDF <- length(which(M_order<M_order[i]))/length(M_order)
		A[k] <- A[k]+((M_order[i]-M_order[i-1])*CDF)
		}
}
Delta <- rep(0,kmax)
Delta[2] <- A[2]
for(i in 3:(kmax-1)){
	Delta[3:kmax] <- (A[i+1]-A[i])/A[i]
	}
k_hat <- which(Delta==max(Delta))[1]
return(k_hat)
}	




#' maxmean function
#' 
#' function to summarize the information of a Z vector into a module
#' @param Z are the Z to summarize, so corresponding to one module
#' @param R is the correlation matrix corresponding to Z.
#' @return the summarized version of z for the module
#' @export

maxmean <- function(Z,R){
	toreturn <- vector("list",3)
	names(toreturn) <- c("Z","n","var")
	if(length(which(is.na(Z)))==length(Z)){toreturn$Z <- NA
		toreturn$n <- length(Z)
		toreturn$var <- NA
		return(toreturn)
		}
	else{
	s_plus <- sum(Z[which(Z>=0)],na.rm=TRUE)/length(Z)
	s_moins <- sum(Z[which(Z<0)],na.rm=TRUE)/length(Z)
	
	if(abs(s_plus)>=abs(s_moins)){
		toreturn$Z <- s_plus
		toreturn$n <- length(which(Z>=0))
		if(toreturn$n!=1){
			R_plus <- R[which(Z>=0),which(Z>=0)]
			diag(R_plus) <- rep(0,toreturn$n)
			toreturn$var <- (1/length(Z)^2)*(toreturn$n+sum(as.vector(R_plus)))
			}
		else toreturn$var <- 1
		return(toreturn)
		} 
	else {
		toreturn$Z <- s_moins
		toreturn$n <- length(which(Z<0))
		if(toreturn$n!=1){
		R_moins <- R[which(Z<0),which(Z<0)]
		diag(R_moins) <- rep(0,toreturn$n)
		toreturn$var <- (1/length(Z)^2)*(toreturn$n+sum(as.vector(R_moins)))
		}
		else toreturn$var <- 1
		return(toreturn)
		}
	}
}

#' rank_mean function
#' 
#' function to summarize ranks into modules
#' @param R the ranks to be summarized into one value
#' @export

rank_mean <- function(R){
	number_na <- length(which(is.na(R)))
	number_obs <- length(R)-length(number_na)
	R_trans <- qnorm(R/(length(R)+1))
	return(pnorm(mean(R_trans,na.rm=TRUE)*(length(R)+1)))
}
