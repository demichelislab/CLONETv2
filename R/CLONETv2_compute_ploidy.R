
findOptimumPloidy <- function(observedLog2Beta,maxCN,minAdm=0.1,maxAdm=0.7,stepCN=0.01,stepAdm = 0.05,library_type="WES",ncores){
  confs <- expand.grid(adm=seq(minAdm,maxAdm,stepAdm),cn=seq(1,maxCN,stepCN), stringsAsFactors=F)

  computeObservedRMSE.ind <- function(i, confs,observedLog2Beta,maxCN){
    Conf <-  confs[i,]
    Conf$RMSE <- computeObservedRMSE(Conf=c(confs$adm[i],confs$cn[i]), observedLog2Beta,maxCN)
    return(Conf    )
  }

  confs.list <-mclapply(seq(1,nrow(confs),1),FUN=computeObservedRMSE.ind,confs=confs,observedLog2Beta=observedLog2Beta,maxCN=maxCN,mc.preschedule = T,mc.cores=ncores)
  confs <- fromListToDF(confs.list)

  if (library_type == "WES"){
	  # select confs with minimum RMSE
	  candidate.confs <- confs[which(confs$RMSE <= stats::quantile(confs$RMSE,probs=0.05)),]
	  # take only cases with minimum adm
	  if (nrow(candidate.confs) > 1){
	    discontinuity_points <- which(candidate.confs$cn[seq(1,nrow(candidate.confs)-1,1)] + 2*stepCN < candidate.confs$cn[seq(2,nrow(candidate.confs),1)])
	    if (length(discontinuity_points) > 0 ){
	      candidate.confs <- candidate.confs[seq(1,  min(discontinuity_points,na.rm = T),1),]
	    }
	  }

	  valToRet <- candidate.confs$cn[which.min(candidate.confs$RMSE)]
	  return(valToRet)

  }else if (library_type == "WGS"){
		# as you have more segments in WS you can use dscan algorithm to find clusters of valid points
  	observedLog2Beta_clusters <- hdbscan(observedLog2Beta, minPts = 5)
  	#plot(x = observedLog2Beta$log2,y = observedLog2Beta$beta, pch=20, xlim=c(-1.2,1), ylim=c(0,1), col=cl$cluster+1)

  	## compute median of each cluster
  	found_clusters <- setdiff(unique(observedLog2Beta_clusters$cluster),0)
  	cl_stats <- data.frame(row.names = found_clusters)
  	cl_stats$cluster_id <- found_clusters
  	cl_stats$log2_median <- NA
  	cl_stats$log2_mean <- NA
  	for(cl in found_clusters ){
  		cl_stats[as.character(cl),"log2_median"] <- round(stats::median(observedLog2Beta$log2[which(observedLog2Beta_clusters$cluster == cl)]),2)
  		#cl_stats[as.character(cl),"log2_mean"] <- round(mean(observedLog2Beta$log2[which(observedLog2Beta_clusters$cluster == cl)]),2)
  	}

  	## use leftmost cluster median
  	valToRet <- round(2*2^(-1*min(cl_stats$log2_median)),2)
	  return(valToRet)
  }else{ stop("Value ",library_type," for parameter library_type not supported")}
}

computeObservedRMSE <- function(Conf,observedLog2Beta,maxCN ){
  nPoints <- maxCN + 1
  adm.local <- Conf[1]
  ploidy <- Conf[2]
  #at("adm=",adm.local," ploid=",ploidy,"\n")
  randomData <- data.frame(row.names=seq(1,nPoints,1))
  randomData$amd.local <- rep(adm.local,nPoints)
  randomData$ploidy <- rep(ploidy,nPoints)
  randomData$cnA <- c(0,seq(1,maxCN,1))
  randomData$cnB  <- c(0,seq(1,maxCN,1))
  syntheticSegs.list <- lapply(seq(1,nrow(randomData),1),perturbDataDF,randomData)
  syntheticSegs<-cbind(randomData,do.call(rbind,syntheticSegs.list))
  #plot(syntheticSegs$log2Val,syntheticSegs$betaVal.apparent,pch=20,ylim=c(0,1))
  #points(syntheticSegs$log2Val,syntheticSegs$betaVal.apparent,pch=20,col="orange2")
  expectedLog2Beta <- syntheticSegs[,c("log2Val","betaVal.apparent")]
  colnames(expectedLog2Beta ) <- c("log2Val","betaVal")
  return(computeRMSE(observedLog2Beta=observedLog2Beta,expectedLog2Beta=expectedLog2Beta))
}

perturbData <- function(cnA,cnB,adm.local,MeanPloidy,coeffVariation=0.2){
  plT <- cnA+cnB
  #pltSd <- plT * coeffVariation
  #plT <- max(0,plT + rnorm(n=1,sd=pltSd))
  plN<-2

  betaVal <- ( adm.local * plN )  / ( plT - adm.local * ( plT - plN ) )

  n.cells.mono.apparent <-  max(cnA,cnB) - min(cnA,cnB) #(cnA + cnB) - min(cnA,cnB)
  n.cells.bi.apparent <- 2*min(cnA,cnB)

  adm.local.apparent <- ( adm.local + (1-adm.local)*(n.cells.bi.apparent)) / ( adm.local + (1-adm.local)*(n.cells.bi.apparent) + (1-adm.local)*(n.cells.mono.apparent)  )
  betaVal.apparent <- ( adm.local.apparent * 2 )  / ( 1 + adm.local.apparent )

  #adm.local.apparent <- adm.local + (1-adm.local)*(min(cnA,cnB)/max(cnA,cnB))
  #betaVal.apparent <- round(( adm.local.apparent * 2 )  / ( 1 + adm.local.apparent ) ,3 )
  #betaVal.apparent <- if (cnA==cnB){ 1 }else{ betaVal + (1-betaVal) * ( 2*min(cnA,cnB)/(cnA+cnB)) }
  #adm.local.apparent <- betaVal.apparent / (2-betaVal.apparent)

  log2Val <- log2(( ( betaVal * plN + ( 1 - betaVal ) * plT ) / plN) / (MeanPloidy / plN))


  log2Val.corr <- log2(plT / plN)

  ## add noise to betaBav
  #betaValSD

  #betaVal.apparent <- if (cnA==cnB){ 1 }else{ betaVal + ( 2*min(cnA,cnB)/(cnA+cnB))}
  #((MeanPloidy/2)))

  outDF <- data.frame(row.names=1)
  outDF$cnA <- cnA
  outDF$cnB <- cnB
  outDF$ploidy <- MeanPloidy
  outDF$adm.local <- adm.local
  outDF$adm.local.apparent <- adm.local.apparent
  outDF$betaVal <- betaVal
  outDF$betaVal.apparent <- betaVal.apparent
  outDF$log2Val <- log2Val
  outDF$log2Val.corr <- log2Val.corr
  #outDF$log2Val.Segmentation <- log2(plT/MeanPloidy)
  return(outDF)
}

perturbDataDF <- function(i,randomData){
  return(perturbData(cnA=randomData$cnA[i],cnB=randomData$cnB[i],adm.local=randomData$amd.local[i],MeanPloidy=randomData$ploidy[i]))
}

computeRMSE <- function(observedLog2Beta,expectedLog2Beta,nCores=1){
  #for each observed compute the minimum distance from all expected
  distances <- unlist(mclapply(seq(1,nrow(observedLog2Beta),1),
                               computeMinDistanceFromExpected,
                               observedLog2Beta=observedLog2Beta,
                               expectedLog2Beta=expectedLog2Beta,mc.preschedule = T,mc.cores=nCores))
  RMSE <-  sqrt(sum(distances^2)/length(distances))
  return(RMSE)
}

computeMinDistanceFromExpected <- function(i,observedLog2Beta,expectedLog2Beta){
  thisLog2 <- observedLog2Beta$log2[i]
  return(min(abs(expectedLog2Beta$log2Val - thisLog2),na.rm=T))
}

#for testing
getPossibleCnComb <- function(maxCN,step=1){
  out<-c()
  for(i in seq(0,maxCN,step)){
    for(j in seq(0,i,step)){
      out <- rbind(out,c(j,(i-j)))
    }
  }
  return(out)
}

#' Function to compute ploidy from a beta table.
#'
#' This function takes the beta table of a tumor sample and returns its ploidy.
#'
#' @param beta_table data.frame formatted as the output of function
#' \code{\link[CLONETv2:compute_beta_table]{compute_beta_table}}
#' @param max_homo_dels_fraction estimated maximum proportion of genomic
#'   segments corresponding to an homozygous deletion (default=0.01)
#' @param beta_limit_for_neutral_reads minimum beta value of a segment valid for
#'   computing ploidy (default=0.90)
#' @param min_coverage minimum coverage of a segment valid for computing ploidy
#'   (default=20)
#' @param min_required_snps minimum number of informative snps in  a segment
#'   valid for computing ploidy (default=10)
#' @param library_type WES, WGS (default=WES)
#' @param n_digits number of digits in the output table (default=3)
#' @param n_cores number of available cores for computation (default=1)
#' @return A data.frame with two columns: sample that corresponds to column
#'   sample of the input beta_table, and ploidy computed
#' @examples
#' \donttest{
#' ## Compute ploidy table with default parameters
#' pl_table_toy <- compute_ploidy(bt_toy)
#' }
#' @author Davide Prandi
#' @export
#' @md
compute_ploidy <- function(beta_table,
													 max_homo_dels_fraction = 0.01,
													 beta_limit_for_neutral_reads = 0.90,
													 min_coverage=20,
													 min_required_snps=10,
													 library_type="WES",
													 n_digits=3,
													 n_cores=1){
  available_library_types <- c("WES","WGS")
	if (!library_type %in% available_library_types){stop("Parameter library_type must be one of ",paste(available_library_types,collapse = ", "))}

	sample_id <- unique(beta_table$sample)

	if (length(sample_id) != 1){stop("Beta table must contains only one sample")}

	## remove potential homo dels
	if (library_type == "WES"){
		beta_table <- beta_table[which(beta_table$log2 > stats::quantile(beta_table$log2,probs=max_homo_dels_fraction)),]
	}else if (library_type == "WGS"){
		beta_table <- beta_table[which(beta_table$log2 > stats::quantile(beta_table$log2[which(!is.na(beta_table$beta))],probs=max_homo_dels_fraction)),]
	}else{
		stop("Value ",library_type," for parameter library_type not supported")
	}

  ## only putative copy number neutral segments
  beta_table <- beta_table[which(beta_table$beta >=  beta_limit_for_neutral_reads &
  															 beta_table$nsnp >= min_required_snps &
  															 beta_table$cov >= min_coverage	),]

  if (nrow(beta_table) > 0 ){
    ## prepare data
    observedLog2Beta <- data.frame(row.names=seq(1,nrow(beta_table),1))
    observedLog2Beta$log2 <-  beta_table$log2
    observedLog2Beta$beta <-  beta_table$beta
    maxCN <- max(3,ceiling( 2/(2^min(observedLog2Beta$log2)) * 10)/10)
    plComputed <- findOptimumPloidy(observedLog2Beta,maxCN=maxCN,library_type = library_type ,ncores=n_cores  )

  }else{
    plComputed <- NA
  }
  dfOut <- data.frame(row.names=1,stringsAsFactors=F)
  dfOut$sample <- sample_id
  dfOut$ploidy <- round(plComputed, n_digits)

  return(dfOut)


}



