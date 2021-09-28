## single sample adm.global estimate in 2D space
## use the cluster defined by localBeta.good to define the clonality
getAdmGlobal.2D <- function(betaT, G.seq, minVal = 0, ncores = 1 ){

  #betaT <- betaTable.this.ext.good
  #admG <- G.seq[1]
  localTransform <- function(admG,betaT,minVal){
    cnAB.list <- lapply(seq(1,nrow(betaT),1),
                        function(i){return( from_BetaLogR_to_cnAB(G=admG, LogR=betaT$log2.plCorr[i], Beta=betaT$beta[i] ))   })
    cnAB <- do.call(rbind,cnAB.list)
    dfOut <- data.frame(row.names = admG)
    dfOut$adm <- admG
    dfOut$dist <- abs(stats::median(cnAB$cnB) - minVal)
    dfOut$mindist <- abs(min(cnAB$cnB) - minVal)
    dfOut$maxdist <- abs(max(cnAB$cnB) - minVal)
    return(dfOut)
  }

  dists.list <- mclapply(G.seq, localTransform, betaT = betaT,  minVal = minVal, mc.preschedule = T, mc.cores = ncores )
  adm.dist <- fromListToDF(dists.list)
  adm.dist <- adm.dist[which(adm.dist$dist <=1),]

  # to add for compiling reports
  #plot(adm.dist$adm,adm.dist$dist,pch=20 )
  #points(adm.dist$adm,adm.dist$mindist ,pch=20,col="green3" )
  #points(adm.dist$adm,adm.dist$maxdist,pch=20,col="blue4" )

  globalAdm <- data.frame(row.names=betaT$sample[1])
  globalAdm$sample <- betaT$sample[1]
  globalAdm$adm <- adm.dist$adm[which.min(adm.dist$dist)]
  globalAdm$adm.min <- adm.dist$adm[which.min(adm.dist$mindist)]
  globalAdm$adm.max <- adm.dist$adm[which.min(adm.dist$maxdist)]
  globalAdm$n.segs <- nrow(betaT)
  globalAdm$n.SNPs <- sum(betaT$nsnp)

  return(globalAdm)

}

## Function for cnA vs cnB space management ----
from_BetaLogR_to_cnAB <- function(G, LogR, Beta ){

  cnB = ( Beta * 2^(LogR) - G) / (1 - G)
  cnA = ( (2 - Beta) * ( Beta * 2^LogR  - G)  + 2 * G * (1-Beta) ) / (  (1-G) * Beta )

  dfOut <- data.frame(cnA=cnA,cnB=cnB,stringsAsFactors=F)
  return(dfOut)

}




#####
## compute DNA admixture----

#' Function to compute DNA admixture of a tumor sample from the associatd beta
#' table and ploidy table
#'
#' This function takes a beta table and the associated ploidy table and computes
#' DNA admixture.
#'
#' @param beta_table data.frame formatted as the output of function
#'   \code{\link[CLONETv2:compute_beta_table]{compute_beta_table}}
#' @param ploidy_table data.frame formatted as the output of function
#'   \code{\link[CLONETv2:compute_ploidy]{compute_ploidy}}
#' @param min_coverage minimum coverage of a segment valid for computing ploidy
#'   (default=20)
#' @param min_required_snps minimum number of informative snps in  a segment
#'   valid for computing ploidy (default=10)
#' @param error_tb data.frame that reports for each combination of coverage and
#'   number informative SNPs the expected estimation error around beta. The
#'   data.frame error_tb must contains 3 columns: \describe{
#'   \item{mean.cov}{mean coverage} \item{n.info.snps}{number of informative
#'   SNPs} \item{adm.estimation.error}{estimated error on computed beta on a
#'   segment with coverage mean.cov and n.info.snps informative SNPs} } Package
#'   CLONETv2 have built in error_tb named error_table (default=error_table)
#' @param library_type WES, WGS (default=WES)
#' @param n_digits number of digits in the output table (default=3)
#' @param n_cores number of available cores for computation (default=1)
#' @param debug return extra columns for debugging (default=F)
#' @return A data.frame with two columns: sample that corresponds to column
#'   sample of the input beta_table, and amd that represent the fraction of
#'   estimated DNA admixture
#' @examples
#'
#' ## Compute admixture table with default parameters
#' adm_table_toy <- compute_dna_admixture(beta_table = bt_toy, ploidy_table = pl_table_toy)
#'
#' @author Davide Prandi
#' @export
#' @md
compute_dna_admixture <- function(beta_table,
															ploidy_table,
															min_required_snps=10,
															min_coverage=20,
															error_tb = error_table,
															library_type="WES",
															n_digits=3,
															n_cores=1,
															debug=F){

	available_library_types <- c("WES","WGS")
	if (!library_type %in% available_library_types){stop("Parameter library_type must be one of ",paste(available_library_types,collapse = ", "))}

	sample_id <- unique(beta_table$sample)
	if (length(sample_id) != 1){stop(paste("[",Sys.time() ,"] beta_table must contain exactly one sample\n",sep=""))}

  ## use only chr 1 - 22
  beta_table <- beta_table[which(suppressWarnings(as.numeric(gsub("chr|Chr|CHR","",beta_table$chr))) %in% seq(1,22,1 )),]

  ## add ploidy information
  ploidy_table$log2shift <- round(-log2(ploidy_table$ploidy/2),3)
  beta_table <- merge(x = beta_table, y = ploidy_table, by.x = "sample", by.y = "sample")
  beta_table$log2.plCorr <- beta_table$log2 - beta_table$log2shift

  ## add error information
  beta_table.list <- by(beta_table, INDICES = beta_table$sample, FUN = addErrorToBetaTable, errorTable = error_tb, ncores = n_cores )
  beta_table <- fromListToDF(beta_table.list)
  rm(beta_table.list)

  ## prepare output data.frame
  adm_tb <- data.frame(row.names = sample_id)
  adm_tb$sample <- sample_id
  adm_tb$adm <- NA
  adm_tb$adm.min <- NA
  adm_tb$adm.max <- NA
  adm_tb$n.segs <- NA
  adm_tb$n.SNPs <- NA

  ## filter on nsnps and coverage
  beta_table <- beta_table[which(beta_table$nsnp >= min_required_snps & beta_table$cov >= min_coverage),]


  if (library_type == "WGS"){
  	## filter beta on normal less than 0.9
  	beta_table <- beta_table[which(beta_table$n_beta > 0.9),]
  }

  ## check if empty beta_table or not defined ploidy
  if (nrow(beta_table) == 0 || is.na(beta_table$ploidy[1])){
    return(adm_tb)
  }

  ## compute allele specific copy number supposing a 100% pure sample
  beta_table <- extendBetaTableWithCopyNumbers(BetaTable = beta_table, G =0, ncores = n_cores)


  ## find maximum cnB value for the low cluster
  if (library_type %in% c("WES")){
  	localBeta <- beta_table[which(beta_table$cnB < 1),]
  }else if (library_type %in% c("WGS")){
  	localBeta <- beta_table[which(beta_table$cnB < 1 & beta_table$cnA > 0.75),]
  }else{
  	stop("Parameter library_type  ",library_type," not fully supported")
  }

  if (nrow(localBeta) < 2){
    return(adm_tb)
  }

  ## sort by increasing cnB
  localBeta <- localBeta[with(localBeta, expr = order(cnB)),]

  ## extract localBeta starting from position2 and adding an empty row at the end
  localBeta.next <- rbind(localBeta[seq(2,nrow(localBeta),1),],rep(NA,ncol(localBeta)))

  ## report on the same line segments and the next one in the genomic order
  colnames(localBeta.next) <- paste0("next_",colnames(localBeta.next))
  localBeta <- cbind(localBeta,localBeta.next)

  ## find segment such that the difference in the cnB with the next segment can't be explained with the error on beta estimate
  localBeta$max_cnB_current <- pmax(localBeta$cnB,localBeta$cnB.betaMin,localBeta$cnB.betaMax, na.rm = T)
  localBeta$min_cnB_next <- pmin(localBeta$next_cnB,localBeta$next_cnB.betaMin,localBeta$next_cnB.betaMax,na.rm = T)
  max.cnB <- suppressWarnings(localBeta$cnB[min(which(localBeta$max_cnB_current < localBeta$min_cnB_next ))])

  ## only segments in the cluster with cnB less that max.cnB
  localBeta.good <- localBeta[which(localBeta$cnB <= max.cnB ),]

  if (nrow(localBeta.good) > 0){
    G.seq <- seq(0,0.99,0.01)
    adm_tb <- getAdmGlobal.2D(betaT = localBeta.good, G.seq = G.seq, minVal = 0, ncores = n_cores)

    adm_tb$adm <- round(adm_tb$adm, n_digits)
    adm_tb$adm.min <- round(adm_tb$adm.min, n_digits)
    adm_tb$adm.max <- round(adm_tb$adm.max, n_digits)

    if (!debug){
    	col_to_save <- c("sample","adm", "adm.min", "adm.max")
    	adm_tb <- adm_tb[,col_to_save]
    }


  }
	return(adm_tb)
}
