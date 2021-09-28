
getClonalityWithError<-function(integerCN,
                                beta,
                                adm.global.thisZone,
                                adm.global.thisZone.min,
                                adm.global.thisZone.max,
                                clonalityThreshold,
                                local.error,
                                roundDec=3){
  if(integerCN == 0){

    ## to fix when beta = 0
    if (beta==0){
      betaCorr <- NA
      zone.ctm.local <- NA
      clonality <- NA
      clonality.int <- interval(NULL,NULL)
      clonality.status <- "not.analysed"
    }else{
      betaCorr <- beta
      zone.ctm.local <- beta / (2 - beta)
      G <- adm.global.thisZone
      L <- zone.ctm.local
      clonality <- setBetween0and1(1- ((G - ( L * G)) / ( L * ( 1 - G))))

      zone.ctm.local.min <- zone.ctm.local - local.error
      zone.ctm.local.max <- zone.ctm.local + local.error

      G.int <- interval(adm.global.thisZone.min,adm.global.thisZone.max)
      L.int <- interval(zone.ctm.local.min,zone.ctm.local.max  )

      clonality.int <-  1 - ((G.int - ( L.int * G.int)) / ( L.int * ( 1 - G.int)))

      clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))
    }

  }

  if (integerCN == 1){
    betaCorr <- beta
    zone.ctm.local <- beta / (2 - beta)
    clonality <- setBetween0and1((1 - zone.ctm.local) / (1 - adm.global.thisZone))

    zone.ctm.local.min <- zone.ctm.local - local.error
    zone.ctm.local.max <- zone.ctm.local + local.error

    zone.ctm.local.int <- interval(zone.ctm.local.min,zone.ctm.local.max)
    adm.global.thisZone.int <- interval(adm.global.thisZone.min,adm.global.thisZone.max)

    clonality.int <- ( 1 - zone.ctm.local.int) / (1 - adm.global.thisZone.int)

    clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))

  }
  if (integerCN == 2){
    betaCorr <- beta
    zone.ctm.local <- betaCorr
    clonality <- zone.ctm.local
    zone.ctm.local.min <- zone.ctm.local - local.error
    zone.ctm.local.max <- zone.ctm.local + local.error
    clonality.int <- interval(zone.ctm.local.min,zone.ctm.local.max)

    clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))

  }
  if (integerCN == 3){

    betaCorr  <-    1 - (3 * (1 - beta) )
    if ( betaCorr > 0){
      zone.ctm.local <- ( 3 * betaCorr ) / (( 3 - 2 ) * betaCorr + 2)
      clonality <- setBetween0and1((1 - zone.ctm.local) / (1 - adm.global.thisZone))

      zone.ctm.local.min <- zone.ctm.local - local.error
      zone.ctm.local.max <- zone.ctm.local + local.error

      zone.ctm.local.int <- interval(zone.ctm.local.min,zone.ctm.local.max)
      adm.global.thisZone.int <- interval(adm.global.thisZone.min,adm.global.thisZone.max)

      clonality.int <- ( 1 - zone.ctm.local.int) / (1 - adm.global.thisZone.int)
      clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))
    }else{
      betaCorr <- NA
      zone.ctm.local <- NA
      clonality <- NA
      clonality.int <- interval(NULL,NULL)
      clonality.status <- "not.analysed"
    }

  }
  if (integerCN == 4){
    betaCorr <- beta
    zone.ctm.local <- ( 4 * betaCorr ) / ((4-2) * betaCorr +2 )
    clonality <- setBetween0and1(zone.ctm.local  )
    zone.ctm.local.min <- zone.ctm.local - local.error
    zone.ctm.local.max <- zone.ctm.local + local.error
    clonality.int <- interval(zone.ctm.local.min,zone.ctm.local.max)
    clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))
  }

  if (integerCN == 5){
    betaCorr  <-    1 - (5 * (1 - beta) )

    if ( betaCorr > 0){
      zone.ctm.local <- ( 5 * betaCorr ) / (( 5 - 2 ) * betaCorr + 2)
      clonality <- setBetween0and1((1 - zone.ctm.local) / (1 - adm.global.thisZone))

      zone.ctm.local.min <- zone.ctm.local - local.error
      zone.ctm.local.max <- zone.ctm.local + local.error

      zone.ctm.local.int <- interval(zone.ctm.local.min,zone.ctm.local.max)
      adm.global.thisZone.int <- interval(adm.global.thisZone.min,adm.global.thisZone.max)

      clonality.int <- ( 1 - zone.ctm.local.int) / (1 - adm.global.thisZone.int)
      clonality.int <- interval(setBetween0and1(min(clonality.int)),setBetween0and1(max(clonality.int)))
    }else{
      betaCorr <- NA
      zone.ctm.local <- NA
      clonality <- NA
      clonality.int <- interval(NULL,NULL)
      clonality.status <- "not.analysed"
    }
  }

  if (integerCN > 5){
    betaCorr <- NA
    zone.ctm.local <- NA
    clonality <- NA
    clonality.int <- interval(NULL,NULL)
    clonality.status <- "not.analysed"
  }

  if (interval_is_empty(clonality.int)){
    clonality.min <- NA
    clonality.max <- NA
  }else{
    clonality <- round(clonality,roundDec)
    clonality.min <- round(min(clonality.int),roundDec)
    clonality.max <- round(max(clonality.int),roundDec)
    if ( clonality.min >= clonalityThreshold  ){ clonality.status <- "clonal"}
    if ( clonality.max <= clonalityThreshold ){ clonality.status <- "subclonal"}
    if ( clonality.min <= clonalityThreshold  &&  clonality.max >= clonalityThreshold && clonality < clonalityThreshold ) { clonality.status <- "uncertain.subclonal" }
    if ( clonality.min <= clonalityThreshold &&  clonality.max >= clonalityThreshold && clonality >= clonalityThreshold ){ clonality.status <- "uncertain.clonal"}
  }

  data<-data.frame(row.names=1,stringsAsFactors=F)
  data$betaCorr <- betaCorr
  data$zone.ctm.local <- zone.ctm.local
  data$clonality <- clonality
  data$clonality.min <- clonality.min
  data$clonality.max <- clonality.max
  data$clonality.status <- clonality.status
  return(data)
}



compute_clonality <- function(betaTable,errorTable,clonalityThreshold,betaThreshold=NULL,roundDec = 3,n_cores=1){

  #cat(betaTable$sample[1],"\n")
  ## backward compatibility
  if (is.null(betaThreshold)){betaThreshold<-clonalityThreshold}

  betaTable$integerCN <- NA
  betaTable$clonality <- NA
  betaTable$clonality.min <- NA
  betaTable$clonality.max <- NA
  betaTable$clonality.status <- "not.analysed"

  #i<-1
  getIndexClonality <- function(i){
    #cat(i,"\n")
    zone <- betaTable[i,,drop=F]
    zone$adm.global <- round(zone$adm,roundDec)
    zone$adm.global.min <- round(zone$adm.min,roundDec)
    zone$adm.global.max <- round(zone$adm.max,roundDec)

    #for(colN in c(1:ncol(zone)))

    #can be analysed
    if (is.na(zone$beta) || is.na(zone$log2.corr) || is.na( zone$adm.global.min) || is.na(zone$adm.global.max)){
      return(zone)
    }

    #determine error from nsnps and cov
    nsnps <- zone$nsnp
    cov <- zone$cov
    idxError <- which(errorTable$n.info.snps <= nsnps & errorTable$mean.cov <= cov)
    if (length(idxError) == 0 ){idxError <- 1}
    local.error <- errorTable$adm.estimation.error[max(idxError)]

    #cn given log2
    cn <- 2*2^zone$log2.corr

    #valid cn are floor(cn) and ceil(cn)
    clonCNlow <- getClonalityWithError(integerCN=floor(cn),
                                       beta=zone$beta,
                                       adm.global.thisZone=zone$adm.global,
                                       adm.global.thisZone.min=zone$adm.global.min,
                                       adm.global.thisZone.max=zone$adm.global.max,
                                       clonalityThreshold,
                                       local.error,
                                       roundDec)

    clonCNhigh <- getClonalityWithError(integerCN=ceiling(cn),
                                        beta=zone$beta,
                                        adm.global.thisZone=zone$adm.global,
                                        adm.global.thisZone.min=zone$adm.global.min,
                                        adm.global.thisZone.max=zone$adm.global.max,
                                        clonalityThreshold,
                                        local.error,
                                        roundDec)

    if (is.null(clonCNlow) | is.null(clonCNhigh) | is.na(clonCNlow$clonality) | is.na(clonCNhigh$clonality) ){
      return(zone)
    }



    if (cn >= 0 & cn < 1){
      #if ( ! clonCNhigh$clonality.status %in% c("clonal","uncertain.clonal")){
      #if ( betaThreshold
      #  zone$integerCN <- floor(cn)
      #  zone$clonality <- clonCNlow$clonality
      #  zone$clonality.min <- clonCNlow$clonality.min
      #  zone$clonality.max <- clonCNlow$clonality.max
      #  zone$clonality.status <- clonCNlow$clonality.status
      #}else{
      #  zone$integerCN <- ceil(cn)
      #  zone$clonality <- clonCNhigh$clonality
      #  zone$clonality.min <- clonCNhigh$clonality.min
      #  zone$clonality.max <- clonCNhigh$clonality.max
      #  zone$clonality.status <- clonCNhigh$clonality.status
      #}
      if ( zone$beta >= betaThreshold ){
        zone$integerCN <- 0
        zone$clonality <- clonCNlow$clonality
        zone$clonality.min <- clonCNlow$clonality.min
        zone$clonality.max <- clonCNlow$clonality.max
        zone$clonality.status <- "clonal"
      }else if (cn <= 0.8) {
        zone$integerCN <- floor(cn)
        zone$clonality <- clonCNlow$clonality
        zone$clonality.min <- clonCNlow$clonality.min
        zone$clonality.max <- clonCNlow$clonality.max
        zone$clonality.status <- clonCNlow$clonality.status
      }else{
        zone$integerCN <- ceiling(cn)
        zone$clonality <- clonCNhigh$clonality
        zone$clonality.min <- clonCNhigh$clonality.min
        zone$clonality.max <- clonCNhigh$clonality.max
        zone$clonality.status <- clonCNhigh$clonality.status
      }
    }
    if (cn >= 1 & cn < 2){
      #if ( cn < 1+clonalityThreshold && ! clonCNhigh$clonality.status %in% c("clonal")){
      if (zone$beta <= betaThreshold){
        zone$integerCN <- floor(cn)
        zone$clonality <- clonCNlow$clonality
        zone$clonality.min <- clonCNlow$clonality.min
        zone$clonality.max <- clonCNlow$clonality.max
        zone$clonality.status <- clonCNlow$clonality.status
      }else{
        zone$integerCN <- ceiling(cn)
        zone$clonality <- clonCNhigh$clonality
        zone$clonality.min <- clonCNhigh$clonality.min
        zone$clonality.max <- clonCNhigh$clonality.max
        zone$clonality.status <- "not.analysed"
      }
    }
    if (cn >= 2 & cn <= 3){
      #if ( cn > 2 + ( 1 - clonalityThreshold) && ! clonCNlow$clonality.status %in% c("clonal")  ){
      if (zone$beta <= betaThreshold){
        zone$integerCN <- ceiling(cn)
        zone$clonality <- clonCNhigh$clonality
        zone$clonality.min <- clonCNhigh$clonality.min
        zone$clonality.max <- clonCNhigh$clonality.max
        zone$clonality.status <- clonCNhigh$clonality.status
      }else{
        zone$integerCN <- floor(cn)
        zone$clonality <- clonCNlow$clonality
        zone$clonality.min <- clonCNlow$clonality.min
        zone$clonality.max <- clonCNlow$clonality.max
        zone$clonality.status <- "not.analysed"
      }
    }


    #if(zone$nsnp < minSNPs || zone$cov < minCov){
    #  zone$clonality.status <- "not.analysed"
    #}

    return(zone)
  }

  #i<-346
  #clonalityTable.list <-lapply(c(1:nrow(betaTable)),getIndexClonality)
  clonalityTable.list <- mclapply(seq(1,nrow(betaTable),1), getIndexClonality, mc.preschedule = T, mc.cores = n_cores)
  clonalityTable <- fromListToDF(clonalityTable.list)

  return(clonalityTable)
  #hist(clonalityTable$ )
}

#'Function to compute clonality of somatic copy number data
#'
#'This function takes the beta table of a tumor sample together with the
#'associated ploidy and admixtures tables and computes the clonality of each
#'segment in the beta table.
#'
#'@param beta_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_beta_table]{compute_beta_table}}
#'@param ploidy_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_ploidy]{compute_ploidy}}
#'@param admixture_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_dna_admixture]{compute_dna_admixture}}
#'@param error_tb  data.frame that reports for each combination of coverage and
#'  number informative SNPs the expected estimation error around beta. The
#'  data.frame error_tb must contains 3 columns: \describe{ \item{mean.cov}{mean
#'  coverage} \item{n.info.snps}{number of informative SNPs}
#'  \item{adm.estimation.error}{estimated error on computed beta on a segment
#'  with coverage mean.cov and n.info.snps informative SNPs} } Package CLONETv2
#'  have built in error_tb named error_table (default=error_table)
#'@param clonality_threshold threshold to discretize continuous clonality value
#'  (default=0.85)
#'@param beta_threshold threshold on beta value to determine clonality direction
#'  (default=0.90)
#'@param n_digits number of digits in the output table  (default=3)
#'@param n_cores number of cores (default=1)
#'@param debug return extra columns for debugging (default=F)
#'@return A data.frame that extends input beta_table with columns \describe{
#'  \item{clonality}{estimated fraction of tumor cell with log2 copy number}
#'  \item{clonality.min}{minum estimated fraction of tumor cell with log2 copy
#'  number} \item{clonality.max}{minum estimated fraction of tumor cell with
#'  log2 copy number} \item{clonality.status}{discretized clonality status into
#'  five values: \emph{clonal}, large majority of the tumor cells has the same
#'  copy number; \emph{subclonal}, not all the tumor cells has the same copy
#'  number; \emph{not.analysed}, is is not possible to determine clonality;
#'  \emph{uncertain.clonal} and \emph{uncertain.subclonal} correspond
#'  respectively to \emph{clonal} and \emph{subclonal} populations but with less
#'  reliable clonality estimate } }
#' @examples
#' \donttest{
#'
#' ## Compute clonality table with default parameters
#' scna_clonality_table_toy <- compute_scna_clonality_table(beta_table = bt_toy,
#'   ploidy_table = pl_table_toy, admixture_table = adm_table_toy)
#' }
#'@author Davide Prandi
#'@export
#'@md
compute_scna_clonality_table<-function(beta_table,
																			 ploidy_table,
																			 admixture_table,
																			 error_tb = error_table,
																			 clonality_threshold = 0.85,
																			 beta_threshold = 0.90,
																			 n_digits=3,
																			 n_cores=1,
																			 debug=F){

	## check consistency of sample names
	sample_id <- unique(beta_table$sample)
	if (length(sample_id) != 1){stop(paste("[",Sys.time() ,"] beta_table must contain exactly one sample\n",sep=""))}

	if (nrow(ploidy_table) !=1 || nrow(admixture_table) !=1 || (pls = ploidy_table$sample) != sample_id || (ads = admixture_table$sample) != sample_id){
		stop(paste("[",Sys.time() ,"] ploidy_table and  admixture_table must contain only sample ",sample_id,"\n",sep=""))
	}

	## merge ploidy table, admixture table and beta_table and adjust log2 by ploidy and admixture
  beta_table <- merge(beta_table, ploidy_table,  by="sample")
	beta_table$log2shift <- round(-log2(beta_table$ploidy/2),n_digits)
	beta_table$log2.plCorr <- round(beta_table$log2 - beta_table$log2shift, n_digits)
	beta_table <- merge(x = beta_table, y = admixture_table, by="sample")
	beta_table$log2.corr <-  suppressWarnings(log2(pmax(( 2 ^ (beta_table$log2.plCorr ) - admixture_table$adm[1] ) / (1 - admixture_table$adm[1]), 0)))

	## compute clonality
	clonality_table <- compute_clonality(betaTable = beta_table,
																			 errorTable = error_tb,
																			 clonalityThreshold = clonality_threshold,
																			 betaThreshold = beta_threshold,
																			 roundDec = n_digits,
																			 n_cores = n_cores )

	if (!debug){
		cols_to_save <- c("sample","chr","start","end","num.mark","log2","beta","nsnp","cov","n_beta","clonality","clonality.min","clonality.max","clonality.status")
		clonality_table <- clonality_table[,intersect(colnames(clonality_table),cols_to_save)]
	}

	return(clonality_table)
}
