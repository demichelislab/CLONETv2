
computeAllelicImbalance <- function(betaTable,
                                    cnTh = 0.5, #the max (min) distance from the integer cn value
                                    Ncores = 4
                                     ){

  #round cn values
  line<-1
  computeAllelicImbalance.line <- function(line,betaTable,cnTh){

    dfOut <- betaTable[line,]
    #cat(line,"\n")
    dfOut$AllelicImbalance <- NA
    dfOut$AllelicImbalance.int <- NA
    dfOut$cnA.int <- NA
    dfOut$cnB.int <- NA
    dfOut$isCNNL <- 0
    dfOut$isLOHg <- 0
    dfOut$isLOHseg <- 0
    dfOut$isUnbalancedGain <- 0

    if ( !is.na(dfOut$cnA) & !is.na(dfOut$cnB) ){
      if ( dfOut$cnA < floor(dfOut$cnA) + cnTh ){
        dfOut$cnA.int <- floor(dfOut$cnA)
      }else if (dfOut$cnA >= ceiling(dfOut$cnA) - cnTh){
        dfOut$cnA.int <- ceiling(dfOut$cnA)
      }

      if ( dfOut$cnB < floor(dfOut$cnB) + cnTh ){
        dfOut$cnB.int <- floor(dfOut$cnB)
      }else if (dfOut$cnB >= ceiling(dfOut$cnB) - cnTh){
        dfOut$cnB.int <- ceiling(dfOut$cnB)
      }

      dfOut$AllelicImbalance <- dfOut$cnA - dfOut$cnB
      dfOut$AllelicImbalance.int <- dfOut$cnA.int - dfOut$cnB.int
      if (!is.na(dfOut$cnA.int) & !is.na(dfOut$cnB.int)){
        if (dfOut$cnA.int == 2 & dfOut$cnB.int == 0){
          dfOut$isCNNL <- 1
        }
        if (dfOut$cnA.int > 2 & dfOut$cnB.int == 0){
          dfOut$isLOHg <- 1
        }
        if (dfOut$cnA.int > 1 & dfOut$cnB.int == 0){
          dfOut$isLOHseg <- 1
        }
        if (dfOut$cnA.int > 2 &  dfOut$cnB.int > 0 & ( dfOut$cnA.int - dfOut$cnB.int > 1)){
          dfOut$isUnbalancedGain <- 1
        }

      }
    }
    return(dfOut)
  }

  betaTableOut.List <- mclapply(seq(1,nrow(betaTable),1), computeAllelicImbalance.line, betaTable, cnTh, mc.preschedule = T, mc.cores = Ncores )
  betaTableOut.List <- lapply(seq(1,nrow(betaTable),1), computeAllelicImbalance.line, betaTable, cnTh)
  betaTableOut <- fromListToDF(betaTableOut.List)

}



#'Function to compute allele specific somatic copy number
#'
#'This function takes the beta table of a tumor sample together with the
#'associated ploidy and admixtures tables and computes the allele specific copy
#'number of each segment in the beta table.
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
#'@param allelic_imbalance_th maximum distance from allele spefici copy number
#'  of a segment to define integer alelle specific copy number value. Value 0.5
#'  corresponds to round cnA and cnB (default=0.5)
#'@param n_digits number of digits in the output table  (default=3)
#'@param n_cores number of cores (default=1)
#'@param debug return extra columns for debugging (default=F)
#'@return A data.frame that extends input beta_table with columns \describe{
#'  \item{log2.corr}{log2 ratio adjusted by ploidy and admixture}
#'  \item{cnA}{copy number of the major allele} \item{cnB}{copy number of the
#'  minor allele} \item{cnA.int}{integet copy number of the major allele}
#'  \item{cnB.int}{integet copy number of the minor allele} }
#' @examples
#'
#' ## Compute clonality table with default parameters
#' allele_specific_cna_table_toy <- compute_allele_specific_scna_table(
#'   beta_table = bt_toy, ploidy_table = pl_table_toy,
#'   admixture_table = adm_table_toy)
#'
#'@author Davide Prandi
#'@export
#'@md
compute_allele_specific_scna_table<-function(beta_table,
																		   			 ploidy_table,
																			       admixture_table,
																			       error_tb = error_table,
																			       allelic_imbalance_th = 0.5,
																			       n_digits=3,
																			       n_cores=1,
																						 debug=F){

  ## check consistency of sample names
	sample_id <- unique(beta_table$sample)
	beta_cols <- colnames(beta_table)
	if (length(sample_id) != 1){stop(paste("[",Sys.time() ,"] beta_table must contain exactly one sample\n",sep=""))}

	if (nrow(ploidy_table) !=1 || nrow(admixture_table) !=1 || (pls = ploidy_table$sample) != sample_id || (ads = admixture_table$sample) != sample_id){
		stop(paste("[",Sys.time() ,"] ploidy_table and  admixture_table must contain only sample ",sample_id,"\n",sep=""))
	}

	## merge ploidy table, admixture table and beta_table and adjust log2 by ploidy and admixture
  beta_table <- merge(beta_table, ploidy_table,  by="sample")
	beta_table$log2shift <- round(-log2(beta_table$ploidy/2),n_digits)
	beta_table$log2.plCorr <- beta_table$log2 - beta_table$log2shift
	beta_table <- merge(x = beta_table, y = admixture_table, by="sample")
  ## note: consider setting to 0 when log2.corr is NA
	beta_table$log2.corr <-  suppressWarnings(log2( ( 2 ^ (beta_table$log2.plCorr ) - admixture_table$adm[1] ) / (1 - admixture_table$adm[1])))

	## add error information
  beta_table.list <- by(beta_table, INDICES = beta_table$sample, FUN = addErrorToBetaTable, errorTable = error_tb, ncores = n_cores )
  beta_table <- fromListToDF(beta_table.list)
  rm(beta_table.list)

	## add error information
  beta_table.list <- by(beta_table, INDICES = beta_table$sample, FUN = addErrorToBetaTable, errorTable = error_tb, ncores = n_cores )
  beta_table <- fromListToDF(beta_table.list)
  rm(beta_table.list)

	beta_table <- extendBetaTableWithCopyNumbers(BetaTable = beta_table, G = admixture_table$adm[1], ncores = n_cores)
  beta_table <- computeAllelicImbalance(betaTable = beta_table, cnTh = allelic_imbalance_th, Ncores = n_cores)



  if (!debug){
		valid_cols <- c(beta_cols,"log2.corr","cnA", "cnB","cnA.int", "cnB.int")
		beta_table <- beta_table[,valid_cols]
  }
  return(beta_table)
}
