## add to vep annotated snvs columns for clonality computations
preprocess_snvtable_vep <- function (snv_table, normal_sample_name, tumor_sample_name) {
  filt_snv_table <- snv_table[snv_table$VARIANT_CLASS == 'SNV', ]
  with(filt_snv_table, {
           get_el    <- function (i) function (ee) ee[i]
           split_loc <- strsplit(Location, ':')
           pos       <- as.numeric(sapply(split_loc, get_el(2)))
           data.frame(
             filt_snv_table,
             Chromosome                  = sapply(split_loc, get_el(1)),
             Start_position              = pos,
             End_position                = pos,
             Tumor_Sample_Barcode        = tumor_sample_name,
             Matched_Norm_Sample_Barcode = normal_sample_name,
             stringsAsFactors            = F,
             check.names                  = F
           )
         })
}



## merge a table reporting snvs and a beta table
## snvs table needs columns: sample, chr, start, end
## beta table needs columns: sample, chr, start, end
extend_SNVt_with_bt <- function(SNVtable,
                                betaTable, # five columns table: chr  start	end	HUGO	cyto	arm
                                Ncores = 4){
  #snvID<-278
  findSNVinBT <- function(snvID, SNVtable, betaTable){
    #cat(snvID,"\n")
    thisBT <- betaTable[which(betaTable$sample == SNVtable$Tumor_Sample_Barcode[snvID]),]
    #betaPoses <- getSegmentsPos(SNVtable$chr[snvID], SNVtable$start[snvID], SNVtable$end[snvID], segments = thisBT[,c("chr","start","end")] )
    betaPoses <- getSegmentsPos(SNVtable$Chromosome[snvID], SNVtable$Start_position[snvID], SNVtable$End_position[snvID], segments = thisBT[,c("chr","start","end")] )

    if (length(betaPoses) > 1){
      stop("Error ",snvID,": many segments intersect one position\n")

    }
    if (length(betaPoses) == 0){
      ## cat("Warning ",snvID,": no segments\n")
      ## create an empty bt
      emptyBT <- thisBT[1,,drop=F]
      # emptyBT$chr <- SNVtable$chr[snvID]
      # emptyBT$start <- SNVtable$start[snvID]
      # emptyBT$end  <- SNVtable$end[snvID]
      emptyBT$chr <- SNVtable$Chromosome[snvID]
      emptyBT$start <- SNVtable$Start_position[snvID]
      emptyBT$end  <- SNVtable$End_position[snvID]
      emptyBT[,setdiff(colnames(emptyBT),c("sample","chr","start","end","adm","ploidy"))] <- NA
      thisBT <- emptyBT
      betaPoses <- 1
    }

    thisGene <- SNVtable[snvID,,drop=F]
    thisGene <- cbind(thisGene,thisBT[betaPoses,])

    return(thisGene)
  }

  geneBT.list <- mclapply(seq(1,nrow(SNVtable),1), findSNVinBT, SNVtable = SNVtable, betaTable = betaTable, mc.preschedule = T, mc.cores = Ncores)
  #lapply(seq(20,nrow(SNVtable),1), findSNVinBT, SNVtable = SNVtable, betaTable = betaTable)
  geneBT <- fromListToDF(geneBT.list)
  return(geneBT)
}


## given and snv table extended with allelic specific copy nubmer data compute the clonality of each snvs
## requried columns
#SNVid <-12
#SNVtable.ext <- snv_read_count_ext
computeSNVclonality <- function(SNVid, SNVtable.ext){

  ## select snvs
  ## cat(SNVid,"\n")
  thisSNV <- SNVtable.ext[SNVid,]


  thisSNV$t_cov <- NA
  thisSNV$t_af <- NA
  thisSNV$n_admReads <- NA
  thisSNV$t_ref_count_corr <- NA
  thisSNV$t_af_corr <- NA
  thisSNV$cn.int <- NA
  thisSNV$CN_SNVmut <- NA
  thisSNV$VAFexp <- NA
  thisSNV$SNV.clonality <- NA
  thisSNV$SNV.clonality.int <- NA
  if (is.na(thisSNV$cnA) || is.na(thisSNV$rc_alt_tumor)){
    return(thisSNV)
  }

  thisSNV$t_cov <- thisSNV$rc_alt_tumor + thisSNV$rc_ref_tumor
  if (thisSNV$t_cov == 0){
    return(thisSNV)
  }



  thisSNV$t_af <- thisSNV$rc_alt_tumor / thisSNV$t_cov
  #thisSNV$n_admReads  <- round((2 * thisSNV$ASEQ_t_cov * thisSNV$adm) / ( 2 * thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm)))
  thisSNV$n_admReads  <- round((2 * thisSNV$t_cov * thisSNV$adm) /
                         ( 2 * thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm)))


  thisSNV$t_ref_count_corr <- thisSNV$rc_ref_tumor - thisSNV$n_admReads
  thisSNV$t_af_corr <-  thisSNV$rc_alt_tumor / (thisSNV$rc_alt_tumor + thisSNV$t_ref_count_corr)
  #thisSNV$ASEQ_t_alt/thisSNV$ASEQ_t_cov
  #### compute clonality (% of tumor cell with the mutation)
  ## firt find the number of allele mutated
  thisSNV$cn.int <- thisSNV$cnA.int + thisSNV$cnB.int
  if (thisSNV$cn.int <= 0){return(thisSNV)}

  possibleAF <- seq(1,thisSNV$cn.int,1)/thisSNV$cn.int
  #possibleAF <- c(1,2,3)/3
  #CN_M <- seq(1,thisSNV$cn.int,1)[which.min(abs(possibleAF - thisSNV$ASEQ_t_AF_corr))]
  #seq(1,3,1)[which.min(abs(possibleAF - thisSNV$ASEQ_t_AF_corr))]
  thisSNV$CN_SNVmut <- seq(1,thisSNV$cn.int,1)[which.min(abs(possibleAF - thisSNV$t_af_corr))]
  thisSNV$VAFexp <- (thisSNV$CN_SNVmut * (1 - thisSNV$adm)) /
    ( 2 * thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm))
  thisSNV$SNV.clonality =  ((thisSNV$t_af - thisSNV$VAFexp) *  ( 2 * thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm))) / (thisSNV$CN_SNVmut * (1 -thisSNV$adm)) + 1
  thisSNV$SNV.clonality.int =  ((thisSNV$t_af - thisSNV$VAFexp) * ( 2 * thisSNV$adm + (thisSNV$cnA.int + thisSNV$cnB.int) * (1 - thisSNV$adm))) / (thisSNV$CN_SNVmut * (1 -thisSNV$adm)) + 1

  return(thisSNV)
}



#'Function to compute clonality of SNVs
#'
#'This function takes as input the genomic position of a SNVs and computes the
#'percentage of genomic homogeneus cells harboring the mutation.
#'
#'
#'@param sample_id the id of the analyzed sample. It must be the same value
#'  reported in column sample of tables beta_table, ploidy_table, and
#'  admixture_table
#'@param snv_read_count data.frame reporting in each row the genomic coordinates
#'  of an SNV together with number of reference and alternative reads covering
#'  the position in columns rc_ref_tumor and rc_alt_tumor, respectively. See
#'  parameter annotation_style for details about column names
#'@param beta_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_beta_table]{compute_beta_table}}
#'@param ploidy_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_beta_table]{compute_ploidy}}
#'@param admixture_table data.frame formatted as the output of function
#'  \code{\link[CLONETv2:compute_beta_table]{compute_dna_admixture}}
#'@param error_tb  data.frame that reports for each combination of coverage and
#'  number informative SNPs the expected estimation error around beta. The
#'  data.frame error_tb must contains 3 columns: \describe{ \item{mean.cov}{mean
#'  coverage} \item{n.info.snps}{number of informative SNPs}
#'  \item{adm.estimation.error}{estimated error on computed beta on a segment
#'  with coverage mean.cov and n.info.snps informative SNPs} } Package CLONETv2
#'  have built in error_tb named error_table (default=error_table)
#'@param error_rate expected fraction of SNV positions with outlier variant
#'  allelic fraction (default=0.05)
#'@param n_digits number of digits in the output table  (default=3)
#'@param n_cores number of cores (default=1)
#'@param annotation_style a string that corresponds to the format of the columns
#'  that describe the genomic coordinates of a SNV. Accepted values are VEP and
#'  MAF. [VEP
#'  annotation](https://www.ensembl.org/info/docs/tools/vep/index.html)
#'  describes genomic coordinates with a single column named Location. [MAF
#'  format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) has
#'  columns Chromosome, Start_position, and End_position for each aberrant
#'  position
#'@param debug return extra columns for debugging (default=F)
#'@return A data.frame that extends input table snv_read_count with columns
#'  sample, cnA, cnB, t_af, t_af_corr, SNV.clonality, and SNV.clonality.status.
#'  Columns cnA and cnB report the allele specific copy number of the genomic
#'  segment containing the SNV position. Columns t_af and t_af_corr are
#'  respectively raw and ploidy/purity adjusted tumor varian allelic fractions.
#'  SNV.clonality reports the percentage of tumor cells harboring the SNV and
#'  with allele specific copy number cnA and cnB. SNV.clonality.status column
#'  lists dicretized SNV.clonality values. Discrete states are clonal,
#'  uncertain.clonal, uncertain.subclonal, and subclonal based in threshold
#'  automatically computed on the SNV.clonality values. Empty
#'  SNV.clonality.status of an SNV indicates that clonality cannot be assessed.
#' @examples
#'
#' ## Compute SNVs clonality
#' snv_clonality_table_toy <- compute_snv_clonality("toy_sample",
#'   snv_reads_toy, bt_toy, pl_table_toy, adm_table_toy)
#'
#'@author Davide Prandi, Tarcisio Fedrizzi
#'@export
#'@md
compute_snv_clonality<-function(sample_id,
																snv_read_count,
  															beta_table,
 		                            ploidy_table,
																admixture_table,
																error_tb = error_table,
																error_rate=0.05,
		                            n_digits=3,
		                            n_cores=1,
																annotation_style = "VEP",
																debug=F
){

	## check consistency of sample names
	bt_id <- unique(beta_table$sample)
  snv_cols <- colnames(snv_read_count)
	if (length(sample_id) != 1){stop(paste("[",Sys.time() ,"] beta_table must contain exactly one sample\n",sep=""))}

	if (nrow(ploidy_table) !=1 || nrow(admixture_table) !=1 || (pls = ploidy_table$sample) != bt_id || (ads = admixture_table$sample) != bt_id){
		stop(paste("[",Sys.time() ,"] ploidy_table and  admixture_table must contain only sample ",sample_id,"\n",sep=""))
	}

	if(sample_id != bt_id){
		stop(paste("[",Sys.time() ,"] sample_id and sample in beta_table do not match",sample_id,"\n",sep=""))
	}

	beta_table <- compute_allele_specific_scna_table(beta_table, ploidy_table, admixture_table, error_tb =  error_tb, n_digits =  n_digits, n_cores = n_cores, debug = T)

	if (annotation_style == "VEP"){
		snv_read_count <- preprocess_snvtable_vep(snv_table = snv_read_count, normal_sample_name = paste0(sample_id,"_normal"), tumor_sample_name = sample_id)
	}else if (annotation_style == "MAF"){

	}else{
		stop(paste("[",Sys.time() ,"] only VEP or MAF annotations are supported",sample_id,"\n",sep=""))
	}

	###################################
	## compute clonality
	snv_read_count_ext <- extend_SNVt_with_bt(SNVtable = snv_read_count, betaTable = beta_table, Ncores = n_cores)

	SNVdata.Cl.list <- mclapply(seq(1,nrow(snv_read_count_ext),1), computeSNVclonality, snv_read_count_ext, mc.preschedule = T, mc.cores = n_cores)
	SNVdata.Cl <- fromListToDF(SNVdata.Cl.list)

	## correct clonality > 1
	limitVal <- stats::quantile(SNVdata.Cl$SNV.clonality,na.rm = T,probs = 1-error_rate)
	SNVdata.Cl$SNV.clonality[which(SNVdata.Cl$SNV.clonality >= limitVal)] <- NA
	SNVdata.Cl$SNV.clonality[which(SNVdata.Cl$SNV.clonality > 1)] <- 2-SNVdata.Cl$SNV.clonality[which(SNVdata.Cl$SNV.clonality > 1)]

	## discretize calls
	## discretize clonality
	SNVdata.Cl$SNV.clonality.status <- ""
	if (nrow(SNVdata.Cl) > 4){
		clonality_intervals <- suppressWarnings(discretize(x = SNVdata.Cl$SNV.clonality, method = "cluster", categories = 4, onlycuts = T, nstart=20, iter.max=500))
		SNVdata.Cl$SNV.clonality.status[which(SNVdata.Cl$SNV.clonality >= clonality_intervals[1]  & SNVdata.Cl$SNV.clonality < clonality_intervals[2] )] <- "subclonal"
		SNVdata.Cl$SNV.clonality.status[which(SNVdata.Cl$SNV.clonality >= clonality_intervals[2]  & SNVdata.Cl$SNV.clonality < clonality_intervals[3] )] <- "uncertain.subclonal"
		SNVdata.Cl$SNV.clonality.status[which(SNVdata.Cl$SNV.clonality >= clonality_intervals[3]  & SNVdata.Cl$SNV.clonality < clonality_intervals[4] )] <- "uncertain.clonal"
		SNVdata.Cl$SNV.clonality.status[which(SNVdata.Cl$SNV.clonality > clonality_intervals[4]   )] <- "clonal"
	}


	if (!debug){
     save_cols <- c("sample","cnA", "cnB", "t_af", "t_af_corr", "SNV.clonality", "SNV.clonality.status")
     SNVdata.Cl <-  SNVdata.Cl[,c(snv_cols,save_cols)]
	}

	return(SNVdata.Cl)


}
