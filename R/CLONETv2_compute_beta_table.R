#source("R/CLONETv2_basic_functions.R")

getValsfromPDF <- function(pdf,n)
{
  return(pdf[sample(1:length(pdf),n,replace=T)])
}

# af <- af.control
# cov <- cov.control
# size <- 1000
# subsample_data <- function(af, cov, size){
#
# 	af_valid_ids <- which(af >= quantile(af)[2] & af <= quantile(af)[4])
#
# 	af_valid <- af[af_valid_ids]
# 	cov_valid <- cov[af_valid_ids]
#
# 	if (length(af_valid) > size){
#
# 		cov_quant_ids <- which(cov_valid >= quantile(cov_valid)[2] & cov_valid <= quantile(cov_valid)[4])
#     cov_valid <- cov_valid[cov_quant_ids]
#     af_valid <- af_valid[cov_quant_ids]
#
#     if (length(af_valid) > size){
#     	set.seed(3451)
#     	random_ids <- sample(x = length(af_valid),size = size, replace = T)
#   		cov_valid <- cov_valid[random_ids]
# 			af_valid <- af_valid[random_ids]
#     }
#
# 	}
#
# 	return(list(cov = cov_valid, af = af_valid))
# }

BetaDistr <- function(af.control,cov.control)
{
  rep = 1000
	#sample_data <- subsample_data(af.control, cov.control, 10000)
  #tumor <- sample_data$cov
  #tumor.af <- sample_data$af
	beta.distr = list()
  for(beta in seq(1,0.01,-0.01))
  {
    tumor = getValsfromPDF(cov.control,min(rep, length(cov.control)))
    tumor[which(tumor<20)] = 20
    tumor.adm = round(tumor*beta)
    tumor.del = round(tumor-tumor.adm)
    tumor.af = getValsfromPDF(af.control,min(rep, length(af.control)))
    tumor.alt = round(tumor.adm*tumor.af)
    tumor.ref = tumor.adm-tumor.alt

    number = sample(1:length(tumor.ref),1)
    dir = sample(1:length(tumor.ref),number)

    if(beta<1)
    {
      tumor.ref[dir] = tumor.ref[dir] + tumor.del[dir]
      if(length(dir)<length(tumor.ref))
        tumor.alt[-dir] = tumor.alt[-dir] + tumor.del[-dir]
    }
    tumor.af = tumor.alt/(tumor.alt+tumor.ref)

    af.mirror=tumor.alt/(tumor.ref+tumor.alt)
    af.mirror[which(af.mirror<0.5)] = 1-af.mirror[which(af.mirror<0.5)]

    beta.distr[[length(beta.distr)+1]] = af.mirror
  }
  return(beta.distr)
}

BetaDistrMonotone <- function(beta.distr)
{
  for(i in 2:length(beta.distr))
  {
    if(stats::median(beta.distr[[i]]) < stats::median(beta.distr[[i-1]]))
    {
      shift  = stats::median(beta.distr[[i-1]]) - stats::median(beta.distr[[i]])
      beta.distr[[i]] = beta.distr[[i]]+shift
    }
  }
  return(beta.distr)
}

computeBeta <- function(snps.distr,beta.distr,p.thr=0.01)
{
  betas = seq(1,0.01,-0.01)
  evidence = 1
  if (stats::wilcox.test(snps.distr,beta.distr[[1]],alternative="greater")$p.value>p.thr)
    evidence = 0
  greater = 1
  for(i in 1:length(beta.distr))
  {
    if(stats::wilcox.test(snps.distr,beta.distr[[i]],alternative="greater")$p.value<p.thr) { greater=i } else { break }
  }
  less = greater = i
  for(j in i:length(beta.distr))
  {
    if(stats::wilcox.test(snps.distr,beta.distr[[j]],alternative="less")$p.value<p.thr) { less=j;break }
  }
  less = j
  if(less>greater)
    less = j-1
  q = 1-(stats::quantile(snps.distr)-min(snps.distr))/(max(snps.distr)-min(snps.distr))
  beta = betas[less] + (betas[greater]-betas[less])*q[3]
  beta.max = betas[less] + (betas[greater]-betas[less])*q[2]
  beta.min = betas[less] + (betas[greater]-betas[less])*q[4]
  return(c(beta,beta.min,beta.max,betas[greater],betas[less],evidence))
}


#' Function to compute beta table
#'
#' This function takes segmented data and per base pileup of tumor and matched
#' normal of a sample as input and associates a beta value to each genomic segment.
#'
#' @param seg_tb data.frame in [SEG
#'   format](http://software.broadinstitute.org/software/igv/SEG). Rows
#'   report per segment log2 ratio
#'   numeric value. CLONETv2 inteprets first column as sample name, columns two
#'   to four as genomic coordinates (chromosome, start location, and end
#'   location), column five is not used, and column six is the log2 ratio
#'   returned by segmentation algorithm.
#' @param pileup_tumor,pileup_normal data.frame reporting pileup of SNPs in
#'   tumor and normal samples respectively. First row contains column names and
#'   subsequent rows report the pileup of a specific genomic positions. Required
#'   information for each genomic position includes
#'   chromosome,	position,	allelic fraction, and	coverage. Required column names
#'   are chr, pos, af, and cov
#' @param min_coverage minimum number of reads for considering a pileup position
#'   valid (default=20)
#' @param min_required_snps minimum number of snps to call beta for a segment
#'   (default=10)
#' @param min_af_het_snps minimum allowed allelic fraction of a SNP genomic
#'   position (default=0.2)
#' @param max_af_het_snps maximum allowed allelic fraction of a SNP genomic
#'   position (default=0.8)
#' @param n_digits number of digits in the output table  (default=3)
#' @param n_cores number of available cores for computation (default=1)
#' @param plot_stats plot summary statistics of the computed beta table (default=F)
#' @param debug return extra columns for debugging (default=F)
#' @return A data.frame that extends input seg_tb with columns beta, nsnp, cov,
#'   n_beta. Moreover, CLONETv2 renames colums of seg_tb as sample, chr,
#'   start, end, XYZ, log2, with XYZ being the original name of column five
#'   As for seg_tb, each raw of the output table represents a genomic
#'   segments. For each raw,  the value of beta is the proportion of neutral
#'   reads in the segment, while nsnp and cov represents respectively the number
#'   of informative SNPs and the mean coverage of the given segment. The value
#'   n_beta is the proportion of neutral reads in the normal sample. The value
#'   of n_beta should be 1 as in normal samples parental chromosomes are equally
#'   represented. Values lower than 1 of n_beta could indicate the presence of germline CNVs
#'   or sequencing errors.
#' @examples
#'
#' ## Compue beta table with default parameters
#' bt_toy <- compute_beta_table(seg_tb_toy, pileup_tumor_toy, pileup_normal_toy)
#' @author Davide Prandi, Alessandro Romanel
#' @export
#' @md
compute_beta_table<-function(seg_tb,
                             pileup_tumor,
                             pileup_normal,
                             min_coverage=20,
                             min_required_snps=10,
														 min_af_het_snps=0.2,
														 max_af_het_snps=0.8,
                             n_digits=3,
                             n_cores=1,
														 plot_stats=F,
														 debug=F
){

	## rename colnames
	colnames(seg_tb)[c(1:4,6)] <- c("sample","chr","start","end","log2")

	## check if sample is uniqe in seg file
	sample_id <- unique(seg_tb$sample)
	if(length(sample_id) > 1){stop("Table seg_tb contains more than 1 sample")}

	## check pileup columns
	if (!all(c("chr","pos","af","cov") %in% colnames(pileup_tumor)) | !all(c("chr","pos","af","cov") %in% colnames(pileup_normal))){
		stop("Tables pileup_tumor and pileup_normal must have columns chrm pos, af and cov")
	}

  ##############
  ### filter pileup
	pileup_normal <- pileup_normal[which( pileup_normal$ref %in% c("A","C","G","T") &
	 																			 pileup_normal$alt %in% c("A","C","G","T")),]
	pileup_tumor <- pileup_tumor[which( pileup_tumor$ref %in% c("A","C","G","T") &
	 																			 pileup_tumor$alt %in% c("A","C","G","T")),]


  pileup_normal <- pileup_normal[which( pileup_normal$af >= min_af_het_snps &
  																			pileup_normal$af <= max_af_het_snps &
  																			pileup_normal$cov >= min_coverage     ),]
  pileup_normal$UID <- paste0(pileup_normal$chr,":",pileup_normal$pos)
  pileup_tumor$UID <- paste0(pileup_tumor$chr,":",pileup_tumor$pos)
  pileup_tumor <- pileup_tumor[which( pileup_tumor$cov >= min_coverage & pileup_tumor$UID %in% pileup_normal$UID     ),]

  # use only chromosomes 1-22
  pileup_tumor <- pileup_tumor[which(suppressWarnings(as.numeric(gsub("chr|CHR|Chr","",pileup_tumor$chr))) %in% seq(1,22,1)),]


  ## some checks before computing beta
  if (nrow(pileup_tumor) == 0){
    stop("No valid heterozygous SNPs identified in tumor_pileup table")
  }

  if (nrow(seg_tb) == 0){
    stop("No segments in seg_tb table")
  }

  ## for reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(utf8ToInt(sample_id)[1]*sum(utf8ToInt(sample_id)))

  #############
  ### compute the distribution of beta on normal to asses noise
  beta.distr.raw = BetaDistr(af.control = pileup_normal$af,cov.control = pileup_normal$cov)
  beta.distr = BetaDistrMonotone(beta.distr.raw)

  # compute beta for each tumor sample raw
  extend_seg_with_beta_value <- function(seg_number, seg_tb, pileup_tumor, pileup_normal,beta.distr){

  	  #cat(seg_number,"\n")

			thisSeg <- seg_tb[seg_number,,drop=F]

      ## init output
      thisSeg$beta = NA
      thisSeg$nref = NA
      thisSeg$nsnp = NA
      thisSeg$cov = NA
      thisSeg$sd = NA
      thisSeg$AFmean = NA
      thisSeg$AFsd = NA
      thisSeg$pval = NA
      thisSeg$delevidence = NA
      thisSeg$beta75 = NA
      thisSeg$beta25 = NA
      thisSeg$betamin = NA
      thisSeg$betamax = NA

      ## add beta on normal
      thisSeg$n_beta = NA
      thisSeg$n_delevidence = NA
      thisSeg$n_beta75 = NA
      thisSeg$n_beta25 = NA
      thisSeg$n_betamin = NA
      thisSeg$n_betamax = NA

  		representativeSNPs <-   pileup_tumor[which(
        pileup_tumor$chr == thisSeg$chr &
        pileup_tumor$pos >= thisSeg$star &
        pileup_tumor$pos <= thisSeg$end) ,]

  		if(nrow(representativeSNPs)>=min_required_snps)
      {
        af.seg = representativeSNPs$af
        af.seg[which(af.seg<0.5)] = 1-af.seg[which(af.seg<0.5)]
        beta = computeBeta(snps.distr=af.seg,beta.distr=beta.distr,p.thr=0.01)
        if(length(beta)>0)
        {
          thisSeg$beta = round(beta[1],digits=n_digits)
          thisSeg$nref = NA
          thisSeg$nsnp = nrow(representativeSNPs)
          thisSeg$cov = mean(representativeSNPs$cov)
          thisSeg$sd = stats::sd(representativeSNPs$cov)
          thisSeg$AFmean = mean(representativeSNPs$af)
          thisSeg$AFsd = stats::sd(representativeSNPs$af)
          thisSeg$pval = 0.01
          thisSeg$delevidence = round(beta[6])
          thisSeg$beta75 = round(beta[2],digits=n_digits)
          thisSeg$beta25 = round(beta[3],digits=n_digits)
          thisSeg$betamin = round(beta[5],digits=n_digits)
          thisSeg$betamax = round(beta[4],digits=n_digits)
        }

	       ## compute on normal
	       n_representativeSNPs <- pileup_normal[which(pileup_normal$UID %in% representativeSNPs$UID),]
	       n_af.seg = n_representativeSNPs$af
	       n_af.seg[which(n_af.seg<0.5)] = 1-n_af.seg[which(n_af.seg<0.5)]
	       n_beta = computeBeta(snps.distr=n_af.seg,beta.distr=beta.distr,p.thr=0.01)
	       if(length(n_beta)>0)
	       {
	        	thisSeg$n_beta = round(n_beta[1],digits=n_digits)
	        	thisSeg$n_delevidence = round(n_beta[6])
	        	thisSeg$n_beta75 = round(n_beta[2],digits=n_digits)
	        	thisSeg$n_beta25 = round(n_beta[3],digits=n_digits)
	        	thisSeg$n_betamin = round(n_beta[5],digits=n_digits)
	        	thisSeg$n_betamax = round(n_beta[4],digits=n_digits)
	       }
	      }else{
	        thisSeg$nsnp = nrow(representativeSNPs)
	        thisSeg$cov = mean(representativeSNPs$cov)
	        thisSeg$sd = stats::sd(representativeSNPs$cov)
	        thisSeg$AFmean = mean(representativeSNPs$af)
	        thisSeg$AFsd = stats::sd(representativeSNPs$af)
	      }
  		return(thisSeg)

  }

  res = mclapply(seq(1,nrow(seg_tb),1),extend_seg_with_beta_value,seg_tb, pileup_tumor, pileup_normal,beta.distr,
  							 mc.preschedule = T, mc.cores = n_cores, mc.set.seed = T)

  outTable <- fromListToDF(res)

  if (plot_stats){
    n_segments <- nrow(outTable)
    n_segments_with_beta <- length(which(!is.na(outTable$beta)))
    fraction_analyzed_segments <- n_segments_with_beta / n_segments
    seg_lenght_distribution <- stats::quantile(outTable$end-outTable$start + 1)
    seg_cov_distribution <- stats::quantile(outTable$cov, na.rm = T)

    n_snps_distr <- stats::quantile(outTable$nsnp)

    out_text <-
    	paste(
    		"Computed beta table of sample \"",sample_id,"\"\n ",
    		"Number of processed segments: ",n_segments,"\n ",
    		"Number of segments with valid beta: ",n_segments_with_beta," (",round(fraction_analyzed_segments*100),"%)\n ",
    		"Quantiles of input segment lenghts:\n  ",
    		paste(format(names(seg_lenght_distribution)), format(seg_lenght_distribution), sep = ":", collapse = "\n  "),
    		"\n ",
    		"Quantiles of input segment coverage:\n  ",
    		paste(format(names(seg_cov_distribution)), format(seg_cov_distribution), sep = ":", collapse = "\n  "),
    		"\n ",
    		"Quantiles of number of informative SNPs per input segment:\n  ",
    	  paste(format(names(n_snps_distr)), format(n_snps_distr), sep = ":", collapse = "\n  "),
    		sep = ""
    	)
    cat(out_text)
  }

  if (!debug){
  	seg_columns <-colnames(outTable)[1:6]
  	extra_cols <- c("beta", "nsnp", "cov","n_beta")
  	outTable <- outTable[,c(seg_columns,extra_cols)]
  }


  return(outTable)

}
