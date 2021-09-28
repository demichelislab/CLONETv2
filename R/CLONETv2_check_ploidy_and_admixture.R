#f given a max cn value return all the possible allele specific combinatios
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
#' @param ploidy_table data.frame formatted as the output of function
#' \code{\link[CLONETv2:compute_ploidy]{compute_ploidy}}
#' @param admixture_table data.frame formatted as the output of function
#' \code{\link[CLONETv2:compute_dna_admixture]{compute_dna_admixture}}
#' @return A ggplot2 plot reporting log2 on the x axis and beta and the y axis.
#'   Each dot represents a segment of the input beta_table. Red transparent
#'   circles corresponds to expected log2 vs beta position for different allele
#'   specific copy number combinations given ploidy and admixture reported in
#'   tables ploidy_table and admixture_table, respectively. Labels in the form
#'   (cnA, cnB) indicate repsectively the major and minor allele copy number
#'   value. Labels above the plot comprises sample name and ploddy/admixture estimates.
#' @examples
#'
#' ## check ploidy and admixture estimates
#' check_plot_toy <- check_ploidy_and_admixture(beta_table = bt_toy, ploidy_table = pl_table_toy,
#'   admixture_table = adm_table_toy)
#'
#' @author Davide Prandi
#' @export
#' @md
check_ploidy_and_admixture <- function(beta_table,
																		   ploidy_table,
																			 admixture_table){

	sample_id <- unique(beta_table$sample)

	if (length(sample_id) != 1){stop("Beta table must contains only one sample")}

  # filter good segs
  beta_table <- beta_table[which(!is.na(beta_table$beta)),]

  # compute expected log2, beta
  n_arms <- 4
  possibleCNs <- getPossibleCnComb(n_arms)
  expected_log2_beta <- as.data.frame(possibleCNs)
  expected_log2_beta <- expected_log2_beta[which(expected_log2_beta[,1] >= expected_log2_beta[,2]),]
  expected_log2_beta <- expected_log2_beta[which(expected_log2_beta[,1] !=0 | expected_log2_beta[,2] != 0),]

  colnames(expected_log2_beta) <- c("cnA","cnB")
  expected_log2_beta$admixture <- rep(admixture_table$adm,nrow(expected_log2_beta))
  expected_log2_beta$ploidy <- rep(ploidy_table$ploidy,nrow(expected_log2_beta))
  expected_log2_beta$plT <- expected_log2_beta$cnA + expected_log2_beta$cnB
	expected_log2_beta$plN <- 2
  expected_log2_beta$betaVal <- ( expected_log2_beta$admixture * expected_log2_beta$plN )  /
  	( expected_log2_beta$plT - expected_log2_beta$admixture * ( expected_log2_beta$plT - expected_log2_beta$plN ) )

  expected_log2_beta$n.cells.mono.apparent <-  pmax(expected_log2_beta$cnA,expected_log2_beta$cnB) - pmin(expected_log2_beta$cnA,expected_log2_beta$cnB) #(cnA + cnB) - min(cnA,cnB)
  expected_log2_beta$n.cells.bi.apparent <- 2*pmin(expected_log2_beta$cnA,expected_log2_beta$cnB)

  expected_log2_beta$adm.local.apparent <- ( expected_log2_beta$admixture + (1-expected_log2_beta$adm)*(expected_log2_beta$n.cells.bi.apparent)) /
  	( expected_log2_beta$admixture + (1-expected_log2_beta$admixture)*(expected_log2_beta$n.cells.bi.apparent) +
  			(1-expected_log2_beta$adm)*(expected_log2_beta$n.cells.mono.apparent)  )
  expected_log2_beta$betaVal.apparent <- ( expected_log2_beta$adm.local.apparent * 2 )  / ( 1 + expected_log2_beta$adm.local.apparent )


  expected_log2_beta$log2Val <-
  	log2(( ( expected_log2_beta$betaVal * expected_log2_beta$plN + ( 1 - expected_log2_beta$betaVal ) * expected_log2_beta$plT ) / expected_log2_beta$plN) /
  			 	 (expected_log2_beta$ploidy / expected_log2_beta$plN))

 xmin <- min(-1, min(beta_table$log2, na.rm = T))
 xmax <- max(1, max(beta_table$log2, na.rm = T))

 check_plot <- ggplot(beta_table, aes(x=log2, y=beta)) +
		geom_point(col="gray40") +
		geom_point(data = data.frame(), aes(x=expected_log2_beta$log2Val, y=expected_log2_beta$betaVal.apparent), col="red", size=10, alpha=0.1) +
		ylim(0,1) +
		xlim(xmin,xmax  ) +
		theme_light() +
		xlab("LogR") +
		ylab("beta") +
		labs(title =  sample_id, subtitle = paste0("Ploidy ",ploidy_table$ploidy,"\nAdmixture ",admixture_table$adm))

 expected_log2_beta$plot_label <- paste0("(",expected_log2_beta$cnA,",",expected_log2_beta$cnB,")")
 check_plot <- check_plot + geom_text_repel(data = data.frame(), aes(x=expected_log2_beta$log2Val, y=expected_log2_beta$betaVal.apparent, label = expected_log2_beta$plot_label))

 return(check_plot)


}



