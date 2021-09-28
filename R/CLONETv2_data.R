#' Beta estimation error.
#'
#' A precomputed table reporting for different combinations of coverage and
#' number of informative SNPs  the expected error of the  beta value computed by
#' function \code{\link[CLONETv2:compute_beta_table]{compute_beta_table}}.
#'
#' @format A data frame column names mean.cov, n.info.snps, and
#'   adm.estimation.error \describe{ \item{mean.cov}{genomic segment coverage}
#'   \item{n.info.snps}{number of informative SNPs}
#'   \item{adm.estimation.error}{expected error on beta estimate } }
"error_table"

#' Toy example of segmetd data.
#'
"seg_tb_toy"

#' Toy example of tumor pileup data.
#'
"pileup_tumor_toy"

#' Toy example of normal pileup data.
#'
"pileup_normal_toy"

#' Toy example of snv data.
#'
"snv_reads_toy"

#' Toy example of beta table.
#'
"bt_toy"

#' Toy example of ploidy table.
#'
"pl_table_toy"

#' Toy example of admixture table.
#'
"adm_table_toy"

#' Toy example of clonality table of somatic copy number.
#'
"scna_clonality_table_toy"

#' Toy example of allele specific table of somatic copy number.
#'
"allele_specific_cna_table_toy"

#' Toy example of snv clonality table.
#'
"snv_clonality_table_toy"


