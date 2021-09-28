\name{CLONETv2-package}
\alias{CLONETv2-package}
\docType{package}
\title{
CLONETv2
}
\description{
This package is designed to analyze data from next-generation sequencing
experiments on genomic samples. It offers a
set of functions to compute allele specific copy number and
clonality from segmented data and SNPs position pileup. The library
also calculated the clonality of single nucleotide variants given
read counts at mutated positions.

The package has been developed at the laboratory of Computational and Functional Oncology,
Department of CIBIO, University of Trento (Italy), under the supervision of prof.
Francesca Demichelis.
}
\author{
Maintainer:
Davide Prandi

Contributors:
Tarcisio Fedrizzi, Alessandro Romanel
}
\references{
Prandi, D., Baca, S. C., Romanel, A., Barbieri, C. E., Mosquera, J. M., Fontugne, J., Beltran,
H., Sboner,
A., Garraway, L. A., Rubin, M. A., and Demichelis, F. (2014). Unraveling the clonal hierarchy of
somatic genomic aberrations. Genome biology 15, 439.

Carreira, S., Romanel, A., Goodall, J., Grist, E., Ferraldeschi, R., Miranda, S., Prandi, D.,
Lorente, D.,
Frenel, J. S., Pezaro, C., et al. (2014). Tumor clone dynamics in lethal prostate cancer. Science
translational medicine 6, 254ra125.

Beltran, H., Eng, K., Mosquera, J. M., Sigaras, A., Romanel, A., Rennert, H., Kossai, M.,
Pauli, C.,
Faltas, B., Fontugne, J., et al. (2015). Whole-Exome Sequencing of Metastatic Cancer and
Biomarkers
of Treatment Response. JAMA Oncol 1, 466-474.

Faltas, B. M., Prandi, D., Tagawa, S. T., Molina, A. M., Nanus, D. M., Sternberg, C.,
Rosenberg, J.,
Mosquera, J. M., Robinson, B., Elemento, O., et al. (2016). Clonal evolution of
chemotherapyresistant
urothelial carcinoma. Nature genetics 48, 1490-1499.
}
\keyword{
Purity, Ploidy, Clonality, Cancer Genomics, Allele specific Analysis
}
\examples{
\donttest{
###############
###############
## Diploid tumor sample

## Load example data
seg_tb <- read.table(system.file("sample1.seg", package = "CLONETv2"),header = TRUE, as.is=TRUE)
pileup_tumor <- read.table(
  gzfile(system.file("sample1_tumor_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
pileup_normal <- read.table(
  gzfile(system.file("sample1_normal_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
snv_reads <- read.table(system.file("sample1_snv_read_count.tsv", package = "CLONETv2"),
  header = TRUE, as.is=TRUE, comment.char = "", check.names = FALSE, na.strings = "-")

## Compute beta table with default parameters
bt <- compute_beta_table(seg_tb, pileup_tumor, pileup_normal)

## Compute ploidy table with default parameters
pl_table <- compute_ploidy(bt)

## Compute admixture table with default parameters (admixture= 1-tumor_purity)
adm_table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl_table)

## Check ploidy and admixture estimates
check_plot <-  check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)
print(check_plot)

## Compute clonality table with default parameters
scna_clonality_table <- compute_scna_clonality_table(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)

## Compute allele specific scna
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
  ploidy_table = pl_table, admixture_table = adm_table)

## Compute snvs colonality
sample_id <- "sample1"
snv_clonality_table <- compute_snv_clonality(sample_id = sample_id, snv_read_count = snv_reads,
  beta_table = bt, ploidy_table = pl_table, admixture_table = adm_table)


###############
###############
## Aneuploid tumor sample

## Load example data
seg_tb <- read.table(system.file("sample2.seg", package = "CLONETv2"),header = TRUE, as.is=TRUE)
pileup_tumor <- read.table(
  gzfile(system.file("sample2_tumor_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
pileup_normal <- read.table(
  gzfile(system.file("sample2_normal_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
snv_reads <- read.table(system.file("sample2_snv_read_count.tsv", package = "CLONETv2"),
  header = TRUE, as.is=TRUE, comment.char = "", check.names = FALSE, na.strings = "-")

## Compute beta table with default parameters
bt <- compute_beta_table(seg_tb, pileup_tumor, pileup_normal)

## Compute ploidy table with default parameters
pl_table <- compute_ploidy(bt)

## Compute admixture table with default parameters (admixture= 1-tumor_purity)
adm_table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl_table)

## Check ploidy and admixture estimates
check_plot <-  check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)
print(check_plot)

## Compute clonality table with default parameters
scna_clonality_table <- compute_scna_clonality_table(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)

## Compute allele specific scna
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
  ploidy_table = pl_table, admixture_table = adm_table)

## Compute snvs colonality
sample_id <- "sample2"
snv_clonality_table <- compute_snv_clonality(sample_id = sample_id, snv_read_count = snv_reads,
  beta_table = bt, ploidy_table = pl_table, admixture_table = adm_table)


###############
###############
## Aneuploidy tumor sample with problematic ploidy estimate

## Load example data
seg_tb <- read.table(system.file("sample3.seg", package = "CLONETv2"),header = TRUE, as.is=TRUE)
pileup_tumor <- read.table(
  gzfile(system.file("sample3_tumor_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
pileup_normal <- read.table(
  gzfile(system.file("sample3_normal_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)

## Compute beta table with default parameters
bt <- compute_beta_table(seg_tb, pileup_tumor, pileup_normal)

## Compute ploidy table with default parameters
pl_table <- compute_ploidy(bt)

## Compute admixture table with default parameters (admixture= 1-tumor_purity)
adm_table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl_table)

## Check ploidy and admixture estimates
check_plot <-  check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)
print(check_plot)
## Observed data (gray points) does not fit with expcted positions (Red circles)


###############
###############
## Tumor sample with problem in the segmented input data

## Load example data
seg_tb <- read.table(system.file("sample4.seg", package = "CLONETv2"),header = TRUE, as.is=TRUE)
pileup_tumor <- read.table(
  gzfile(system.file("sample4_tumor_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
pileup_normal <- read.table(
  gzfile(system.file("sample4_normal_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)

## Compute beta table with default parameters
bt <- compute_beta_table(seg_tb, pileup_tumor, pileup_normal)

## Compute ploidy table with default parameters
pl_table <- compute_ploidy(bt)

## Compute admixture table with default parameters (admixture= 1-tumor_purity)
adm_table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl_table)

## Check ploidy and admixture estimates
check_plot <-  check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)
print(check_plot)
## CLONETv2 does not provide an estimate of the DNA admixture because
## (LogR, beta) data does not fit any CLONETv2 model


###############
###############
## Diploid tumor sample with subclonal hemizygous and homozygous deletions

## Load example data
seg_tb <- read.table(system.file("sample5.seg", package = "CLONETv2"),header = TRUE, as.is=TRUE)
pileup_tumor <- read.table(
  gzfile(system.file("sample5_tumor_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
pileup_normal <- read.table(
  gzfile(system.file("sample5_normal_pileup.tsv.gz", package = "CLONETv2")),
  header = TRUE, as.is=TRUE)
snv_reads <- read.table(system.file("sample5_snv_read_count.tsv", package = "CLONETv2"),
  header = TRUE, as.is=TRUE, comment.char = "", check.names = FALSE, na.strings = "-")

## Compute beta table with default parameters
bt <- compute_beta_table(seg_tb, pileup_tumor, pileup_normal)

## Compute ploidy table with default parameters
pl_table <- compute_ploidy(bt)

## Compute admixture table with default parameters (admixture= 1-tumor_purity)
adm_table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl_table)

## Check ploidy and admixture estimates
check_plot <-  check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)
print(check_plot)

## Compute clonality table with default parameters
scna_clonality_table <- compute_scna_clonality_table(beta_table = bt, ploidy_table = pl_table,
  admixture_table = adm_table)

## Compute allele specific scna
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
  ploidy_table = pl_table, admixture_table = adm_table)

## Compute snvs colonality
sample_id <- "sample5"
snv_clonality_table <- compute_snv_clonality(sample_id = sample_id, snv_read_count = snv_reads,
  beta_table = bt, ploidy_table = pl_table, admixture_table = adm_table)

}
}
