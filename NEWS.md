# CLONETv2 Version: 2.2.1

* Set lower bound for corrected CN ratio (this solves issues for homozygous calls in noisy samples).

# CLONETv2 Version: 2.2.0

* better support for WGS data. Functions `compute_ploidy` and `compute_dna_admixture` have paramter `library_type` to specify if sequencing data
is targeted (WES) or whole genome (WGS).

# CLONETv2 Version: 2.1.0

* function `compute_beta_table` has now optional parameter `plot_stats`. 
If `plot_stats` is set to `TRUE`, the function prints statistics about the computation of beta table.

# CLONETv2 Version: 2.0.0

* Package version of CLONET. 
  Previous versions of CLONET were released as stand-alone tools. 
  Details about older CLONET versions are available at bitbucket.org/deid00/clonet/.
