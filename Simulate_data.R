#### set up, just run this once  ####
source("install_packages.R") ## this first step takes some time!! If you're running this through Rstudio cloud, you only need to do this once.
source("functions.R")
save.image() ## save the image so that you don't have to start from scratch

## run simulations from here
## choose parameters ##
parameters <- list(
  folder = "LD_sim",    # files are saved to this folder, change if you want to but not necessary, files will be overwritten, with a warning that you can ignore
  Ne = 2000,            # effective population size
  r = 0.0001,           # recombination rate
  sampleSize = 500,     # number of haploid genomes to simulate
  chrom_sizes = rep(1000, 3),     # chromosome number and size distribution. three chromosomes of 1000 bp each. can also use "human/10"
  maf = 0.05,           # minor allele frequency for main data set
  subsample = TRUE,     # use only nLoci in final data set, do not change if you don't want huge data sets
  nLoci = 500,          # number of loci to subsample
  m = NA,               # number migrants per generation, NA means no population structure. Assumes two populations in migration drif equilibrium
  admixture_nGen = NA,  # how many generations of admixture?, NA means no admixture. At least one of 'admixture_nGen' or 'm' must be NA
  LD = TRUE,            # estimate LD
  div_time = 5000       # if admixture_nGen != NA, div_time gives time of divergence. Start with div_time=Ne
)

## generate simulated data ##
sim_data <- get_LD_sim_data(parameters)


## Principal component analyses
plotPCA(sim_data$geotypes)    ## only makes sense when you have admixture

par(mfcol=c(2,2))
## plot LD vs physical distance, recombination rate is uniform so scales directly with cM ##
plot_dist_vs_r2(sim_data)

## plot LD clustering tree, no need to change parameters ##
ldna <- LDnaRaw(sim_data$LDmat)
extractClusters(ldna, min.edges = 10, plot.tree = TRUE, plot.graph = F, extract = F) ## This also exctracts clusters, but we only use this to plot the LD-clustering tree

## plot networks at given LD thresholds, change threshold to see different graphs
## Color indicates physical position on a chromosome, going from green (one chromosome end) to red (the other chromosome end)
plotLDnetwork(LDmat=sim_data$LDmat, option=1, threshold=0.3, col = rgb(sim_data$pos/max(sim_data$pos), 1-sim_data$pos/max(sim_data$pos),0))

