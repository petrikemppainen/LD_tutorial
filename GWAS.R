### get simulated data from Simulate_data.R
####### GWAS ########

## simulation parameters for phenotypes
parameters$mu <- "equal" ## effect size distribution, "normal" or "gamma"
parameters$nCausal <- 5 ## number of causal loci
parameters$h <- 0.5 ## heritability
parameters$clustered <- "N" ## wheather oci are randomly distributed or clustered in the genome 
## the simulated data
data_gwas <- sim_data$gwaa

## get pheotypic data
Phen <-  getPhenotypicData(data_gwas)
Phen$weights
data_gwas@phdata$phenotype <-  Phen$phenotype_new

# GWAS and plotting 
gkin <- ibs(data_gwas,w="freq") ## typically a relatedness matrix is used to correct for relatedness
gkin[] <- 0 ## we assume all individuals are unrelated, we will check P-value inflation below
# here's the actual analysis, no need to know the details
h2ht <- polygenic(phenotype,kin=gkin, data_gwas)
gr <- qtscore(h2ht$pgres,data=data_gwas,clam=FALSE)
# now plot
plot.scan.gwaa(gr, main = "GWAS", type="h") # plot results
#points(as.numeric(rownames(gr@results))[Phen$causal_loci], -log10(gr@results$P1df)[Phen$causal_loci], col="red")

# genome wide significance level
abline(h=-log10(0.05/data_gwas@gtdata@nsnps), lty=2)

## these are the causal loci
results <- cbind(gr[Phen$causal_loci,c(1,2)], Pval = round(gr[Phen$causal_loci,c(15)], 3), effects = round(Phen$weights, 2))
results[order(results$Pval),]
