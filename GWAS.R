### get simulated data from Simulate_data.R
####### GWAS ########

## simulation parameters for phenotypes
parameters$mu <- "equal"       ## effect size distribution, "normal" or "gamma"
parameters$nCausal <- 5        ## number of causal loci
parameters$h <- 0.5            ## heritability
parameters$clustered <- "N"    ## wheather oci are randomly distributed or clustered in the genome 

## the simulated data on which the phenotypes are based on
data_gwas <- sim_data$gwaa

## get pheotypic data, chooses QTL among all loci in data set and simulated phenotypes based on the parameters above
Phen <-  getPhenotypicData(data_gwas)   
Phen$weights                            ## distribution of effect sizes
data_gwas@phdata$phenotype <-  Phen$phenotype_new    ## add phenotype to the simulated data

# GWAS and plotting 
gkin <- ibs(data_gwas,w="freq")                         ## typically a relatedness matrix is used to correct for relatedness
gkin[] <- 0                                             ## we assume all individuals are unrelated, we will check P-value inflation below

### here's the actual analysis, no need to know the details
h2ht <- polygenic(phenotype,kin=gkin, data_gwas)
gr <- qtscore(h2ht$pgres,data=data_gwas,clam=FALSE)

### now plot
plot.scan.gwaa(gr, main = "GWAS", type="h") # plot results

### genome wide significance level
abline(h=-log10(0.05/data_gwas@gtdata@nsnps), lty=2) ## conservative, but works as a ballpark

## these are the causal loci, are all singificant?
Phen$causal_loci                         ## the causal loci index i.e. the position when loci are numbered from 1 to number of loci

## extract results of causal loci from the GWAS
results_causal <- cbind(gr[Phen$causal_loci,c(1,2)], Pval = round(gr[Phen$causal_loci,c(15)], 3), effects = round(Phen$weights, 2))
results_causal <- results_causal[order(results_causal$Pval),]
results_causal

## all results, can you spot any "false positives"??
results_all <- data.table(gr[,])
results_all[order(P2df),][P2df<(0.05/data_gwas@gtdata@nsnps),]      ## showing only significant SNPs

results_causal$Position ## positions of the causal loci in cM


### P-value inlflation?
hist(results_all$P2df)                   ## are p-values uniformly distributed between 0 and 1?

## QQ-plot of -log10(P) aka PP-plot
obs_P <- -log10(results_all$P1df)        ## observed p-values
exp_P <- -log10(runif(results_all$P1df)) ## expected p-values
COL <- rep("black", length(obs_P))       ## set color for causal loci
COL[Phen$causal_loci] <- "red"
Ord_obs <- order(obs_P)
Ord_exp <- order(exp_P)
plot(exp_P[Ord_exp],obs_P[Ord_obs], col=COL[Ord_obs], main="PP-plot", pch=20) ## this is the QQ-plot

abline(0,1, col="red")                                             ## under the null-hypothesis all points are expected to be on this 1:1 line, and only causal loci should be above the line                           
abline(h=-log10(0.05/data_gwas@gtdata@nsnps), lty=2)               ## points above horisontal line are significant QTL after (conservative) correction for multiple testing


### go back and run analyses without "gkin[] <- 0 " to see what happens

### gennomic control ##
## p-value inflation is estiamted by lambda
lambda(gr)
estlambda(gr[, "P1df"], plot=TRUE)          ## lambda is the slope of the red line


obs_P <- -log10(results_all$Pc1df)          ## same as above but with corrected P-values
exp_P <- -log10(runif(results_all$Pc1df)) 
COL <- rep("black", length(obs_P))       
COL[Phen$causal_loci] <- "red"
Ord_obs <- order(obs_P)
Ord_exp <- order(exp_P)
plot(exp_P[Ord_exp],obs_P[Ord_obs], col=COL[Ord_obs], main="PP-plot", pch=20) ## this is the QQ-plot

abline(0,1, col="red")                                           
abline(h=-log10(0.05/data_gwas@gtdata@nsnps), lty=2)   

### Demonstration of what genomic control is??
Lambda <- lambda(gr)$estimate
P_GC1 <- -log10(results_all$Pc1df)
P_GC2 <- -log10(results_all$P1df)/Lambda ## our own correction
P_GC1==P_GC1


