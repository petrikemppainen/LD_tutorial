sortmap.internal <- function (chrom, map, delta = 1) 
{
  chnum <- chrom.char2num(chrom)
  ix <- order(chnum, map)
  map <- map[ix]
  off <- c(0, map[1:(length(map) - 1)])
  off <- map - off
  off[which(off <= 0)] <- delta
  cummap <- cumsum(off)
  out <- list()
  out$ix <- ix
  out$cummap <- cummap
  out$chnum <- chnum
  out
}

chicken = c(4415676, 3365997, 2504917, 2054347, 1346534, 798288, 831597, 
            674280, 542169, 460023, 455071, 448995, 414335, 351094, 287177, 
            246666, 248691, 224610, 317560, 154391, 106453, 130310, 
            141338, 65492, 119507, 127384, 111855, 37135) ## exluding chromosomes 16 and 29:32
human <- c(1731870, 1684824, 1379447, 1323231, 1262874, 1188215, 1108493, 
           1009660, 962745, 930764, 939733, 927132, 795578, 744652, 709504, 
           628440, 579182, 559118, 407775, 448307, 324939, 353520)

tit <- c(579, 415, 700, 596, 356, 103, 346, 177, 176, 134, 130, 148,
         135, 152, 117, 126, 173, 96, 93, 97, 155, 308)



## scale chromosomes to same size
chicken <- round(chicken/(sum(chicken)/20e+06))
human <- round(rev(human)/(sum(human)/1e+06))
tit <- round(rev(tit)/(sum(tit)/1e+04))


mat2el <- function(LDmat){
  
  ############################
  # transforms your upper diagonal pairwise r-squared matrix to an edge list
  ############################
  
  from <- rownames(LDmat)[row(LDmat)[upper.tri(LDmat)]]
  to <- colnames(LDmat)[col(LDmat)[upper.tri(LDmat)]]
  weight <- LDmat[upper.tri(LDmat)]
  el <- cbind(from,to,weight)
}

get_LD_sim_data <- function(parameters){
  ## update chromosome sizes 
  parameters$u <- ((4*parameters$nLoci/(2*parameters$Ne))/sum(1/(1:(parameters$sampleSize))))/sum(parameters$chrom_sizes)
  #parameters$r <- ((4*parameters$nLoci/(2*parameters$Ne))/sum(1/(1:(parameters$sampleSize))))/sum(parameters$chrom_sizes)*10
  
  path_name <- paste0(parameters$folder, '/', parameters$folder)
  
  ## set folders
  
  system(paste0("mkdir ", parameters$folder))
  
  
  ## create par file
  if(!is.na(parameters$m)){
    system(paste0("echo '//Number of population samples (demes)
                  2
                  //Population effective sizes (number of genes)
                  ", parameters$Ne/2,
                  "
                  ", parameters$Ne/2,
                  "
                  //Samples sizes
                  ", parameters$sampleSize/2,
                  "
                  ", parameters$sampleSize/2,
                  "
                  //Growth rates
                  0
                  0
                  //Number of migration matrices : 0 implies no migration between demes
                  1
                  // migration matrix
                  0 ", parameters$m/(parameters$Ne/2),
                  "
                  ", parameters$m/(parameters$Ne/2), " 0
                  //historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
                  0 historical event
                  //Number of independent loci [chromosome]
                  ", length(parameters$chrom_sizes), " 1
                  ",
                  
                  paste(lapply(parameters$chrom_sizes, function(i){
                    paste0("//Per chromosome: Number of linkage blocks
                           1
                           //per Block: data type, num loci, rec. rate and mut rate + optional parameters
                           DNA ", i, " ", parameters$r, " ", parameters$u, " 0.33")
                  }), collapse="\n"), "' > ", path_name, ".par")
                    )
  }
  if(!is.na(parameters$admixture_nGen)){
    system(paste0("echo '//Number of population samples (demes)
                  2
                  //Population effective sizes (number of genes)
                  ", parameters$Ne/2,
                  "
                  ", parameters$Ne/2,
                  "
                  //Samples sizes
                  ", parameters$sampleSize/2,
                  "
                  ", parameters$sampleSize/2,
                  "
                  //Growth rates
                  0
                  0
                  //Number of migration matrices : 0 implies no migration between demes
                  2
                  //migration matrix
                  0	0.5
                  0.5	0	
                  //migration matrix
                  0	0
                  0	0	
                  //historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
                  3 historical event
                  ", parameters$div_time," 0 1 1 2 0 0
                  ", parameters$admixture_nGen,
                  " 0 1 0.5 1 0 1
                  ", parameters$admixture_nGen,
                  " 1 0 0.5 1 0 1
                  //Number of independent loci [chromosome]
                  ", length(parameters$chrom_sizes), " 1
                  ",
                  
                  paste(lapply(parameters$chrom_sizes, function(i){
                    paste0("//Per chromosome: Number of linkage blocks
                           1
                           //per Block: data type, num loci, rec. rate and mut rate + optional parameters
                           DNA ", i, " ", parameters$r, " ", parameters$u, " 0.33")
                  }), collapse="\n"), "' > ", path_name, ".par")
                    )
  }
  
  if(all(is.na(parameters$admixture_nGen), is.na(parameters$m))){
    system(paste0("echo '//Number of population samples (demes)
                1
                  //Population effective sizes (number of genes)
                  ", parameters$Ne,
                  "
                  //Samples sizes
                  ", parameters$sampleSize,
                  "
                  //Growth rates
                  0
                  //Number of migration matrices : 0 implies no migration between demes
                  0
                  //historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
                  0 historical event
                  //Number of independent loci [chromosome]
                  ", length(parameters$chrom_sizes), " 1
                  ",
                  
                  paste(lapply(parameters$chrom_sizes, function(i){
                    paste0("//Per chromosome: Number of linkage blocks
                           1
                           //per Block: data type, num loci, rec. rate and mut rate + optional parameters
                           DNA ", i, " ", parameters$r, " ", parameters$u, " 0.33")
                  }), collapse="\n"), "' > ", path_name, ".par")
                    )  
  }
    
  
  
  
  system(paste0("./fsc_unix -i ./", parameters$folder, "/", parameters$folder, ".par -n 1"))
  
  
  system(paste0('grep -n "Total number of polymorphi" ', parameters$folder, '/', parameters$folder, '_1_1.arp | cut -d ":" -f3 > ', parameters$folder, '/', 'tot_pol_sites.txt'))
  
  temp <- scan(paste0(parameters$folder, "/", "tot_pol_sites.txt"), quiet = TRUE)
  print(paste0("Total number of polymorphic sites are: ", temp))
  
  system(paste0('grep -n "{" ', path_name, '_1_1.arp | cut -d ":" -f1  > ', parameters$folder, '/start2.txt'))
  system(paste0('grep -n "}" ', path_name, '_1_1.arp | cut -d ":" -f1  > ', parameters$folder, '/end2.txt'))
  
  a2 <- scan(paste0(parameters$folder, '/start2.txt'), quiet = TRUE)
  a2 <- a2[-length(a2)]
  b2 <- scan(paste0(parameters$folder, '/end2.txt'), quiet = TRUE)
  b2 <- b2[-length(b2)]
  nSamples <- b2-a2-2
  
  ## get polymorphic positions
  SNPs <- do.call(rbind, lapply(1:length(nSamples), function(i){
    system(paste0("sed -n '", a2[i]+1, ",", b2[i]-2, "p' ", path_name, "_1_1.arp > ", parameters$folder, "/SNPs.txt"))
    test <- read.table(paste0(parameters$folder, "/SNPs.txt"))[,-(1:2)]
    
    test <- apply(as.matrix(test), 1, function(x) strsplit(as.vector(x), ""))
    do.call(rbind, lapply(split(1:nSamples[i], rep(1:(nSamples[i]/2), each=2)), function(x){
      
      temp <- unlist(test[x])
      temp2 <- 1:length(temp)
      temp2[(1:(length(temp)/2))*2-1] <- (1:(length(temp)/2))
      temp2[(1:(length(temp)/2))*2] <- ((length(temp)/2)+1):length(temp)
      temp[temp2]
      
    }))
  }))

  system(paste0('grep -A1 "polymorphic positions on chromosome" ', path_name, '_1_1.arp | grep -v "polymorphic positions on chromosome" | grep -v "-" | cut -d "#" -f2 | paste -d, -s - > ', parameters$folder, '/polyorphic_positions.txt'))
  
  positions <- unlist(read.csv(paste0(parameters$folder, '/polyorphic_positions.txt'), header = F))
  
  system(paste0('grep "polymorphic positions on" ', path_name, '_1_1.arp | cut -d "#" -f2 > ',  parameters$folder, '/chromosomes.txt'))
  
  chromosomes <- fread(paste0(parameters$folder, "/chromosomes.txt"))
  
  chromosomes <- rep(chromosomes[,(V6)],chromosomes[,(V1)])
  
  
  loci <- split(1:ncol(SNPs), rep(1:(ncol(SNPs)/2), each=2))
  
  temp <- sapply(loci, function(x) length(unique(as.vector(SNPs[,x]))))
  keep <- which(temp==2)
  SNPs <- SNPs[,unlist(loci[keep])]
  
  chromosomes <- chromosomes[keep]
  positions <- positions[keep]
  
  nrowSNPS <- nrow(SNPs)
  SNPs <- cbind(rep(1, nrow(SNPs)), 
               1:nrow(SNPs),
               matrix(0, nrow(SNPs), 2),
               rep(1, nrow(SNPs)),
               rep(-9, nrow(SNPs)),
               SNPs)
  write.table(SNPs, file=paste0(path_name, '.ped'), quote = F, col.names = F, row.names = F)
  rm(SNPs)
  ## map file ##
  temp <- split(positions, as.factor(chromosomes))
  #cM <- unlist(lapply(temp, function(x) x/chromosome_length))
  
  map <- cbind(chromosomes, 1:length(chromosomes), positions, positions)
  write.table(map, file=paste0(path_name, '.map'), quote = F, col.names = F, row.names = F)
  
  
  ## pheno file ##
  
  morph <- cbind(rep(1, nrowSNPS),1:nrowSNPS, 1:nrowSNPS)
  colnames(morph) <- c("sex", "id", "phenotype")
  write.table(morph, file=paste0(path_name, '_morph.txt'), quote = F, col.names = T, row.names = F)
  
  
  ## do some filtering in plink_unix already here ##
  
  # filter by maf
  
  if(length(parameters$chrom_sizes>22)){
    system(paste0("./plink_unix --file ", path_name, " --make-bed  --maf ", parameters$maf, " --out ", path_name, " --noweb"), ignore.stdout=FALSE)
    system(paste0("cd ", parameters$folder, "\n", "../plink_unix  --bfile ", parameters$folder, "  --write-snplist --noweb"), ignore.stdout=FALSE)
    if(parameters$subsample){
      # read snp names
      snplist <- scan(paste0('./', parameters$folder, "/plink.snplist"), quiet = TRUE)
      if(length(snplist)>parameters$nLoci){
        loci <- sort(sample(snplist, parameters$nLoci))
        write.table(loci, file=paste0(path_name, ".snplist"), quote = FALSE, col.names = F, row.names = F)
        system(paste0("./plink_unix --bfile ", path_name, "    --make-bed --extract ", path_name, ".snplist --out ", path_name, " --noweb"), ignore.stdout=TRUE)
      }
    }
    system(paste0("./plink_unix --bfile ", path_name, "  --recode --tab --out ", path_name, " --noweb"), ignore.stdout=TRUE)  
    
  }else{
    system(paste0("./plink_unix --file ", path_name, " --chr-set 32 --make-bed  --maf ", parameters$maf, " --out ", path_name, " --noweb"), ignore.stdout=TRUE)
    system(paste0("cd ", parameters$folder, "\n", "../plink_unix   --chr-set 32 --bfile ", parameters$folder, "  --write-snplist --noweb"), ignore.stdout=TRUE)
    if(parameters$subsample){
      # read snp names
      snplist <- scan(paste0('./', parameters$folder, "/plink.snplist"), quiet = TRUE)
      loci <- sort(sample(snplist, parameters$nLoci))
      write.table(loci, file=paste0(path_name, ".snplist"), quote = FALSE, col.names = F, row.names = F)
      system(paste0("./plink_unix --bfile ", path_name, "   --chr-set 32 --make-bed --extract ", path_name, ".snplist --out ", path_name, " --noweb"), ignore.stdout=TRUE)
    }
    system(paste0("./plink_unix --bfile ", path_name, " --chr-set 32 --recode --tab --out ", path_name, " --noweb"), ignore.stdout=TRUE)  
    
  }
  
  
  
  # subsample loci
  
  ### prepare data ###
  convert.snp.ped(paste0(path_name, '.ped'), paste0(path_name, '.map'), paste0(path_name, ".raw"), mapHasHeaderLine=FALSE)
  
  geno <- paste0(path_name, ".raw") #assign the .raw file to genotypes
  
  #load phenotype file
  pheno <- paste0(path_name, '_morph.txt')
  
  #load genotype and phenotye files into a gwaa.data object useable in genable. makemap =T orders the SNP markers.
  gwaa_data <- load.gwaa.data(pheno,geno,makemap=F)
  
  if(parameters$LD){
    r2 <- t(r2fast(gwaa_data))
    r2[upper.tri(r2)] <- NA
    
    r2_el <- mat2el(t(r2))
    temp <- apply(r2_el[,1:2], 2, as.numeric)
    
    Dists <- abs(positions[temp[,1]]-positions[temp[,2]])
    Chr <- data.table(chromosomes[temp[,1]], chromosomes[temp[,2]])
    
    genind <- df2genind(fread(paste0(path_name, '.ped'))[,-c(1:6)], ncode = 1, ploidy = 1)
    
    list(LDmat=r2, LD_dist=data.table(dist=Dists, r2=as.numeric(r2_el[,3]))[apply(Chr, 1, function(x) x[1]==x[2]),], pos=positions[loci], genind=genind, gwaa=gwaa_data)
  }else{
    gwaa_data
  }
}


getPhenotypicData <- function(data_gwas, nCausal=parameters$nCausal, h=parameters$h, mu=parameters$mu, clustered=parameters$clustered){
  
  ### get_phenotype 
  if(clustered=="Y"){
    n <- data_gwas@gtdata@nsnps
    dist1 <- 10
    p0 <- 0.75
    causal_loci <- sample(1:n, 1)
    for(i in 1:(nCausal-1)){
      temp <- ifelse(rbinom(1,1, p0)==1, sample(((causal_loci[1]-dist1):(causal_loci[1]+dist1))[-(dist1+1)], 1), sample(1:n, 1))
      while(temp<1 | temp>n | temp %in% causal_loci){
        temp <- ifelse(rbinom(1,1, p0)==1, sample(((causal_loci[1]-dist1):(causal_loci[1]+dist1))[-(dist1+1)], 1), sample(1:n, 1))
      }
      causal_loci <- c(temp, causal_loci)
    }
  }else{
    causal_loci <- sample(1:data_gwas@gtdata@nsnps, nCausal)
  }
  temp <- as.genotype(data_gwas@gtdata[1:data_gwas@gtdata@nids, causal_loci])
  
  data <- lapply(1:ncol(temp), function(x) do.call(rbind, strsplit(as.vector(temp[,x]), "/", fixed=T)))
  
  causal_alleles <- sapply(data, function(y){
    sample(as.vector(unique(as.vector(y))), 1)
  })
  
  # get p's for causal alleles
  
  if(h==0){
    phenotype <- phenotype_new <- rnorm(nrow(data[[1]]))
    u <- NA
  }else{
    p <- sapply(1:length(data), function(y){
      as.vector((table(as.vector(data[[y]]) == causal_alleles[y])/(nrow(data[[y]])*2))[1])
    })
    
    if(mu=="gamma"){
      a <- qgamma(ppoints(nCausal), shape=0.5, rate=0.2)
      temp <- a-min(a)
      u <- temp/max(temp)
      u <- u/sum(u)
    }
    
    if(mu=="normal"){
      a <- qnorm(ppoints(nCausal))
      temp <- a-min(a)
      u <- temp/max(temp)
      u <- u/sum(u)
    }
    
    #nCausal <- 125
    if(mu=="equal"){
      a <- rep(1, nCausal)  
      u <- a/sum(a)
    }
    
    
    # get allele counts for causal alleles for each locus
    #z <- 2
    phenotype <- sapply(1:nrow(data[[1]]), function(z){
      
      x <- sapply(1:length(causal_loci), function(y){
        length(which(data[[y]][z,]==causal_alleles[y]))
      })
      
      w = (x - 2*p)/sqrt(2*p*(1-p))
      sum(w*u)
    })
    
    phenotype_new <-  phenotype + rnorm(phenotype, sd=(abs(var(phenotype)*(1-1/(h+1e-6))))^0.5)
    
  }
  
  temp <- phenotype_new-min(phenotype_new)
  phenotype_new <- temp/max(temp)
  
  out <- list(phenotype, phenotype_new,  causal_loci, u)
  names(out) <- c('phenotype', 'phenotype_new',  'causal_loci', "weights")
  return(out)
  
}


plot_PCA <- function(sim_data){
  X <- scaleGen(sim_data$genind, NA.method="mean")
  pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
  
  
  plot(pca1$li[,1:2], cex=1, pch=20, col=ifelse(find.clusters(sim_data$genind, n.clust =  2, n.pca = 10)$grp==2, 'cornflowerblue', 'gold3'),
       ylab=paste("Axis 2", paste("(", round(pca1$eig[2]/sum(pca1$eig)*100), " %)", sep="" )), 
       xlab=paste("Axis 1", paste("(", round(pca1$eig[1]/sum(pca1$eig)*100), " %)", sep="" )),
       main='Principal component analysis')
}

plot_dist_vs_r2 <- function(sim_data){
  # suppressWarnings(ggplot(sim_data$LD_dist, aes(dist, r2)) +
  #                    geom_point() +
  #                    theme_bw() +
  #                    theme(panel.grid.major = element_blank(),
  #                          panel.grid.minor  = element_blank(),
  #                          strip.background = element_blank(),
  #                          legend.position = c(.01, .99),
  #                          legend.justification = c(0,1),
  #                          legend.text.align = 0
  #                    ))
  plot(sim_data$LD_dist)
     
}

LDnaRaw <- function(LDmat, mc.cores=NULL){
  if(is.na(LDmat[2,1])) LDmat <- t(LDmat)
  LDmat1 <- round(LDmat, 2)
  temp2 <- hclust(as.dist(1-LDmat1), method="single")
  tree <- as.phylo(temp2)
  tree <- di2multi(tree)
  Ntips <- length(tree$tip.label)
  
  stats.fun <- function(x){
    tree.temp <- extract.clade(tree, x)
    loci <- colnames(LDmat1) %in% tree.temp$tip.label
    LDmat.temp <- LDmat1[loci, loci]
    
    tot.d <- tree$edge.length[which(tree$edge[,1]==x)[1]]
    temp <- x
    while(tree$edge[,2][tree$edge[,1]==temp][1]>x){
      temp <- tree$edge[,2][tree$edge[,1]==temp][1]
      tot.d <- tot.d+tree$edge.length[which(tree$edge[,1]==temp)[1]]
    }
    
    out <- list(nV=length(which(loci)), 
                nE=length(which(LDmat.temp >= 1-2*tot.d)),
                tot.d=tot.d,
                loci=loci)
  }
  
  if(!is.null(mc.cores)){out <- parallel::mclapply(Ntips:length(tree$edge.length)+1, stats.fun, mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{out <- lapply(Ntips:length(tree$edge.length)+1, stats.fun)}
  
  
  nV <- sapply(out, function(x){x$nV})
  nE <- sapply(out, function(x){x$nE})
  tot.d <- sapply(out, function(x){x$tot.d})
  clusterfile <- do.call('cbind', lapply(out, function(x){x$loci}))
  
  tot.d2 <- 1-2*tot.d
  temp <- paste(tree$Nnode+1-c(1:tree$Nnode), round(tot.d2, 3)[order(round(tot.d2, 3))], sep="_")
  c.names <- temp[order(order(round(tot.d2, 3)))]
  colnames(clusterfile) <- c.names
  rownames(clusterfile) <- tree$tip.label
  
  LDmat3 <- LDmat
  LDmat3[upper.tri(LDmat3)] <- t(LDmat3)[upper.tri(LDmat3)]
  
  lambda.fun <- function(x){
    p <- tree$edge[,1][tree$edge[,2] == Ntips+x]-Ntips
    temp <- LDmat3[clusterfile[,x],clusterfile[,x]]
    temp <- temp[lower.tri(temp)]
    LD1 <- median(temp, na.rm=TRUE)
    temp <- c(temp, LDmat3[clusterfile[,p] + clusterfile[,x]==1,clusterfile[,x]])
    LD2 <- median(temp, na.rm=TRUE)
    list(lambda=(LD1 - LD2)*length(which(clusterfile[,x])))
  }
  
  if(!is.null(mc.cores)){out <- parallel::mclapply(2:ncol(clusterfile), lambda.fun, mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{out <- lapply(2:ncol(clusterfile), lambda.fun)}
  
  lambda <- sapply(out, function(x){x$lambda})
  names <- c(tree$tip.label, c.names)
  stats <- data.frame(nV, nE, lambda=c(0,lambda))
  rownames(stats) <- c.names
  
  d <- tree$edge
  d[,2] <- names[tree$edge[,1]]
  d[,1] <- names[tree$edge[,2]]
  d <- cbind(d, c(round(tree$edge.length*2, 3)))
  d <- rbind(c(d[1,2], "root", min(round(tot.d2, 3))), d)
  d <- data.frame(d , stats[match(d[,1], c.names), ])
  names(d) <- c("cluster", "parent_cluster", "distance", "nV", "nE", "lambda")
  d$nV[is.na(d$nV)] <- 1
  d$nE[is.na(d$nE)] <- 0
  d$lambda[is.na(d$lambda)] <- 0
  d <- d[rev(order(d$nV)),]
  rownames(d) <- 1:nrow(d)
  out <- list(clusterfile, d)
  names(out) <- c("clusterfile", "stats")
  return(out)
}

extractClusters <- function(ldna, min.edges=20, phi=2, lambda.lim=NULL, rm.COCs=TRUE, extract=TRUE, plot.tree=TRUE, plot.graph=TRUE){
  # Get file for tree and clusters above min.edges and their lambda values
  tree <- clusterPhylo(ldna, min.edges)
  if(extract){
    clusters <- tree$tip.label
    lambda_new <- ldna$stats$lambda[ldna$stats$nE >= min.edges][-1]
    names(lambda_new) <- as.vector(ldna$stats$cluster[ldna$stats$nE >= min.edges][-1])
    if(min.edges==0) lambda_new <- lambda_new[ldna$stats$nV!=1][-length(lambda_new[ldna$stats$nV!=1])]
    
    # Get thresholds. If lambda.lim is given, this takes precedence.
    if(!is.null(lambda.lim)){threshold <- lambda.lim
    }else{
      if(!is.null(phi)){
        threshold <- median(lambda_new)+mad(lambda_new, constant=phi)
      }  
    }
    
    #get outlier clusters
    clusters.out <- names(lambda_new)[which(lambda_new >= threshold)]
    if(identical(clusters.out, character(0))) stop("No outlier clusters, please decrease phi or lambda.lim")
    
    
    # get SOCs and COCs
    temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% clusters.out]
    if(is.matrix(temp)){
      nested <- matrix(NA, ncol(temp), ncol(temp))
      for(i in 1:ncol(temp)){
        for(j in 1:ncol(temp)){
          if(i!=j & any(apply(cbind(temp[,i], temp[,j]),1, function(x) x[1]==TRUE & x[2]==TRUE))){
            nested[i,j] <- "COC"
          }
        }
      }
      nested[is.na(nested)] <- "SOC"
      nested[lower.tri(nested)] <- NA    
      COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
      SOCs <- colnames(temp)[!colnames(temp) %in% COCs]
    }else{
      SOCs <- clusters.out
      COCs <- NA
    }
    
    if(plot.graph){
      col <- lambda_ord <- lambda_new[order(lambda_new)]
      #if(length(clusters.out)>1) 
      col[which(names(lambda_ord) %in% COCs)] <- "blue"
      col[which(names(lambda_ord) %in% SOCs)] <- "red"
      col[which(!col %in% c("red","blue"))] <- "black"
      plot(lambda_ord, col=col, ylab=expression(lambda))
      lines(c(1,length(lambda_ord)),c(threshold,threshold), col="red", lty=2) 
      text(length(lambda_ord)/3, y=threshold, as.expression(bquote(lambda[lim]*"="*.(signif(threshold,3)))), pos=3, adj=c(0,0))
      col.text <- lambda_ord
      if(rm.COCs==FALSE){
        if(length(clusters.out)>1) col.text[which(names(lambda_ord) %in% COCs)] <- "blue"
      } 
      col.text[which(names(lambda_ord) %in% SOCs)] <- "red"
      col.text[!col.text %in% c("blue","red")] <- "#00000000"
      text(lambda_ord, names(lambda_ord), pos=2, cex=0.75, col=col.text)
      if(!is.null(lambda.lim)){
        title(main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }else{
        title(main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }
    }
    
    # plot tree
    if(plot.tree){
      col <- rep("grey", length(tree$edge))
      if(rm.COCs==FALSE){
        distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% clusters.out)]
        clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% clusters.out)]]
        keep.col <- clusters.temp[distances > 0]
        col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "blue"
        col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% clusters.out[!clusters.out %in% keep.col])]] <- "blue"
      }
      tree$edge[tree$edge[,2] %in% which(tree$tip.label %in% SOCs),]
      distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]
      clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]]
      keep.col <- clusters.temp[distances > 0]
      col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "red"
      col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% SOCs[!SOCs %in% keep.col])]] <- "red"
      col.tip <- rep("#00000000", length(tree$tip.label))
      if(rm.COCs==FALSE){
        if(length(clusters.out)>1) col.tip[tree$tip.label %in% clusters.out] <- "blue"
      }
      col.tip[tree$tip.label %in% SOCs] <- "black"
      plot(tree, show.tip.label=T, edge.width=3, edge.color=col, cex=1, tip.color=col.tip, root.edge=TRUE, underscore=T,x.lim=1)
      #if(min.edges==0){
      axis(1, at=c(0,(1:10)*0.1))
      #}else{
      #  temp <- as.vector(ldna$stats[ldna$stats$nE > min.edges,2][-1])
      #  x <- round(10*max(as.numeric(do.call('rbind', strsplit(temp, "_", fixed=TRUE))[,2])),0)+1
      #  if(x>10) x <- 10
      #  axis(1, at=c(0,(1:x)*0.1))
      #}
      if(!is.null(lambda.lim)){
        title(xlab="LD threshold", main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }else{
        title(xlab="LD threshold", main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }
    }
    
    if(length(clusters.out)>1){
      
      if(rm.COCs==FALSE){out <- clusters.out[order(-as.numeric(do.call('rbind', strsplit(clusters.out, "_"))[,2]))]
      }else{out <- SOCs[order(-as.numeric(do.call('rbind', strsplit(SOCs, "_"))[,2]))]}
      
      temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% out]
      if(is.matrix(temp)){
        temp <- temp[,order(as.numeric(do.call('rbind', strsplit(colnames(temp), "_", fixed=T))[,1]))]
        loci <- apply(temp, 2, function(x) rownames(temp)[x])          
      }else{
        loci <- list(names(temp)[temp])
      }
    }else{
      temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% clusters.out]
      loci <- list(names(temp)[temp])
    }
    
    #rm(out)
    #rm(loci)
    return(loci)
    
  }else{
    plot(tree, edge.width=3,show.tip.label=F,edge.color="grey", cex=1,  root.edge=TRUE, x.lim=1)
    axis(1, at=c(0,(1:10)*0.1))
    title(xlab="LD threshold", main=as.expression(bquote("|E|"[min]*plain("=")* .(min.edges))))
  }
}

clusterPhylo <-  function(ldna, min.edges=0){
  d <- ldna$stats[ldna$stats$nE>=min.edges,]
  
  d$times<-rep(0,dim(d)[1])
  d$offspring<-as.list(rep(NA,dim(d)[1]))
  d$terminal<-rep(0,dim(d)[1])
  d$ancestor<-rep(0,dim(d)[1])
  
  row<-1
  rowOld<-0
  newick<-character()
  
  repeat{
    cluster<-as.character(d$cluster[row])
    d$offspring[[row]]<-which(as.character(d$parent_cluster)==cluster)
    offspNo<-sum(!is.na(d$offspring[[row]]))
    
    if(offspNo==0){d$terminal[row]<-1}
    
    if(d$terminal[row]==1){
      d$ancestor[row]<-rowOld
      if(d$ancestor[row]==0){
        newick<-paste("(",cluster,":0,:0);", sep="")
        break}
      newick<-paste(newick,cluster,":",d$distance[row], sep="")
      row<-d$ancestor[row]
      next}
    
    if(d$times[row]==0){
      newick<-paste(newick,"(", sep="")
      if(offspNo>1){newick<-paste(newick,"(", sep="")}
      d$ancestor[row]<-rowOld
      rowOld<-row
      d$times[row]<-d$times[row]+1
      row<-d$offspring[[row]][1]
      next}
    
    if(d$times[row]>0 & d$times[row]<offspNo){
      newick<-paste(newick,",", sep="")
      if((offspNo-d$times[row])>1){newick<-paste(newick,"(", sep="")}
      rowOld<-row
      d$times[row]<-d$times[row]+1
      row<-d$offspring[[row]][(d$times[row])]
      next}
    
    if(d$times[row]>0 & d$times[row]==offspNo){
      if(offspNo==1){
        newick<-paste(newick,",",cluster,":0):",d$distance[row], sep="")
      }else{
        newick<-paste(newick,paste(rep("):0",(offspNo-1)),sep="", collapse=""),",",cluster,":0):",d$distance[row], sep="", collapse="")            
      }
      if(d$ancestor[row]==0){
        newick<-paste(newick,";", sep="")
        break}
      row<-d$ancestor[row]
      next}
  }  
  tree <- read.tree(text=newick)
}

plotLDnetwork <- function(ldna, LDmat, option, threshold, clusters, summary,
                          exl=NULL, full.network=TRUE, include.parent=FALSE, after.merger=FALSE, graph.object=FALSE, col="grey", pos=NULL){
  if(is.na(LDmat[2,1])) LDmat <- t(LDmat)
  if(option==1) g <- option1(LDmat, threshold, exl, pos, col);
  if(option==2) option2(ldna, LDmat, clusters, summary, exl, full.network, include.parent, after.merger)
  if(option==1 && graph.object) return(g)
}

option1 <- function(LDmat, threshold, exl, pos, col){
  
  LDmat <- LDmat[!(rownames(LDmat) %in% exl),!(rownames(LDmat) %in% exl)]
  g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=T)
  if(!is.null(pos)){
    V(g)$color <- rgb(pos, max(pos)-pos, 0, maxColorValue = max(pos))  
    V(g)$frame.color <- V(g)$color
  }else{
    V(g)$color <- col
    V(g)$frame.color <- V(g)$color
  }
  
  E(g)$weight <- round(E(g)$weight, 2)
  g <- delete.edges(g, which(E(g)$weight<=threshold))
  g <- delete.vertices(g, which(degree(g) == 0))
  
  plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL, edge.width=1, vertex.label=NA)
  title(main=paste(" @", threshold, sep=""))
  return(g)
}

option2 <- function(ldna, LDmat, clusters, summary, exl, full.network, include.parent, after.merger){
  col <- as.vector(summary$Type)
  col[col=="COC"] <- "blue"
  col[col=="SOC"] <- "red"
  
  for(i in 1:nrow(summary)){
    threshold <- as.numeric(as.vector(summary$Merge.at[rownames(summary) == names(clusters)[[i]]]))
    option2raw(ldna, LDmat, exl, clusters[[i]], col=col[i],
               full.network, threshold, include.parent, after.merger,
               names(clusters)[[i]])
  }
}

option2raw <- function(ldna, LDmat, exl, loci, col, full.network, threshold, include.parent, after.merger, cluster.name){
  if(!is.null(exl)){
    LDmat <- LDmat[!(rownames(LDmat) %in% exl),!(rownames(LDmat) %in% exl)]
  }
  
  p <- as.vector(ldna$stats$parent_cluster[ldna$stats$cluster %in%  cluster.name])
  loci_p <- rownames(ldna$clusterfile)[ldna$clusterfile[,colnames(ldna$clusterfile) == p]]
  
  if(full.network==FALSE){
    if(include.parent==FALSE){
      LDmat <- LDmat[(rownames(LDmat) %in% loci), (rownames(LDmat) %in% loci)]
    }else{
      LDmat <- LDmat[(rownames(LDmat) %in% loci_p), (rownames(LDmat) %in% loci_p)]
    }
  }
  
  if(after.merger==TRUE) {
    p2 <- as.vector(ldna$stats$parent_cluster[ldna$stats$cluster %in%  p])
    if(p2=="root") {threshold <- threshold-0.01
    }else{threshold <- as.numeric(strsplit(p2, "_", fixed=T)[[1]][2])}
  }
  
  g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=T)
  E(g)$weight <- round(E(g)$weight, 2)
  g <- delete.edges(g, which(E(g)$weight<=threshold))
  g <- delete.vertices(g, which(degree(g) == 0))
  
  col.frame <- V(g)$name
  col.frame[which(col.frame %in% loci)] <- col
  col.frame[which(!col.frame %in% col)] <- "grey"
  
  plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL, vertex.color=col.frame, edge.width=1, vertex.label=NA, vertex.frame.color=col.frame)
  title(main=paste(cluster.name, " @", threshold, sep=""))
}

