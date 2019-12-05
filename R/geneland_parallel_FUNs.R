library(Geneland)
library(coda)
library(ggplot2)
library(data.table)
library(parallel)

#filter.null.alleles
#miss.loc
#pop.size

#----------------------------------------------------------------------------#
#### Helper functions ####
#----------------------------------------------------------------------------#
#' Run Geneland in parallel
#' nrun Number of parallel runs
#' ncores Number of cores to use. If default (NULL) are automatically selected
#' model Whether "Correlated" or "Uncorrelated"
#' path Path to root directory where to store results
#' main.dir Directory name where to store output (within path)
#' spatial Whether to use the spatial model
#' GEN Input genotypes (for now only codominant are available)
#' COORD Coordinates
#' jitter Spatial uncertainty
#' nxdom,nydom Number of pixel for discretization of the spatial domain in the 
#'         horizontal and vertical direction respectively. If "auto" is automatically set
#'         otherwise a numeric vector of length=1
#' burnin Burnin
#' POPSMAX Maximum number of pops
#' POPINIT Initial number of pops
#' POPMIN Minimum number of pops
#' ITR Number of iteration
#' THN Thinning
#' nullMatrix A matrix with nindiv lines and nloc columns of 0 or 1. For each 
#'          individual, at each locus it says if the locus is genuinely missing 
#'          (no attempt to measure it). This info is used under the option 
#'          filterNA=TRUE do decide how a double missing value should be treated 
#'          (genuine missing data or double null allele).
run_paral_geneland <- function(nrun, ncores=NULL, model="Correlated", 
                               main.dir="output_corr", spatial=TRUE, path,
                               GEN, COORD, jitter=0,  nxdom="auto", nydom="auto",
                               burnin, POPSMAX, POPINIT, POPMIN, ITR, THN, nullMatrix=NULL) {
  library(parallel)
  if(is.null(ncores)){
    ncores <- detectCores() 
    if(ncores>nrun) ncores <- nrun
  }
  folders <- make.folders(file.path(path, main.dir),  nrun=nrun)
  #list2env(list(folders=folders), envir=.GlobalEnv)
  cl <- makeCluster(ncores)
  clusterSetRNGStream(cl, iseed=NULL)
  clusterEvalQ(cl, library("Geneland"))
  clusterExport(cl, 
                varlist=c("GEN", "COORD", "nrun", "folders", "burnin", "ITR", "THN",
                          "POPSMAX", "POPINIT", "POPMIN", "nxdom", "nydom",
                          "model", "spatial", "jitter", "nullMatrix"), 
                envir= #.GlobalEnv
                  environment()) 
  on.exit(stopCluster(cl))
  st <-  system.time(
    l <- parLapply(cl, seq_len(nrun), paral_geneland, folders=folders,
                   GEN=GEN, COORD=COORD, 
                   varnpop=TRUE, 
                   npopmin=POPMIN,
                   npopinit=POPINIT,
                   npopmax=POPSMAX,
                   spatial=spatial, 
                   freq.model=model, 
                   nxdom=nxdom, nydom=nydom,
                   burnin=burnin,
                   ITR=ITR,
                   THN=THN,
                   jitter=jitter,
                   nullMatrix=nullMatrix)   # parallel execution
  )
  
  print(st)
  #### Post-process ####
  logs <- lapply(folders, read.logs)
  logs.mcmc <- as.mcmc.list(make.mcmc(logs))
  
  BurnIn <- burnin/THN
  logs.mcmc.bin.rm <- rm.burn(logs.mcmc, burn=BurnIn)
  
  # ESS
  ESSlist <- lapply(logs.mcmc.bin.rm, effectiveSize)
  ESS <- as.data.frame(ESSlist)
  names(ESS) <- paste("Run", 1:length(ESSlist))
  message("ESS summary")
  print(ESS[3:5,])
  
  write.csv(ESS[3:5,], file=file.path(path, paste(main.dir, "ESS.csv", sep="_")))
  
  # Gelman's diagnostic
  g.diag <- gelman.diag(logs.mcmc, autoburnin=FALSE, multivariate = F)
  message("Gelman snd Rubin's convergence diagnostic")
  print(g.diag$psrf[3:5,])
  write.csv(g.diag$psrf[3:5,], file=file.path(path, paste(main.dir, "Gelman.diag.csv", sep="_")))
  
  pdf(file.path(path, paste(main.dir, "plot_trace.pdf", sep="_")))
  plot(as.mcmc.list(logs.mcmc.bin.rm)[, c("likelihood", "posterior", "K"), drop=TRUE])
  dev.off()
  
  # Combine and write
  logs.dt <- rbindlist(logs)
  logs.dt.bin <- logs.dt[iteration>BurnIn,]
  SumReps <- logs.dt.bin[, .(Likelihood=mean(likelihood), Posterior=mean(posterior), K=MODE(K)),
                         by=run]
  SumReps <- SumReps[order(Posterior, decreasing = TRUE),]
  message("Runs summary: Max mean Posterior; Mode K")
  print(SumReps)
  write.csv(SumReps, file.path(path, paste(main.dir, "Run_summary.csv", sep="_")))
  
  lapply(paste0(folders, "/"), make.npop.plots, bur=burnin)
  
  # Null frequencies
  filter.null.alleles <- vector(length = length(folders))
  for(i in seq_along(folders)) {
    param <- as.matrix(read.table(file = file.path(folders[i], "parameters.txt")))
    filter.null.alleles[i] <- as.logical(param[param[, 1] == "filter.null.alleles", 3])
  }
  
  if(all(filter.null.alleles)) {
    nullFreq <- lapply(paste0(folders, "/"), EstimateFreqNA)
    nullFreq_df <- do.call(rbind.data.frame, nullFreq)
    names(nullFreq_df) <- names(nullFreq[[1]])
    nullFreq_df <- rbind(nullFreq_df, sapply(nullFreq_df, mean))
    row.names(nullFreq_df)[nrow(nullFreq_df)] <- "Mean"
    nullFreq_df <- rbind(nullFreq_df, sapply(nullFreq_df, sd))
    row.names(nullFreq_df)[nrow(nullFreq_df)] <- "SD"
    message("Mean null allele frequencies:")
    print(nullFreq_df["Mean",])
    write.csv(nullFreq_df, file.path(path, paste(main.dir, "NullFreq.csv", sep="_")))
  }
  
  return(list(Output_paths=folders, ESS=ESS, GelmansDiagnostic=g.diag, SummaryReps=SumReps,
              Null = if(all(filter.null.alleles)) nullFreq_df else NULL))
}
  


# Execute Geneland and collate results in log
paral_geneland <- function(irun, folders, GEN, COORD, varnpop=TRUE, 
                           spatial=TRUE, freq.model, jitter,
                           nxdom=nxdom, nydom=nydom,
                           burnin=100000,
                           npopinit=1,
                           npopmin=1,
                           npopmax=5,
                           ITR=1000000,
                           THN=100,
                           nullMatrix) {
  if(nxdom=="auto"&nydom=="auto") {
    nxdom <- 400
    nydom <- round(nxdom * (max(COORD[,2]) - min(COORD[,2])) /
      (max(COORD[,1]) - min(COORD[,1])), digits=0)
  }
  path.mcmc<-paste0(folders[irun], "/")
  MCMC(geno.dip.codom=GEN,
       coordinates=COORD, 
       varnpop=varnpop,
       npopmin=npopmin,
       npopinit=npopinit,
       npopmax=npopmax,
       spatial=spatial, 
       freq.model=freq.model,
       nit=ITR,
       thinning=THN,
       path.mcmc=path.mcmc,
       delta.coord=jitter,
       miss.loc=nullMatrix,
       write.size.pop=TRUE)
  
  PostProcessChain(coordinates=COORD, 
                   path.mcmc=path.mcmc, 
                   nxdom=nxdom, 
                   nydom=nydom, 
                   burnin=burnin/THN)
  
  log.likel <- read.table(paste0(path.mcmc, "log.likelihood.txt"))
  iteration <- 1:nrow(log.likel)
  
  log.post <- read.table(paste0(path.mcmc, "log.posterior.density.txt"))
  log.K <- read.table(paste0(path.mcmc, "populations.numbers.txt"))
  log <- cbind(run=irun, 
               iteration, 
               likelihood=log.likel[,1], 
               posterior=log.post[,1], 
               K=log.K[,1])
  write.table(log, file=paste0(path.mcmc, "log.txt"), row.names=FALSE)
}

#### Parallel admix geneland ####
paral_adm_geneland <- function(nrun, ncores=NULL, adm="adm", 
                               GEN, COORD, path.mcmc.noadm, 
                               ITR, THN) {
  library(parallel)
  if(is.null(ncores)){
    ncores <- detectCores() 
    if(ncores>nrun) ncores <- nrun
  }
  folders <- make.folders(file.path(path.mcmc.noadm, adm),  nrun=nrun)
  #list2env(list(folders=folders), envir=.GlobalEnv)
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl))
  clusterSetRNGStream(cl, iseed=NULL)
  clusterEvalQ(cl, library("Geneland"))
  clusterExport(cl, 
                varlist=c("GEN",  "folders", "ITR", "THN", "path.mcmc.noadm", "COORD"), 
                envir= #.GlobalEnv
                  environment()) 
  st <-  system.time(
    l <- parLapply(cl, paste0(folders, "/"), HZ, 
                   coordinates=COORD, geno.dip.codom=GEN, 
                   path.mcmc.noadm=paste0(path.mcmc.noadm, "/"),
                   nit=ITR, thinning=THN, geno.dip.dom=NULL,
                   geno.hap=NULL,
                   dist.IC=NULL,
                   allele.freq=NULL,
                   ncluster=NULL,
                   cluster.indiv=NULL,
                   a.init=NULL,
                   b.init=NULL,
                   c.init=1,
                   a.max=10,
                   b.max=NULL,
                   c.max=1,
                   estimate.a=TRUE,
                   estimate.b=TRUE,
                   estimate.c=FALSE,
                   common.param=TRUE)   # parallel execution
  )
  
  print(st)
  
  return(folders)
}


#### Create folders for each run in output ####
make.folders <- function(main, nrun) {
  folders<-paste0(main, "R", 1:nrun)
  lapply(folders, dir.create, showWarnings=FALSE, recursive=FALSE, mode="0777") 
  return(folders)
}


# Convert object in mcmc
make.mcmc <- function(x) {
  x.mcmc <- lapply(x, as.mcmc)
  return(x.mcmc)
}

read.logs <- function(path) {
  log <- read.table(file.path(path, "log.txt"), colClasses="numeric", header=TRUE)
  return(log)
}

# Remove burn-in
rm.burn <- function(x, burn) {
  x.burn <- lapply(x, window, start=burn)
  return(x.burn)
}


mode.param <- function(df, param="K") {
    ux <- unique(df[, param])
    ux[which.max(tabulate(match(df[, param], ux)))]
  }

MODE <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

count.vals<-function(i) {
  length(which(Ind.clust.list[[i]]==clust.mode.list[[i]]))/(ITR/THN)
  }

# Make npopplots
make.npop.plots <- function(d, bur=burnin, th=THN) {
  on.exit(dev.off())
  pdf(paste0(d, "Plotnpop", ".pdf"), width=15, height=7)
  Plotnpop(path.mcmc=d, burnin=bur/th)
}


write.Fst <- function(GEN,dirIn, digs=3) {
  Fst <- Fstat.output(genotypes=GEN,path.mcmc=dirIn)
  fst <- do.call(rbind, Fst)
  fst <- round(fst, digits=digs)
  row.names(fst) <- c("Fis", rep("Fst", length.out=nrow(fst)-1))
  colnames(fst) <- paste0("C", seq_len(ncol(fst)))
  write.csv(file.path(dirIn, "Fst-stats.csv"), x=fst)
  return(fst)
}
#' burning as burning/thinning
adm_process <- function(path.mcmc.noadm, adm, nrun, burnin, coord, pal="Set1") {
  library(data.table)
  library(ggplot2)
  dirIn <- make.folders(file.path(path.mcmc.noadm, adm),  nrun=nrun)
  a <- fread(file.path(dirIn, "a.txt"))
  nK <- ncol(a)
  iteration <- 1:nrow(a)
  Knames <- paste0("C", seq_len(nK))
  setnames(a, paste0("a_", Knames))
  a[, Iter:=iteration]
  a_long <- melt.data.table(a, id.vars="Iter")  
  a_ptrace <- ggplot(a_long, aes(Iter, value)) +
    geom_line() +
    facet_grid(variable~.)
  #a_ptrace
  
  ggsave(file.path(dirIn, "a_ptrace.pdf"), plot=a_ptrace, 
         height=5*length(a_long[, unique(variable)]), units="cm")
  
  
  b <- fread(file.path(dirIn, "b.txt"))
  setnames(b, paste0("b_", Knames))
  b[, Iter:=iteration]
  b_long <- melt.data.table(b, id.vars="Iter")  
  b_ptrace <- ggplot(b_long, aes(Iter, value)) +
    geom_line() +
    facet_grid(variable~.)
  #b_ptrace
  
  ggsave(file.path(dirIn, "b_ptrace.pdf"), plot=b_ptrace, 
         height=5*length(b_long[, unique(variable)]), units="cm")
  
  #set key and subset rm burning
  log_long <- rbindlist(list(a_long, b_long))
  setkey(log_long, Iter)
  desc_stats <- log_long[J(burnin:nrow(a)), .(Mean=mean(value), SD=sd(value), median(value),
                   Lower=quantile(value, 0.0275), Upper=quantile(value, 0.975)),
           by=variable]
  
  write.csv(x=desc_stats, file.path(dirIn, "a_b_summary.csv"), row.names=FALSE)
  
  # Q
  q <- fread(file.path(dirIn, "q.txt"))
  setnames(q, Knames)
  nsamples <- nrow(q)/nrow(a)
  q[, Iter:=rep(iteration, each=nsamples)]
  q[, Sample:=rep(seq_len(nsamples), nrow(a))]
  setkey(q, Iter)
  q <- q[J(burnin:nrow(a)), ]
  q_long <- melt.data.table(q, id.vars=c("Iter", "Sample"), 
                            variable.name="Cluster", value.name="Ancestry")
  q_mean <- q[, lapply(.SD, mean), .SDcols=Knames, by=Sample]
  q_SD <- q[, lapply(.SD, sd), .SDcols=Knames, by=Sample]
  
  write.csv(x=q_mean, file.path(dirIn, "q_mean.csv"), row.names=FALSE)
  write.csv(x=q_SD, file.path(dirIn, "q_SD.csv"), row.names=FALSE)
  
  q_mean <- cbind(q_mean, coord)
  q_mean <- q_mean[order(q_mean$Long),]
  q_mean[, Sample:=factor(Sample, levels = unique(Sample))]
  q_mean_long <- melt.data.table(q_mean, id.vars="Sample", 
                                 measure.vars=patterns("C"),
                                 variable.name="Cluster",
                                 value.name="Ancestry")
  p <- ggplot(q_mean_long, aes(Sample, Ancestry, fill=Cluster)) + 
    #geom_col() +
    geom_bar(stat="identity", width=1) +
    scale_fill_brewer(palette=pal) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          #axis.text.x=element_text(colour="black", angle=90), 
          legend.position="none") +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(gsub(pattern="\\.|/", "",  dirIn)) #+
    #scale_x_discrete(breaks=qind$Sample[nperPop$Bks], labels=nperPop$Label) +
    #geom_vline(xintercept=nperPop[, PopLim], size=0.1)
  
  #p
  
  ggsave(file.path(dirIn, paste0("BarPlot_K", nK, gsub(pattern="\\.|/", "", dirIn), ".pdf")), 
  plot=p, width=24, height=8, units="cm", dpi=250)
 return(p) 
}


