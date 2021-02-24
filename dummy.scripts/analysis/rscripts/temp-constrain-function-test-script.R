# Terrence Sylvester
# 27 October 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(chromePlus)
library(evobiR)
# library(coda)
library(diversitree)
library(viridis)

# load helper functions
source("helper.functions.R")

# define the number of simulations that will be performed
nsim <- 100
prior <- make.prior.exponential(r = .5)

# make a place holder for results
results <- vector(mode = "list", length = 100)

for(i in 1:100){
  
  # get data
  dat <- GetData(trees = "../data/Trees/posterior.trees.nex",
                 data = "../data/chrom.data/chroms.csv")
  
  # read in trees
  tree <- read.nexus("../data/Trees/posterior.trees.nex")[[i]]
  
  # lets try this to Orthoptera
  dat <- dat[dat$order == "Orthoptera",]
  
  # rename some of the scs for clarity'=
  dat$SCS[dat$notes == "X1X1X2X2/X1X2Y1"] <- "XXY"
  dat$SCS[dat$notes == "X1X1/X1Y1Y2"] <- "XYY"
  dat$SCS[dat$notes == "X1X1X2X2/X1X2Y"] <- "XXY"
  dat$SCS[dat$SCS == "XY|homomorphic"] <- "XY"
  
  # isolate those data that we need
  dat <- dat[dat$SCS %in% c("XO", "XY", "XXY", "XYY"),]
  # remove species that have no chromosome number data
  dat <- dat[!(is.na(dat$haploid)),]
  
  # get the haploid chromosome number value
  dat$hap.auto[dat$SCS == "XO"] <- dat$haploid[dat$SCS == "XO"] - 1
  dat$hap.auto[dat$SCS == "XY"] <- dat$haploid[dat$SCS == "XY"] - 1
  dat$hap.auto[dat$SCS == "XXY"] <- dat$haploid[dat$SCS == "XXY"] - 2
  dat$hap.auto[dat$SCS == "XYY"] <- dat$haploid[dat$SCS == "XYY"] - 1
  
  # keep these tips only
  phy <- keep.tip(tree, dat$species)
  
  # make tree unit length
  tree.depth <-  max(branching.times(phy))
  phy$edge.length <- phy$edge.length / tree.depth
  
  # make a data table to hold the species names, chromosome number and
  # sex chromosome system
  dat.new <-  as.data.frame(matrix(data=NA, nrow = Ntip(phy), ncol = 3))
  colnames(dat.new) <- c("SpeciesName", "chroms", "scs")
  
  # fill in the data table
  dat.new$SpeciesName <- dat$species
  dat.new$chroms <- dat$hap.auto
  dat.new$scs <- dat$SCS
  
  # now we get the qmatrix and pmatrix
  inputs <- get.matrixes(haploid.scs = T,
                         autosome.as.input = T,
                         Neo.sex = F,
                         complex = T,
                         chrom.range.expansion = 0,
                         dat = dat.new,
                         trees = phy,
                         def.rates = list(r01 = 1,  # r1 = 1,    # AA fusion  XO
                                          r02 = 2,  # r2 = 2,    # AA fission XO
                                          r03 = 1,  # r3 = 3,    # AA fusion  XY
                                          r04 = 2,  # r4 = 4,    # AA fission XY
                                          r05 = 1,  # r5 = 5,    # AA fusion  Neo.XY
                                          r06 = 2,  # r6 = 6,    # AA fission Neo.XY
                                          r07 = 1,  # r7 = 7,    # AA fusion  XXY
                                          r08 = 2,  # r8 = 8,    # AA fission XXY
                                          r09 = 1,  # r9 = 9,    # AA fusion  XYY
                                          r10 = 2,  # r10 = 10,  # AA fission XYY
                                          r11 = 3,  # r11 = 11,  # SA fusion  XO -> XY
                                          r12 = 3,  # r12 = 12,  # SA fusion  XY -> Neo.XY
                                          r13 = 3,  # r13 = 13,  # SA fusion  XY -> XXY
                                          r14 = 3,  # r14 = 14,  # SA fusion  XY -> XYY
                                          r15 = 4,  # r15 = 15,  # transision Neo.XY -> XY
                                          r16 = 5,  # r16 = 16,  # X fission  XY -> XXY
                                          r17 = 6,  # r17 = 17,  # Y fission  XY -> XYY
                                          r18 = 7,  # r18 = 18,  # Y loss     XY -> XO
                                          r19 = 8,  # r19 = 19,  # Y loss     XYY -> XY
                                          r20 = 9,  # r20 = 20,  # X fusion   XXY -> XY
                                          r21 = 10,  # r21 = 21,  # Y capture  XO -> XY
                                          r22 = 0,  # r22 = 22,  # polyploidy XO
                                          r23 = 0,  # r23 = 23,  # polyploidy XY
                                          r24 = 0,  # r24 = 24,  # polyploidy Neo.XY
                                          r25 = 0,  # r25 = 25,  # polyploidy XXY
                                          r26 = 0,  # r26 = 26,  # polyploidy XYY
                                          r27 = 0)) # r27 = 27,  # translocation XO -> XXY (White(1973), Animal cytology and evolution)
  
  # get the relavent inputs
  qmat <- inputs$qmat
  pmat <- inputs$pmat
  # trees <- inputs$trees
  states <- inputs$states
  karyotypes <- inputs$karyotypes
  
  # make a new column to prepare data to make the likelihood function in 
  # diversitree
  dat.new$Karyotype <- paste(dat.new$chroms, dat.new$scs, sep = "")
  dat.new$Karyotype.states <- NA
  # fill in the states names
  for(j in 1:nrow(dat.new)){
    dat.new$Karyotype.states[j] <- which(karyotypes %in% dat.new$Karyotype[j])
  }
  # make a vector which includes states names and species names
  MuSSE.states <- dat.new$Karyotype.states
  names(MuSSE.states) <- dat.new$SpeciesName
  # make the likelihood function
  lik <- make.musse(tree = phy,
                    states = MuSSE.states,
                    strict = F,
                    k = length(karyotypes),
                    control = list(method = "ode"))
  
  # constrain the likelihood function
  con.lik <- constrainQmat(qmat = qmat, lik = lik)
  # run a temp mcmc to get a value to w
  temp <-  mcmc(con.lik,
                x.init = rep(1,length(argnames(con.lik))),
                nsteps = 20,
                w = 1,
                prior = prior,
                upper = 20)
  w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
  
  # run the mcmc
  results[[i]] <- mcmc(con.lik,
                       x.init = rep(1,length(argnames(con.lik))),
                       nsteps = 100,
                       w = w, 
                       prior = prior,
                       upper = 10)
  results[[i]][2:(length(argnames(con.lik))+1)] <- results[[i]][2:(length(argnames(con.lik))+1)] / tree.depth 
}
phy$edge.length <- phy$edge.length * tree.depth
post.burnin <- results[[i]][51:100,]

colnames(post.burnin)
qmat[qmat == 1] <- mean(post.burnin$par01)
qmat[qmat == 2] <- mean(post.burnin$par02)
qmat[qmat == 3] <- mean(post.burnin$par03)
qmat[qmat == 5] <- mean(post.burnin$par05)
qmat[qmat == 7] <- mean(post.burnin$par07)
qmat[qmat == 9] <- mean(post.burnin$par09)
qmat[qmat == 10] <- mean(post.burnin$par10)

hists <- make.simmap(tree = phy,
                     x = pmat,
                     nsim = 1,
                     pi = "equal",
                     Q = qmat)
