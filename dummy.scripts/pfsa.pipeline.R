library(phytools)
library(evobiR)
library(chromePlus)
library(diversitree)

# helper functions
source("pfsa.pipeline.functions.R")

# Load in trees
# Here we will make 100 phylogenies with three tips. this resembles the
# posterior distribution from a bayesian phylogenetic inference

starttime <- Sys.time()

trees <- vector(mode = "list", length = 100)
names(trees) <- paste("tree", 1:100, sep = "")

for(i in 1:100){
  set.seed(1)
  trees[[i]] <- rcoal(n = 3,
                      br = runif(n = 1000))
}

class(trees) <- "multiPhylo"


# this is to store results from the mcmc runs
results.trees <- vector("list", length = 100)
names(results.trees) <- names(trees)

# this is to store stochastic mapps
simmaps <- simmap.summaries <- vector("list", length = 100)
names(simmaps) <- names(simmap.summaries) <- names(trees)

# this is to store the observed pfsa values and expected pfsa values from stochastic mapping
pfsaObserved <- pfsaExpected <- matrix(data = NA,
                                       nrow = 100,
                                       ncol = 100)

colnames(pfsaObserved) <- colnames(pfsaExpected) <- names(trees)
rownames(pfsaObserved) <- rownames(pfsaExpected) <- paste("simulation", 1:100, sep = "")

# get the rates of chromosome number evolution #

for(i in 1:100){
  # sample a single tree from the posterior distribution
  tree <- trees[[i]]
  
  # make thre tree unit length
  tree.depth <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth
  
  # load chromosome data
  # Here we will make a data table with four columns
  
  # col 1 - species names
  # col 2 - chromosome number (haploid)
  # col 3 - sex chromosome system 
  dat <- as.data.frame(matrix(data = NA,
                              nrow = length(tree$tip.label),
                              ncol = 3))
  
  # lable columns
  colnames(dat) <- c("species", "chroms", "scs")
  
  # fill the data frame
  dat$species <- tree$tip.label
  dat$chroms <- c(10,9,11)
  dat$scs <- c("XY","XY","XO")
  
  # get the qmatrix
  qmat <- qmatGen(chrom.range = range(dat$chroms))
  
  # get the probability matrix
  probMat <- getProbMatrix(dat)
  
  # make a copy of the probmat that will be used in mcmc
  probMatMCMC <- probMat
  
  colnames(probMatMCMC) <-  gsub(pattern = "XO",
                                 replacement = "",
                                 x = colnames(probMatMCMC))
  
  colnames(probMatMCMC) <-  gsub(pattern = "XY",
                                 replacement = "h",
                                 x = colnames(probMatMCMC)) 
  
  # make the likelihood function
  lik <- make.musse(tree = tree,
                    states = probMatMCMC,
                    k = ncol(probMat),
                    strict = F,
                    control = list(method="ode"))
  
  # constrain the likelihood functoin
  con.lik <- constrainMuSSE(data = probMatMCMC,
                          lik = lik,
                          hyper = T,
                          s.lambda = T,
                          s.mu = T, 
                          polyploidy = F,
                          constrain =  list(drop.poly=T,
                                            drop.demi=T))
  
  # mcmc
  print(paste("mcmc of tree:", i))
  
  results.trees[[i]] <-  mcmc(lik = con.lik,
                              x.init = rep(1, length(argnames(con.lik))),
                              nsteps = 100,
                              w = 1,
                              prior = make.prior.exponential(.5))
  
  # scale units to millions of years
  # results.trees[[i]][1:length(argnames(con.lik)) + 1] <- results.trees[[1]][1:length(argnames(con.lik)) + 1] / tree.depth
  
  # stochastic mapping
  print(paste("stockastic mapping of tree:", i))
  # apply rates of chromosome number evolution taken from the above mcmc to the
  # qmatrix
  qmat.stoch <- qmat
  
  # define burnim
  burnin <- .5
  
  # define post burnin range
  post.burnin.range <- ((nrow(results.trees[[i]]) * burnin) + 1 ): nrow(results.trees[[i]])
  
  # autosomal fissions at XO
  qmat.stoch[qmat.stoch == 1] <- mean(results.trees[[i]]$asc1[post.burnin.range])
  # autosomal fusions at XO
  qmat.stoch[qmat.stoch == 2] <- mean(results.trees[[i]]$desc1[post.burnin.range])
  # autosomal fissions at XY
  qmat.stoch[qmat.stoch == 3] <- mean(results.trees[[i]]$asc2[post.burnin.range])
  # autosomal fusions at XY
  qmat.stoch[qmat.stoch == 4] <- mean(results.trees[[i]]$desc2[post.burnin.range])
  # sex chromosome autosome fusions 
  qmat.stoch[qmat.stoch == 5] <- mean(results.trees[[i]]$tran12[post.burnin.range])
  # ay chromosome loss
  qmat.stoch[qmat.stoch == 6] <- mean(results.trees[[i]]$tran21[post.burnin.range])
  
  # make diagonal the negative sum
  for(j in 1:nrow(qmat.stoch)){
    qmat.stoch[j,j] <- 0 - sum(qmat.stoch[j,])
  }
  
  # lets assume 11XO as the root state
  # following will get the probability of each karyotype ar the root
  rootChrom <- which(colnames(probMat) == "11XO")
  rootProbs <- rep(0, length(colnames(probMat)))
  rootProbs[rootChrom] <- 1
  
  simmaps[[i]] <- make.simmap(tree = tree,
                              x = probMat,
                              pi = rootProbs,
                              nsim = 100,
                              message = TRUE,
                              Q = qmat.stoch)
  
  simmap.summaries[[i]] <- describe.simmap(simmaps[[i]])
  
  ## following will give the observed distribution of the proportion of
  ## sex chromosome autosome fusions
  print(paste("calculatoin of the observed proportion of sex chromosome autosome fusions of tree:", i))
  
  # get the proportion of sex chromosome autosome fusions
  # using the qmatrix get the rownames and colnamaes of events of sex chromosome
  # autosome fusion. (according to the above qmatrix it is rate parameter 3)
  
  sa.hit <- vector(mode = "list", length = length(qmat[qmat == 5]))
  counter <- 1
  for(j in 1:nrow(qmat)){
    for(k in 1:ncol(qmat)){
      if(qmat[j,k] == 5){
        sa.hit[[counter]] <- c(row.names(qmat)[j],colnames(qmat)[k])
        counter <- counter + 1
      }
    }
  }
  
  
  # get the number of autosome autosome fusions
  # (according to the above qmatrix it is rate parameter 2)
  aa.hit <- vector(mode = "list", length = length(qmat[qmat %in% c(2,4)]))
  counter <- 1
  for(j in 1:nrow(qmat)){
    for(k in 1:ncol(qmat)){
      if(qmat[j,k] %in% c(2,4)){
        aa.hit[[counter]] <- c(row.names(qmat)[j],colnames(qmat)[k])
        counter <- counter + 1
      }
    }
  }
  
  #get the columns that represent sex chromosome autosome fusions
  sa.col.hit <- vector(mode = "numeric", length = length(sa.hit))
  
  for(j in 1:length(sa.hit)){
    if(length(which(colnames(simmap.summaries[[i]]$count) == paste(unlist(sa.hit[[j]])[1],",", unlist(sa.hit[[j]])[2], sep = ""))) != 0){
      sa.col.hit[j] <- which(colnames(simmap.summaries[[i]]$count) == paste(unlist(sa.hit[[j]])[1],",", unlist(sa.hit[[j]])[2], sep = ""))
    }
  }
  
  #get the columns that represent autosome autosome fusions
  aa.col.hit <- vector(mode = "numeric", length = length(aa.hit))
  
  for(j in 1:length(aa.hit)){
    if(length(which(colnames(simmap.summaries[[i]]$count) == paste(unlist(aa.hit[[j]])[1],",", unlist(aa.hit[[j]])[2], sep = ""))) != 0){
      aa.col.hit[j] <- which(colnames(simmap.summaries[[i]]$count) == paste(unlist(aa.hit[[j]])[1],",", unlist(aa.hit[[j]])[2], sep = ""))
    }
  }
  
  # if zeros are present remove them
  sa.col.hit <- sa.col.hit[sa.col.hit != 0]
  aa.col.hit <- aa.col.hit[aa.col.hit != 0]
  
  # for each tree calculate the frequency of sex chromosome autosome fusions
  for (j in 1:100) {
    pfsaObserved[j,i] <- sum(simmap.summaries[[i]]$count[j,sa.col.hit]) / (sum(simmap.summaries[[i]]$count[j,sa.col.hit]) + sum(simmap.summaries[[i]]$count[j,aa.col.hit]))
  }
  
  ## following will give the expected distribution of the proportion of
  ## sex chromosome autosome fusions
  print(paste("calculatoin of the expected proportion of sex chromosome autosome fusions of tree:", i))
  
  for(j in 1:100){
    
    totalColumn <- which(colnames(simmap.summaries[[i]]$times) == "total")
    
    # get the weighted pfsa for each tree
    weightedPfsa <- matrix(data = NA,
                           nrow = ncol(simmap.summaries[[i]]$times) - 1,
                           ncol = 5)
    
    colnames(weightedPfsa) <- c("karyotype",
                                "DiploidCount",
                                "SCS",
                                "timeSpendinTree",
                                "weightedPfsa")
    
    weightedPfsa[,1] <- colnames(simmap.summaries[[i]]$times)[1:(ncol(simmap.summaries[[i]]$times) - 1)]
    weightedPfsa[,2] <- as.numeric(gsub(pattern = "[aA-zZ]",replacement = "",x = colnames(simmap.summaries[[i]]$times)[1:(ncol(simmap.summaries[[i]]$times) - 1)])) * 2
    weightedPfsa[,3] <- gsub(pattern = "[0-9]",replacement = "",x = colnames(simmap.summaries[[i]]$times)[1:(ncol(simmap.summaries[[i]]$times) - 1)])
    weightedPfsa[,4] <- simmap.summaries[[i]]$times[j,1:(ncol(simmap.summaries[[i]]$times) - 1)] / simmap.summaries[[i]]$times[j,totalColumn]
    
    for(k in 1:nrow(weightedPfsa)){
      weightedPfsa[k,5] <-  Pfsa(Da = as.numeric(weightedPfsa[k,2]),
                                 scs = weightedPfsa[k,3]) * as.numeric(weightedPfsa[k,4])
    }
    
    pfsaExpected[j,i] <- sum(as.numeric(weightedPfsa[,5]))
  }
  
}
# plot the simmap
# plotSimmap(map)

# plot the observed and expected distributions

pfsa.observed.density <- density(pfsaObserved[!is.na(pfsaObserved)])
pfsa.espected.density <- density(pfsaExpected)

# get the HPD intervals
pfsa.espected.HPD <- coda::HPDinterval(coda::as.mcmc(as.numeric(pfsaExpected)))
pfsa.observed.HPD <- coda::HPDinterval(coda::as.mcmc(as.numeric(pfsaObserved[!is.na(pfsaObserved)])))

# get y max xmax amd y min
ymax <- max(c(pfsa.espected.density$y, pfsa.observed.density$y))
ymin <- 0 - (ymax * .05)
xmax <- max(c(pfsa.espected.density$x, pfsa.observed.density$x))

plot(pfsa.observed.density,
     ylim = c(ymin, ymax),
     main = "Distributions of proportion of observed \nexpected sex chromosome autosome fusions",
     xlab = "Proportion",
     col = "red")

polygon(pfsa.observed.density,col = rgb(1,0,0,.5),border = NA)
polygon(pfsa.espected.density,col = rgb(0,0,1,.5),border = rgb(0,0,1,1))

# plot the HPD intervals
segments(x0 = pfsa.espected.HPD[,1],
         x1 = pfsa.espected.HPD[,2],
         y0 = ymin * .5,
         y1 = ymin * .5,
         lwd = 3,
         col = rgb(0,0,1,.5))


segments(x0 = pfsa.observed.HPD[,1],
         x1 = pfsa.observed.HPD[,2],
         y0 = ymin * 1,
         y1 = ymin * 1,
         lwd = 3,
         col = rgb(1,0,0,.5))

# plot the legent
points(x = rep(0.8, 2),
       y = c(ymax,(ymax - (ymax * .05))),
       pch = 16,
       col = c("red", "blue"),
       cex = 1.5)

text(x = rep(0.8, 2),
     y = c(ymax,(ymax - (ymax * .05))),
     labels = c("Observed", "Expected"),
     pos = 4,
     cex = .97)

endtime <- Sys.time()

print(paste("it took ", 
            endtime - starttime,
            " minutes to complete this analysis",sep = ""))

# Time difference of 7.187194 mins

