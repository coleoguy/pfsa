library(phytools)
library(evobiR)

# helper functions
source("functions.R")

# make a random tree
# set.seed(1)
tree <- rcoal(n = 3,
              br = runif(n = 100))

# plot the tree without tip lables
# plot(tree, show.tip.label = F)

# make thre tree unit length
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# make a data table with four columns
# col 1 - species names
# col 2 - chromosome number (haploid)
# col 3 - sex chromosome system state 0 - XO, 1 - XY
dat <- as.data.frame(matrix(data = NA,
                            nrow = length(tree$tip.label),
                            ncol = 3))

# lable columns
colnames(dat) <- c("species", "chroms", "scs")

# fill the data frame
dat$species <- tree$tip.label
dat$chroms <- c(10,9,11)
dat$scs <- c("XY","XY","XO")

# make the qmatrix
qmat <- qmatGen(chrom.range = range(dat$chroms))

# get the probability matrix
probMat <- getProbMatrix(dat)

# stochastic mapping
# lets assume 11XO as the root state
# following will get the probability of each karyotype ar the root
rootChrom <- which(colnames(probMat) == "11XO")
rootProbs <- rep(0, length(colnames(probMat)))
rootProbs[rootChrom] <- 1

map <- make.simmap(tree = tree,
                   x = probMat,
                   model = qmat,
                   pi = rootProbs,
                   nsim = 100,
                   message = TRUE,
                   Q = "mcmc")

# get the number of transisions
simmapSummary <- describe.simmap(map)

# plot the simmap
# plotSimmap(map)

## following will give the observed distribution of the proportion of
## sex chromosome autosome fusions

# get the proportion of sex chromosome autosome fusions
# using the qmatrix get the rownames and colnamaes of events of sex chromosome
# autosome fusion. (according to the above qmatrix it is rate parameter 3)

sa.hit <- vector(mode = "list", length = length(qmat.new[qmat.new == 3]))
counter <- 1
for(i in 1:ncol(qmat.new)){
  for(j in 1:nrow(qmat.new)){
    if(qmat.new[j,i] == 3){
      sa.hit[[counter]] <- c(row.names(qmat.new)[j],colnames(qmat.new)[i])
      counter <- counter + 1
    }
  }
}


# get the number of autosome autosome fusions
# (according to the above qmatrix it is rate parameter 2)
aa.hit <- vector(mode = "list", length = length(qmat.new[qmat.new == 2]))
counter <- 1
for(i in 1:ncol(qmat.new)){
  for(j in 1:nrow(qmat.new)){
    if(qmat.new[j,i] == 2){
      aa.hit[[counter]] <- c(row.names(qmat.new)[j],colnames(qmat.new)[i])
      counter <- counter + 1
    }
  }
}


#get the columns that represent sex chromosome autosome fusions
sa.col.hit <- vector(mode = "numeric", length = length(sa.hit))

for(i in 1:length(sa.hit)){
  if(length(which(colnames(simmapSummary$count) == paste(unlist(sa.hit[[i]])[1],",", unlist(sa.hit[[i]])[2], sep = ""))) != 0){
    sa.col.hit[i] <- which(colnames(simmapSummary$count) == paste(unlist(sa.hit[[i]])[1],",", unlist(sa.hit[[i]])[2], sep = ""))
  }
}

#get the columns that represent autosome autosome fusions
aa.col.hit <- vector(mode = "numeric", length = length(aa.hit))

for(i in 1:length(aa.hit)){
  if(length(which(colnames(simmapSummary$count) == paste(unlist(aa.hit[[i]])[1],",", unlist(aa.hit[[i]])[2], sep = ""))) != 0){
    aa.col.hit[i] <- which(colnames(simmapSummary$count) == paste(unlist(aa.hit[[i]])[1],",", unlist(aa.hit[[i]])[2], sep = ""))
  }
}

# if zeros are present remove them
sa.col.hit <- sa.col.hit[sa.col.hit != 0]
aa.col.hit <- aa.col.hit[aa.col.hit != 0]

# for each tree calculate the frequency of sex chromosome autosome fusions
pfsaObserved <- vector(mode = "numeric", length = 100)

for (i in 1:100) {
  pfsaObserved[i] <- sum(simmapSummary$count[i,sa.col.hit]) / (sum(simmapSummary$count[i,sa.col.hit]) + sum(simmapSummary$count[i,aa.col.hit]))
}


## following will give the expected distribution of the proportion of
## sex chromosome autosome fusions

# for each tree get the proportion of time each state spent in the tree
# get the column number which has the total branch length
# make a table to store weigthed pfsa values
pfsaExpected <- matrix(data = NA,
               nrow = 100,
               ncol = 2)

pfsaExpected[,1] <- 1:100

for(j in 1:100){

  totalColumn <- which(colnames(simmapSummary$times) == "total")

  # get the weighted pfsa for each tree
  weightedPfsa <- matrix(data = NA,
                         nrow = ncol(simmapSummary$times) - 1,
                         ncol = 5)

  colnames(weightedPfsa) <- c("qmatState",
                              "DiploidCount",
                              "SCS",
                              "timeSpendinTree",
                              "weightedPfsa")

  weightedPfsa[,1] <- colnames(simmapSummary$times)[1:(ncol(simmapSummary$times) - 1)]
  weightedPfsa[,2] <- as.numeric(qmatStates[c(as.numeric(weightedPfsa[,1])),2]) * 2
  weightedPfsa[,3] <- qmatStates[c(as.numeric(weightedPfsa[,1])),3]
  weightedPfsa[,4] <- simmapSummary$times[j,1:(ncol(simmapSummary$times) - 1)] / simmapSummary$times[j,totalColumn]

  for(i in 1:nrow(weightedPfsa)){
    weightedPfsa[i,5] <-  Pfsa(Da = as.numeric(weightedPfsa[i,2]),
                               scs = weightedPfsa[i,3]) * as.numeric(weightedPfsa[i,4])
  }

  pfsaExpected[j,2] <- sum(as.numeric(weightedPfsa[,5]))
}

# plot the observed and expected distributions

plot(density(pfsaObserved),
     ylim = c(0, 50),
     main = "Distributions of proportion of observed \nexpected sex chromosome autosome fusions",
     xlab = "Proportion")

polygon(density(pfsaExpected[,2]))








############ snippets ###############


# haploid autosome count
maxChromValue <- max(dat$chroms)
nChroms <- (maxChromValue - 1) * 2

# get the max chrom number at XO state
nChromsXO <- length(dat$chroms[dat$scs==1])

# fill in the states of the qmat
# XO
dat$qmatState[dat$scs==0] <- dat$chroms[dat$scs==0] - 1

#XY
dat$qmatState[dat$scs==1] <- dat$chroms[dat$scs==1] - 1 + (nChroms/2)

# get the range of chromosome number at each sex chromosome state
rngXO <- range(dat$qmatState[dat$scs == 0])
rngXY <- range(dat$qmatState[dat$scs == 1])

# define the miargin
# this is basically by how much you change the min and max chromosome number
margin <- 1

# get the full range of chromosome numbers icluding all sex chromosome systems
rng <- c(c((rngXO[1] - margin):(rngXO[2] + margin)),
         c((rngXY[1] - margin):(rngXY[2] + margin)))


# get the probability matrix
probMat <- matrix(data = 0,
                  nrow = length(dat$species),
                  ncol = length(rng))

colnames(probMat) <- rng
rownames(probMat) <- dat$species

# assign tip probabilities
for(i in 1:nrow(dat)){
  probMat[i, which(colnames(probMat) == dat$qmatState[i])] <- 1
}

# get the named tip states
tip.states <- dat$qmatState
names(tip.states) <- dat$species

# plot tipl tables
tip.labs <- paste(dat$chroms, c("XO","XY")[dat$scs + 1])

tiplabels(tip.labs,
          frame = "n",
          offset = .02,
          cex = .7)

# make the qmatrix
qmat <- qmatGen(maxChromValue)

# make a table that gives the qmat state respective diploid count and
# the respective sex chromosome system
qmatStates <- matrix(data = NA,
                     nrow = nChroms,
                     ncol = 3)
# fill in qmat states
qmatStates[,1] <- colnames(qmat)

# fill in chromosome numners
qmatStates[(1:(nChroms/2)) ,2] <- (1:(nChroms/2)) + 1
qmatStates[((nChroms/2) + 1):(nChroms) ,2] <- (1:(nChroms/2)) + 1

# fill in sex chromosome systems
qmatStates[(1:(nChroms/2)) ,3] <- "XO"
qmatStates[((nChroms/2) + 1):(nChroms) ,3] <- "XY"


# get the subset of qmatrix
qmat.new <- qmat[rng, rng]





