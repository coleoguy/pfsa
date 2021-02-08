# Terrence Sylvester
# 27 October 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(chromePlus)
library(evobiR)
library(coda)
library(diversitree)
library(viridis)

# load helper functions
source("helper.functions.R")
# defininf the number of simulations that will be performed
nsim <- 100

# get data
dat <- GetData(trees = "../data/Trees/posterior.trees.nex",
               data = "../data/chrom.data/chroms.csv")
# read in trees
trees <- read.nexus("../data/Trees/posterior.trees.nex")
# lets try this to orthoptera
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
dat$hap.auto[dat$SCS == "XO"] <- dat$haploid[dat$SCS == "XO"] - 1
dat$hap.auto[dat$SCS == "XY"] <- dat$haploid[dat$SCS == "XY"] - 1
dat$hap.auto[dat$SCS == "XXY"] <- dat$haploid[dat$SCS == "XXY"] - 2
dat$hap.auto[dat$SCS == "XYY"] <- dat$haploid[dat$SCS == "XYY"] - 1
# keep these tips only
phy <- vector(mode = "list", length = length(trees))
for(i in 1:100){
  phy[[i]] <- keep.tip(trees[[i]], dat$species)
}
class(phy) <- "multiPhylo"
# make a data table to hold the species names, chromosome number and
# sex chromosome system
dat.new <-  as.data.frame(matrix(data=NA, nrow = Ntip(phy[[1]]), ncol = 3))
colnames(dat.new) <- c("SpeciesName", "chroms", "scs")
# fill in the data table
dat.new$SpeciesName <- dat$species
dat.new$chroms <- dat$hap.auto
dat.new$scs <- dat$SCS
# now we get the qmatrix and pmatrix
inputs <- get.matrixes.v2(haploid.scs = T,
                          autosome.as.input = T,
                          Neo.sex = F,
                          complex = T,
                          chrom.range.expansion = 0,
                          dat = dat.new,
                          trees = phy,
                          def.rates = list(r01 = 1,
                                           r02 = 2,
                                           r03 = 1,
                                           r04 = 2,
                                           r05 = 1,
                                           r06 = 2,
                                           r07 = 1,
                                           r08 = 2,
                                           r09 = 1,
                                           r10 = 2,
                                           r11 = 3,
                                           r13 = 3,
                                           r14 = 3,
                                           r16 = 0,
                                           r17 = 0,
                                           r19 = 0,
                                           r20 = 0))

# get the relavent inputs
qmat <- inputs$qmat
pmat <- inputs$pmat
# trees <- inputs$trees
states <- inputs$states
karyotypes <- inputs$karyotypes

####################

dat.new$Karyotype <- paste(dat.new$chroms, dat.new$scs, sep = "")
dat.new$Karyotype.states <- NA

for(i in 1:nrow(dat.new)){
  dat.new$Karyotype.states[i] <- which(karyotypes %in% dat.new$Karyotype[i])
}

MuSSE.states <- dat.new$Karyotype.states
names(MuSSE.states) <- dat.new$SpeciesName

lik <- make.musse(tree = phy[[1]],
           states = MuSSE.states,
           strict = F,
           k = length(karyotypes),
           control = list(method = "ode"))

argnames(lik)

con.lik <- constrainQmat(qmat = qmat, lik = lik)
argnames(con.lik)
diversitree::mcmc(con.lik,
                  x.init = rep(1,length(argnames(con.lik))),
                  nsteps = 10,
                  w = 1, upper = 1)
