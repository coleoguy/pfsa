# Terrence Sylvester
# 27 October 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(chromePlus)
library(evobiR)
library(coda)
# load helper functions
source("helper.functions.R")
# defininf the number of simulations that will be performed
nsim <- 1000
# simulate a single tree with n number of tips
tree <- rcoal(n = 15)
# simulate a trait dataset
trait <-  simChrom(tree, pars=c(2, 1, 1, 2, 0, 0, 0, 0, 1,  1, 20, 0), 
         limits = c(10, 30), model = "ChromPlus")
# lets make sure we have enough transitions from 1 sex chromosome system
# to the other
while (sum(trait$binary.state) < 3) {
  trait <-  simChrom(tree, pars=c(.3, .6, .3, .6, 0, 0, 0, 0, .2,  .3, 20, 1), 
                     limits = c(10, 30), model = "ChromPlus")
}
# look at the trait distributution
trait$binary.state
trait$chrom.num
# make a data table to hold the species names, chromosome number and 
# sex chromosome system
dat <-  as.data.frame(matrix(data=NA, nrow = Ntip(tree), ncol = 3))
colnames(dat) <- c("SpecisName", "chroms", "scs")
# fill in the data table
dat$SpecisName <- tree$tip.label
dat$chroms <- trait$chrom.num
dat$scs <- c("XO", "XY")[trait$binary.state + 1]
# now we get the qmatrix and pmatrix
inputs <- get.matrixes(chrom.range = range(dat$chroms),
                       Haplodiploidy = T,
                       Neo.sex = F,
                       complex = F,
                       chrom.range.expansion = 1,
                       dat = dat,
                       trees = tree)
# get the relavent inputs
qmat <- inputs$qmat
pmat <- inputs$pmat
dat <- inputs$dat
trees <- inputs$trees
states <- inputs$states
# perform stochastic mappings
hists <- make.simmap(tree = trees,
            x = pmat,
            model = qmat,
            nsim = nsim,
            pi = "estimated")
# lets look at a randtom stochastic map
plot(hists[[sample(1:nsim, 1)]])
# get the number of times each transision has occured
counts <- describe.simmap(hists)$count
# get colnames that represent SA and AA fusions
SAfusionColnames <- vector(mode = "character", length = nrow(qmat))
AAfusionColnames <- vector(mode = "character", length = nrow(qmat))
for(i in 1:nrow(qmat)){
  if(length(which(qmat[i,] == 3)) != 0){
    SAfusionColnames[i] <- paste(i,",",which(qmat[i,] == 3),sep = "")
  }
  if(length(which(qmat[i,] == 2)) != 0){
    AAfusionColnames[i] <- paste(i,",",which(qmat[i,] == 2),sep = "")
  }
}
# remove those that are empty 
SAfusionColnames <- SAfusionColnames[SAfusionColnames != ""]
AAfusionColnames <- AAfusionColnames[AAfusionColnames != ""]
# get the col number that represent each fusion type
SAfusionCol <- which(colnames(counts) %in% SAfusionColnames)
AAfusionCol <- which(colnames(counts) %in% AAfusionColnames)
# get the number of times each fusion has occured
SAfusioncounts <- rowSums(counts[, c(SAfusionCol)])
AAfusioncounts <- rowSums(counts[, c(AAfusionCol)])
# get the observed pfSA
obspropSA <- SAfusioncounts/AAfusioncounts
# remove entries that are either inf or NaN
obspropSA <- obspropSA[!(is.infinite(obspropSA) | is.nan(obspropSA))]
# now we get the expected pfsa
# make a table to hold the pfsa given scs and chrom number
pfSA.tab <- as.data.frame(matrix(data = NA, nrow = nrow(states), ncol = 2))
colnames(pfSA.tab) <- c("state", "pfsa")
for(i in 1:nrow(states)){
  # for XO / ZO
  if(gsub(pattern = "[0-9]",x= states$karyotype[i], replacement = "") %in% c("XO", "ZO")){
    pfSA.tab$state[i] <- states$sim.state[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= states$karyotype[i] ,
                                                  replacement = "")),
                             scs = "XO")
  }
  # for XY / ZW
  if(gsub(pattern = "[0-9]",x= states$karyotype[i], replacement = "") %in% c("XY", "ZW")){
    pfSA.tab$state[i] <- states$sim.state[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= states$karyotype[i] ,
                                                  replacement = "")),
                             scs = "XY")
  }
  # for Neo.XY / Neo.ZW
  if(gsub(pattern = "[0-9]",x= states$karyotype[i], replacement = "") %in% c("Neo.XY", "Neo.ZW")){
    pfSA.tab$state[i] <- states$sim.state[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= states$karyotype[i] ,
                                                  replacement = "")),
                             scs = "XY")
  }
  # for XXY, XYY / ZZW, ZWW
  if(gsub(pattern = ".*/",x= states$karyotype[i], replacement = "") %in% c("XXY", "XYY", "ZWW", "ZZW")){
    pfSA.tab$state[i] <- states$sim.state[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z./]",
                                                  x= states$karyotype[i] ,
                                                  replacement = "")),
                             scs = "XXY")
  }
  
}
# get the expeveted pSA
expSA <- c()
for(i in 1:nsim){
  times <- describe.simmap(hists[[i]])$times[2, -(nrow(states)+1)]
  expSA[i] <- sum(times * pfSA.tab$pfsa)
}

# lets plot
# clear any previous plots
dev.off()
plot(density(expSA, bw = .009),
     xlim = c(-0.5, 1.5), 
     ylim = c(-2,40),
     main = "",
     xlab = "Proportion sex-autosome fusion",
     cex.axis = 1, cex.lab = 1)
polygon(density(expSA, bw = .009),
        col = rgb(1, 0, 0, .3))
lines(density(obspropSA))
polygon(density(obspropSA),
        col = rgb(0, 0, 1, .3))
points(x = c(1, 1),
       y = c(40, 38),
       pch = 15,
       col = c(rgb(1, 0, 0, .5),
               rgb(0, 0, 1, .5)),
       cex = 1.5)
text(x = c(1, 1),
     y = c(40, 38),
     labels = c("Expected", "Inferred"),
     pos = 4,cex = .8)

mean(obspropSA)
HPDobsPsa <- HPDinterval(as.mcmc(obspropSA))
HPDexpPsa <- HPDinterval(as.mcmc(expSA))

# plot the HPD intervals
segments(x0 = HPDexpPsa[,1],
         x1 = HPDexpPsa[,2],
         y0 = -1,
         y1 = -1,
         col = rgb(1, 0, 0, .5),
         lwd = 3)

segments(x0 = HPDobsPsa[,1],
         x1 = HPDobsPsa[,2],
         y0 = -2,
         y1 = -2,
         col = rgb(0, 0, 1, .5),
         lwd = 3)

