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
dat$diploid <- dat$haploid * 2
# remove species that have no chromosome number data
dat <- dat[!(is.na(dat$diploid)),]
# keep these tips only
tree <- keep.tip(trees[[1]], dat$species)
# # simulate a single tree with n number of tips
# tree <- tree.bd(pars=c(3,1), max.taxa=150)
# tree$edge.length <- tree$edge.length/max(branching.times(tree))
# # simulate a trait dataset
# trait <-  simChrom(tree, pars=c(2, 1, 2, 1, 0, 0, 0, 0, .2,  .2, 20, 0),
#          limits = c(10, 30), model = "ChromPlus")
# # lets make sure we have enough transitions from 1 sex chromosome system
# # to the other
# while (sum(trait$binary.state) < 3) {
#   trait <-  simChrom(tree, pars=c(.3, .6, .3, .6, 0, 0, 0, 0, .2,  .3, 20, 1),
#                      limits = c(10, 30), model = "ChromPlus")
# }
# # look at the trait distributution
# trait$binary.state
# trait$chrom.num
# make a data table to hold the species names, chromosome number and
# sex chromosome system
dat.new <-  as.data.frame(matrix(data=NA, nrow = Ntip(tree), ncol = 3))
colnames(dat.new) <- c("SpecisName", "chroms", "scs")
# fill in the data table
dat.new$SpecisName <- dat$species
dat.new$chroms <- dat$diploid
dat.new$scs <- dat$SCS
# now we get the qmatrix and pmatrix
inputs <- get.matrixes(chrom.range = range(dat.new$chroms),
                       Haplodiploidy = T,
                       Neo.sex = F,
                       complex = T,
                       chrom.range.expansion = 1,
                       dat = dat.new,
                       trees = tree)
# this will be the new function to get the matrixes
inputs <- get.matrixes.new(chrom.range = range(dat.new$chroms),
                           Haplodiploidy = T,
                           Neo.sex = F,
                           complex = T,
                           chrom.range.expansion = 1,
                           dat = dat.new,
                           trees = tree,
                           def.rates = list(r1 = 1,
                                            r2 = 2,
                                            r3 = 3,
                                            r4 = 4,
                                            r5 = 5,
                                            r6 = 6,
                                            r7 = 7,
                                            r8 = 8,
                                            r9 = 9,
                                            r10 = 10,
                                            r11 = 11,
                                            r12 = 12))
# get the relavent inputs
qmat <- inputs$qmat
pmat <- inputs$pmat
dat <- inputs$dat
trees <- inputs$trees
states <- inputs$states
karyotypes <- inputs$karyotypes
# perform stochastic mappings
hists <- make.simmap(tree = trees,
                     x = pmat,
                     model = qmat,
                     nsim = nsim,
                     pi = "estimated",
                     Q = "mcmc")
# lets look at a randtom stochastic map
# add colors
cols <- setNames(object = viridis(ncol(pmat)),nm = colnames(pmat))
plot(hists[[sample(1:nsim, 1)]], col = cols)
# get the number of times each transision has occured
counts <- describe.simmap(hists)$count
# get colnames that represent SA, SS and AA fusions
SAfusionColnames.r3 <- vector(mode = "character", length = nrow(qmat))
SAfusionColnames.r4 <- vector(mode = "character", length = nrow(qmat))
SAfusionColnames.r6 <- vector(mode = "character", length = nrow(qmat))
SAfusionColnames.r7 <- vector(mode = "character", length = nrow(qmat))
AAfusionColnames <- vector(mode = "character", length = nrow(qmat))
SSfusionColnames <- vector(mode = "character", length = nrow(qmat))
for(i in 1:nrow(qmat)){
  if(length(which(qmat[i,] == 3)) != 0){
    SAfusionColnames.r3[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 3)),sep = "")
  }
  if(length(which(qmat[i,] == 4)) != 0){
    SAfusionColnames.r4[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 4)),sep = "")
  }
  if(length(which(qmat[i,] == 6)) != 0){
    SAfusionColnames.r6[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 6)),sep = "")
  }
  if(length(which(qmat[i,] == 7)) != 0){
    SAfusionColnames.r7[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 7)),sep = "")
  }
  if(length(which(qmat[i,] == 2)) != 0){
    AAfusionColnames[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 2)),sep = "")
  }
  if(length(which(qmat[i,] == 11)) != 0){
    SSfusionColnames[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == 11)),sep = "")
  }
}
# remove those that are empty
SAfusionColnames.r3 <- SAfusionColnames.r3[SAfusionColnames.r3 != ""]
SAfusionColnames.r4 <- SAfusionColnames.r4[SAfusionColnames.r4 != ""]
SAfusionColnames.r6 <- SAfusionColnames.r6[SAfusionColnames.r6 != ""]
SAfusionColnames.r7 <- SAfusionColnames.r7[SAfusionColnames.r7 != ""]
AAfusionColnames <- AAfusionColnames[AAfusionColnames != ""]
SSfusionColnames <- SSfusionColnames[SSfusionColnames != ""]
# get the col number that represent each fusion type
SAfusionCol.r3 <- which(colnames(counts) %in% SAfusionColnames.r3)
SAfusionCol.r4 <- which(colnames(counts) %in% SAfusionColnames.r4)
SAfusionCol.r6 <- which(colnames(counts) %in% SAfusionColnames.r6)
SAfusionCol.r7 <- which(colnames(counts) %in% SAfusionColnames.r7)
AAfusionCol <- which(colnames(counts) %in% AAfusionColnames)
SSfusionCol <- which(colnames(counts) %in% SSfusionColnames)
# get the number of times each fusion has occured
SAfusioncounts.r3 <- rowSums(counts[, c(SAfusionCol.r3)])
SAfusioncounts.r4 <- rowSums(counts[, c(SAfusionCol.r4)])
SAfusioncounts.r6 <- rowSums(counts[, c(SAfusionCol.r6)])
SAfusioncounts.r7 <- rowSums(counts[, c(SAfusionCol.r7)])
AAfusioncounts <- rowSums(counts[, c(AAfusionCol)])
SSfusioncounts <- rowSums(counts[, c(SSfusionCol)])
# get all SA fusion counts
SAfusioncounts <- SAfusioncounts.r3 + SAfusioncounts.r4 + SAfusioncounts.r6 + SAfusioncounts.r7
# get the observed pfSA
obspropSA <- SAfusioncounts/(AAfusioncounts+SAfusioncounts+SSfusioncounts)
# obspropSA <- SAfusioncounts/(AAfusioncounts+SAfusioncounts)
# obspropSA <- SAfusioncounts/AAfusioncounts
# remove entries that are either inf or NaN
# obspropSA <- obspropSA[!(is.infinite(obspropSA) | is.nan(obspropSA))]
# now we get the expected pfsa
# make a table to hold the pfsa given scs and chrom number
pfSA.tab <- as.data.frame(matrix(data = NA, nrow = length(karyotypes), ncol = 2))
colnames(pfSA.tab) <- c("state", "pfsa")
for(i in 1:length(karyotypes)){
  # for XO / ZO
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XO", "ZO")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")),
                             scs = "XO")
  }
  # for XY / ZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XY", "ZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")),
                             scs = "XY")
  }
  # for Neo.XY / Neo.ZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("Neo.XY", "Neo.ZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")),
                             scs = "XY")
  }
  # for XXY / ZZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XXY","ZZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z./]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")),
                             scs = "XXY")
  }
  # for XYY / ZWW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XYY", "ZWW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z./]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")),
                             scs = "XYY")
  }
}
# get the expeveted pSA
expSA <- vector(mode = "numeric", length = length(karyotypes))
for(i in 1:nsim){
  times <- describe.simmap(hists[[i]])$times[2, -(nrow(pfSA.tab)+1)]
  expSA[i] <- sum(times * pfSA.tab$pfsa)
}
# lets plot
# clear any previous plots
dev.off()
# get the densities of expected and observed pfsa
den.exp <- density(expSA)
den.obs <- density(obspropSA)
# define the limits of the plot region
xlim <- c(min(den.exp$x, den.obs$x), max(den.exp$x, den.obs$x))
ymax <- max(den.exp$y, den.obs$y)
ymin <- 0 - (ymax *.1)
ylim <- c(ymin, ymax)
plot(den.exp,
     xlim = xlim,
     ylim = ylim,
     main = "Orthoptera \nnot accounting for SS fusions",
     xlab = "Proportion sex-autosome fusion",
     cex.axis = 1, cex.lab = 1)
polygon(den.exp,
        col = rgb(1, 0, 0, .3))
lines(den.obs)
polygon(den.obs,
        col = rgb(0, 0, 1, .3))
# make legend
xpoint <- rep(0.05,2)
points(x = xpoint,
       y = c((ymax - ymax * .05), (ymax - ymax * .1)),
       pch = 15,
       col = c(rgb(1, 0, 0, .5),
               rgb(0, 0, 1, .5)),
       cex = 1.5)
text(x = xpoint,
     y = c((ymax - ymax * .05), (ymax - ymax * .1)),
     labels = c("Expected", "Inferred"),
     pos = 4,cex = .8)
# get the HPD intervals
HPDobsPsa <- HPDinterval(as.mcmc(obspropSA))
HPDexpPsa <- HPDinterval(as.mcmc(expSA))
# plot the HPD intervals
# expected
segments(x0 = HPDexpPsa[,1],
         x1 = HPDexpPsa[,2],
         y0 = ymin,
         y1 = ymin,
         col = rgb(1, 0, 0, .5),
         lwd = 3)
# observed
segments(x0 = HPDobsPsa[,1],
         x1 = HPDobsPsa[,2],
         y0 = ymin/2,
         y1 = ymin/2,
         col = rgb(0, 0, 1, .5),
         lwd = 3)

# save results
save.image("../results/pfsa.complex.model.RData")
