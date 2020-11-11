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
dat <- dat[dat$SCS %in% c("XY", "XYY", "XXY"),]
# remove species that have no chromosome number data
dat <- dat[!(is.na(dat$haploid)),]
dat$hap.auto[dat$SCS == "XY"] <- dat$haploid[dat$SCS == "XY"] - 1
dat$hap.auto[dat$SCS == "XXY"] <- dat$haploid[dat$SCS == "XXY"] - 2
# keep these tips only
phy <- vector(mode = "list", length = length(trees))
for(i in 1:100){
  phy[[i]] <- keep.tip(trees[[i]], dat$species)
}
class(phy) <- "multiPhylo"
# tree <- keep.tip(trees[[1]], dat$species)
# make a data table to hold the species names, chromosome number and
# sex chromosome system
dat.new <-  as.data.frame(matrix(data=NA, nrow = Ntip(phy[[1]]), ncol = 3))
colnames(dat.new) <- c("SpeciesName", "chroms", "scs")
# fill in the data table
dat.new$SpeciesName <- dat$species
dat.new$chroms <- dat$hap.auto
dat.new$scs <- dat$SCS
# now we get the qmatrix and pmatrix
# inputs <- get.matrixes(chrom.range = range(dat.new$chroms),
#                        Haplodiploidy = T,
#                        Neo.sex = F,
#                        complex = T,
#                        chrom.range.expansion = 1,
#                        dat = dat.new,
#                        trees = tree)
# this will be the new function to get the matrixes
# inputs <- get.matrixes.new(chrom.range = range(dat.new$chroms),
#                            Haplodiploidy = F,
#                            Neo.sex = F,
#                            complex = T,
#                            chrom.range.expansion = 1,
#                            dat = dat.new,
#                            trees = phy,
#                            def.rates = list(r1 = 1,
#                                             r2 = 2,
#                                             r3 = 3,
#                                             r4 = 4,
#                                             r5 = 5,
#                                             r6 = 6,
#                                             r7 = 7,
#                                             r8 = 0,
#                                             r9 = 9,
#                                             r10 = 10,
#                                             r11 = 0,
#                                             r12 = 12,
#                                             r13 = 13))

inputs <- get.matrixes.v2(haploid.scs = F,
                          autosome.as.input = T,
                          Neo.sex = F,
                          complex = T,
                          chrom.range.expansion = 0,
                          dat = dat.new,
                          trees = phy,
                          def.rates = list(r01 = 1,    # AA fusion  XO
                                           r02 = 2,    # AA fission XO
                                           r03 = 1,    # AA fusion  XY
                                           r04 = 2,    # AA fission XY
                                           r05 = 1,    # AA fusion  Neo.XY
                                           r06 = 2,    # AA fission Neo.XY 
                                           r07 = 1,    # AA fusion  XXY
                                           r08 = 2,    # AA fission XXY
                                           r09 = 1,    # AA fusion  XYY
                                           r10 = 2,  # AA fission XYY
                                           r11 = 11,  # SA fusion  XO -> XY
                                           r12 = 11,  # SA fusion  XY -> Neo.XY
                                           r13 = 11,  # SA fusion  XY -> XXY
                                           r14 = 11,  # SA fusion  XY -> XYY
                                           r15 = 15,  # transision Neo.XY -> XY
                                           r16 = 16,  # X fission  XY -> XXY
                                           r17 = 17,  # Y fission  XY -> XYY
                                           r18 = 0,  # Y loss     XY -> XO
                                           r19 = 0,  # Y loss     XYY -> XY
                                           r20 = 0,  # X fusion   XXY -> XY
                                           r21 = 0,  # Y capture  XO -> XY
                                           r22 = 22,  # polyploidy XO
                                           r23 = 23,  # polyploidy XY
                                           r24 = 24,  # polyploidy Neo.XY
                                           r25 = 25,  # polyploidy XXY
                                           r26 = 26)) # polyploidy XYY)

# get the relavent inputs
qmat <- inputs$qmat
pmat <- inputs$pmat
# trees <- inputs$trees
states <- inputs$states
karyotypes <- inputs$karyotypes
#plot tree and states
plot(phy[[1]], show.tip.label = F)
tiplabels(col = c("red", "blue")[as.factor(inputs$dat$scs)], pch = 16, offset = 1)
tiplabels(inputs$dat$chroms)
# perform stochastic mappings
x <- Sys.time()
hists <- make.simmap(tree = phy,
                     x = pmat,
                     model = qmat,
                     nsim = nsim,
                     pi = "estimated")
time.took <- Sys.time() - x
time.took
# Time difference of 55.65726 mins
# "Time difference of 53.06742 mins"
#lets look at a randtom stochastic map
# add colors
cols <- setNames(object = viridis(ncol(pmat)),nm = colnames(pmat))
plot(hists[[sample(1:nsim, 1)]], col = cols)
# get the number of times each transision has occured
counts <- describe.simmap(hists)$count
# get the obspSA value
obspSA <- obspSA(counts = counts,qmat = qmat,states = states)
obspropSA <- obspSA$obspropSA
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
                                                  replacement = "")) * 2,
                             scs = "XO")
  }
  # for XY / ZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XY", "ZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")) * 2,
                             scs = "XY")
  }
  # for Neo.XY / Neo.ZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("Neo.XY", "Neo.ZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")) * 2,
                             scs = "XY")
  }
  # for XXY / ZZW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XXY","ZZW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z./]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")) * 2,
                             scs = "XXY")
  }
  # for XYY / ZWW
  if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XYY", "ZWW")){
    pfSA.tab$state[i] <- karyotypes[i]
    pfSA.tab$pfsa[i] <- Pfsa(Da = as.numeric(gsub(pattern = "[A-z./]",
                                                  x= karyotypes[i] ,
                                                  replacement = "")) * 2,
                             scs = "XYY")
  }
}
# get the expeveted pSA
expSA <- vector(mode = "numeric", length = length(hists))
for(i in 1:length(hists)){
  times <- describe.simmap(hists[[i]])$times[2, -(nrow(pfSA.tab)+1)]
  expSA[i] <- sum(times * pfSA.tab$pfsa)
}
# lets plot
# clear any previous plots
dev.off()
# get the densities of expected and observed pfsa
den.exp <- density(expSA, bw = 0.009)
den.obs <- density(obspropSA, bw = .009)
# define the limits of the plot region
xlim <- c(min(den.exp$x, den.obs$x), max(den.exp$x, den.obs$x))
ymax <- max(den.exp$y, den.obs$y)
ymin <- 0 - (ymax *.1)
ylim <- c(ymin, ymax)
plot(den.exp,
     xlim = xlim,
     ylim = ylim,
     main = "Orthoptera \nTransisions between XY and Complex sex chromosome systems",
     xlab = "Proportion sex-autosome fusion",
     cex.axis = 1, cex.lab = 1)
polygon(den.exp,
        col = rgb(1, 0, 0, .3))
lines(den.obs)
polygon(den.obs,
        col = rgb(0, 0, 1, .3))
# make legend
xpoint <- rep(0.45,2)
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
save.image("../results/pfsa-orthoptera-xy-complex-getmatv2-hap-auto-simple-model.RData")
