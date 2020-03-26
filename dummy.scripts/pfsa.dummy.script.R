library(phytools)
library(chromePlus)
library(geiger)
library(diversitree)
library(viridis)


# make a random tree
tree <- rcoal(n = 3,
              br = runif(n = 20))

# plot the tree without tip lables
plot(tree, show.tip.label = F)


# tree$edge.length <- tree$edge.length * 100

# chroms <- simChrom(tree = tree,
#                    pars = c(.001, .01, 0, 0, 20),
#                    limits = c(10,30),
#                    model = "2010")

# scs <- sample(x = c(1,0), size = length(tree$tip.label), replace = T)

# names(scs) <- names(chroms)

# make thre tree unit length
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# make a data table with three columns
# col 1 - species names
# col 2 - chromosome number (haploid)
# col 3 - sex chromosome system state 1 - XO, 0 - XY
dat <- matrix(data = NA,
              nrow = length(tree$tip.label),
              ncol = 3)

colnames(dat) <- c("species", "chroms", "scs")

dat <- as.data.frame(dat)

dat$species <- tree$tip.label
dat$chroms <- c(9,10,11)
dat$scs <- c(0,1,1)

# plot tipl tables
tip.labs <- paste(dat$chroms, c("XO","XY")[dat$scs + 1])

tiplabels(tip.labs,
          frame = "n",
          offset = .02,
          cex = .7)

# get the range of the chromosome numbers
rng <- range(dat$chroms)

# covert the data to matrix which will be used in making the likelihood function
dat.mat <- datatoMatrix(dat, 
                        range = c(rng[1] - 2,rng[2] + 2),  
                        hyper = T)

# make the likelihood function
lik <- make.mkn(tree = tree, 
                states = dat.mat, 
                k = ncol(dat.mat), 
                strict = F, 
                control = list(method = "ode"))


#  length(argnames(lik))

# constrain the likelihood function. Here we set verbose = T to get the
# rate matrix
con.lik <- constrainMkn(data = dat.mat,
                        lik = lik,
                        hyper = T,
                        polyploidy = F,
                        verbose = T,
                        oneway = F,
                        constrain = list(drop.demi = T,
                                         drop.poly = F))

#  length(argnames(con.lik$`likelihood function`))

# get the rate matrix
qmat <- con.lik$`parameter matrix`

# these are the rates in the rate identity matrix

# rate1 ascending aneuploidy - state1 
# rate2 descending aneuploidy - state1 
# rate3 ascending aneuploidy - state2 
# rate4 descending aneuploidy - state2 
# rate5 polyploidization of state1 
# rate6 polploidization of state2 
# rate7 rediploidization of a state2 
# rate8 transitions from 1 to 2 for hyperstate 
# rate9 transitions from 2 to 1 for hyperstate 
# rate10 demipolyploidy for state1 - even 
# rate11 demipolyploidy for state1 - odd 
# rate12 demipolyploidy for state2 - even 
# rate13 demipolyploidy for state2 - odd 

# adjust the rate matrix so that we have the neccessary rates

# make demiploidy rates to zero
qmat[qmat == 10] <- 0 # rate10 demipolyploidy for state1 - even 
qmat[qmat == 11] <- 0 # rate11 demipolyploidy for state1 - odd
qmat[qmat == 12] <- 0 # rate12 demipolyploidy for state2 - even 
qmat[qmat == 13] <- 0 # rate13 demipolyploidy for state2 - odd 

# correct the possitoin of the transision rates between states
qmat[qmat == 8] <- 0 # rate8 transitions from 1 to 2 for hyperstate 
qmat[qmat == 9] <- 0 # rate9 transitions from 2 to 1 for hyperstate

# if present remove redeploidisations as well
qmat[qmat == 7] <- 0 # rate7 rediploidization of a state2

# put the transision rates in there correct places
# from state 1 > 2
for(i in 2:(nrow(qmat)/2)){
  qmat[i,(nrow(qmat)/2) + i -1] <- 8
}

# from state 2 < 1
for(i in ((nrow(qmat)/2)+2):nrow(qmat)){
  qmat[i, i-(nrow(qmat)/2)-1] <- 9
}

# put rate parameter for y chromosome loss - rate14
for(i in ((nrow(qmat)/2)+1):nrow(qmat)){
  qmat[i, i-(nrow(qmat)/2)] <- 14
}

# put rate parameter for y chromosome capture - rate15
for(i in 1:(nrow(qmat)/2)){
  qmat[i, (nrow(qmat)/2) + i] <- 15
}

# see the qmatrix
# qmat

# stochastic mapping
map <- make.simmap(tree = tree,
                   x = dat.mat,
                   model = qmat,
                   pi = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                   nsim = 1)

# visualize stochastic map
plotSimmap(map, colors = setNames(viridis(n = ncol(dat.mat)), colnames(dat.mat)),
           lwd = 4)

# get the proportion of time each trait has spent in the tree
# also get the number of transisions from XO to Xy
simmap.summary <- describe.simmap(map)

simmap.summary$Tr

