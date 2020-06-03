# Terrence Sylvester
# 26th March 2020
# pradakshanas@gmail.com
# make transition matrix

# make a matrix with 100 chromosomes on each state. Here we will be using
# sex chromosome systems as the states. These states are XO and XY.

# we will go from 2:100 chromosomes. So we will have a total of 198 rows and cols

qmatGen <- function(chrom.range = NULL){

  # define total number of chromosome states
  #chrom.range <- c(10, 19)
  maxChroms <- chrom.range[2] + 1
  minChroms <- chrom.range[1] - 1
  nChroms <- (maxChroms - 1) * 2

  # make the initial q matrix
  qmat <- matrix(data = 0,
                 nrow = nChroms,
                 ncol = nChroms)

  # define rownames and colnames
  # although row names and col names starts from 1, chromsome counts starts from
  # 2
  row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                          paste(2:maxChroms,"XY",sep=""))

  # following are the rate parameters that we define
  # rate1 <- autosomal fusions in XO and XY sex chromosome sytems
  # rate2 <- autosomal fissions in XO anf XY sex chromosome systems
  # rate3 <- sex chromosome autosome fusions
  # rate4 <- Y chromosome loss

  # autosomal fissions and fusions in XO systems
  for(i in 1){
    while (i < (nrow(qmat)/2)) {
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
      i <- i + 1
    }
  }

  # autosomal fissions and fusions in XY systems
  for(i in ((nrow(qmat)/2)+1):nrow(qmat)){
    if(i < nrow(qmat)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # sex chromosome autosome fusion
  for(i in 2:(nrow(qmat)/2)){
    qmat[i,(nrow(qmat)/2) + i -1] <- 3
  }
  # Y chromosome loss
  for(i in ((nrow(qmat)/2)+1):nrow(qmat)){
    qmat[i, i-(nrow(qmat)/2)] <- 4
  }

  if(chrom.range[1] > 2){

  }
  hp <- ncol(qmat)/2 + 1
  drops <- c(1:(minChroms-2), hp:(hp + minChroms-3))
  qmat <- qmat[-drops, -drops]
  return(qmat)
}


getProbMatrix <- function(dat){
    chrom.range <- range(dat$chroms)
    maxChroms <- chrom.range[2] + 1
    minChroms <- chrom.range[1] - 1
    nChroms <- (maxChroms - 1) * 2

    # make the initial q matrix
    qmat <- matrix(data = 0,
                   nrow = nChroms,
                   ncol = nChroms)

    # define rownames and colnames
    # although row names and col names starts from 1, chromsome counts starts from
    # 2
    row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                            paste(2:maxChroms,"XY",sep=""))
    pmat <- qmat[NULL,]
    hp <- ncol(qmat)/2 + 1
    
    # add rows 
    # make a new matrix. number of rows in thes new matrix is equal to the
    # number of species in our data table. 
    temp.mat <- matrix(data = 0,
                       nrow = nrow(dat),
                       ncol = ncol(pmat))
    
    # name the row names according to the species names in our data table
    rownames(temp.mat) <- dat$species

    # combine the new matrix with the pmat
    pmat <- rbind(pmat, temp.mat)
    
        # fill out the pmat
    for(i in 1:nrow(dat)){
      if(dat$scs[i] == "XO"){
        pmat[i, dat$chroms[i] - 1] <- 1
      }
      if(dat$scs[i] == "XY"){
        pmat[i, (hp + dat$chroms[i] - 2)] <- 1
      }
    }
    
    drops <- c(1:(minChroms-2), hp:(hp + minChroms-3))
    pmat <- pmat[, -drops]
    
    return(pmat)
}



