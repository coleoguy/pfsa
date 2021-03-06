# Terrence Sylvester
# 27 October 2020
# pradakshanas@gmail.com
# Helper functions

# this function will create a Q matrix and a P matrix (probability) given the
# chromosome numbers and sex chromosome systems.
get.matrixes <- function(chrom.range = NULL,
                         Haplodiploidy = NULL,
                         Neo.sex = NULL,
                         complex = NULL,
                         sex.system = "XY",
                         chrom.range.expansion = NULL,
                         dat = NULL,
                         trees = NULL){
  # parameter definitions
  # chrome.range <- minimum and maximum value of the given chromosome number
  # copmlex <- are complex sex chromosome systems included
  # sex.system <- XY or ZW
  
  # error checkinig
  # chrom.range should be present and numeric
  if(is.null(chrom.range)){
    stop("chromosome range is empty. Please provide the upper and lover limit of
         the chromosome number")
  }
  if(!is.numeric(chrom.range)){
    stop("Chromosome range should be numeric")
  }
  # when present Haplodiploidy should be logical
  if(!is.null(Haplodiploidy) & !is.logical(Haplodiploidy)){
    stop("Haplodiploidy should be logical")
  }
  # when present Neo.sex should be logical
  if(!is.null(Neo.sex) & !is.logical(Neo.sex)){
    stop("Neo.sex should be logical")
  }
  # when present complex should be logical
  if(!is.null(complex) & !is.logical(complex)){
    stop("complex should be logical")
  }
  # sex.system should be present and numeric
  if(is.null(sex.system)){
    stop("sex.system  is empty. Please provide which is the sex limited chromosome
         in the form of 'XY' or 'ZW'")
  }
  if(!is.character(sex.system)){
    stop("sex.system should be character string")
  }
  # chrom.range.expansion should be numeric when present
  if(!is.null(chrom.range.expansion) & !is.numeric(chrom.range.expansion)){
    stop("chrom.range.expansion should be a single positive numeric value which 
         describes by how much the chromosome range is extended")
  }
  # make sure data table is present to make the pmatrix
  if(is.null(dat)){
    stop("data table with species names chromosome numer and karyotype should be
         included in order to make the probability matrix")
  }
  # make sure trees are included to re order the data frame so that it matches 
  # with the order of the tree
  if(is.null(trees)){
    stop("include the list of trees to re order the datatable so that the species
         order of the datatable matches with the soecies order of the tree")
  }
  # set default parameters for those that are not provided
  if(is.null(Haplodiploidy)){
    Haplodiploidy <- F
  }
  if(is.null(complex)){
    complex <- F
  }
  if(is.null(Neo.sex)){
    Neo.sex <- F
  }
  # make it so that the letter case does not matter for sex.system
  sex.system <- toupper(sex.system)
  # if the chrom.range.expansion is not provided then it is set to 1
  if(is.null(chrom.range.expansion)){
    chrom.range.expansion <- 1
  }
  # also if the minimum chromosome number is 4 or less then the chrom.range.expansion
  # is set to 1
  if(chrom.range[1] <= 4){
    chrom.range.expansion <- 1
  }
  # define total number of chromosome states
  maxChroms <- chrom.range[2] + chrom.range.expansion
  minChroms <- chrom.range[1] - chrom.range.expansion
  
  ## devider ##
  # this is the parameter that defines the number of parts in the qmatrix.
  # we have  4 parts on our q matrix. 
  # a part for XO / ZO
  # a part for XY / ZW
  # a part for Neo.XY / Neo.ZW
  # a part for complex XY / complex ZW
  devider <- 4
  nChroms <- (maxChroms -1) * devider
  
  # here we will make the complete q matrix wichi will include all 4 parts and
  # depending on the imput data we will remove those parts that are not needed
  
  # make the initial q matrix
  qmat <- matrix(data = 0,
                 nrow = nChroms,
                 ncol = nChroms)
  # define rownames and colnames
  # although row names and col names starts from 1, chromsome counts starts from
  # 2
  if(sex.system == "XY"){
    row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                            paste(2:maxChroms,"XY",sep=""),
                                            paste(2:maxChroms,"Neo.XY",sep=""),
                                            paste(2:maxChroms,"XXY/XYY",sep=""))
  }
  if(sex.system == "ZO"){
    row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                            paste(2:maxChroms,"ZW",sep=""),
                                            paste(2:maxChroms,"Neo.ZW",sep=""),
                                            paste(2:maxChroms,"ZZW/ZWW",sep=""))
  }
  # following are the rate parameters that we define
  # rate1 <- autosomal fusion 
  # rate2 <- autosomal fission 
  # rate3 <- sex chromosome autosome fusions
  # rate4 <- transition back to XY from Neo-XY
  # rate5 <- Y chromosome loss in XO
  # rate6 <- Y/X chromosome loss in complex systems
  
  # here we define the number of rows and columns assign for each scs
  limiter <- (nrow(qmat)/devider)
  
  # autosomal fissions and fusions in XO systems
  for(i in 1){
    while (i < limiter) {
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
      i <- i + 1
    }
  }
  # autosomal fissions and fusions in XY systems
  for(i in ((limiter+1):(limiter*2))){
    if(i < (limiter*2)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # autosomal fisions and fusions in Neo.XY systems
  for(i in ((limiter*2)+1):(limiter*3)){
    if(i < (limiter*3)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # autosomal fissions and fusions in complex systems
  for(i in ((limiter*3)+1):(limiter*4)){
    if(i < (limiter*4)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # sex chromosome autosome fusion from XO to XY
  for(i in 2:limiter){
    qmat[i,limiter + i -1] <- 3
  }
  # sex chromosome autosome fusion from XY to Neo-XY
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter)+i-1)] <- 3
  }
  # transision back to XY from Neo-XY
  for(i in (limiter*2 + 1):(limiter*3)){
    qmat[i,(i-limiter)] <- 4
  }
  # sex chromosome autosome fusion from XY to complex-XY
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter*2)+i-1)] <- 3
  }
  # sex chromosome fissoin from XY to complex-XY
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter*2)+i)] <- 7
  }
  # y chromosome loss in XY systems
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,(i-limiter)] <- 5
  }
  # sex chromosome loss in complex systems (transision back to XY)
  for(i in (limiter*3 + 1):(limiter*4)){
    qmat[i,(i-limiter*2)] <- 6
  }
  
  # make the probability matrix
  pmat <- qmat[NULL,]
  # make a temporary matrix
  temp.mat <- matrix(data = 0,
                     nrow = nrow(dat),
                     ncol = ncol(pmat))
  # combine the new matrix with the pmat
  pmat <- rbind(pmat, temp.mat)
  
  # sort the order of taxa in the data table to match with that of the tree
  temp.dat <- dat
  
  if(class(trees) == "multiPhylo"){
    tree <- trees[[sample((1:length(trees)),1)]]
  }else{
    tree <- trees
  }
  
  for(i in 1:nrow(dat)){
    temp.dat[i,] <- dat[dat$SpecisName == tree$tip.label[i],]
  }
  dat <- temp.dat
  rownames(pmat) <- dat$SpecisName
  
  # fill the p matrix 
  if(sex.system == "XY"){
    for(i in 1:nrow(dat)){
      if(dat$scs[i] == "XO"){
        pmat[i, dat$chroms[i] - 1] <- 1
      }
      if(dat$scs[i] == "XY"){
        pmat[i, (limiter+1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "Neo.XY"){
        pmat[i, (limiter*2 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] %in% c("XYY","XXY")){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 3)] <- 1
      }
    }
  }
  # fill out the pmat when the sex limited chromosome is W
  if(sex.system == "ZW"){
    for(i in 1:nrow(dat)){
      if(dat$scs[i] == "ZO"){
        pmat[i, dat$chroms[i] - 1] <- 1
      }
      if(dat$scs[i] == "ZW"){
        pmat[i, (limiter+1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "Neo.ZW"){
        pmat[i, (limiter*2 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] %in% c("ZZW","ZWW")){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 3)] <- 1
      }
    }
  }
  # now we remove unwanted parts of the Q matrix and P matrix based on the 
  # input parameters
  if(Haplodiploidy == T & Neo.sex == T & complex == T){
    qmat <- qmat
    pmat <- pmat
  }
  if(Haplodiploidy == T & Neo.sex == T & complex == F){
    qmat <- qmat[-c((limiter*3 + 1):(limiter*4)),-c((limiter*3 + 1):(limiter*4))]
    pmat <- pmat[,-c((limiter*3 + 1):(limiter*4))]
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == T){
    qmat <- qmat[-c((limiter*2 + 1):(limiter*3)),-c((limiter*2 + 1):(limiter*3))]
    pmat <- pmat[,-c((limiter*2 + 1):(limiter*3))]
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == F){
    qmat <- qmat[-c((limiter*2 + 1):(limiter*4)),-c((limiter*2 + 1):(limiter*4))]
    pmat <- pmat[,-c((limiter*2 + 1):(limiter*4))]
  }
  if(Haplodiploidy == F & Neo.sex == F & complex == T){
    qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3)),-c(1:(limiter),(limiter*2 + 1):(limiter*3))]
    pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3))]
  }
  if(Haplodiploidy == F & Neo.sex == T & complex == F){
    qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*4)),-c(1:(limiter),(limiter*2 + 1):(limiter*4))]
    pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*4))]
  }
  
  # Here we will remove the first portion of the qmatrices for each
  # type of sex chromosome system to match with the given range of chromosome number
  if(Haplodiploidy == T & Neo.sex == T & complex == T){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
  }
  if(Haplodiploidy == T & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))  
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == T){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  if(Haplodiploidy == F & Neo.sex == F & complex == T){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  if(Haplodiploidy == F & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  qmat <- qmat[-drops, -drops]
  pmat <- pmat[,-drops]
  
  # assign numerical value for each state
  # states <-as.data.frame(matrix(data = NA,
                                # nrow = nrow(qmat),
                                # ncol = 2))
  # colnames(states) <- c("sim.state", "karyotype")
  # make the state names so that they reflect 1:number of states
  # states$sim.state <- 1:nrow(qmat)
  # states$karyotype <- rownames(qmat)
  # rename the row and colnames for qmat
  # rownames(qmat) <- colnames(qmat) <- colnames(pmat) <- 1:nrow(qmat)
  
  # make a new table to define all the parameters in the qmatrix
  states <- as.data.frame(matrix(, nrow = 7, ncol = 2))
  colnames(states) <- c("par", "meaning")
  states$par <- 1:7
  states$meaning <- c("AA fission",
                      "AA fusion",
                      "SA fusion",
                      "transision from Neo-XY to XY",
                      "Y loss in XY",
                      "transision from complex-XY to XY",
                      "SS fissions")
  
  # make a vector to hold the karyotypes
  kar <- colnames(qmat)
  # store results and give them appropriate names
  results <- list(qmat,
                  pmat,
                  dat,
                  trees,
                  states,
                  kar)
  names(results) <- c("qmat",
                      "pmat",
                      "dat",
                      "trees",
                      "states",
                      "karyotypes")
  return(results)
}
