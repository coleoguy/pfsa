# Terrence Sylvester
# 27 October 2020
# pradakshanas@gmail.com
# Helper functions

GetData <- function(trees="../data/trees/posterior.trees.nex", 
                    data="../data/chrom.data/chroms.csv"){
  # get packages
  library(chromePlus)
  library(diversitree)
  
  # read trees
  trees <- read.nexus(file = trees)
  
  # read data
  chroms <- read.csv(file = data, 
                     as.is = T)
  
  
  # get exact matches
  new.dat <- as.data.frame(matrix(,1,6))
  counter <- 1
  hit.genera <- c()
  hit.family <- c()
  colnames(new.dat) <- c("order", "species", "haploid", "sex.system", "SCS", "notes")
  for(i in 1:length(trees[[1]]$tip.label)){
    current <- trees[[1]]$tip.label[i]
    if(current %in% chroms$species){
      hit <- which(chroms$species == current)
      if(length(hit) > 1) hit <- sample(hit, 1)
      new.dat[counter, ] <- chroms[hit, c(1,4,10,5,6,7)]
      hit.genera <- c(hit.genera, chroms[hit,3])
      hit.family <- c(hit.family, chroms[hit, 2])
      counter <- counter + 1
    }
  }
  # find genera found
  hit.genera <- unique(hit.genera)
  # find genera not found
  unhit.genera <- unique(chroms$genus[!chroms$genus %in% hit.genera])
  # split names at underscores
  tree.taxa <- strsplit(trees[[1]]$tip.label, "_")
  # making a table with genus and species epethet split out to columns
  tree.taxa <- matrix(unlist(tree.taxa), 
                      length(unlist(tree.taxa)), 
                      2, byrow=T)
  for(i in 1:length(unhit.genera)){
    if(unhit.genera[i] %in% tree.taxa[,1]){
      hit <- which(tree.taxa[,1] == unhit.genera[i])
      if(length(hit) > 1) hit <- sample(hit, 1)
      hit <- which(chroms$genus == tree.taxa[hit, 1])
      if(length(hit) > 1) hit <- sample(hit, 1)
      new.dat[counter, ] <- chroms[hit, c(1,4,10,5,6,7)]
      new.dat[counter, 2] <- paste(strsplit(new.dat[counter,2], "_")[[1]][1], "_sp",sep="")
      counter <- counter + 1
      hit.family <- c(hit.family, chroms[hit, 2])
    }
  }
  return(new.dat)
}

# this function will create a Q matrix and a P matrix (probability) given the
# chromosome numbers and sex chromosome systems.
get.matrixes.v2 <- function(haploid.scs = NULL,
                             Neo.sex = NULL,
                             complex = NULL,
                             sex.system = "XY",
                             chrom.range.expansion = NULL,
                             dat = NULL,
                             trees = NULL,
                             def.rates = list(NULL), # polyploidy XYY
                             autosome.as.input = T){
  # parameter definitions
  # chrome.range <- minimum and maximum value of the given chromosome number
  # haploid.scs <- whether males are XO or ZO
  # copmlex <- are complex sex chromosome systems included
  # sex.system <- XY or ZW
  # error checkinig
  # # chrom.range should be present and numeric
  # if(is.null(chrom.range)){
  #   stop("chromosome range is empty. Please provide the upper and lover limit of
  #        the chromosome number")
  # }
  # lets make sure that the correct column names are present in the data frame
  if(sum(c("SpeciesName", "chroms", "scs") %in% colnames(dat)) !=3){
    stop("input data is either missing following columns or they are labled in a
         different way. Please rename your columns to following: 'SpeciesName', 'chroms', 'scs'")
  }
  # make sure data table is present to make the pmatrix
  if(is.null(dat)){
    stop("data table with species names chromosome numer and karyotype should be
         included in order to make the probability matrix")
  }
  # when present haploid.scs should be logical
  if(!is.null(haploid.scs) & !is.logical(haploid.scs)){
    stop("haploid.scs should be logical")
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
    stop("chrom.range.expansion should be either Zero or a single positive numeric value which 
         describes by how much the chromosome range is extended")
  }
  # make sure trees are included to re order the data frame so that it matches 
  # with the order of the tree
  if(is.null(trees)){
    stop("include the list of trees to re order the datatable so that the species
         order of the datatable matches with the soecies order of the tree")
  }
  # make sure input chromosome type is present
  if(!is.logical(autosome.as.input)){
    stop("Please make sure autosome.as.input is logical. T is for autosome counts F for chromosome counts")
  }
  # set default parameters for those that are not provided
  if(is.null(haploid.scs)){
    haploid.scs <- F
  }
  if(is.null(complex)){
    complex <- F
  }
  if(is.null(Neo.sex)){
    Neo.sex <- F
  }
  #set chrom.range
  chrom.range <- range(dat$chroms)
  # define rate parameters
  if(is.null(def.rates$r01)){
    def.rates$r01 <- 1 
  }
  if(is.null(def.rates$r02)){
    def.rates$r02 <- 2
  }
  if(is.null(def.rates$r03)){
    def.rates$r03 <- 3
  }
  if(is.null(def.rates$r04)){
    def.rates$r04 <- 4
  }
  if(is.null(def.rates$r05)){
    def.rates$r05 <- 5
  }
  if(is.null(def.rates$r06)){
    def.rates$r06 <- 6
  }
  if(is.null(def.rates$r07)){
    def.rates$r07 <- 7
  }
  if(is.null(def.rates$r08)){
    def.rates$r08 <- 8
  }
  if(is.null(def.rates$r09)){
    def.rates$r09 <- 9
  }
  if(is.null(def.rates$r10)){
    def.rates$r10 <- 10
  }
  if(is.null(def.rates$r11)){
    def.rates$r11 <- 11
  }
  if(is.null(def.rates$r12)){
    def.rates$r12 <- 12
  }
  if(is.null(def.rates$r13)){
    def.rates$r13 <- 13
  }
  if(is.null(def.rates$r14)){
    def.rates$r14 <- 14
  }
  if(is.null(def.rates$r15)){
    def.rates$r15 <- 15
  }
  if(is.null(def.rates$r16)){
    def.rates$r16 <- 16
  }
  if(is.null(def.rates$r17)){
    def.rates$r17 <- 17
  }
  if(is.null(def.rates$r18)){
    def.rates$r18 <- 18
  }
  if(is.null(def.rates$r19)){
    def.rates$r19 <- 19
  }
  if(is.null(def.rates$r20)){
    def.rates$r20 <- 20
  }
  if(is.null(def.rates$r21)){
    def.rates$r21 <- 21
  }
  if(is.null(def.rates$r22)){
    def.rates$r22 <- 0
  }
  if(is.null(def.rates$r23)){
    def.rates$r23 <- 0
  }
  if(is.null(def.rates$r24)){
    def.rates$r24 <- 0
  }
  if(is.null(def.rates$r25)){
    def.rates$r25 <- 0
  }
  if(is.null(def.rates$r26)){
    def.rates$r26 <- 0
  }
  print("currently polyploidy is not set up in this function. therefore these rate parameters are set to Zero")
  # make it so that the letter case does not matter for sex.system
  sex.system <- toupper(sex.system)
  # also if the minimum chromosome number is 4 or less then the chrom.range.expansion
  # is set to 1
  if(chrom.range[1] <= 3){
    if(chrom.range.expansion > 1){
      print("provided chrom.range.expansion will cause the minimum chromosome number \nto extend beyond 2. Therefore chrom.range.expansion will be set to 1")
      chrom.range.expansion <- 1
    }
  }
  # define total number of chromosome states
  maxChroms <- chrom.range[2] + chrom.range.expansion
  minChroms <- chrom.range[1] - chrom.range.expansion
  ## devider ##
  # this is the parameter that defines the number of parts in the qmatrix.
  # we have 5 parts on our q matrix. 
  # a part for XO / ZO
  # a part for XY / ZW
  # a part for Neo.XY / Neo.ZW
  # a part for XXY / ZZW
  # a part for XYY / ZWW
  devider <- 5
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
                                            paste(2:maxChroms,"XXY",sep=""),
                                            paste(2:maxChroms,"XYY",sep=""))
  }
  if(sex.system == "ZW"){
    row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                            paste(2:maxChroms,"ZW",sep=""),
                                            paste(2:maxChroms,"Neo.ZW",sep=""),
                                            paste(2:maxChroms,"ZWW",sep=""),
                                            paste(2:maxChroms,"ZZW",sep=""))
  }
  # here we define the number of rows and columns assign for each scs
  limiter <- (nrow(qmat)/devider)
  # autosomal fissions and fusions in XO systems
  for(i in 1){
    while (i < limiter) {
      qmat[i, i+1] <- def.rates$r01 # fissions
      qmat[i+1, i] <- def.rates$r02 # fusions
      i <- i + 1
    }
  }
  # autosomal fissions and fusions in XY systems
  for(i in ((limiter+1):(limiter*2))){
    if(i < (limiter*2)){
      qmat[i, i+1] <- def.rates$r03 # fissions
      qmat[i+1, i] <- def.rates$r04 # fusions
    }
  }
  # autosomal fisions and fusions in Neo.XY systems
  for(i in ((limiter*2)+1):(limiter*3)){
    if(i < (limiter*3)){
      qmat[i, i+1] <- def.rates$r05 # fissions
      qmat[i+1, i] <- def.rates$r06 # fusions
    }
  }
  # autosomal fissions and fusions in XXY systems
  for(i in ((limiter*3)+1):(limiter*4)){
    if(i < (limiter*4)){
      qmat[i, i+1] <- def.rates$r07 # fissions
      qmat[i+1, i] <- def.rates$r08 # fusions
    }
  }
  # autosomal fissions and fusions in XYY systems
  for(i in ((limiter*4)+1):(limiter*5)){
    if(i < (limiter*5)){
      qmat[i, i+1] <- def.rates$r09  # fissions
      qmat[i+1, i] <- def.rates$r10 # fusions
    }
  }
  # sex chromosome autosome fusion from XO to XY
  for(i in 2:limiter){
    qmat[i,limiter + i -1] <- def.rates$r11
  }
  # sex chromosome autosome fusion from XY to Neo-XY
  for(i in (limiter + 1):(limiter*2)){
    if(((limiter)+i-1) == limiter*2){
      qmat[i,((limiter)+i-1)] <- 0
    }else{
      qmat[i,((limiter)+i-1)] <- def.rates$r12
    }
  }
  # sex chromosome autosome fusion from XY to XXY when the input is haploid autosome count
  for(i in (limiter + 1):(limiter*2)){
    if(((limiter*2)+(i)-1) == limiter*3){
      qmat[i,((limiter*2)+i-1)] <- 0
    }else{
      qmat[i,((limiter*2)+i-1)] <- def.rates$r13
    }
  }
  # sex chromosome autosome fusion from XY to XYY when the input is haploid autosome count
  if(autosome.as.input == T){
    for(i in (limiter + 1):(limiter*2)){
      if(((limiter*3)+(i)-1) == limiter*4){
        qmat[i,((limiter*3)+i-1)] <- 0
      }else{
        qmat[i,((limiter*3)+i-1)] <- def.rates$r14
      }
    }  
  }
  # sex chromosome autosome fusion from XY to XYY when the input is haploid chromosome count
  if(autosome.as.input == F){
    for(i in (limiter + 1):(limiter*2)){
      if(((limiter*3)+(i)-1) == limiter*4){
        qmat[i,((limiter*3)+i-1)] <- 0
      }else{
        qmat[i,((limiter*3)+i-1)] <- def.rates$r14
      }
    }  
  }
  # transision back to XY from Neo-XY
  for(i in (limiter*2 + 1):(limiter*3)){
    qmat[i,(i-limiter)] <- def.rates$r15
  }
  # sex (X or Z) chromosome fissoin from XY to complex-XY when the input is haploid autosome count
  if(autosome.as.input == T){
    for(i in (limiter + 1):(limiter*2)){
      qmat[i,((limiter*2)+i)] <- def.rates$r16
    }  
  }
  # sex (X or Z) chromosome fissoin from XY to complex-XY when the input is haploid chromosome count
  if(autosome.as.input == F){
    for(i in (limiter + 1):(limiter*2)){
      if(((limiter*2)+i+1) == (limiter * 4 + 1)){
        qmat[i,((limiter*2)+i+1)] <- 0
      }else{
        qmat[i,((limiter*2)+i+1)] <- def.rates$r16 
      }
    }
  }
  # sex chromosome (Y or W) fissoin from XY to complex-XY 
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter*3)+i)] <- def.rates$r17
  }
  # y chromosome loss in XY systems
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,(i-limiter)] <- def.rates$r18
  }
  # Y loss in XYY
  for(i in (limiter*4 + 1):(limiter*5)){
    qmat[i,(i-limiter*3)] <- def.rates$r19
  }
  # X fusion in XXY systems (transision back to XY) when the input is haploid autosome count
  if(autosome.as.input == T){
    for(i in (limiter*3 + 1):(limiter*4)){
      qmat[i,(i-limiter*2)] <- def.rates$r20
    }
  }
  # X fusion in XXY systems (transision back to XY) when the input is haploid chromosome count
  if(autosome.as.input == F){
    for(i in (limiter*3 + 1):(limiter*4)){
      if(((i-limiter*2)-1) == limiter){
        qmat[i,(i-limiter*2)-1] <- 0
      }else{
        qmat[i, (i-limiter*2)-1] <- def.rates$r20
      }
    }
  }
  # transition from XO to XY through capture
  for(i in 1:limiter){
    qmat[i,limiter + i] <- def.rates$r21
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
    print("Tree is a class of multiphylo. Randomly sampling a single tree to make the Pmat")
    tree <- trees[[sample((1:length(trees)),1)]]
  }else{
    tree <- trees
  }
  for(i in 1:nrow(dat)){
    temp.dat[i,] <- dat[dat$SpeciesName == tree$tip.label[i],]
  }
  dat <- temp.dat
  rownames(pmat) <- dat$SpeciesName
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
      if(dat$scs[i] == "XXY"){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "XYY"){
        pmat[i, (limiter*4 + 1 + dat$chroms[i] - 2)] <- 1
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
      if(dat$scs[i] == "ZWW"){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "ZZW"){
        pmat[i, (limiter*4 + 1 + dat$chroms[i] - 2)] <- 1
      }
    }
  }
  # now we remove unwanted parts of the Q matrix and P matrix based on the 
  # input parameters
  comp.scs.xy <- c("XXY", "XYY")
  comp.scs.zw <- c("ZWW", "ZZW")
  dat.scs <-toupper(unique(dat$scs))
  # lets look at each scenario
  # for XY systems
  if(sum(comp.scs.xy %in% dat.scs) == 2){
    type <- 1
  }
  if(sum(comp.scs.xy %in% dat.scs) == 1){
    if(comp.scs.xy[which(comp.scs.xy %in% dat.scs)] == "XYY"){
      type <- 2
    }
    if(comp.scs.xy[which(comp.scs.xy %in% dat.scs)] == "XXY"){
      type <- 3
    }
  }
  # for ZW systems
  if(sum(comp.scs.zw %in% dat.scs) == 2){
    type <- 1
  }
  if(sum(comp.scs.zw %in% dat.scs) == 1){
    if(comp.scs.zw[which(comp.scs.zw %in% dat.scs)] == "ZZW"){
      type <- 2
    }
    if(comp.scs.zw[which(comp.scs.zw %in% dat.scs)] == "ZWW"){
      type <- 3
    }
  }
  # filter out unwanted parts based on each scenario
  if(haploid.scs == T & Neo.sex == T & complex == T){
    if(type == 1){
      qmat <- qmat
      pmat <- pmat  
    }
    if(type == 2){
      qmat <- qmat[-c((limiter*3 + 1):(limiter*4)),-c((limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c((limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c((limiter*4 + 1):(limiter*5)),-c((limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c((limiter*4 + 1):(limiter*5))] 
    }
  }
  if(haploid.scs == T & Neo.sex == T & complex == F){
    qmat <- qmat[-c((limiter*3 + 1):(limiter*5)),-c((limiter*3 + 1):(limiter*5))]
    pmat <- pmat[,-c((limiter*3 + 1):(limiter*5))]
  }
  if(haploid.scs == T & Neo.sex == F & complex == T){
    if(type == 1){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3)),-c((limiter*2 + 1):(limiter*3))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3))]
    }
    if(type == 2){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4)),-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5)),-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
    }
  }
  if(haploid.scs == T & Neo.sex == F & complex == F){
    qmat <- qmat[-c((limiter*2 + 1):(limiter*5)),-c((limiter*2 + 1):(limiter*5))]
    pmat <- pmat[,-c((limiter*2 + 1):(limiter*5))]
  }
  if(haploid.scs == F & Neo.sex == F & complex == T){
    if(type == 1){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3)),-c(1:(limiter),(limiter*2 + 1):(limiter*3))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3))] 
    }
    if(type == 2){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3)(limiter*3 + 1):(limiter*4)),-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5)),-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
    }
  }
  if(haploid.scs == F & Neo.sex == T & complex == F){
    qmat <- qmat[-c(1:(limiter),(limiter*3 + 1):(limiter*5)),-c(1:(limiter),(limiter*3 + 1):(limiter*5))]
    pmat <- pmat[,-c(1:(limiter),(limiter*3 + 1):(limiter*5))]
  }
  if(haploid.scs == F & Neo.sex == F & complex == F){
    qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*5)),-c(1:(limiter),(limiter*2 + 1):(limiter*5))]
    pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*5))]
  }
  # Here we will remove the first portion of the qmatrices for each
  # type of sex chromosome system to match with the given range of chromosome number
  if(haploid.scs == T & Neo.sex == T & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2),(1+limiter*4):(limiter*4+ minChroms-2)) 
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
  }
  if(haploid.scs == T & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))  
  }
  if(haploid.scs == T & Neo.sex == F & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
  }
  if(haploid.scs == T & Neo.sex == F & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  if(haploid.scs == F & Neo.sex == F & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
    }
  }
  if(haploid.scs == F & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  if(haploid.scs == F & Neo.sex == F & complex == F){
    drops <- c(1:(minChroms-2))
  }
  # if the minChrom value is equal to 2 then we shall not reduce the qmatrix
  # if the minChrim value is greater than 2 then qmatrix is reduced to the desierd 
  # limit
  if(minChroms == 2){
    qmat <- qmat
    pmat <- pmat
  }
  if(minChroms > 2){
    qmat <- qmat[-drops, -drops]
    pmat <- pmat[,-drops]
  }
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
  states <- as.data.frame(matrix(data = NA, nrow = length(def.rates), ncol = 3))
  colnames(states) <- c("rate","par", "meaning")
  states$rate <- names(unlist(def.rates))[order(names(unlist(def.rates)))]
  states$par <- unlist(def.rates)[order(names(unlist(def.rates)))]
  if(sex.system == "XY"){
    states$meaning <- c(
      "AA fusion in XO",
      "AA fission in XO",
      "AA fusion in XY",
      "AA fission in XY",
      "AA fusion in Neo.XY",
      "AA fission in Neo.XY",
      "AA fusion in XXY",
      "AA fission in XXY",
      "AA fusion in XYY",
      "AA fission in XYY",
      "SA fusion in XO to XY",
      "SA fusion in XY to Neo.XY",
      "SA fusion in XY to XXY",
      "SA fusion in XY to XYY",
      "transision in Neo.XY to XY",
      "X fission in XY to XXY",
      "Y fission in XY to XYY",
      "Y loss in XY to XO",
      "Y loss in XYY to XY",
      "X fusion in XXY to XY",
      "Y capture in XO to XY",
      "polyploidy in XO",
      "polyploidy in XY",
      "polyploidy in Neo.XY",
      "polyploidy in XXY",
      "polyploidy in XYY")
  }
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

obspSA <- function(counts = NULL,
                   qmat = NULL,
                   states = NULL){
  # get colnames that represent SA, SS and AA fusions
  # AA fusions
  AAfusionColnames.r01 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r03 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r05 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r07 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r09 <- vector(mode = "character", length = nrow(qmat))
  # SA fusions
  SAfusionColnames.r11 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r12 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r13 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r14 <- vector(mode = "character", length = nrow(qmat))
  # SS fusions
  SSfusionColnames <- vector(mode = "character", length = nrow(qmat))
  # get the row and column names of the qmat corrosponding to each type of fission
  for(i in 1:nrow(qmat)){
    if(length(which(qmat[i,] == states$par[1])) != 0){
      AAfusionColnames.r01[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[1])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[3])) != 0){
      AAfusionColnames.r03[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[3])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[5])) != 0){
      AAfusionColnames.r05[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[5])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[7])) != 0){
      AAfusionColnames.r07[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[7])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[9])) != 0){
      AAfusionColnames.r09[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[9])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[11])) != 0){
      SAfusionColnames.r11[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[11])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[12])) != 0){
      SAfusionColnames.r12[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[12])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[13])) != 0){
      SAfusionColnames.r13[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[13])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[14])) != 0){
      SAfusionColnames.r14[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[14])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[20])) != 0){
      SSfusionColnames[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[20])),sep = "")
    }
  }
  
  # for AA and SA fusions all these fusions are devided in multiple objects.
  # combine them to a sinlge object and keep only the unique records
  # AA fusions
  ALLAAfusionColnames <- c(AAfusionColnames.r01,
                           AAfusionColnames.r03,
                           AAfusionColnames.r05,
                           AAfusionColnames.r09)
  ALLAAfusionColnames <- unique(ALLAAfusionColnames)
  # SA fusions
  ALLSAfusionColsnames <- c(SAfusionColnames.r11,
                            SAfusionColnames.r12,
                            SAfusionColnames.r13,
                            SAfusionColnames.r14)
  ALLSAfusionColsnames <- unique(ALLSAfusionColsnames)
  # remove those that are empty
  ALLAAfusionColnames <- ALLAAfusionColnames[ALLAAfusionColnames != ""]
  ALLSAfusionColsnames <- ALLSAfusionColsnames[ALLSAfusionColsnames != ""]
  SSfusionColnames <- SSfusionColnames[SSfusionColnames != ""]
  # get the col number that represent each fusion type
  AAfusionCol <- which(colnames(counts) %in% ALLAAfusionColnames)
  SAfusionCol <- which(colnames(counts) %in% ALLSAfusionColsnames)
  SSfusionCol <- which(colnames(counts) %in% SSfusionColnames)
  # get the number of times each fusion has occured
  AAfusioncounts <- rowSums(counts[, c(AAfusionCol)])
  SAfusioncounts <- rowSums(counts[, c(SAfusionCol)])
  SSfusioncounts <- rowSums(counts[, c(SSfusionCol)])
  # get the observed pfSA
  obspropSA <- SAfusioncounts/(AAfusioncounts+SAfusioncounts+SSfusioncounts)
  # store these in results
  results <- list(obspropSA,
                  AAfusioncounts,
                  SAfusioncounts,
                  SSfusioncounts)
  # give these results appropriate names
  names(results) <- c("obspropSA",
                      "AAfusioncounts",
                      "SAfusioncounts",
                      "SSfusioncounts")
  
  return(results)
}


#### depricated functions ####
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

# this function will create a Q matrix and a P matrix (probability) given the
# chromosome numbers and sex chromosome systems.
get.matrixes.new <- function(chrom.range = NULL,
                             Haplodiploidy = NULL,
                             Neo.sex = NULL,
                             complex = NULL,
                             sex.system = "XY",
                             chrom.range.expansion = NULL,
                             dat = NULL,
                             trees = NULL,
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
                                              r12 = 12,
                                              r12 = 13)){
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
  # define rate parameters
  if(is.null(def.rates$r1)){
    def.rates$r1 <- 1 
  }
  if(is.null(def.rates$r2)){
    def.rates$r2 <- 2
  }
  if(is.null(def.rates$r3)){
    def.rates$r3 <- 3
  }
  if(is.null(def.rates$r4)){
    def.rates$r4 <- 4
  }
  if(is.null(def.rates$r5)){
    def.rates$r5 <- 5
  }
  if(is.null(def.rates$r6)){
    def.rates$r6 <- 6
  }
  if(is.null(def.rates$r7)){
    def.rates$r7 <- 7
  }
  if(is.null(def.rates$r8)){
    def.rates$r8 <- 8
  }
  if(is.null(def.rates$r9)){
    def.rates$r9 <- 9
  }
  if(is.null(def.rates$r10)){
    def.rates$r10 <- 10
  }
  if(is.null(def.rates$r11)){
    def.rates$r11 <- 11
  }
  if(is.null(def.rates$r12)){
    def.rates$r12 <- 12
  }
  if(is.null(def.rates$r13)){
    def.rates$r13 <- 13
  }
  # make it so that the letter case does not matter for sex.system
  sex.system <- toupper(sex.system)
  # if the chrom.range.expansion is not provided then it is set to 1
  if(is.null(chrom.range.expansion)){
    chrom.range.expansion <- 1
  }
  # also if the minimum chromosome number is 4 or less then the chrom.range.expansion
  # is set to 1
  if(chrom.range[1] <= 3){
    if(!(is.null(chrom.range.expansion))){
      if(chrom.range.expansion > 1){
        print("provided chrom.range.expansion will cause the minimum chromosome number
            to extend beyond 2. Therefore chrom.range.expansion will be set to 1")
        chrom.range.expansion <- 1
      }else{
        chrom.range.expansion <- 0  
      }
      
    }
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
  devider <- 5
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
                                            paste(2:maxChroms,"XXY",sep=""),
                                            paste(2:maxChroms,"XYY",sep=""))
  }
  if(sex.system == "ZO"){
    row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                            paste(2:maxChroms,"ZW",sep=""),
                                            paste(2:maxChroms,"Neo.ZW",sep=""),
                                            paste(2:maxChroms,"ZZW",sep=""),
                                            paste(2:maxChroms,"ZWW",sep=""))
  }
  # following are the rate parameters that we define
  # rate1 <- "AA fission",
  # rate2 <- "AA fusion",
  # rate3 <- "SA fusion - XO to XY",
  # rate4 <- "SA fusion - XY to Neo.XY",
  # rate5 <- "transision from Neo.XY to XY",
  # rate6 <- "SA fusion - XY to XXY",
  # rate7 <- "SA fusion - XY to XYY",
  # rate8 <- "X fission - XY to XXY",
  # rate9 <- "Y fission - XY to XYY",
  # rate10 <- "Y loss - XY to XO",
  # rate11 <- "X fission - XXY to XY",
  # rate12 <- "Y loss - XYY to XY"
  
  # here we define the number of rows and columns assign for each scs
  limiter <- (nrow(qmat)/devider)
  # autosomal fissions and fusions in XO systems
  for(i in 1){
    while (i < limiter) {
      qmat[i, i+1] <- def.rates$r1 # fissions
      qmat[i+1, i] <- def.rates$r2 # fusions
      i <- i + 1
    }
  }
  # autosomal fissions and fusions in XY systems
  for(i in ((limiter+1):(limiter*2))){
    if(i < (limiter*2)){
      qmat[i, i+1] <- def.rates$r1 # fissions
      qmat[i+1, i] <- def.rates$r2 # fusions
    }
  }
  # autosomal fisions and fusions in Neo.XY systems
  for(i in ((limiter*2)+1):(limiter*3)){
    if(i < (limiter*3)){
      qmat[i, i+1] <- def.rates$r1 # fissions
      qmat[i+1, i] <- def.rates$r2 # fusions
    }
  }
  # autosomal fissions and fusions in XXY systems
  for(i in ((limiter*3)+1):(limiter*4)){
    if(i < (limiter*4)){
      qmat[i, i+1] <- def.rates$r1 # fissions
      qmat[i+1, i] <- def.rates$r2 # fusions
    }
  }
  # autosomal fissions and fusions in XYY systems
  for(i in ((limiter*4)+1):(limiter*5)){
    if(i < (limiter*5)){
      qmat[i, i+1] <- def.rates$r1 # fissions
      qmat[i+1, i] <- def.rates$r2 # fusions
    }
  }
  # sex chromosome autosome fusion from XO to XY
  for(i in 2:limiter){
    qmat[i,limiter + i -1] <- def.rates$r3
  }
  # sex chromosome autosome fusion from XY to Neo-XY
  for(i in (limiter + 1):(limiter*2)){
    if(((limiter)+i-1) == limiter*2){
      qmat[i,((limiter)+i-1)] <- 0
    }else{
      qmat[i,((limiter)+i-1)] <- def.rates$r4
    }
  }
  # transision back to XY from Neo-XY
  for(i in (limiter*2 + 1):(limiter*3)){
    qmat[i,(i-limiter)] <- def.rates$r5
  }
  # sex chromosome autosome fusion from XY to XXY
  for(i in (limiter + 1):(limiter*2)){
    if(((limiter*2)+(i)-1) == limiter*3){
      qmat[i,((limiter*2)+i-1)] <- 0
    }else{
      qmat[i,((limiter*2)+i-1)] <- def.rates$r6
    }
  }
  # sex chromosome autosome fusion from XY to XYY
  for(i in (limiter + 1):(limiter*2)){
    if(((limiter*3)+(i)-1) == limiter*4){
      qmat[i,((limiter*3)+i-1)] <- 0
    }else{
      qmat[i,((limiter*3)+i-1)] <- def.rates$r7
    }
  }
  # sex chromosome fissoin from XY to complex-XY
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter*2)+i)] <- def.rates$r8
  }
  # sex chromosome fissoin from XY to complex-XY
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,((limiter*3)+i)] <- def.rates$r9
  }
  # y chromosome loss in XY systems
  for(i in (limiter + 1):(limiter*2)){
    qmat[i,(i-limiter)] <- def.rates$r10
  }
  # sex chromosome loss in complex systems (transision back to XY)
  for(i in (limiter*3 + 1):(limiter*4)){
    qmat[i,(i-limiter*2)] <- def.rates$r11
  }
  for(i in (limiter*4 + 1):(limiter*5)){
    qmat[i,(i-limiter*3)] <- def.rates$r12
  }
  # transition from XO to XY through capture
  for(i in 1:limiter){
    qmat[i,limiter + i] <- def.rates$r13
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
    print("Tree is a class of multiphylo. Randomly sampling a single tree to make the Pmat")
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
      if(dat$scs[i] == "XXY"){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "XYY"){
        pmat[i, (limiter*4 + 1 + dat$chroms[i] - 2)] <- 1
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
      if(dat$scs[i] == "ZZW"){
        pmat[i, (limiter*3 + 1 + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] == "ZWW"){
        pmat[i, (limiter*4 + 1 + dat$chroms[i] - 2)] <- 1
      }
    }
  }
  # now we remove unwanted parts of the Q matrix and P matrix based on the 
  # input parameters
  comp.scs <- c("XXY", "XYY")
  dat.scs <-toupper(unique(dat$scs))
  # lets look at each scenario
  if(sum(comp.scs %in% dat.scs) == 2){
    type <- 1
  }
  if(sum(comp.scs %in% dat.scs) == 1){
    if(comp.scs[which(comp.scs %in% dat.scs)] == "XYY"){
      type <- 2
    }
    if(comp.scs[which(comp.scs %in% dat.scs)] == "XXY"){
      type <- 3
    }
  }
  # filter out unwanted parts based on each scenario
  if(Haplodiploidy == T & Neo.sex == T & complex == T){
    if(type == 1){
      qmat <- qmat
      pmat <- pmat  
    }
    if(type == 2){
      qmat <- qmat[-c((limiter*3 + 1):(limiter*4)),-c((limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c((limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c((limiter*4 + 1):(limiter*5)),-c((limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c((limiter*4 + 1):(limiter*5))] 
    }
  }
  if(Haplodiploidy == T & Neo.sex == T & complex == F){
    qmat <- qmat[-c((limiter*3 + 1):(limiter*4)),-c((limiter*3 + 1):(limiter*4))]
    pmat <- pmat[,-c((limiter*3 + 1):(limiter*4))]
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == T){
    if(type == 1){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3)),-c((limiter*2 + 1):(limiter*3))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3))]
    }
    if(type == 2){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4)),-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5)),-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c((limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
    }
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == F){
    qmat <- qmat[-c((limiter*2 + 1):(limiter*5)),-c((limiter*2 + 1):(limiter*5))]
    pmat <- pmat[,-c((limiter*2 + 1):(limiter*5))]
  }
  if(Haplodiploidy == F & Neo.sex == F & complex == T){
    if(type == 1){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3)),-c(1:(limiter),(limiter*2 + 1):(limiter*3))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3))] 
    }
    if(type == 2){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3)(limiter*3 + 1):(limiter*4)),-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*3 + 1):(limiter*4))]
    }
    if(type == 3){
      qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5)),-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
      pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*3),(limiter*4 + 1):(limiter*5))]
    }
  }
  if(Haplodiploidy == F & Neo.sex == T & complex == F){
    qmat <- qmat[-c(1:(limiter),(limiter*2 + 1):(limiter*4)),-c(1:(limiter),(limiter*2 + 1):(limiter*4))]
    pmat <- pmat[,-c(1:(limiter),(limiter*2 + 1):(limiter*4))]
  }
  # Here we will remove the first portion of the qmatrices for each
  # type of sex chromosome system to match with the given range of chromosome number
  if(Haplodiploidy == T & Neo.sex == T & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2),(1+limiter*4):(limiter*4+ minChroms-2)) 
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
  }
  if(Haplodiploidy == T & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))  
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2),(1+limiter*3):(limiter*3+ minChroms-2))
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
  }
  if(Haplodiploidy == T & Neo.sex == F & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  if(Haplodiploidy == F & Neo.sex == F & complex == T){
    if(type == 1){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2),(1+limiter*2):(limiter*2+ minChroms-2))
    }
    if(type == 2){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
    }
    if(type == 3){
      drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
    }
  }
  if(Haplodiploidy == F & Neo.sex == T & complex == F){
    drops <- c(1:(minChroms-2), (1+limiter):(limiter+ minChroms-2))
  }
  # if the minChrom value is equal to 2 then we shall not reduce the qmatrix
  # if the minChrim value is greater than 2 then qmatrix is reduced to the desierd 
  # limit
  if(minChroms == 2){
    qmat <- qmat
    pmat <- pmat
  }
  if(minChroms > 2){
    qmat <- qmat[-drops, -drops]
    pmat <- pmat[,-drops]
  }
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
  states <- as.data.frame(matrix(data = NA, nrow = length(def.rates), ncol = 3))
  colnames(states) <- c("rate","par", "meaning")
  states$rate <- names(unlist(def.rates))
  states$par <- unlist(def.rates)
  states$meaning <- c("AA fission",
                      "AA fusion",
                      "SA fusion - XO to XY",
                      "SA fusion - XY to Neo.XY",
                      "transision from Neo.XY to XY",
                      "SA fusion - XY to XXY",
                      "SA fusion - XY to XYY",
                      "X fission - XY to XXY",
                      "Y fission - XY to XYY",
                      "Y loss - XY to XO",
                      "X fusion - XXY to XY",
                      "Y loss - XYY to XY",
                      "Y capture - XO to XY")
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


