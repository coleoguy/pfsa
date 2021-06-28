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
getMatrixes <- function(haploid.scs = NULL,
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
  # complex <- are complex sex chromosome systems included
  # sex.system <- XY or ZW
  # error checking
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
  
  # following are the rate parameters and what they mean
  # r1 = 1,    # AA fission  XO
  # r2 = 2,    # AA fusion XO
  # r3 = 3,    # AA fission  XY
  # r4 = 4,    # AA fusion XY
  # r5 = 5,    # AA fission  Neo.XY
  # r6 = 6,    # AA fusion Neo.XY 
  # r7 = 7,    # AA fission  XXY
  # r8 = 8,    # AA fusion XXY
  # r9 = 9,    # AA fission  XYY
  # r10 = 10,  # AA fusoin XYY
  # r11 = 11,  # SA fusion  XO -> XY
  # r12 = 12,  # SA fusion  XY -> Neo.XY
  # r13 = 13,  # SA fusion  XY -> XXY
  # r14 = 14,  # SA fusion  XY -> XYY
  # r15 = 15,  # transision Neo.XY -> XY
  # r16 = 16,  # X fission  XY -> XXY
  # r17 = 17,  # Y fission  XY -> XYY
  # r18 = 18,  # Y loss     XY -> XO
  # r19 = 19,  # Y loss     XYY -> XY
  # r20 = 20,  # X fusion   XXY -> XY
  # r21 = 21,  # Y capture  XO -> XY
  # r22 = 22,  # polyploidy XO
  # r23 = 23,  # polyploidy XY
  # r24 = 24,  # polyploidy Neo.XY
  # r25 = 25,  # polyploidy XXY
  # r26 = 26,  # polyploidy XYY
  # r27 = 27,  # translocation XO -> XXY (White(1973), Animal cytology and evolution)
  
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
  if(is.null(def.rates$r27)){
    def.rates$r27 <- 27
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
  # translation from XO to XXY when the input is haploid autosome count 
  if(autosome.as.input == T){
    for(i in 1:limiter){
      if((limiter*3 + i - 1) == limiter*3 ){
        qmat[i,limiter*3 + i - 1] <- 0
      }else{
        qmat[i,limiter*3 + i - 1] <- def.rates$r27
      }
    }
  }
  # translation from XO to XXY  when the input is haploid chromosome count 
  if(autosome.as.input == F){
    for(i in 1:limiter){
      qmat[i,limiter*3 + i] <- def.rates$r27
    }
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
  # make a new table to define all the parameters in the qmatrix
  states <- as.data.frame(matrix(data = NA, nrow = length(def.rates), ncol = 3))
  colnames(states) <- c("rate","par", "meaning")
  states$rate <- names(unlist(def.rates))[order(names(unlist(def.rates)))]
  states$par <- unlist(def.rates)[order(names(unlist(def.rates)))]
  if(sex.system == "XY"){
    states$meaning <- c(
      "AA fission in XO",
      "AA fusion in XO",
      "AA fissio in XY",
      "AA fusion in XY",
      "AA fission in Neo.XY",
      "AA fusion in Neo.XY",
      "AA fission in XXY",
      "AA fusion in XXY",
      "AA fission in XYY",
      "AA fusion in XYY",
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
      "polyploidy in XYY",
      "translocation in XO to XXY")
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

# This cunction will calculate the observed proportion of sex chromosome 
# autosome fusions
obspSA <- function(counts = NULL,
                   qmat = NULL,
                   states = NULL){
  # get colnames that represent SA, SS and AA fusions
  # AA fusions
  AAfusionColnames.r02 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r04 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r06 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r08 <- vector(mode = "character", length = nrow(qmat))
  AAfusionColnames.r10 <- vector(mode = "character", length = nrow(qmat))
  # SA fusions
  SAfusionColnames.r11 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r12 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r13 <- vector(mode = "character", length = nrow(qmat))
  SAfusionColnames.r14 <- vector(mode = "character", length = nrow(qmat))
  # SS fusions
  SSfusionColnames <- vector(mode = "character", length = nrow(qmat))
  # get the row and column names of the qmat corrosponding to each type of fission
  for(i in 1:nrow(qmat)){
    if(length(which(qmat[i,] == states$par[2])) != 0){
      AAfusionColnames.r02[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[2])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[4])) != 0){
      AAfusionColnames.r04[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[4])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[6])) != 0){
      AAfusionColnames.r06[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[6])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[8])) != 0){
      AAfusionColnames.r08[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[8])),sep = "")
    }
    if(length(which(qmat[i,] == states$par[10])) != 0){
      AAfusionColnames.r10[i] <- paste(rownames(qmat)[i],",",names(which(qmat[i,] == states$par[10])),sep = "")
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
  ALLAAfusionColnames <- c(AAfusionColnames.r02,
                           AAfusionColnames.r04,
                           AAfusionColnames.r06,
                           AAfusionColnames.r08,
                           AAfusionColnames.r10)
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

constrainQmat <-  function (qmat, lik) {
  # pad determines how many digits there are in the states. for example if there
  # are 100 states then fist state will be written as 001 instead of just 1. 
  # this is done because diversitree adds the zero digits in front of the states
  # depending in the number of states
  if (ncol(qmat) < 100) 
    pad <- 2
  if (ncol(qmat) >= 100) 
    pad <- 3
  if (ncol(qmat) < 10) 
    pad <- 1
  # lets store the input qmatrix in a new object called parMat
  parMat <- qmat  
  # rename the colnames and row names of the parMat to numerals.
  colnames(parMat) <- sprintf(paste("%0", pad, "d", sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  # make a new table which have three columns. col 1 and col 2 have the state
  # names.
  rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat),3))
  rate.table[, 1] <- rep(as.character(row.names(parMat)), each = ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[2],]
  # now fill the third column of the rate table. this will define which parameters
  # we will use to constrain the initial likelihood function
  rate.table[rate.table[, 3] == 1, 3] <- "par01"
  rate.table[rate.table[, 3] == 2, 3] <- "par02"
  rate.table[rate.table[, 3] == 3, 3] <- "par03"
  rate.table[rate.table[, 3] == 4, 3] <- "par04"
  rate.table[rate.table[, 3] == 5, 3] <- "par05"
  rate.table[rate.table[, 3] == 6, 3] <- "par06"
  rate.table[rate.table[, 3] == 7, 3] <- "par07"
  rate.table[rate.table[, 3] == 8, 3] <- "par08"
  rate.table[rate.table[, 3] == 9, 3] <- "par09"
  rate.table[rate.table[, 3] == 10, 3] <- "par10"
  rate.table[rate.table[, 3] == 11, 3] <- "par11"
  rate.table[rate.table[, 3] == 12, 3] <- "par12"
  rate.table[rate.table[, 3] == 13, 3] <- "par13"
  rate.table[rate.table[, 3] == 14, 3] <- "par14"
  rate.table[rate.table[, 3] == 15, 3] <- "par15"
  rate.table[rate.table[, 3] == 16, 3] <- "par16"
  rate.table[rate.table[, 3] == 17, 3] <- "par17"
  rate.table[rate.table[, 3] == 18, 3] <- "par18"
  rate.table[rate.table[, 3] == 19, 3] <- "par19"
  rate.table[rate.table[, 3] == 20, 3] <- "par20"
  rate.table[rate.table[, 3] == 21, 3] <- "par21"
  rate.table[rate.table[, 3] == 22, 3] <- "par22"
  rate.table[rate.table[, 3] == 23, 3] <- "par23"
  rate.table[rate.table[, 3] == 24, 3] <- "par24"
  rate.table[rate.table[, 3] == 25, 3] <- "par25"
  rate.table[rate.table[, 3] == 26, 3] <- "par26"
  rate.table[rate.table[, 3] == 27, 3] <- "par27"
  # make the formulae
  formulae <- vector(mode = "character", length = nrow(rate.table))
  for (i in 1:nrow(rate.table)) {
    formulae[i] <- paste("q",
                         rate.table[i, 1], 
                         rate.table[i,2],
                         " ~ ", 
                         rate.table[i, 3], 
                         collapse = "", 
                         sep = "")
    lambda <- mu <- vector()
    for(i in 1:ncol(parMat)){
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], 
                                " ~ lambda", sep = ""))
      mu <- c(mu, paste("mu", colnames(parMat)[i], 
                        " ~ mu", sep = ""))
    }
  }
  extras <- c("par01",
              "par02",
              "par03",
              "par04",
              "par05",
              "par06",
              "par07",
              "par08",
              "par09",
              "par10",
              "par11",
              "par12",
              "par13",
              "par14",
              "par15",
              "par16",
              "par17",
              "par18",
              "par19",
              "par20",
              "par21",
              "par22",
              "par23",
              "par24",
              "par25",
              "par26",
              "par27",
              "lambda",
              "mu")
  lik.con <- constrain(lik, formulae = c(formulae, lambda, mu), extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(qmat)
  return(lik.con)
}

fillQmat <- function(post.burnin, qmat){
  # get the parameters of interest
  pars <- colnames(post.burnin)[-c(1, (length(colnames(post.burnin))-2):length(colnames(post.burnin)))]
  # get the numerical value of the parameters
  pars.values <- gsub(pattern = "par", replacement = "",x = pars)
  # fill in the rate values at the necessary places in the qmatrix
  for(i in 1:length(pars.values)){
    qmat[qmat == as.numeric(pars.values)[i]] <- mean(post.burnin[[which(colnames(post.burnin) == pars[i])]])
  }
  # make row sums to zero
  # rows of the qmatrix should sum to 0
  for(i in 1:nrow(qmat)){
    qmat[i,i] <- -(sum(qmat[i,]))
  }
  return(qmat)
}

# this function is under development
getPfsaTab <- function(karyotypes, autosomes.as.input){
  pfSA.tab <- as.data.frame(matrix(data = NA, nrow = length(karyotypes), ncol = 2))
  colnames(pfSA.tab) <- c("state", "pfsa")
  # if autosomes are the inputs then make necessary adjustments to calculate the
  # expected pfsa
  if(autosomes.as.input == T){
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
    # if chromosomes are the inputs then make necessary adjustments to calculate the
    # expected pfsa
    # assuming female data are given
  }else{
    for(i in 1:length(karyotypes)){
      # for XO / ZO
      if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XO", "ZO")){
        pfSA.tab$state[i] <- karyotypes[i]
        pfSA.tab$pfsa[i] <- Pfsa(Da = (as.numeric(gsub(pattern = "[A-z]",
                                                       x= karyotypes[i] ,
                                                       replacement = "")) * 2) - 2,
                                 scs = "XO")
      }
      # for XY / ZW
      if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XY", "ZW")){
        pfSA.tab$state[i] <- karyotypes[i]
        pfSA.tab$pfsa[i] <- Pfsa(Da = (as.numeric(gsub(pattern = "[A-z]",
                                                       x= karyotypes[i] ,
                                                       replacement = "")) * 2) - 2,
                                 scs = "XY")
      }
      # for Neo.XY / Neo.ZW
      if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("Neo.XY", "Neo.ZW")){
        pfSA.tab$state[i] <- karyotypes[i]
        pfSA.tab$pfsa[i] <- Pfsa(Da = (as.numeric(gsub(pattern = "[A-z]",
                                                       x= karyotypes[i] ,
                                                       replacement = "")) * 2)-2,
                                 scs = "XY")
      }
      # for XXY / ZZW
      if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XXY","ZZW")){
        pfSA.tab$state[i] <- karyotypes[i]
        pfSA.tab$pfsa[i] <- Pfsa(Da = (as.numeric(gsub(pattern = "[A-z./]",
                                                       x= karyotypes[i] ,
                                                       replacement = "")) * 2) - 4,
                                 scs = "XXY")
      }
      # for XYY / ZWW
      if(gsub(pattern = "[0-9]",x= karyotypes[i], replacement = "") %in% c("XYY", "ZWW")){
        pfSA.tab$state[i] <- karyotypes[i]
        pfSA.tab$pfsa[i] <- Pfsa(Da = (as.numeric(gsub(pattern = "[A-z./]",
                                                       x= karyotypes[i] ,
                                                       replacement = "")) * 2) - 2,
                                 scs = "XYY")
      }
    } 
  }
  return(pfSA.tab)
}

sp.matches <- function(dat, trees, use.sub.species = F){
  dat.new <- dat
  if(use.sub.species == F){
    for(i in 1:nrow(dat)){
      dat.new$binomial[i] <- paste(unlist(strsplit(dat.new$binomial[i], split = "_"))[1],
                                   unlist(strsplit(dat.new$binomial[i], split = "_"))[2], 
                                   sep = "_")
    }
  }
  if(class(trees) == "multiPhylo"){
  tree <- trees[[1]]
  }
  
  if(class(trees) == "phylo"){
    tree <- trees
  }
  # get the species level matches
  dat.new.s.m <- dat.new[dat.new$binomial %in%  tree$tip.label,]
  
  # from these species level matches if there are multiple hits for a single
  # species get only 1 hit
  unique.sp.matches <- unique(dat.new.s.m$binomial)
  # make an empty data frame to store values unique species matches
  dat.new.s.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.s.m$binomial)),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.s.m.unique) <- colnames(dat.new)
  # fill in the data frame
  for(i in 1:length(unique.sp.matches)){
    dat.new.temp <- dat.new.s.m[dat.new.s.m$binomial %in% unique.sp.matches[i],]
    dat.new.s.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
  }
  # subset the original dataset by the species which had no hit to the tree 
  dat.new.unhit <- dat.new[!(dat.new$binomial %in% tree$tip.label),]
  # get the genera names of the species had matches to the tree dataset
  hit.genera <- c(unique(dat.new.s.m$Genus))
  # get the genera names that have no species level matches to the tree dataset
  unhit.gen <- unique(dat.new$Genus[!(dat.new$Genus %in% hit.genera)])
  # from the dataset that have species which have no species level matches remove
  # species that share the genus with the species that have a match with the tree
  dat.new.unhit <- dat.new.unhit[dat.new.unhit$Genus %in% unhit.gen,]
  # get the names of the species from the tree that have no species level match
  tree.unhit <- tree$tip.label[!(tree$tip.label %in% unique.sp.matches)]
  # get the genus names of the species hat have no species level match with the
  # trait dataset
  tree.unhit.gen <- c()
  for(i in 1:length(tree.unhit)){
    tree.unhit.gen[i] <- unlist(strsplit(tree.unhit[i], split = "_"))[[1]]
  }
  # get the index of the genera which have species level matches
  tree.hit.gen.num <- which(tree.unhit.gen %in% hit.genera)
  # get the names of the genera that appear on the phylogeny that have no species level match
  tree.true.unhit.gen <- tree.unhit.gen[!(tree.unhit.gen %in% hit.genera)]
  # get genera level matches
  dat.new.g.m <- dat.new.unhit[dat.new.unhit$Genus %in%  tree.true.unhit.gen,]
  # get the names of unique genera level matches
  unique.gn.matches <- unique(dat.new.g.m$Genus)
  # from these genera level matches if there are multiple hits for a single
  # genera get only 1 hit
  dat.new.g.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.g.m$Genus)),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.g.m.unique) <- colnames(dat.new)
  # fill in the data frame
  for(i in 1:length(unique.gn.matches)){
    dat.new.temp <- dat.new.g.m[dat.new.g.m$Genus %in% unique.gn.matches[i],]
    dat.new.g.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
    dat.new.g.m.unique$binomial[i] <- paste(dat.new.g.m.unique$Genus[i], "_sp.", sep = "")
  }
  # get the index of genera level matches
  tree..unhit.hit.gen.num <- which(tree.unhit.gen %in% dat.new.g.m.unique$Genus)
  # get the names of the species as they appear on the phylogeny which have a
  # genera level match
  tree..unhit.hit.sp.names <- tree.unhit[tree..unhit.hit.gen.num]
  # this will be used to rename the names in the tree
  # make a data frame to store the species names of those species that have a
  # genera level match
  unhit.temp.dat.new <- as.data.frame(matrix(data = NA,
                                             nrow = length(tree..unhit.hit.sp.names),
                                             ncol = 2))
  # give column names
  colnames(unhit.temp.dat.new) <- c("gen", "sp")
  # fill in the species names 
  unhit.temp.dat.new$sp <- tree..unhit.hit.sp.names
  # give them the proper naming for genera level match
  for(i in 1:length(tree..unhit.hit.sp.names)){
    unhit.temp.dat.new$gen[i] <- unlist(strsplit(tree..unhit.hit.sp.names[i], split = "_"))[[1]]
  }
  # some genera are duplicated. get the unique genera name and randomly sample
  # a single species to represent these genera
  tree.gen.matches <- as.data.frame(matrix(data = NA,
                                           nrow = length(unique(unhit.temp.dat.new$gen)),
                                           ncol = 2))
  # give column names
  colnames(tree.gen.matches) <- c("gen", "sp")
  # fill in the data table  
  tree.gen.matches$gen <- unique(unhit.temp.dat.new$gen)
  for(i in 1:nrow(tree.gen.matches)){
    tree.gen.matches$sp[i] <- sample(unhit.temp.dat.new$sp[unhit.temp.dat.new$gen == tree.gen.matches$gen[i]],1)
  }
  # rename the genera
  tree.gen.matches$gen <- paste(tree.gen.matches$gen, "_sp.", sep = "")
  # for tree name correction
  tree.names <- as.data.frame(matrix(data = NA,
                                     nrow = (nrow(dat.new.s.m.unique) + nrow(dat.new.g.m.unique)),
                                     ncol = 2))
  colnames(tree.names) <- c("species_name", "name_on_tree")
  # fill in
  tree.names$species_name[1:nrow(dat.new.s.m.unique)] <- tree.names$name_on_tree[1:nrow(dat.new.s.m.unique)] <- dat.new.s.m.unique$binomial
  tree.names$species_name[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$gen
  tree.names$name_on_tree[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$sp
  # get the finalized dataset
  finaldat.new <- rbind(dat.new.s.m.unique, dat.new.g.m.unique)
  # get the results
  results <- list(finaldat.new, tree.names)
  names(results) <- c("chroms", "name_corrections")
  return(results)
}

# this function will combine XXY and XYY systems into one group
# must work on this function inprder to apply it to ZY systems

combineComplexSCS <- function(qmat,pmat){
  scs <- gsub(pattern = "[0-9]", x = colnames(qmat), replacement = "")
  scs <- unique(scs)
  
  limiter <- ncol(qmat)/length(scs)
  
  
  hit.qmat.drop <- which(gsub(pattern = "[0-9]", x = colnames(qmat), replacement = "") == "XYY")
  hit.pmat.drop <- which(gsub(pattern = "[0-9]", x = colnames(pmat), replacement = "") == "XYY")
  qmat <- qmat[-hit.qmat.drop, -hit.qmat.drop]
  colnames(qmat) <- gsub(pattern = "XXY", x = colnames(qmat), replacement = "Complex.XY")
  rownames(qmat) <- gsub(pattern = "XXY", x = rownames(qmat), replacement = "Complex.XY")
  
  for(i in 1:nrow(pmat)){
    if(gsub(pattern = "[0-9]", x = names(pmat[i,][pmat[i,] == 1]), replacement = "") == "XYY"){
      hit.pmat.drop.XYY <- which(pmat[i,] == 1)
      pmat[i, (hit.pmat.drop.XYY - limiter)] <- 1
    }
  }
  pmat <- pmat[,-hit.pmat.drop]
  colnames(pmat) <- gsub(pattern = "XXY", x = colnames(pmat), replacement = "Complex.XY")
  
  results <- list(qmat, pmat)
  names(results) <- c("qmat", "pmat")
  return(results)
}



