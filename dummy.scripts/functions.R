# Terrence Sylvester
# 26th March 2020
# pradakshanas@gmail.com
# make transition matrix

# make a matrix with 100 chromosomes on each state. Here we will be using
# sex chromosome systems as the states. These states are XO and XY.

get.Qmatrix <- function(chrom.range = NULL, 
                        complex = F, sex.system = "XY"){
  # parameter definitions
  # chrome.range <- minimum and maximum value of the given chromosome number
  # copmlex <- are complex sex chromosome systems included
  # sex.system <- XY or ZW
  
  # make it so that the letter case does not matter for sex.system
  sex.system <- toupper(sex.system)
  
  # define total number of chromosome states
  maxChroms <- chrom.range[2] + 1
  minChroms <- chrom.range[1] - 1
  if(complex == T){
    nChroms <- (maxChroms -1) * 3
  }else{
    nChroms <- (maxChroms -1) * 2
  }
  
  # make the initial q matrix
  qmat <- matrix(data = 0,
                 nrow = nChroms,
                 ncol = nChroms)
  
  # define rownames and colnames
  # although row names and col names starts from 1, chromsome counts starts from
  # 2
  if(complex == T){
    if(sex.system == "XY"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                              paste(2:maxChroms,"XY",sep=""),
                                              paste(2:maxChroms,"XXY/XYY",sep=""))
    }
    if(sex.system == "ZW"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                              paste(2:maxChroms,"ZW",sep=""),
                                              paste(2:maxChroms,"ZZW/ZWW",sep=""))
    }
  }else{
    if(sex.system == "XY"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                              paste(2:maxChroms,"XY",sep=""))
    }
    if(sex.system == "ZW"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                              paste(2:maxChroms,"ZW",sep=""))
    }
  }
  # following are the rate parameters that we define
  # rate1 <- autosomal fusion 
  # rate2 <- autosomal fission 
  # rate3 <- sex chromosome autosome fusions
  # rate4 <- Y chromosome loss
  # rate5 <- Y chromosome loss in complex systems
  
  # this is the parameter that defines the number of parts in the qmatrix. if
  # we have complex sex chromosome systems then the q matrix have 3 parts. 
  # a part for XO / ZO
  # a part for XY / ZW
  # a part for complex XY / complex ZW
  if(complex == T){
    devider <- 3
  }else{
    devider <- 2
  }
  # from here onwards we fill the q matrix
  # autosomal fissions and fusions in XO systems
  for(i in 1){
    while (i < (nrow(qmat)/devider)) {
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
      i <- i + 1
    }
  }
  # autosomal fissions and fusions in XY systems
  for(i in ((nrow(qmat)/devider)+1):((nrow(qmat)/devider)*2)){
    if(i < ((nrow(qmat)/devider)*2)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # autosomal fissions and fusions in complex systems
  for(i in (((nrow(qmat)/devider)*2)+1):nrow(qmat)){
    if(i < nrow(qmat)){
      qmat[i, i+1] <- 1 # fissions
      qmat[i+1, i] <- 2 # fusions
    }
  }
  # sex chromosome autosome fusion from XO to XY
  for(i in 2:(nrow(qmat)/devider)){
    qmat[i,(nrow(qmat)/devider) + i -1] <- 3
  }
  # sex chromosome autosome fusion from XY to multi-XY
  if(complex == T){
    for(i in ((nrow(qmat)/devider)+1):((nrow(qmat)/devider)*2)){
      qmat[i,(((nrow(qmat)/devider)) + i-1)] <- 3
    }
  }
  # Y chromosome loss in XY
  for(i in ((nrow(qmat)/devider)+1):((nrow(qmat)/devider)*2)){
    qmat[i, i-(nrow(qmat)/devider)] <- 4
  }
  # Y chromosome loss in multi-xy
  if(complex == T){
    for(i in (((nrow(qmat)/devider)*2)+1):(nrow(qmat))){
      qmat[i, i-(nrow(qmat)/devider)] <- 5
    }
  }
  # this defines number of columns each sex chromosome systems have in the
  # qmat
  hp <- ncol(qmat)/devider + 1
  # this will remove unwanted columnts making the qmat +- 1 the given chrom
  # range
  if(complex == T){
    drops <- c(1:(minChroms-2), hp:(hp + minChroms-3),(hp*2 - 1):((hp*2 - 1) + minChroms-3))  
  }else{
    drops <- c(1:(minChroms-2), hp:(hp + minChroms-3))
  }
  qmat <- qmat[-drops, -drops]
  return(qmat)
}

# this function will generate a probability matrix given chromosome number and
# sex chromosome systems
get.Pmatrix <- function(dat, complex = F, sex.system = "XY"){
  
  # make it so that the letter case does not matter when selecting columns and
  # for sex.system
  colnames(dat) <- tolower(colnames(dat))
  sex.system <- toupper(sex.system)
  # define total number of chromosome states
  chrom.range <- range(dat$chroms)
  maxChroms <- chrom.range[2] + 1
  minChroms <- chrom.range[1] - 1
  if(complex == T){
    nChroms <- (maxChroms -1) * 3
  }else{
    nChroms <- (maxChroms -1) * 2
  }
  # make the initial q matrix
  qmat <- matrix(data = 0,
                 nrow = nChroms,
                 ncol = nChroms)
  # define rownames and colnames
  # although row names and col names starts from 1, chromsome counts starts from
  # 2
  if(complex == T){
    if(sex.system == "XY"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                              paste(2:maxChroms,"XY",sep=""),
                                              paste(2:maxChroms,"XXY/XYY",sep=""))
    }
    if(sex.system == "ZW"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                              paste(2:maxChroms,"ZW",sep=""),
                                              paste(2:maxChroms,"ZZW/ZWW",sep=""))
    }
  }else{
    if(sex.system == "XY"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"XO",sep=""),
                                              paste(2:maxChroms,"XY",sep=""))
    }
    if(sex.system == "ZW"){
      row.names(qmat)  <- colnames(qmat) <- c(paste(2:maxChroms,"ZO",sep=""),
                                              paste(2:maxChroms,"ZW",sep=""))
    }
  }
  # this is the parameter that defines the number of parts in the qmatrix. if
  # we have complex sex chromosome systems then the q matrix have 3 parts. 
  # a part for XO / ZO
  # a part for XY / ZW
  # a part for complex XY / complex ZW
  if(complex == T){
    devider <- 3
  }else{
    devider <- 2
  }
  # we start making the p matrix
  pmat <- qmat[NULL,]
  # this defines number of columns each sex chromosome systems have in the
  # qmat
  hp <- ncol(qmat)/devider + 1
  # add rows 
  # make a new matrix. number of rows in thes new matrix is equal to the
  # number of species in our data table. 
  temp.mat <- matrix(data = 0,
                     nrow = nrow(dat),
                     ncol = ncol(pmat))
  # name the row names according to the species names in our data table
  rownames(temp.mat) <- dat$specisname
  # combine the new matrix with the pmat
  pmat <- rbind(pmat, temp.mat)
  # fill out the pmat when the sex limited chromosome is Y
  if(sex.system == "XY"){
    for(i in 1:nrow(dat)){
      if(dat$scs[i] == "XO"){
        pmat[i, dat$chroms[i] - 1] <- 1
      }
      if(dat$scs[i] == "XY"){
        pmat[i, (hp + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] %in% c("XYY","XXY")){
        pmat[i, ((hp*2) + dat$chroms[i] - 3)] <- 1
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
        pmat[i, (hp + dat$chroms[i] - 2)] <- 1
      }
      if(dat$scs[i] %in% c("ZWW","ZZW")){
        pmat[i, ((hp*2) + dat$chroms[i] - 3)] <- 1
      }
    }
  }
  # this will remove unwanted columnts making the qmat +- 1 the given chrom
  # range
  if(complex == T){
    drops <- c(1:(minChroms-2), hp:(hp + minChroms-3),(hp*2 - 1):((hp*2 - 1) + minChroms-3))  
  }else{
    drops <- c(1:(minChroms-2), hp:(hp + minChroms-3))
  }
  pmat <- pmat[, -drops]
  
  return(pmat)
}

# this function will look for species and genus level matches given a datafram 
# with species names and a set of trees
get.SpeciesMatches <- function(dat, trees){
  dat.new <- dat
  colnames(dat.new) <- tolower(colnames(dat.new))
  tree <- trees[[1]]
  # get the species level matches
  dat.new.s.m <- dat.new[dat.new$specisname %in%  tree$tip.label,]
  
  # from these species level matches if there are multiple hits for a sinlge
  # species get only 1 hit
  unique.sp.matches <- unique(dat.new.s.m$specisname)
  dat.new.s.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.s.m$specisname)),
                                             ncol = ncol(dat.new)))
  colnames(dat.new.s.m.unique) <- colnames(dat.new)
  
  for(i in 1:length(unique.sp.matches)){
    dat.new.temp <- dat.new.s.m[dat.new.s.m$specisname %in% unique.sp.matches[i],]
    dat.new.s.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
  }
  dat.new.unhit <- dat.new[!(dat.new$specisname %in% tree$tip.label),]
  
  hit.genera <- c(unique(dat.new.s.m$Genus))
  unhit.gen <- unique(dat.new$Genus[!(dat.new$Genus %in% hit.genera)])
  dat.new.unhit <- dat.new.unhit[dat.new.unhit$Genus %in% unhit.gen,]
  tree.unhit <- tree$tip.label[!(tree$tip.label %in% unique.sp.matches)]
  tree.unhit.gen <- c()
  
  for(i in 1:length(tree.unhit)){
    tree.unhit.gen[i] <- unlist(strsplit(tree.unhit[i], split = "_"))[[1]]
  }
  
  tree.hit.gen.num <- which(tree.unhit.gen %in% hit.genera)
  tree.true.unhit.gen <- tree.unhit.gen[!(tree.unhit.gen %in% hit.genera)]
  
  # get genera level matches
  dat.new.g.m <- dat.new.unhit[dat.new.unhit$Genus %in%  tree.true.unhit.gen,]
  # sum(tolower(unique(tree.true.unhit.gen)) %in% unique(dat.new.g.m$Genus))
  
  # from these species level matches if there are multiple hits for a sinlge
  # species get only 1 hit
  unique.gn.matches <- unique(dat.new.g.m$Genus)
  dat.new.g.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.g.m$Genus)),
                                             ncol = ncol(dat.new)))
  colnames(dat.new.g.m.unique) <- colnames(dat.new)
  
  for(i in 1:length(unique.gn.matches)){
    dat.new.temp <- dat.new.g.m[dat.new.g.m$Genus %in% unique.gn.matches[i],]
    dat.new.g.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
    dat.new.g.m.unique$specisname[i] <- paste(dat.new.g.m.unique$Genus[i], "_sp.", sep = "")
  }
  
  tree..unhit.hit.gen.num <- which(tree.unhit.gen %in% dat.new.g.m.unique$Genus)
  tree..unhit.hit.sp.names <- tree.unhit[tree..unhit.hit.gen.num]
  unhit.temp.dat.new <- as.data.frame(matrix(data = NA,
                                             nrow = length(tree..unhit.hit.sp.names),
                                             ncol = 2))
  
  colnames(unhit.temp.dat.new) <- c("gen", "sp")
  
  unhit.temp.dat.new$sp <- tree..unhit.hit.sp.names
  
  for(i in 1:length(tree..unhit.hit.sp.names)){
    unhit.temp.dat.new$gen[i] <- unlist(strsplit(tree..unhit.hit.sp.names[i], split = "_"))[[1]]
  }
  
  tree.gen.matches <- unhit.temp.dat.new[!(duplicated(unhit.temp.dat.new$gen)),]
  tree.gen.matches$gen <- paste(tree.gen.matches$gen, "_sp.", sep = "")
  
  # for tree name correction
  tree.names <- as.data.frame(matrix(data = NA,
                                     nrow = (nrow(dat.new.s.m.unique) + nrow(dat.new.g.m.unique)),
                                     ncol = 2))
  colnames(tree.names) <- c("species_name", "name_on_tree")
  
  tree.names$species_name[1:nrow(dat.new.s.m.unique)] <- tree.names$name_on_tree[1:nrow(dat.new.s.m.unique)] <- dat.new.s.m.unique$specisname
  tree.names$species_name[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$gen
  tree.names$name_on_tree[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$sp
  
  finaldat.new <- rbind(dat.new.s.m.unique, dat.new.g.m.unique)
  results <- list(finaldat.new, tree.names)
  
  names(results) <- c("chroms", "name_corrections")
  
  
  return(results)
}