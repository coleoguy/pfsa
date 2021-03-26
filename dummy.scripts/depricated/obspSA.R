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