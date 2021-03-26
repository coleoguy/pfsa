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
  ## here is a description of the different parameters that will be used in the
  ## qmatrix
  
  # r1 = 1,    # AA fusion  XO
  # r2 = 2,    # AA fission XO
  # r3 = 3,    # AA fusion  XY
  # r4 = 4,    # AA fission XY
  # r5 = 5,    # AA fusion  Neo.XY
  # r6 = 6,    # AA fission Neo.XY 
  # r7 = 7,    # AA fusion  XXY
  # r8 = 8,    # AA fission XXY
  # r9 = 9,    # AA fusion  XYY
  # r10 = 10,  # AA fission XYY
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
  rate.table[rate.table[, 3] == 11, 3] <- "parq1"
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
                                " ~ lambda1", sep = ""))
      mu <- c(mu, paste("mu", colnames(parMat)[i], 
                        " ~ mu1", sep = ""))
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
              "lambda1",
              "mu1")
  lik.con <- constrain(lik, formulae = c(formulae, lambda, mu), extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(qmat)
    return(lik.con)
}
