constrainQmat <-  function (qmat, lik, verbose = F) {
  if (ncol(qmat) < 100) 
    pad <- 2
  if (ncol(qmat) >= 100) 
    pad <- 3
  if (ncol(qmat) < 10) 
    pad <- 1
  parMat <- qmat  
  colnames(parMat) <- sprintf(paste("%0", pad, "d", sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat),3))
  rate.table[, 1] <- rep(as.character(row.names(parMat)), each = ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[2],]
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
  rate.table[rate.table[, 3] == 1, 3] <- "r01"
  rate.table[rate.table[, 3] == 2, 3] <- "r02"
  rate.table[rate.table[, 3] == 3, 3] <- "r03"
  rate.table[rate.table[, 3] == 4, 3] <- "r04"
  rate.table[rate.table[, 3] == 5, 3] <- "r05"
  rate.table[rate.table[, 3] == 6, 3] <- "r06"
  rate.table[rate.table[, 3] == 7, 3] <- "r07"
  rate.table[rate.table[, 3] == 8, 3] <- "r08"
  rate.table[rate.table[, 3] == 9, 3] <- "r09"
  rate.table[rate.table[, 3] == 10, 3] <- "r10"
  rate.table[rate.table[, 3] == 11, 3] <- "rq1"
  rate.table[rate.table[, 3] == 12, 3] <- "r12"
  rate.table[rate.table[, 3] == 13, 3] <- "r13"
  rate.table[rate.table[, 3] == 14, 3] <- "r14"
  rate.table[rate.table[, 3] == 15, 3] <- "r15"
  rate.table[rate.table[, 3] == 16, 3] <- "r16"
  rate.table[rate.table[, 3] == 17, 3] <- "r17"
  rate.table[rate.table[, 3] == 18, 3] <- "r18"
  rate.table[rate.table[, 3] == 19, 3] <- "r19"
  rate.table[rate.table[, 3] == 20, 3] <- "r20"
  rate.table[rate.table[, 3] == 21, 3] <- "r21"
  rate.table[rate.table[, 3] == 22, 3] <- "r22"
  rate.table[rate.table[, 3] == 23, 3] <- "r23"
  rate.table[rate.table[, 3] == 24, 3] <- "r24"
  rate.table[rate.table[, 3] == 25, 3] <- "r25"
  rate.table[rate.table[, 3] == 26, 3] <- "r26"
  rate.table[rate.table[, 3] == 27, 3] <- "r27"
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
  extras <- c("r01",
              "r02",
              "r03",
              "r04",
              "r05",
              "r06",
              "r07",
              "r08",
              "r09",
              "r10",
              "r11",
              "r12",
              "r13",
              "r14",
              "r15",
              "r16",
              "r17",
              "r18",
              "r19",
              "r20",
              "r21",
              "r22",
              "r23",
              "r24",
              "r25",
              "r26",
              "r27",
              "lambda1",
              "mu1")
  lik.con <- constrain(lik, formulae = c(formulae, lambda, mu), extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(qmat)
  if (verbose == T) 
    return(list(lik.con, parMat))
  if (verbose == F) 
    return(lik.con)
}
