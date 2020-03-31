# Terrence Sylvester
# 26th March 2020
# pradakshanas@gmail.com
# make transision matrix

# make a matrix with 100 chromosomes on each state. Here we will be using 
# sex chromosome systems as the states. These states are XO and XY.

# we will go from 2:100 chromosomes. So we will have a total of 198 rows and cols

qmatGen <- function(maxChroms = NULL){

# define total number of chromosome states
#maxChroms <- 15
nChroms <- (maxChroms - 1) * 2

# make the initial q matrix
qmat <- matrix(data = 0,
               nrow = nChroms,
               ncol = nChroms) 

# define rownames and colnames
# although row names and col names starts from 1, chromsome counts starts from
# 2
row.names(qmat)  <- colnames(qmat) <- 1:nChroms 

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
  if(i + 1 < nrow(qmat)){
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

# # put rate parameter for y chromosome capture - rate15
# for(i in 1:(nrow(qmat)/2)){
#   qmat[i, (nrow(qmat)/2) + i] <- 15
# }
return(qmat)} 

qmatGen(maxChroms = 4)
