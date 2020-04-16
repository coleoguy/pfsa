Pfsa <- function(Da, scs){
  if(scs=="XO"){
    Xs <- 1
    Y <- 0
    Ds <- Da + 1
    Dd <- Da + 2
  }
  if(scs=="XY"){
    Xs <- 1
    Y <- 1
    Ds <- Da + 2
    Dd <- Da + 2
  }
  if(scs=="XYY"){
    Xs <- 1
    Y <- 2
    Ds <- Da + 3
    Dd <- Da + 2
  }
  if(scs=="XXY"){
    Xs <- 2
    Y <- 1
    Ds <- Da + 3
    Dd <- Da + 4
  }
  res <- 1 - ((Da*(Da-2)+2*Xs*(2*Xs-2))/(2*Dd*(Dd-2))) - 
    ((Da*(Da-2)+max(c(Xs,Y))*(max(c(Xs,Y))-1))/(2*Ds*(Ds-2)))
  return(res)
}