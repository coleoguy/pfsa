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

maxnum <- 60
XO <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XO")
XY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XY")
XYY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XYY")
XXY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XXY")
rates <- c(XO,XY,XYY,XXY)
types <- rep(c("XO", "XY", "XYY","XXY"), each=length(XO))
autosomes <- rep(seq(from=2, to=maxnum, by=2), times=4)
res <- data.frame(rates, types, autosomes)

ggplot(res, aes(y=rates, x=autosomes)) + geom_point(aes(colour=types), stat="identity", position="identity", alpha=0.5, size=3) + geom_line(aes(colour=types), stat="identity", position="identity", alpha=0.5) + theme_bw() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + xlab("Diploid autosome count") + ylab("Proportion of fusions joining autosome and gonosome")
