# load libraries
library(ggplot2)
library(ape)
library(stringr)
# read in data
dat <- read.csv("mammal_chroms.csv", as.is=T)
tree <- read.tree("mammals.newick")
# remove species with no chromosome number data
dat <- dat[dat$female2n != "",]
# process data to get a single chromosome number when there is a
# range of chromosome number given
for(i in 1:nrow(dat)){
  dat$female2n[i] <- as.numeric(strsplit(dat$female2n[i], split = "-")[[1]][1])
}
# get only the data that we need
dat <- dat[,c(1,5)]
dat$female2n <- as.numeric(dat$female2n)
# get the haploid chromosome number
dat$femalen <- (dat$female2n) / 2
# round the decimal values
dat$femalen <- round(dat$femalen)
# plot
## here mammal families are plotted in the alphabetical order
ggplot(dat, aes(x= femalen, y=Order)) + 
  geom_tile(stat="bin2d", position="identity", alpha=1, bins = 50) + 
  theme_grey() + 
  theme(text=element_text(family="sans",
                          face="plain", 
                          color="#000000",
                          size=7, 
                          hjust=0.5,
                          vjust=0.5)) + 
  xlab("haploid chromosome number") + ylab("Order")+
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  scale_x_continuous(breaks = seq(0,60, by=1))

### families in phylogenetic order ###
## process the tree to get the phylogenetic order of Mammal families
order <- tree$tip.label
for(i in 1:length(order)){
  order[i] <- strsplit(order[i], split = "_")[[1]][4]
}
# remove the last entry because that entry does not have family information
order <- order[-4099]
# remove duplicated entries
order <- order[!duplicated(order)]
# make them sentence case
order <- str_to_sentence(order, locale = "en")
# this will make the order of the levels to be in phylogenetic order
dat$labs <- factor(dat$Order,
                   levels = order)
# in this phylogeny few families are missing. So here I am removing those missing
# families (7 families to be exact and a total of 14 taxa)
dat <- dat[!is.na(dat$labs),]
#plot
ggplot(dat, aes(x= femalen, y=labs)) + 
  geom_tile(stat="bin2d", position="identity", alpha=1, bins = 75) + 
  theme_grey() + 
  theme(text=element_text(family="sans",
                          face="plain", 
                          color="#000000",
                          size=10, 
                          hjust=0.5,
                          vjust=0.5)) + 
  xlab("haploid chromosome number") + ylab("Order")+
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  scale_x_continuous(breaks = seq(0,60, by=2.5))

