#create the figure.1 from the paper

gamma <- seq(0,1, by = 0.001)
Fst.3 <- matrix(0, nrow = length(gamma),
                ncol = length(allele.frequencies.no.extra))
rownames(Fst.3) <- gamma
colnames(Fst.3) <- names(allele.frequencies.no.extra)
for(allele in names(allele.frequencies.no.extra))
{
  print(allele)
  
  ##get the allele frequencies
  p2 <- allele.frequencies.no.extra[[allele]]["EUROPE",]
  p3 <- allele.frequencies.no.extra[[allele]]["AMERICA",]
  p1 <- allele.frequencies.no.extra[[allele]]["LATINO",]
  
  for(a in gamma)
  {
    ##first get admixture of populations 2 and 3
    admix.2.and.3 <- a*p2+(1-a)*p3
    ##now get Fst between this admixed population and population 1
    Fst.3[as.character(a), allele] <- get.Fst(p1, admix.2.and.3)
  }
}
set.seed(81048)
sample.loci <- sample(names(allele.frequencies.no.extra), 20)
sample.loci[1:5]
ymin <- min(Fst.3[,sample.loci])
ymax <- max(Fst.3[,sample.loci])

par(cex.axis = 1.2, cex.lab = 1.2)
plot(log10(Fst.3[,sample.loci[1]]) ~ gamma,
     ylim = c(log10(ymin), log10(ymax)), type = "l",
     xlab = expression(paste("Admixture coefficient (", gamma, ")")),
     ylab = expression(paste(log[10],"[",G(gamma),"]")))
##get maximum
mm <- which.max(Fst.3[,sample.loci[1]])
points(log10(Fst.3[mm,sample.loci[1]])~gamma[mm], type = "p",
       pch = 20, col="black", lwd = 5)    
for(allele in sample.loci[2:20])
{
  points(log10(Fst.3[,allele])~gamma, type = "l")
  ##get maximum
  mm <- which.max(Fst.3[,allele])
  points(log10(Fst.3[mm,allele])~gamma[mm], type = "p",
         pch = 20, col="black", lwd = 5)    
}

#create the figure.3 from the paper.

gamma <- seq(0,1, by = 0.001)
Fst <- matrix(0, nrow = length(gamma),
              ncol = length(allele.frequencies.no.extra))
rownames(Fst) <- gamma
colnames(Fst) <- names(allele.frequencies.no.extra)
for(allele in names(allele.frequencies.no.extra))
{
  print(allele)
  
  ##get the allele frequencies
  p1 <- allele.frequencies.no.extra[[allele]]["EUROPE",]
  p2 <- allele.frequencies.no.extra[[allele]]["AMERICA",]
  
  for(a in gamma)
  {
    ##first get admixture of populations 1 and 2
    admix.1.and.2 <- a*p1+(1-a)*p2
    ##now get Fst between this admixed population and population 1
    Fst[as.character(a), allele] <- get.Fst(p1, admix.1.and.2)
  }
}
set.seed(81048)
locus1 <- which(apply(Fst, 2, which.max)==1)[1]
sample.loci <- sample(names(allele.frequencies.no.extra), 20)

ymin <- min(Fst[,sample.loci])
ymax <- max(Fst[,sample.loci])

par(cex.axis = 1.2, cex.lab = 1.2)
plot(Fst[,sample.loci[1]] ~ gamma,
     ylim = c(0, 0.125), type = "l",
     xlab = expression(paste("Admixture coefficient (", gamma, ")")),
     ylab = expression(F[st]))
for(allele in sample.loci[1:20])
{
  points(Fst[,allele]~gamma, type = "l")
}