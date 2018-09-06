load("D:/sdxjw/Documents/R working directory/allele_matrices_no_extra_allele.RData")

plot_HorF <- function (x,# x should be a list including the allele frequencies of three populations(EUROPE, AMERICA, AFRICA)
                       y, # y should be either "Fst" or "Hadm"
                       breaks, # breaks are numeric vectors to categorize  Fst or Hamd.
                       label # the label of the legend.
){
  require(Fst)
  require(RColorBrewer)
  require(ggplot2)
  require(ggtern)
  
  if (!(y %in% c("Fst", "Hadm"))) stop("spacing must be either 'Fst' or 'Hadm'")
  alpha1 <- alpha2 <- seq(0,100,by=1)
  alpha3 = 100-alpha1-alpha2
  d <- expand.grid(alpha1, alpha2)
  
  row_sums <- rowSums(d)
  var3 <- 100-row_sums
  
  var3[var3<0] <- NA
  d$var3 <- var3
  d <- na.omit(d)
  
  gamma <- d/100
  
  sum_rows <- rowSums(gamma)
  
  
  colnames(gamma) <- paste("gamma", 1:3)
  
  allele.frequencies.no.extra$x <-  as.matrix(allele.frequencies.no.extra[x])[[1]][1:3,]
  position <- allele.frequencies.no.extra["x"]
  
  for(locus in names(position)){
    print(locus)
    
    ## get the allele frequencies
    
    p1 <- position[[locus]]["EUROPE",]
    p2 <- position[[locus]]["AMERICA",]
    p3 <- position[[locus]]["AFRICA",]
    
    print(p1)
    print(p2)
    print(p3)
  }
  
  ## get the admixed population
  gamma_1 <- gamma$`gamma 1`
  gamma_2 <- gamma$`gamma 2`
  gamma_3 <- gamma$`gamma 3`
  
  pop1 <- list() ## get p1*gamma1
  t <- 1
  for(a in gamma_1){
    admix.pop1 <- a*p1
    pop1[[t]] <- c(admix.pop1)
    t = t+1
  }
  
  
  pop2 <- list() ## get p2*gamma2
  t <- 1
  for(i in gamma_2){
    admix.pop2 <- i*p2
    pop2[[t]] <- c(admix.pop2)
    t = t+1
  }
  
  
  pop3 <- list() ## get p3*gamma3
  t <- 1
  for(n in gamma_3){
    admix.pop3 <- n*p3
    pop3[[t]] <- c(admix.pop3)
    t = t+1
  }
  
  
  pop1 <- t(as.data.frame(pop1))
  pop2 <- t(as.data.frame(pop2))
  pop3 <- t(as.data.frame(pop3))
  
  admix.pop <- pop1 + pop2 + pop3 ## get admix.pop
  row_sums <- data.frame(rowSums(admix.pop)) ## check if p1*gamma1 + p2*gamma2 + p3*gamma3 equals to 1
  
  
  ## get the hadm or Fst of the admixed population and create thire plots.
  if(y == "Hadm"){
    newlist <- list()
    t <- 1
    for(i in 1:nrow(admix.pop)){
      a <- admix.pop[i, ]
      H.1 <- get.het(a)
      newlist[[t]] <- c(H.1)
      t = t+1
    }
    
    hadm <- gamma$hadm <- newlist
    colnames(gamma)[1:3] <- c("gamma.1", "gamma.2", "gamma.3")
    
    gamma$hadm <- unlist(gamma$hadm)
    
    gamma$Hadm = cut(gamma$hadm, breaks)
    levels(gamma$Hadm) <- label
    
    ## create the plot
    
    
    mycol <- brewer.pal(9, 'YlOrRd')
    names(mycol) <- levels(gamma$Hadm)
    colScale <- scale_colour_manual(name = expression(H[adm]),values = mycol)
    Max <- gamma[which(gamma$hadm == max(gamma$hadm)), ]
    Max <- matrix(Max)
    g1 <- Max[1,1]
    g2 <- Max[2,1]
    g3 <- Max[3,1] 
    g1 <- unlist(g1)
    g2 <- unlist(g2)
    g3 <- unlist(g3)
    
    pl <- ggtern(gamma, aes(gamma.2, gamma.1, gamma.3, value=hadm)) + geom_point(aes(color=Hadm)) + colScale 
    pl1 <- pl + theme_hideticks() + theme_hidelabels() + geom_mask() + annotate("point", x=g2, y=g1, z=g3, size=0.5) + 
      annotate(geom = "text", 
               x = c(0.5, 0.41, 0.0),
               y = c(0.0, 0.5, 0.59),
               z = c(0.5, 0.0, 0.5),
               angle = c(0, 0, 0),
               vjust = c(2.0, 0.5, 0.5),
               hjust = c(0.5, 1.3, -0.3),
               label = c("(0, 0.5, 0.5)", "(0.5, 0.5, 0)", "(0.5, 0 ,0.5)"),
               color = c("black", "black", "black")) 
    pl2 <- pl1 + theme(legend.title = element_text(colour="black", size=10,  face="bold")) + guides(col = guide_legend(reverse = TRUE)) + 
      annotate(geom = "text",
               x = c(0, 1, 0),
               y = c(1, 0, 0),
               z = c(0, 0, 1),
               label = c("(1, 0, 0)", "(0, 1, 0)", "(0, 0, 1)"),
               vjust = c(-1, 3, 3))
    temp1 <- expression("Native American"~gamma[2])
    temp2 <- expression("European"~gamma[1])
    temp3 <- expression("African"~gamma[3])
    temp4 <- "H[adm]"
    
    pl3 <- pl2 + theme_hidetitles() + annotate(geom = "text", 
                                               x = c(1, 0, 0), 
                                               y = c(0, 1, 0),
                                               z = c(0, 0, 1),
                                               label = c(as.character(temp1), as.character(temp2), as.character(temp3)),
                                               vjust = c(1, -2, 1),
                                               parse = TRUE)
    pl4 <- pl3 + annotate("text", x = 0, y = 0.75, z = 0.3, label = temp4, hjust = 4.0, vjust = -1.5, fontface = "bold", parse=TRUE) + 
      annotate("segment", x=0.5, xend=0.15, y=0.5, yend=0.15, z=0, zend=0.05) + annotate("segment", x=0, xend=0.05, y=0.5, yend=0.15, z=0.5, zend=0.15) +
      annotate("segment", x=0.5, xend=0.15, y=0, yend=0.05, z=0.5, zend=0.15) 
    
    
    pl4 
  } else {
    newlist <- list()
    t <- 1
    for(i in 1:nrow(admix.pop)){
      a <- admix.pop[i, ]
      Fst.1 <- get.Fst(p1, a)
      newlist[[t]] <- c(Fst.1)
      t = t+1
    }
    
    F_st <- gamma$F_st <- newlist
    colnames(gamma)[1:3] <- c("gamma.1", "gamma.2", "gamma.3")
    
    gamma$F_st <- unlist(gamma$F_st)
    
    gamma$Fst = cut(gamma$F_st, breaks)
    levels(gamma$Fst) <- label
    
    ## create the plot
    
    mycol <- brewer.pal(9, 'YlOrRd')
    names(mycol) <- levels(gamma$Fst)
    colScale <- scale_colour_manual(name = expression(F[st]),values = mycol)
    Max <- gamma[which(gamma$F_st == max(gamma$F_st)), ]
    Max <- matrix(Max)
    g1 <- Max[1,1]
    g2 <- Max[2,1]
    g3 <- Max[3,1] 
    g1 <- unlist(g1)
    g2 <- unlist(g2)
    g3 <- unlist(g3)
    
    p <- ggtern(gamma, aes(gamma.2, gamma.1, gamma.3, value=F_st)) + geom_point(aes(color=Fst)) + colScale 
    p1 <- p + theme_hideticks() + theme_hidelabels() + geom_mask() + annotate("point", x = g2, y = g1, z = g3, size = 0.5) + 
      annotate(geom = "text", 
               x = c(0.5, 0.41, 0.0),
               y = c(0.0, 0.5, 0.59),
               z = c(0.5, 0.0, 0.5),
               angle = c(0, 0, 0),
               vjust = c(2.0, 0.5, 0.5),
               hjust = c(0.5, 1.3, -0.3),
               label = c("(0, 0.5, 0.5)", "(0.5, 0.5, 0)", "(0.5, 0 ,0.5)"),
               color = c("black", "black", "black")) 
    p2 <- p1 + theme(legend.title = element_text(colour="black", size=10,  face="bold")) + guides(col = guide_legend(reverse = TRUE)) + 
      annotate(geom = "text",
               x = c(0, 1, 0),
               y = c(1, 0, 0),
               z = c(0, 0, 1),
               label = c("(1, 0, 0)", "(0, 1, 0)", "(0, 0, 1)"),
               vjust = c(-1, 3, 3))
    temp1 <- expression("Native American"~gamma[2])
    temp2 <- expression("European"~gamma[1])
    temp3 <- expression("African"~gamma[3])
    p3 <- p2 + theme_hidetitles() + annotate(geom = "text", 
                                             x = c(1, 0, 0), 
                                             y = c(0, 1, 0),
                                             z = c(0, 0, 1),
                                             label = c(as.character(temp1), as.character(temp2), as.character(temp3)),
                                             vjust = c(1, -2, 1),
                                             parse = TRUE)
    
    p3
  }
}





plot_HorF("D2S1399", y="Hadm", breaks = c(0, 0.8371, 0.8504, 0.8578, 0.8620, 0.8661, 0.8702, 0.8741, 0.8781, 1),
          label = c(0, 0.8371, 0.8504, 0.8578, 0.8620, 0.8661, 0.8702, 0.8741, 0.8781))

plot_HorF("GATA101G01_12", y="Fst", breaks = c(-1, 0.00018, 0.00086, 0.00242, 0.00617, 0.01318, 0.02512, 0.04713, 0.13772, 1),
          label = c(0, 0.00018, 0.00086, 0.00242, 0.00617, 0.01318, 0.02512, 0.04713, 0.13772))









