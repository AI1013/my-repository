plot_HorF <- function (x,# x should be a list including the allele frequencies of three populations(EUROPE, AMERICA, AFRICA)
                       y, # y should be either "Fst" or "Hadm"
                       breaks = default_break, # breaks are numeric vectors to categorize  Fst or Hamd.
                       label = default_label # the label of the legend.
){
  
  
  
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
  
  allele.frequencies$x <-  as.matrix(allele.frequencies[[x]])[1:3,]
  position <- allele.frequencies["x"]
  
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
    
    the_break <- quantile(as.numeric(hadm), probs = seq(0,1,0.125))
    the_break[9] <- 1
    the_break1 <- c(-1, the_break)
    default_break <- as.numeric(as.matrix(the_break1))
    
    
    the_break2 <- as.list(the_break1)
    the_break2[[10]] <- NULL
    the_break2[[1]] <- 0
    default_label <- as.numeric(the_break2)
    
    gamma$Hadm = cut(gamma$hadm, breaks)
    levels(gamma$Hadm) <- label
    
    ## create the plot
    
    
    mycol <- RColorBrewer::brewer.pal(9, 'YlOrRd')
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
    
    pl <- ggtern::ggtern(gamma, aes(gamma.2, gamma.1, gamma.3, value=hadm)) + geom_point(aes(color=Hadm)) + colScale 
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
    
    the_break <- quantile(as.numeric(F_st), probs = seq(0,1,0.125))
    the_break[9] <- 1
    the_break1 <- c(-1, the_break)
    default_break <- as.numeric(as.matrix(the_break1))
    
    
    the_break2 <- as.list(the_break1)
    the_break2[[10]] <- NULL
    the_break2[[1]] <- 0
    default_label <- as.numeric(the_break2)
    
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

get.het <- function(allele.freq) {
  
  if (any(is.na(allele.freq) == TRUE)) warning("There should not be NA in allele.freq")
  1 - sum(allele.freq^2)
}

get.Hadm.admix.general <- function(pops, gamma){
  if (!(rowSums(pops) == 1)) stop("The sum of the allele frequency for each population should add up to 1")
  if (!(sum(gamma) == 1)) stop("The sum of gamma should add up to 1")
  if (any(is.na(pops) == TRUE)) warning("There should not be NA in pops")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  ##get admixed population
  admix.pop <- t(pops) %*% gamma
  ##calculate the heterozygosity of the admixed population
  Hadm <- get.het(admix.pop)
  Hadm
}


get.Fst <- function(pop.1, pop.2) {
  if (any(is.na(pop.1) == TRUE)) warning("There should not be NA in the pop.1")
  if (any(is.na(pop.2) == TRUE)) warning("There should not be NA in the pop.2")
  meta.pop <- (pop.1 + pop.2)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.2)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

get.Fst.admix <- function(pop.1, pop.2, gamma) {
  if (any(is.na(pop.1) == TRUE)) warning("There should not be NA in pop.1")
  if (any(is.na(pop.2) == TRUE)) warning("There should not be NA in pop.2")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  pop.3 <- gamma * pop.1 + (1 - gamma) * pop.2
  meta.pop <- (pop.1 + pop.3)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.3)))
  Fst <- (Ht - Hs)/Ht
  Fst
}




get.Fst.admix.general <- function(pops, gamma) {
  if (!(rowSums(pops) == 1)) stop("The sum of the allele frequency for each population should add up to 1")
  if (!(sum(gamma) == 1)) stop("The sum of gamma should add up to 1")
  if (any(is.na(pops) == TRUE)) warning("There should not be NA in pops")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  ## get admixed population
  admix.pop <- t(pops) %*% gamma
  
  ## get allele frequencies of 'meta-population'
  meta.pop <- (pops[1, ] + admix.pop)/2
  
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pops[1, ]), get.het(admix.pop)))
  Fst <- (Ht - Hs)/Ht
  Fst
}