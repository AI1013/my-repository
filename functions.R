#Calculate the heterozygosity.

get.het <- function(allele.freq) {
  
  if (any(is.na(allele.freq) == TRUE)) warning("There should not be NA in allele.freq")
  1 - sum(allele.freq^2)
}

#Calculate the heterozygosity of the admixed population

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

#Calculate Fst

get.Fst <- function(pop.1, pop.2) {
  if (any(is.na(pop.1) == TRUE)) warning("There should not be NA in the pop.1")
  if (any(is.na(pop.2) == TRUE)) warning("There should not be NA in the pop.2")
  meta.pop <- (pop.1 + pop.2)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.2)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

#Calculate Fst: between a founder population and the admixed population.

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

#Calculate Fst for pop.1 and admixed population.

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