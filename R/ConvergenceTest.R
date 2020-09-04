###Convergence tests for GMM analyses###
## The functions implemented here are based on the equations and analyses described by Stayton (2015) ##

## Sub-functions ##
# Calculate the change in trait values across the phylogeny #
calcchangephy = function(phy,allvals){
  #phy = phylogeny of taxa in sample
  #allvals = trait dataset for all nodes and tips (result of multianc)
  pdims <- dim(phy$edge)
  changes <- rep(0, pdims[1])
  j <- 1
  for (j in 1:pdims[1]) {
    ancnode <- phy$edge[j, 1]
    desnode <- phy$edge[j, 2]
    ancvals <- allvals[ancnode, ]
    desvals <- allvals[desnode, ]
    change <- (sum((ancvals - desvals)^2))^0.5
    changes[j] <- change
  }
  answer <- changes
}
#
# Calculation for C metrics for a pair of taxa #
calc.conv.rates = function(distmat,phy,tax1,tax2,allphendata){
  #distmat = distance matrix of all tips and all nodes
  #phy = phylogeny of taxa in sample
  #tax1 = name of first taxon of interest
  #tax2 = name of second taxon of interest
  #allphendata = trait dataset for all nodes and tips (result of multianc)
  mrca.t1t2 <- getMRCA(phy,c(tax1,tax2))
  tip.nums = labelstonumbers(phy,c(tax1,tax2))
  anc.t1t2 = unlist(Ancestors(phy,tip.nums,type="all"))
  anc.t1t2 = anc.t1t2[anc.t1t2 >= mrca.t1t2]

  t.combo = expand.grid(c(tax1,tax2),anc.t1t2)
  t.combo = unite(data=t.combo,col="pair",Var1:Var2,sep="_",remove=TRUE)

  distmat.paired = unite(data=distmat,col="pair",n1:n2,sep="_",remove=TRUE)
  tip.dist = subset(distmat.paired, pair == paste(tax1,sep="_",tax2))
  node.dists = rbind(subset(distmat.paired, pair %in% t.combo$pair),tip.dist)

  tdist = tip.dist$value
  mxdist = max(node.dists$value)

  C1 = 1-(tdist/mxdist)
  C2 = mxdist - tdist
  wholephylchanges <- sum(calcchangephy(phy, allphendata))
  C3 <- C2/wholephylchanges

  subtree <- extract.clade(phy, mrca.t1t2)
  subtreephen <- allphendata[sort(c(mrca.t1t2,Descendants(phy,mrca.t1t2,type="all"))),]
  subtreechanges <- sum(calcchangephy(subtree, subtreephen))

  C4 <- C2/subtreechanges
  cbind(C1,C2,C3,C4)
}
#
# Simulation of trait data based on the phylogeny, for statistical comparison #
make.simdata = function(nsim,phy,trait){
  # Create simulated data for significance testing convergence
  # Simulate phenotyic data, using variance-covariance matrix and phyloeny specified by ratematrix()
  phendata.sim = sim.char(phy, ratematrix(phy, trait), nsim=nsim, model = c("BM"), root=0)
  # Run multianc on each set of simulated data to get simulated data for internal nodes and tips
  allphendata.sim = c()

  message("Simulating Phenotypic Evolution")
  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style = 3)
  step_n = 1

  for (i in 1:nsim){
    x = phendata.sim[,,i]
    y = multianc(phy,x) %>%
      "rownames<-"(c(phy$tip.label,seq(from  = length(phy$tip.label)+1, to = length(phy$tip.label)+phy$Nnode, by = 1)))
    allphendata.sim = abind(allphendata.sim,y,along=3)

    setTxtProgressBar(pb,step_n)
    step_n <- step_n + 1
  }
  # Returns a 3D array of simulated data, with ancestral states calculated
  allphendata.sim
}
#
##

## Final complete function ##
calc.conv.rates.sig = function(distmat,phy,tax1,tax2,allphendata,allphendata.sim){
  #distmat = distance matrix of all tips and all nodes
  #phy = phylogeny of taxa in sample
  #tax1 = name of first taxon of interest
  #tax2 = name of second taxon of interest
  #allphendata = trait dataset for all nodes and tips (result of multianc)
  #allphendata.sim = result of make.simdata
  # Observed levels of convergence
  obs <- calc.conv.rates(distmat,phy,tax1,tax2,allphendata)
  # Empty vectors to store results of simulation tests
  C1s <- c(); C2s <- c(); C3s <- c(); C4s <- c()
  # identify nsim from simulated data
  nsim <- dim(allphendata.sim)[[3]]
  # Loop through simulated data and run tests for convergence on each
  for (k in 1:nsim) {
    rsimphendata <- allphendata.sim[, , k]
    simdist = as.matrix(dist(rsimphendata))
    simdist = reshape2::melt(as.matrix(simdist), varnames = c("n1", "n2"))
    simresults = calc.conv.rates(simdist,phy,tax1,tax2,rsimphendata)
    C1s <- c(C1s, simresults[1])
    C2s <- c(C2s, simresults[2])
    C3s <- c(C3s, simresults[3])
    C4s <- c(C4s, simresults[4])
  }
  # Check number of times simulated converegence was greater than observed convergence
  C1greater <- 0
  C2greater <- 0
  C3greater <- 0
  C4greater <- 0
  for (i in 1:nsim) {
    if (C1s[i] >= obs[1]) {
      C1greater <- C1greater + 1
    }
    if (C2s[i] >= obs[2]) {
      C2greater <- C2greater + 1
    }
    if (C3s[i] >= obs[3]) {
      C3greater <- C3greater + 1
    }
    if (C4s[i] >= obs[4]) {
      C4greater <- C4greater + 1
    }
  }
  # Generate p-value
  sig <- cbind(C1greater/(nsim + 1), C2greater/(nsim + 1), C3greater/(nsim + 1), C4greater/(nsim + 1))
  colnames(sig) <- c("C1sig", "C2sig", "C3sig", "C4sig")
  # Return both observed values, and p-values to denote significance
  answer = cbind(obs,sig)
  answer
}
##

###
