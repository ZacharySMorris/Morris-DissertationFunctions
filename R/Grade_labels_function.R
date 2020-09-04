###Grade_labels Function###
Gradelabels <- function (tree = NULL, text, node, node_exclude, offset = NULL, wing.length = NULL,
                         cex = 1, orientation = "vertical")
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (is.null(tree)) {
    wing.length <- 1
    if (is.null(offset))
      offset <- 8
    tree <- list(edge = lastPP$edge, tip.label = 1:lastPP$Ntip,
                 Nnode = lastPP$Nnode)
    H <- matrix(lastPP$xx[tree$edge], nrow(tree$edge), 2)
    tree$edge.length <- H[, 2] - H[, 1]
    class(tree) <- "phylo"
  }
  if (is.null(offset))
    offset <- 0.5
  xx <- mapply(GradelabelSubTree, node, node_exclude, text, MoreArgs = list(tree = tree,
                pp = lastPP, offset = offset, wl = wing.length, cex = cex,
                orientation = orientation))
}


## internal function used by cladelabels
## written by Liam J. Revell 2014, 2015
GradelabelSubTree<-function(tree,nn,nn_exclude,label,pp,offset,wl,cex,orientation){
  if(is.null(wl)) wl<-1
  tree<-reorder(tree)
  tips<-getDescendants(tree,nn)
  tips_exclude<-getDescendants(tree,nn_exclude)
  tips<-tips[-match(tips_exclude,tips)]
  tips<-tips[tips<=Ntip(tree)]
  ec<-0.7 ## expansion constant
  sw<-pp$cex*max(strwidth(tree$tip.label[tips]))
  sh<-pp$cex*max(strheight(tree$tip.label))
  cw<-mean(strwidth(LETTERS)*cex)
  h<-max(sapply(tips,function(x,tree)
    nodeHeights(tree)[which(tree$edge[,2]==x),2],
    tree=tree))+sw+offset*cw
  y<-range(pp$yy[tips])
  lines(c(h,h),y+ec*c(-sh,sh))
  lines(c(h-wl*cw,h),
        c(y[1]-ec*sh,y[1]-ec*sh))
  lines(c(h-wl*cw,h),
        c(y[2]+ec*sh,y[2]+ec*sh))
  text(h+cw,mean(y),
       label,srt=if(orientation=="horizontal") 0 else 90,
       adj=if(orientation=="horizontal") 0 else 0.5,cex=cex)
}
