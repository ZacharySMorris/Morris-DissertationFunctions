###Major axis analysis using pca###
PmaxGMM <- function(X,PCData,PCs=c(1,2),PrintPoints=FALSE)
{
for (i in 1:length(X$taxa)){
  if(length(X$CSize[[i]]) > 4){
    Pmatrix <- princomp(X$PCvalues[[i]][,1:4])

    ###Create Max and Min newPC values lists
    Min_List <- apply(Pmatrix$scores,2, "min")
    Max_List <- apply(Pmatrix$scores,2, "max")

    Pmax_Matrix <- matrix(0, byrow=FALSE, nrow = length(c(Min_List,Max_List)), ncol=4)
    for(j in 1:length(Min_List)){
      Pmax_Matrix[j,j] <- Min_List[j]
      Pmax_Matrix[4+j,j] <- Max_List[j]
    }

    Pmax_Col_names <- colnames(Pmatrix$loadings)
    Pmax_Row_names <- c(paste(Pmax_Col_names, "min", sep = "_"),paste(Pmax_Col_names, "max", sep = "_"))

    colnames(Pmax_Matrix) <- Pmax_Col_names
    rownames(Pmax_Matrix) <- Pmax_Row_names

    Transformed_Pmax_Matrix <- t(t(Pmax_Matrix %*% t(Pmatrix$loadings)) + (Pmatrix$center))

    Pmax <- lm(Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max"),PCs[2]]
               ~ Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max"),PCs[1]])
    P2 <- lm(Transformed_Pmax_Matrix[c("Comp.2_min","Comp.2_max"),PCs[2]]
             ~ Transformed_Pmax_Matrix[c("Comp.2_min","Comp.2_max"),PCs[1]])
    P3 <- lm(Transformed_Pmax_Matrix[c("Comp.3_min","Comp.3_max"),PCs[2]]
             ~ Transformed_Pmax_Matrix[c("Comp.3_min","Comp.3_max"),PCs[1]])
    P4 <- lm(Transformed_Pmax_Matrix[c("Comp.4_min","Comp.4_max"),PCs[2]]
             ~ Transformed_Pmax_Matrix[c("Comp.4_min","Comp.4_max"),PCs[1]])

    LineColor <- unique(X$color[[i]])

###Plotting
    Xlim<-c(floor(min(PCData$pc.scores[,1])*10)/10,ceiling(max(PCData$pc.scores[,1])*10)/10)
    Ylim<-c(floor(min(PCData$pc.scores[,2])*10)/10,ceiling(max(PCData$pc.scores[,2])*10)/10)

    plot(0, 0, type = "n",
         xlim = Xlim,
         ylim = Ylim,
         xlab = paste("Principal Componenet ", PCs[1], " (", round(100*PCData$pc.summary$importance[2,PCs[1]], digits = 1), "%)", sep = ""),
         ylab = paste("Principal Componenet ", PCs[2], " (", round(100*PCData$pc.summary$importance[2,PCs[2]], digits = 1), "%)", sep = ""),
         axes = FALSE,
         frame.plot = FALSE,
         asp=F)

    axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
    axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1])
    abline(h=0, lty=3)
    abline(v=0, lty=3)

    clip(-0.3,0.3,-0.3,0.3)
    if (PrintPoints == TRUE){
      points(X$PCvalues[[i]][,PCs[1]], X$PCvalues[[i]][,PCs[2]], pch=X$shape[[i]], cex=X$size[[i]], bg=alpha(X$color[[i]], 0.75), asp=F)
      mtext(X$taxa[[i]],, at = -0.25)
    }else

    clip(-0.3,0.3,-0.3,0.3)

    dataEllipse(X$PCvalues[[i]][,PCs[1]],
                X$PCvalues[[i]][,PCs[2]],
                add = TRUE, plot.points = FALSE, levels = c(0.75),
                col = "grey", fill = TRUE)

    dataEllipse(X$PCvalues[[i]][,PCs[1]],
                X$PCvalues[[i]][,PCs[2]],
                add = TRUE, plot.points = FALSE, levels = c(0.75),
                col = LineColor, fill = FALSE)
    # ###
    # dataEllipse(Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max","Comp.2_min","Comp.2_max"),1],
    #             Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max","Comp.2_min","Comp.2_max"),2],
    #             add = TRUE, plot.points = FALSE, levels = c(0.75),
    #             col = "grey", fill = TRUE)
    #
    # dataEllipse(Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max","Comp.2_min","Comp.2_max"),1],
    #             Transformed_Pmax_Matrix[c("Comp.1_min","Comp.1_max","Comp.2_min","Comp.2_max"),2],
    #             add = TRUE, plot.points = TRUE, levels = c(0.75),
    #             col = CrocDorsal.PlottingValues$legendcolor[i], fill = FALSE)
    # ###

    clip(
      min(X$PCvalues[[i]][,PCs[1]]),
      max(X$PCvalues[[i]][,PCs[1]]),
      min(X$PCvalues[[i]][,PCs[2]]),
      max(X$PCvalues[[i]][,PCs[2]])
    )

    abline(coef(Pmax), col = LineColor, lwd=3, lty=1 )
    abline(coef(P2), col = LineColor, lwd=3, lty=2)
    abline(coef(P3), col = LineColor, lwd=3, lty=3)
    abline(coef(P4), col = LineColor, lwd=3, lty=4)

    clip(-0.3,max(Xlim),-0.3,max(Ylim))

    text(c(0.2,0.2,0.2,0.2), c(-0.15,-0.16,-0.17,-0.18),
         labels = c("Pmax","P2","P3","P4"), pos=4)
    segments(0.18,-0.15,0.2,-0.15, lty=1)
    segments(0.18,-0.16,0.2,-0.16, lty=2)
    segments(0.18,-0.17,0.2,-0.17, lty=3)
    segments(0.18,-0.18,0.2,-0.18, lty=4)

  }

  else{next}

  }

}
