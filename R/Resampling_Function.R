


###PCvCS Regression Confidence Interval Function###
Major_Axis_Regression_ConfInt <- function(X, PCs, MA_number, method="bootstrap", iter = 999, CI_values=c(0.975,0.025)){

  if(length(CI_values) > 2){
    warning("More than two values supplied for confidence interval. Only maximum and minimum values used.")
  }
  CI_upper <- max(CI_values)
  CI_lower <- min(CI_values)

  Groups <- X$taxa
  Groups_to_remove <- c()

  #Specify method of resampling
  method <- method

  ##assign PCs to be assessed##
  if(length(PCs) > 2){
    warning("More than two values supplied as PCs to assess. Only the first two values were used.")
    }
  X_PC <- PCs[[1]]
  Y_PC <- PCs[[2]]
  ###

  NEW_Xs <- list()
  Pmatrix_list <- list()
  Model_list <- list()
  Loadings_list <- list()
  Slope_list <- list()
  Intercept_list <- list()
  Elevation_list <- list()
  TransformedPmax_list <- list()
  Mean_X_list <- list()
  resampled_loadings <- list()
  resampled_slope_list <- list()
  resampled_intercept_list <- list()
  resampled_elevation_list <- list()
  resampled_transformedPmax_list <- list()
  preds <- list()
  preds_min <- list()
  preds_max <- list()
  preds_CIupper <- list()
  preds_CIlower <- list()
  slope_CI <- list()
  elevation_CI <- list()
  intercept_CI <- list()
  loadings_CI <- list()

# i=4
# X_PC <- 1
# Y_PC <- 2
  for (i in 1:length(Groups)){
    if (nrow(X$PCvalues[[i]]) < 4) {
      warning(Groups[i], " has too few individuals for resampling CI and was dropped from the analysis")
      Groups_to_remove <- append(Groups_to_remove, i)
      next}else{
      Group_Pmatrix <- Pmatrix_pca(X$PCvalues[[i]]) #persform Pmatrix calculation using PCA for current group
      Group_LM <- Pmax_lms(Group_Pmatrix, PC_Axes = c(X_PC,Y_PC)) #perform Major Axis regression for current group

      Group_resampled <- resampled.pcs(pc.scores=X$PCvalues[[i]], method=method, iter=iter)
      sequence_length = dim(X$PCvalues[[i]])[1]
      NEW_X <- seq(min(X$PCvalues[[i]][,X_PC]),max(X$PCvalues[[i]][,X_PC]), length.out = sequence_length)
      Mean_X <- mean(X$PCvalues[[i]][,X_PC])

      Group_Elevation <- ((Mean_X * Group_LM$Coefficients_Matrix[,"slope"])  + Group_LM$Coefficients_Matrix[,"intercept"])


      Pred_Matrix <- matrix(nrow=Group_resampled$iter, ncol=length(NEW_X))


      resampled_slope <- matrix(nrow=Group_resampled$iter, ncol=Group_Pmatrix$N_comp)
      resampled_intercept <-  matrix(nrow=Group_resampled$iter, ncol=Group_Pmatrix$N_comp)
      resampled_elevation <- matrix(nrow=Group_resampled$iter, ncol=Group_Pmatrix$N_comp)
      resampled_loadings_matrix <- matrix(nrow=Group_resampled$iter, ncol=ncol(Group_resampled$original_PCs))
      resampled_L <- list()
      resampled_T <- list()

      for (j in 1:Group_resampled$iter){
        current_Pmatrix <- Pmatrix_pca(Group_resampled$resampled_PCs[[j]]) #persform Pmatrix calculation using PCA for current iteration

        #create and name loadings matrix for current iteration
        current_loadings <- current_Pmatrix$Pmatrix$loadings
        colnames(current_loadings) <- colnames(current_Pmatrix$major.axes)
        resampled_L[[j]] <- current_loadings

        current_LM <- Pmax_lms(current_Pmatrix, PC_Axes = c(X_PC,Y_PC)) #perform Major Axis regression for current iteration
        resampled_T[[j]] <- current_LM$transformed.major.axe

        Pred_Matrix[j,] <- (current_LM$Coefficients_Matrix[MA_number,"slope"]*NEW_X) + current_LM$Coefficients_Matrix[MA_number,"intercept"]

        resampled_slope[j,] <- t(current_LM$Coefficients_Matrix[,"slope"])
        resampled_intercept[j,] <- t(current_LM$Coefficients_Matrix[,"intercept"])
        resampled_elevation[j,] <- ((Mean_X * current_LM$Coefficients_Matrix[,"slope"])  + current_LM$Coefficients_Matrix[,"intercept"])
        resampled_loadings_matrix[j,1:nrow(current_loadings)] <- current_loadings[,MA_number]
      }
      message("Completed modeling for resampled versions of the ", Groups[i], " dataset.")
      colnames(resampled_slope) <- rownames(current_LM$Coefficients_Matrix)
      colnames(resampled_intercept) <- rownames(current_LM$Coefficients_Matrix)
      colnames(resampled_elevation) <- rownames(current_LM$Coefficients_Matrix)



      if(anyNA(Pred_Matrix)==TRUE){
        warning("Some error in the matrix of predicted values based on the resampled versions of the ", Groups[i], " dataset created NA values. Use Pmax_lms function to investigate model error.")
        print(Pred_Matrix)

      }else{
        #assign values to variables
        NEW_Xs[[X$taxa[i]]] <- NEW_X
        Pmatrix_list[[X$taxa[i]]] <- Group_Pmatrix$Pmatrix
        Model_list[[X$taxa[i]]] <- Group_LM$MajorAxes_lms
        Loadings_list[[X$taxa[i]]] <- Group_Pmatrix$Pmatrix$loadings
        Slope_list[[X$taxa[i]]] <- Group_LM$Coefficients_Matrix[,"slope"]
        Intercept_list[[X$taxa[i]]] <- Group_LM$Coefficients_Matrix[,"intercept"]
        Elevation_list[[X$taxa[i]]] <- Group_Elevation
        TransformedPmax_list[[X$taxa[i]]] <- Group_LM$transformed.major.axe
        Mean_X_list[[X$taxa[i]]] <- Mean_X

        resampled_loadings[[X$taxa[i]]] <- resampled_L
        resampled_transformedPmax_list[[X$taxa[i]]] <- resampled_T
        resampled_slope_list[[X$taxa[i]]] <- resampled_slope
        resampled_intercept_list[[X$taxa[i]]] <- resampled_intercept
        resampled_elevation_list[[X$taxa[i]]] <- resampled_elevation
        preds[[X$taxa[i]]] <- Pred_Matrix
        preds_min[[X$taxa[i]]] <- apply(Pred_Matrix, 2, min)
        preds_max[[X$taxa[i]]] <- apply(Pred_Matrix, 2, max)
        preds_CIupper[[X$taxa[i]]] <- apply(Pred_Matrix, 2, quantile, probs = CI_upper)
        preds_CIlower[[X$taxa[i]]] <- apply(Pred_Matrix, 2, quantile, probs = CI_lower)
        slope_CI[[X$taxa[i]]] <- apply(resampled_slope, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE)
        elevation_CI[[X$taxa[i]]] <- apply(resampled_elevation, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE)
        intercept_CI[[X$taxa[i]]] <- apply(resampled_intercept, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE)
        loadings_CI[[X$taxa[i]]] <- apply(resampled_loadings_matrix, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE)

      }
    }
  }
if (length(Groups_to_remove) > 0){
  Groups <- Groups[-Groups_to_remove]
}

    out <- list(Groups = Groups, PCs = PCs, iter = Group_resampled$iter, MA_number = MA_number,
                Pmatrix_list = Pmatrix_list, Model_list = Model_list,
                Slope_list = Slope_list, Intercept_list = Intercept_list, Elevation_list = Elevation_list,
                Loadings_list = Loadings_list, transformedPmax = TransformedPmax_list,
                resampled_loadings = resampled_loadings, resampled_transformedPmax = resampled_transformedPmax_list,
                resampled_slopes = resampled_slope_list, resampled_intercepts = resampled_intercept_list, resampled_elevation_list = resampled_elevation_list,
                NEW_Xs = NEW_Xs, Mean_X_list = Mean_X_list, preds_min = preds_min, preds_max = preds_max,
                preds_CIupper = preds_CIupper, preds_CIlower = preds_CIlower,
                slope_CI = slope_CI, elevation_CI = elevation_CI, intercept_CI = intercept_CI, loadings_CI = loadings_CI)
    class(out) = "ConfIntList"
    out

}
###


###Trajectory CI Polygon plotting###
CI_Poly_Plot <- function(X, PCData, ConfIntList, Colors, OnePlot=FALSE){

  MA_number <- ConfIntList$MA_number
  ConfIntList <- ConfIntList
  Slopes <- ConfIntList$Slope_list
  Intercepts <- ConfIntList$Intercept_list
  Groups <- ConfIntList$Groups
  PCs <- ConfIntList$PCs
  LineColor <- Colors[match(ConfIntList$Groups,names(Colors))]
  PCData <- PCData

  if(OnePlot == TRUE){
    Xlim<-c(floor(min(PCData$pc.scores[,PCs[[1]]])*10)/10,ceiling(max(PCData$pc.scores[,PCs[[1]]])*10)/10)
    Ylim<-c(floor(min(PCData$pc.scores[,PCs[[2]]])*10)/10,ceiling(max(PCData$pc.scores[,PCs[[2]]])*10)/10)

    plot(0, 0, type = "n",
         xlim = Xlim,
         ylim = Ylim,
         xlab = paste("Principal Component ", PCs[[1]], " (", round(100*PCData$pc.summary$importance[2,PCs[[1]]], digits = 1), "%)", sep = ""),
         ylab = paste("Principal Component ", PCs[[2]], " (", round(100*PCData$pc.summary$importance[2,PCs[[2]]], digits = 1), "%)", sep = ""),
         axes = FALSE,
         frame.plot = FALSE,
         asp=F)

    axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
    axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
    clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
    abline(h=0, lty=3)
    abline(v=0, lty=3)

    for (i in 1:length(Groups)){
      clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
      polygon(c(rev(ConfIntList$NEW_Xs[[i]]), ConfIntList$NEW_Xs[[i]]),
              c(rev(ConfIntList$preds_CIupper[[i]]), ConfIntList$preds_CIlower[[i]]),
              col = alpha('grey',0.6), border = NA)

      clip(
        min(ConfIntList$NEW_Xs[[i]]),
        max(ConfIntList$NEW_Xs[[i]]),
        min(ConfIntList$preds_CIlower[[i]]),
        max(ConfIntList$preds_CIupper[[i]])
      )

      lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIlower[[i]], lty = 'dashed', col = LineColor[i])
      lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIupper[[i]], lty = 'dashed', col = LineColor[i])
      abline(Intercepts[[i]][MA_number],Slopes[[i]][MA_number], col = LineColor[i], lwd=3, lty=1 )

      # clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
      # dataEllipse(X$PCvalues[[i]][,1],
      #             X$PCvalues[[i]][,2],
      #             add = TRUE, plot.points = FALSE, levels = c(0.75),
      #             col = LineColor[i], fill = FALSE)
    }

  }else{
    for (i in 1:length(Groups)){
      Xlim<-c(floor(min(PCData$pc.scores[,PCs[[1]]])*10)/10,ceiling(max(PCData$pc.scores[,PCs[[1]]])*10)/10)
      Ylim<-c(floor(min(PCData$pc.scores[,PCs[[2]]])*10)/10,ceiling(max(PCData$pc.scores[,PCs[[2]]])*10)/10)

      plot(0, 0, type = "n",
           xlim = Xlim,
           ylim = Ylim,
           xlab = paste("Principal Component ", PCs[[1]], " (", round(100*PCData$pc.summary$importance[2,PCs[[1]]], digits = 1), "%)", sep = ""),
           ylab = paste("Principal Component ", PCs[[2]], " (", round(100*PCData$pc.summary$importance[2,PCs[[2]]], digits = 1), "%)", sep = ""),
           axes = FALSE,
           frame.plot = FALSE,
           asp=F)

      axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
      axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
      clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
      abline(h=0, lty=3)
      abline(v=0, lty=3)

      mtext(Groups[[i]], col=LineColor[i], at = -0.25)

      polygon(c(rev(ConfIntList$NEW_Xs[[i]]), ConfIntList$NEW_Xs[[i]]),
              c(rev(ConfIntList$preds_CIupper[[i]]), ConfIntList$preds_CIlower[[i]]),
              col = alpha('grey',0.6), border = NA)

      clip(
        min(ConfIntList$NEW_Xs[[i]]),
        max(ConfIntList$NEW_Xs[[i]]),
        min(ConfIntList$preds_CIlower[[i]]),
        max(ConfIntList$preds_CIupper[[i]])
      )

      lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIlower[[i]], lty = 'dashed', col = LineColor[i])
      lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIupper[[i]], lty = 'dashed', col = LineColor[i])
      abline(Intercepts[[i]][MA_number],Slopes[[i]][MA_number], col = LineColor[i], lwd=3, lty=1 )


      # clip(-0.3,3,-0.3,3)
      # dataEllipse(X$PCvalues[[i]][,1],
      #             X$PCvalues[[i]][,2],
      #             add = TRUE, plot.points = FALSE, levels = c(0.75),
      #             col = LineColor[i], fill = FALSE)
    }
  }

}
###

###Function to calculate CI interval for slope and elevation###
MA_CI_calculation <- function(ConfIntList, Comp_ConfIntList){

  ConfIntList <- ConfIntList
  Comp_ConfIntList <- Comp_ConfIntList

  MA_number <- ConfIntList$MA_number
  Slopes <- ConfIntList$Slope_list
  Loadings <- ConfIntList$Loadings_list
  Intercepts <- ConfIntList$Intercept_list
  Slopes_CI <- ConfIntList$slope_CI
  Loadings_CI <- ConfIntList$loadings_CI
  Intercepts_CI <- ConfIntList$intercept_CI
  Groups <- ConfIntList$Groups

  Comp_Slopes <- Comp_ConfIntList$Slope_list
  Comp_Loadings <- Comp_ConfIntList$Loadings_list
  Comp_Intercepts <- Comp_ConfIntList$Intercept_list
  Comp_Slopes_CI <- Comp_ConfIntList$slope_CI
  Comp_Loadings_CI <- Comp_ConfIntList$Loadings_CI
  Comp_Intercepts_CI <- Comp_ConfIntList$intercept_CI
  Comp_Groups <- Comp_ConfIntList$Groups

  Slope_Diff_Table <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  Int_Diff_Table <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))

  Table1s <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  # Table1l <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups)*4, dimnames=(list(Groups,rep(Comp_Groups,times=4))))
  Table1i <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))

  for (i in 1:length(Groups)){
    temp_group <- Groups[[i]]
    temp_slope <- Slopes[[i]][MA_number]
    # temp_loadings <- Loadings[[i]][,MA_number]
    temp_intercept <- Intercepts[[i]][MA_number]

    for (j in 1:length(Comp_Groups)){
      temp_group <- Comp_Groups[[j]]
      temp_comp_slope <- Comp_Slopes[[j]][MA_number]
      temp_comp_int <- Comp_Intercepts[[j]][MA_number]
      temp_slopeCI <- Comp_Slopes_CI[[j]][,MA_number]
      # temp_loadingsCI <- Comp_Loadings_CI[[j]]
      temp_interceptCI <- Comp_Intercepts_CI[[j]][,MA_number]

      temp_ans <- temp_slope >= min(temp_slopeCI) & temp_slope <= max(temp_slopeCI)
      # temp_loadings_ans <- temp_loadings[1:4] >= temp_loadingsCI[1,1:4] & temp_loadings[1:4] <= temp_loadingsCI[2,1:4]
      temp_intercept_ans <- temp_intercept >= min(temp_interceptCI) & temp_intercept <= max(temp_interceptCI)
      temp_slope_diff <- abs(temp_slope - temp_comp_slope)
      temp_int_diff <- abs(temp_intercept-temp_comp_int)

      Slope_Diff_Table[i,j] <- temp_slope_diff
      Int_Diff_Table[i,j] <- temp_int_diff
      Table1s[i,j] <- temp_ans
      # Table1l[i,j] <- temp_loadings_ans
      Table1i[i,j] <- temp_intercept_ans
    }
  }

  Table2s <- matrix(data = NA , nrow = length(Comp_Groups), ncol = length(Groups), dimnames=(list(Comp_Groups,Groups)))
  # Table2l <- matrix(data = NA , nrow = length(Comp_Groups), ncol = length(Groups)*4, dimnames=(list(Comp_Groups,rep(Groups,times=4))))
  Table2i <- matrix(data = NA , nrow = length(Comp_Groups), ncol = length(Groups), dimnames=(list(Comp_Groups,Groups)))

  for (i in 1:length(Comp_Groups)){
    temp_group <- Comp_Groups[[i]]
    temp_slope <- Comp_Slopes[[i]][MA_number]
    # temp_loadings <- Comp_Loadings[[i]][MA_number]
    temp_intercept <- Comp_Intercepts[[i]][MA_number]

    for (j in 1:length(Groups)){
      temp_slopeCI <- Slopes_CI[[j]][,MA_number]
      # temp_loadingsCI <- Loadings_CI[[j]][,MA_number]
      temp_interceptCI <- Intercepts_CI[[j]][,MA_number]

      temp_ans <- temp_slope >= min(temp_slopeCI) & temp_slope <= max(temp_slopeCI)
      # temp_loadings_ans <- temp_loadings[1:4] >= temp_loadingsCI[1,1:4] & temp_loadings[1:4] <= temp_loadingsCI[2,1:4]
      temp_intercept_ans <- temp_intercept >= min(temp_interceptCI) & temp_intercept <= max(temp_interceptCI)


      Table2s[i,j] <- temp_ans
      # Table2l[i,j] <- temp_loadings_ans
      Table2i[i,j] <- temp_intercept_ans
    }
  }

  out <- list(Slope_Differences_Table = Slope_Diff_Table,
              Intercept_Differences_Table = Int_Diff_Table,
              Groups_1v2_slope = Table1s, Groups_1v2_intercept = Table1i,
              # Groups_1v2_loadings = Table1l, Groups_2v1_loadings = Table2l,s
              Groups_2v1_slope = Table2s, Groups_2v1_intercept = Table2i)

  out

}
###

###Function to calculate CI interval for slope and elevation###
MA_P_calculation <- function(ConfIntList, Comp_ConfIntList){

  ConfIntList <- ConfIntList
  Comp_ConfIntList <- Comp_ConfIntList

  MA_number <- ConfIntList$MA_number

  Groups <- ConfIntList$Groups
  Slopes <- ConfIntList$Slope_list
  Intercepts <- ConfIntList$Intercept_list
  Loadings <- ConfIntList$Loadings_list
  Resampled_Slopes <- ConfIntList$resampled_slopes
  Resampled_Intercepts <- ConfIntList$resampled_intercepts
  Resampled_Loadings <- ConfIntList$resampled_loadings

  Comp_Groups <- Comp_ConfIntList$Groups
  Comp_Slopes <- Comp_ConfIntList$Slope_list
  Comp_Intercepts <- Comp_ConfIntList$Intercept_list
  Comp_Loadings <- Comp_ConfIntList$Loadings_list
  Comp_Resampled_Slopes <- Comp_ConfIntList$resampled_slopes
  Comp_Resampled_Intercepts <- Comp_ConfIntList$resampled_intercepts
  Comp_Resampled_Loadings <- Comp_ConfIntList$resampled_loadings

  Slope_Diff_Table <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  Intercept_Diff_Table <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))

  Table1s <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  # Table1l <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups)*4, dimnames=(list(Groups,rep(Comp_Groups,times=4))))
  Table1i <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))

  normality_tests <- list()

  for (i in 1:length(Groups)){
    temp_group <- Groups[[i]]
    temp_slope <- Slopes[[i]][MA_number]
    temp_intercept <- Intercepts[[i]][MA_number]
    # temp_loadings <- Loadings[[i]][,MA_number]

    temp_resampled_slopes <- Resampled_Slopes[[i]][,MA_number]
    temp_resampled_intercepts <- Resampled_Intercepts[[i]][,MA_number]
    # temp_loadings_dist <- Resampled_Loadings[[i]][[1]][,MA_number]

    temp_slope_dist <- sample(temp_resampled_slopes) - temp_resampled_slopes
    temp_intercept_dist <- sample(temp_resampled_intercepts) - temp_resampled_intercepts


    for (j in 1:length(Comp_Groups)){
      temp_comp_group <- Comp_Groups[[j]]
      temp_comp_slope <- Comp_Slopes[[j]][MA_number]
      temp_comp_intercept <- Comp_Intercepts[[j]][MA_number]
      # temp_loadings <- Comp_Loadings_CI[[i]][,MA_number]

      temp_comp_resampled_slopes <- Comp_Resampled_Slopes[[j]][,MA_number]
      temp_comp_resampled_intercepts <- Comp_Resampled_Intercepts[[j]][,MA_number]
      # temp_comp_loadings_dist <- Comp_Resampled_Loadings[[i]][[1]][,MA_number]

      temp_comp_slope_dist <- sample(temp_comp_resampled_slopes) - temp_comp_resampled_slopes
      temp_comp_intercept_dist <- sample(temp_comp_resampled_intercepts) - temp_comp_resampled_intercepts

      Pooled_slope_dist <- c(temp_slope_dist,temp_comp_slope_dist)
      Pooled_intercept_dist <- c(temp_intercept_dist,temp_comp_intercept_dist)

      normality_tests[paste(temp_group, "&", temp_comp_group, "pooled variance of slopes", sep = " ")] <- shapiro.test(Pooled_slope_dist)[[2]]
      normality_tests[paste(temp_group, "&", temp_comp_group, "pooled variance of intercepts", sep = " ")] <- shapiro.test(Pooled_intercept_dist)[[2]]

      temp_slope_diff <- abs(temp_slope - temp_comp_slope)
      temp_slope_P_val <- sum(abs(temp_slope_diff) < abs(Pooled_slope_dist)) / length(Pooled_slope_dist)

      temp_intercept_diff <- abs(temp_intercept - temp_comp_intercept)
      temp_intercept_P_val <- sum(abs(temp_intercept_diff) < abs(Pooled_intercept_dist)) / length(Pooled_intercept_dist)


      Slope_Diff_Table[i,j] <- temp_slope_diff
      Intercept_Diff_Table[i,j] <- temp_intercept_diff
      Table1s[i,j] <- temp_slope_P_val
      Table1i[i,j] <- temp_intercept_P_val
    }
  }

  out <- list(Slope_Differences_Table = Slope_Diff_Table,
              Intercept_Differences_Table = Intercept_Diff_Table,
              Slope_P_values = Table1s,
              Intercept_P_values = Table1i,
              normality_tests = normality_tests)

  out

}
###

###Function to calculate CI interval for slope and elevation###
MA_CI_comp <- function(ConfIntList, Comp_ConfIntList){

  ConfIntList <- ConfIntList
  Comp_ConfIntList <- Comp_ConfIntList

  MA_number <- ConfIntList$MA_number

  Slopes_CI <- ConfIntList$slope_CI
  Intercepts_CI <- ConfIntList$intercept_CI
  Loadings_CI <- ConfIntList$loadings_CI
  Groups <- ConfIntList$Groups

  Comp_Slopes_CI <- Comp_ConfIntList$slope_CI
  Comp_Intercepts_CI <- Comp_ConfIntList$intercept_CI
  Comp_Loadings_CI <- Comp_ConfIntList$loadings_CI
  Comp_Groups <- Comp_ConfIntList$Groups

  Table1s <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  Table1i <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))
  Table1l <- matrix(data = NA , nrow = length(Groups), ncol = length(Comp_Groups), dimnames=(list(Groups,Comp_Groups)))

  for (i in 1:length(Groups)){
    temp_group <- Groups[[i]]
    temp_slopes <- seq(max(Slopes_CI[[i]][,MA_number])*100,
                       min(Slopes_CI[[i]][,MA_number])*100) / 100
    temp_intercepts <- seq(max(Intercepts_CI[[i]][,MA_number])*100,
                           min(Intercepts_CI[[i]][,MA_number])*100) / 100
    temp_loadings <- seq(max(Loadings_CI[[i]][,MA_number])*100,
                         min(Loadings_CI[[i]][,MA_number])*100) / 100


    for (j in 1:length(Comp_Groups)){
      temp_group <- Comp_Groups[[j]]
      temp_slopeCI <- Comp_Slopes_CI[[j]][,MA_number]
      temp_interceptCI <- Comp_Intercepts_CI[[j]][,MA_number]
      temp_loadingCI <- Comp_Loadings_CI[[j]][,MA_number]

      temp_ans <- temp_slopes >= min(temp_slopeCI) & temp_slopes <= max(temp_slopeCI)
      temp_intercept_ans <- temp_intercepts >= min(temp_interceptCI) & temp_intercepts <= max(temp_interceptCI)
      temp_loading_ans <- temp_loadings >= min(temp_loadingCI) & temp_loadings <= max(temp_loadingCI)

      if (any(temp_ans) == TRUE){
        temp_ans <- TRUE
      }else {temp_ans <- FALSE}

      if (any(temp_intercept_ans) == TRUE){
        temp_intercept_ans <- TRUE
      }else {temp_intercept_ans <- FALSE}

      if (any(temp_loading_ans) == TRUE){
        temp_loading_ans <- TRUE
      }else {temp_loading_ans <- FALSE}


      Table1s[i,j] <- temp_ans
      Table1i[i,j] <- temp_intercept_ans
      Table1l[i,j] <- temp_loading_ans
    }
  }

  out <- list(Groups_1v2_slope = Table1s, Groups_1v2_intercept = Table1i, Groups_1v2_loading = Table1l)

  out

}
###

# ###Major Axis Analysis Wrapper Function###
# Major_Axis <- function(X, PCs, W, MajorAxes_N, plot="simple", method = "bootstrap"){
#
#   #load in datasets and variables
#   {X <- X
#   PCs <- PCs
#   W <- W
#   MajorAxes_N <- MajorAxes_N
#   plot <- plot
#   method = method
#
#   Groups <- X$taxa
#   N_groups <- length(Groups)
#
#   Pre_N_comp <- sapply(X$PCvalues, length)
#   N_comp <- min(Pre_N_comp)}
#
#   #Create blank lists
#   {Pmatrix_list <- list()
#   MA_loadings <- list()
#   MA_coefficients <- list()}
#
#   #Run loop to calculate Pmatrix and Pmax models for each subgroup and save to lists
#   for (i in 1:length(Groups)){
#     current_Pmatrix <- Pmatrix_pca(X$PCvalues[[i]],N_comp)
#     current_lm <- Pmax_lms(current_Pmatrix)
#
#     Pmatrix_list[[X$taxa[i]]] <- current_Pmatrix
#     MA_loadings[[X$taxa[i]]] <- current_Pmatrix$Pmatrix$loadings
#     MA_coefficients[[X$taxa[i]]] <- current_lm$Coefficients_Matrix
#
#   }
#
#   #Calculate the modified alpha level based on the number of pairwise comparisons to be made
#   CI_correction <- p.adjust(0.05, method = "bonferroni", n = (N_groups*(N_groups-1)))
#   corrected_CI_levels <- c(100,0) + c(-CI_correction,CI_correction)
#
#   #Calculate CI around major axes
#   MA_ConfIntList <- list()
#   for (i in 1:MajorAxes_N){
#     MA_ConfIntList[[paste("MajorAxis #", i, sep = "")]] <- Major_Axis_Regression_ConfInt(X = X, PCs = PCs, N_comp = N_comp, MA_number = i, method = method, CI_values=corrected_CI_levels)
#
#
#   }
#
#   for (i in 1:length(Fossils_on_Ontogeny.MajorAxis$Groups)){
#     Try_One[[Fossils_on_Ontogeny.MajorAxis$Groups[i]]] <- apply(Fossils_on_Ontogeny.MajorAxis$resampled_slopes[[i]], 2, quantile, probs = probs, na.rm=TRUE)
#   }
#
#   MA_ConfIntList <- Major_Axis_Regression_ConfInt(X = X, PCs = PCs, N_comp = N_comp, method = method, CI_values=corrected_CI_levels)
#
#   MA_SlopeCI <- list()
#   for (i in 1:length(MA_ConfIntList$Groups)){
#     MA_SlopeCI[MA_ConfIntList$$Groups[i]] <-  apply(MA_ConfIntList$resampled_slopes, 2, quantile, probs = probs, na.rm=TRUE)
#   }
#
#   MA_ConfIntList$Slope_list
#
#   matrix(data = unlist(Fossils_on_Ontogeny.MajorAxis$Slope_list), byrow = TRUE,
#          nrow = length(Fossils_on_Ontogeny.MajorAxis$Groups),
#          ncol = length(Fossils_on_Ontogeny.MajorAxis$Slope_list[[1]]),
#          dimnames = list(c(Fossils_on_Ontogeny.MajorAxis$Groups),names(Fossils_on_Ontogeny.MajorAxis$Slope_list[[1]])))
#
#   unlist(Fossils_on_Ontogeny.MajorAxis$resampled_slopes)
#
#   MA_SlopeCI[[i]]
#
#   for (i in 1:length(Fossils_on_Ontogeny.MajorAxis$Groups)){
#     Current_CIs <- apply(Fossils_on_Ontogeny.MajorAxis$resampled_slopes[[i]], 2, quantile, probs = probs, na.rm=TRUE)
#
#   }
#
#
#   MA_coefficients[[j]]
#
#   sapply(MA_ConfIntList[[i]]$Slope_list)
#   Fossils_on_Ontogeny.MajorAxis$Groups
#          Fossils_on_Ontogeny.MajorAxis$Slope_list
#          na.exclude(Fossils_on_Ontogeny.MajorAxis$resampled_slopes[[i]])
#          dim(Fossils_on_Ontogeny.MajorAxis$resampled_slopes[[i]])
#
#          Fossils_on_Ontogeny.MajorAxis$slope_CI
#
#   Fossils_on_Ontogeny.MajorAxis$Pmatrix_list[[i]]
#   Fossils_on_Ontogeny.MajorAxis$Slope_list[[i]]
#   Fossils_on_Ontogeny.MajorAxis$Intercept_list[[i]]
#   Fossils_on_Ontogeny.MajorAxis$Loadings_list[[i]]
#
#
#   #Plotting
#   if (plot="CI Plot"){
#     CI_Poly_Plot((X = X, PCs = PCs, ConfIntList = MA_ConfIntList, Colors = W$legendcolor, OnePlot=FALSE)
#   }
#   if (plot="Simple"){
#     CI_Poly_Plot((X = X, PCs = PCs, ConfIntList = MA_ConfIntList, Colors = W$legendcolor, OnePlot=TRUE)
#   }
#
#
#
# }
#
#
#
#
# ###



# ###Function incorporating everything to perform resampling to calculate CI
# PmaxGMM <- function(X,iter=999){
#
#   group_names <- X$taxa
#   resample_n <- iter
#
#   for (i in 1:length(group_names)){
#     group_to_resample <-X$PCvalues[[i]]
#     specimen_length <- c(1:dim(group_to_resample)[1])
#   }
#
# }

# ### Function for calculating Major Axes trajectories in the original morphospace ###
# Pmax_lms  <- function(Pmatrix_out,
#                       major.axes=Pmatrix_out$major.axes,
#                       N_comp=Pmatrix_out$N_comp,
#                       PC_Axes=c(1,2)){
#
#   if (class(Pmatrix_out)=="Pmatrix"){
#     major.axes <- Pmatrix_out$major.axes
#     N_comp <- Pmatrix_out$N_comp
#     Pmatrix <- Pmatrix_out$Pmatrix
#   }
#
#   # if (class(Pmatrix)=="princomp"){
#   #   Pmatrix <- Pmatrix
#   #   major.axes <- major.axes
#   #   N_comp <- N_comp
#   # }
#
#   ###
#   # Uses the loadings and center from the major axis PCA to transform Max
#   # and Min values along each major axis into PC scores in the original
#   # morphospace. This can then be used in lm to calculate coefficients
#   # of each major axis when translated into the original morphospace.
#   ###
#
#   transformed.major.axe <- t(t(major.axes %*% t(Pmatrix$loadings)) + (Pmatrix$center))
#   names <- rownames(transformed.major.axe)
#   axes_names <- colnames(major.axes)
#
#   if (length(PC_Axes)>2){
#     message("Multiple PCs were input, therefore the first entry was used an independent variable and others were used as dependent variables in multivariate linear model.")
#     First_PC <- PC_Axes[1]
#     Multi_PC <- PC_Axes[-1]
#     PC_Coefficient_dimnames <- sort(apply(expand.grid(paste("PC", Multi_PC, sep = ""), c("intercept", "slope")), 1, paste, collapse="."))
#     Coefficients_Matrix <- matrix(nrow=N_comp, ncol=2*length(Multi_PC), dimnames = list(axes_names,PC_Coefficient_dimnames))
#   }else{
#     First_PC <- PC_Axes[1]
#     Second_PC <- PC_Axes[2]
#     Coefficients_Matrix <- matrix(nrow=N_comp, ncol=2, dimnames = list(axes_names,c("intercept", "slope")))
#     }
#
#   #cycles through each major axis and calculates regression coefficients
#   MajorAxes_lms <- list()
#   for (i in 1:N_comp){
#
#     if (length(PC_Axes)==2){
#       Model_temp <- lm(transformed.major.axe[c(names[[i]],names[[i+N_comp]]),Second_PC]
#                        ~ transformed.major.axe[c(names[[i]],names[[i+N_comp]]),First_PC])
#       MajorAxes_lms[[i]] <- Model_temp
#       Coefficients_Matrix[i,] <- Model_temp$coefficients
#     }
#     if (length(PC_Axes)>2){
#       Model_temp <- lm(transformed.major.axe[c(names[[i]],names[[i+N_comp]]),Multi_PC]
#                        ~ transformed.major.axe[c(names[[i]],names[[i+N_comp]]),First_PC])
#       MajorAxes_lms[[i]] <- Model_temp
#       Coefficients_Matrix[i,] <- Model_temp$coefficients
#     }
#
#   }
#
#   out <- list(transformed.major.axe = transformed.major.axe,
#               Coefficients_Matrix = Coefficients_Matrix,
#               MajorAxes_lms = MajorAxes_lms)
#   class(out) <- "Pmax_lms"
#   out
# }
# ###
