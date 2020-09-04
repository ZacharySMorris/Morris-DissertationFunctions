### Functions to perform resampling and calculate CI for major axes ###

## R-squared calculations for major axis analysis ##
MA.r.squared <- function(Pmatrix, MA_number = 1){
  # Pmatrix = a matrix of MASA subgroup specific PC scores (output from princomp on original PC scores)
  # MA_number = the subgroup major axis number for calculation (defaul is 1, for the first major axis)
  pc.scores <- Pmatrix$scores

  temp_residuals <- apply(pc.scores[,-MA_number], 1, FUN = function(x) sum((x)^2)) ##Calculate the specimen specific residuals

  temp_SS_tot <- sum((pc.scores)^2)                 ##Calculate the multivariate sum of squared differences of PC scores from the mean PC scores
  temp_SS_res <- sum(temp_residuals)                ##Calculate the multivariate sum of squared differences of PC scores from the fitted value for that specimens' PC1 score
  temp_r.squared <- 1 - (temp_SS_res/temp_SS_tot)   ##Calculate r.squared by dividing top from bottom and subtracting this from 1

  out <- list(r.squared = temp_r.squared, residuals = temp_residuals)
  out
}

### Function to perform Pmatrix calculation from PC scores ###
major.axis.pca  <- function(pc.scores){

  pc.scores <- pc.scores
  N_comp <- min(c(nrow(pc.scores)-1), ncol(pc.scores))

  #uses PCA to identify a subgroup's major axes of variation (scores) and
  #calculate the contribution of each original PC axes towards each Major Axis (loadings)
  Pmatrix <- princomp(pc.scores[,1:N_comp])

  #calculate the coefficient of determination
  r.squared <- MA.r.squared(Pmatrix)

  #Calculates the miniumum and maximum values along each major axis among all specimens
  Min_List <- apply(Pmatrix$scores, 2, "min")
  Max_List <- apply(Pmatrix$scores, 2, "max")

  #constructs a matrix with extreme values for each major axis
  X <- matrix(0, byrow=FALSE, nrow = length(c(Min_List,Max_List)), ncol=N_comp)
  for(j in 1:length(Min_List)){
    X[j,j]  <- Min_List[j]
    X[N_comp+j,j]  <- Max_List[j]
  }
  major.axes <- X

  #Generate and assign unique column and r0w names for major.axes
  Pmax_Col_names <- gsub("PC", "MA", rownames(Pmatrix$loadings))
  Pmax_Row_names <- c(paste(Pmax_Col_names, "min", sep = "_"),paste(Pmax_Col_names, "max", sep = "_"))
  colnames(major.axes) <- Pmax_Col_names
  rownames(major.axes) <- Pmax_Row_names

  ##
  # Uses the loadings and center from the major axis PCA to transform Max
  # and Min values along each major axis into PC scores in the original
  # morphospace. This can then be used in MA_angular_comp to calculate angular
  # differences and perform post-hoc statistic tests among subgroups or used
  # to plot major axes when translated into the original morphospace.

  transformed.major.axes <- t(t(major.axes %*% t(Pmatrix$loadings)) + (Pmatrix$center))

  ##

  out <- list(Pmatrix = Pmatrix, major.axes = major.axes, transformed.major.axes = transformed.major.axes,
              loadings = Pmatrix$loadings, r.squared = r.squared$r.squared,
              N_comp = N_comp, pc.scores = pc.scores)
  class(out) <- "Pmatrix"

  out

  # print.Pmatrix_pca <- function(x, ...)
  # {
  #   cat("Call:\n"); dput(x$Pmatrix$call, control=NULL)
  #   cat("\nStandard deviations:\n")
  #   print(x$loadings)
  #   summary(x$Pmattix)
  #   print(x$r.squared, ...)
  #   cat("\n", length(x$Pmatrix$scale), " variables and ", x$Pmatrix$n.obs,
  #       "observations.\n")
  #   invisible(x)
  # }

  # out$r.squared
  # out$loadings
  # out$Pmatrix[c(1,N_comp+1),]

}
###

### Slope calculation for resampled major axes ###
resample.matrix <- function(resampled_MA, MA_number = 1){

  MA_number <- MA_number
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }

  MA_rows.list <- lapply(lapply(resampled_MA$transformed.major.axis, FUN = "rownames"), FUN = "grep", pattern = MA_name)

  resampled_lm_slopes <- matrix(nrow = length(resampled_MA), ncol = ncol(resampled_MA[[1]])-1,
                                dimnames = list(c(1:length(resampled_MA)),c(colnames(resampled_MA[[1]])[-1])))

  for(i in 1:length(resampled_MA)){
    current_LM <- lm(resampled_MA[[i]][,-1] ~ resampled_MA[[i]][,1])
    resampled_lm_slopes[i,] <- current_LM$coefficients[2,]
  }



  out <- list(resampled_lm_slopes = resampled_lm_slopes)
  class(out) <- "resampled_lm_slopes"

  out
}
###

### Slope calculation for resampled major axes ###
resample.lm <- function(resampled_MA,PC_comp){
  resampled_lm_slopes <- matrix(nrow = length(resampled_MA), ncol = ncol(resampled_MA[[1]])-1,
                                dimnames = list(c(1:length(resampled_MA)),c(colnames(resampled_MA[[1]])[-1])))
  resampled_lm_ints <- matrix(nrow = length(resampled_MA), ncol = ncol(resampled_MA[[1]])-1,
                                dimnames = list(c(1:length(resampled_MA)),c(colnames(resampled_MA[[1]])[-1])))
  for(i in 1:length(resampled_MA)){
    current_LM <- lm(resampled_MA[[i]][,-1] ~ resampled_MA[[i]][,1])
    resampled_lm_slopes[i,] <- current_LM$coefficients[2,]
    resampled_lm_ints[i,] <- current_LM$coefficients[1,]
  }
  out <- list(resampled_lm_slopes = resampled_lm_slopes,
              resampled_lm_ints = resampled_lm_ints)
  out
}
###

### Major Axis Calculation and Confidence Interval Function###
resample.major.axis <- function(X, PCs = c(1:4), MA_number = 1, method = c("bootstrap","jack-knife"), iter = 999, alpha = 0.05){

    # X = a formatted list of shape data and covariates (output by subesttingGMM)
    # MA_number = a single value to identify which axis of variation to study (default is 1, for the first major axis)
    # min_n = a single numberic value representing the minimum number of specimens required for a subgroup analysis (default is 4, arbitrarily)
    # method = method of resampling to be used by resampled.pcs; can be either bootstrap or jack-knife
    # iter = number of resampling interations used to generate variance distrubtion (default is 999)
    # alpha = statistical threshold (default if 0.05)

    # Specify groups
    groups <- X$taxa # grab the list of groups in the dataset
    groups_to_remove <- c() # create a list to fill with groups to drop due to undersampling

    # Specify method of resampling
    method <- match.arg(method, c("bootstrap","jack-knife"))

    # Specify confidence interval upper and lower bounds
    CI_lower <- alpha/2
    CI_upper <- 1 - (alpha/2)

    # create lists to capture values
    Major.Axis_list <- list()
    Loadings_list <- list()
    Transformed.MA_list <- list()
    R.squared_list <- list() # loadings of MA
    Importance_list <- list() # loadings of MA

    resampled_loadings <- list()
    resampled_transformed.ma_list <- list()

    loadings_CI <- list()

    n_comp <- length(groups)*iter

    min_n <- length(PCs) + 1

    message("Calculating Major Axes of Variation")
    pb = txtProgressBar(min = 0, max = n_comp, initial = 0, style = 3)
    step_n = 1

    for (i in 1:length(groups)){
      if (nrow(X$PCvalues[[i]]) < min_n) {
        warning(groups[i], " has too few individuals for resampling CI and was dropped from the analysis")
        groups_to_remove <- append(groups_to_remove, i)
        step_n <- step_n + iter
        next} else{
          group_major.axis <- major.axis.pca(X$PCvalues[[i]][,PCs]) #perform subgroup principal componenet rotation to identify major and minor axes
          group_transformed.major.axis <- group_major.axis$transformed.major.axes
          group_r.squared <- group_major.axis$r.squared

          group_resampled <- suppressMessages(resample.pcs(pc.scores=X$PCvalues[[i]][,PCs], method=method, iter=iter)) #resample current subgroup dataset

          resampled_loadings_matrix <- matrix(nrow=group_resampled$iter, ncol=ncol(group_resampled$original_PCs))
          resampled_L <- list()
          resampled_T <- list()

          for (j in 1:group_resampled$iter){
            setTxtProgressBar(pb,step_n)
            step_n <- step_n + 1

            current_major.axis <- major.axis.pca(group_resampled$resampled_PCs[[j]])               #perform major axis pca for current iteration
            current_loadings <- current_major.axis$loadings                                        #save current loadings

            resampled_L[[j]] <- current_loadings                                                   #save PC loadings of the major axis for current iteration
            resampled_loadings_matrix[j,1:nrow(current_loadings)] <- current_loadings[,MA_number]  #save loadings into matrix to caclulate loading CIs
            resampled_T[[j]] <- current_major.axis$transformed.major.axes                          #save transformed major axis for current iteration

          }

            #assign values to variables
            # Pmatrix_list[[X$taxa[i]]] <- group_major.axis$Pmatrix
            Major.Axis_list[[X$taxa[i]]] <- group_major.axis$major.axes
            Loadings_list[[X$taxa[i]]] <- group_major.axis$Pmatrix$loadings
            Transformed.MA_list[[X$taxa[i]]] <- group_transformed.major.axis
            R.squared_list[[X$taxa[i]]] <- group_r.squared
            Importance_list[[X$taxa[i]]] <- summary(group_major.axis)

            resampled_loadings[[X$taxa[i]]] <- resampled_L
            resampled_transformed.ma_list[[X$taxa[i]]] <- resampled_T


            loadings_CI[[X$taxa[i]]] <- apply(resampled_loadings_matrix, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE)

          }
    }
    if (length(groups_to_remove) > 0){
      groups <- groups[-groups_to_remove]
    }

    ##Perform pairwise comparisons of the major axes among subgroups##
    # MA_angular_comp()


    out <- list(groups = groups, iter = iter, MA_number = MA_number, method = method, alpha = alpha,
                Major.Axes = Major.Axis_list,
                Loadings = Loadings_list,
                Transformed.MA = Transformed.MA_list,
                R.squared = R.squared_list,
                Importance = Importance_list,
                resampled_loadings = resampled_loadings,
                resampled_transformed.MA = resampled_transformed.ma_list,
                loadings_CI = loadings_CI)
    class(out) = "ConfIntList"
    out

  }
###

### Major Axis Slope Calculation ###
major.axis.slopes <- function(resampled.major.axis){
  temp_slopes <- list()
  temp_ints <- list()
  temp_upper_preds <- list()
  temp_lower_preds <- list()
  for (i in 1:length(resampled.major.axis$groups)){
    temp_group <- resampled.major.axis$groups[[i]]
    temp_out <- major.axis.lm(temp_group,
                              resampled.major.axis$Transformed.MA[[i]],
                              resampled.major.axis$resampled_transformed.MA[[i]],
                              MA_number = 1, PC_comp = 1)
    temp_slopes[[temp_group]] <- temp_out$Slope_Table
    temp_ints[[temp_group]] <- temp_out$Intercepts
    temp_upper_preds[[temp_group]] <- temp_out$UpperPreds
    temp_lower_preds[[temp_group]] <- temp_out$LowerPreds
  }
  out <- list(Slopes = temp_slopes, Intercepts = temp_ints,
              Upper_CI_Values = temp_upper_preds, Lower_CI_Values = temp_lower_preds)
  out
}
###

###
plot.major.axis <- function(Slopes_obj,PCData,Axes_data,PCs){

Slopes_obj <- Slopes_obj
PCData <- PCData
Axes_data <- Axes_data
PCs <- PCs
PCn <- paste("PC", PCs[[2]], sep = "")

#create X and Y limits
Xlim<-c(floor(min(Axes_data$pc.scores[,PCs[[1]]])*10)/10,ceiling(max(Axes_data$pc.scores[,PCs[[1]]])*10)/10)
Ylim<-c(floor(min(Axes_data$pc.scores[,PCs[[2]]])*10)/10,ceiling(max(Axes_data$pc.scores[,PCs[[2]]])*10)/10)


for (i in 1:length(Slopes_obj)){
  temp_group <- names(Slopes_obj)[i]
  temp_PCs <- PCData[[temp_group]][,PCs]
  group_mean <- apply(temp_PCs,2,"mean")

  temp_obj <- Slopes_obj[[temp_group]]
  temp_slope <- temp_obj$Slope_Table[PCn,"Slope"]
  temp_int <- temp_obj$Intercepts[[PCn]]
  temp_upperpreds <- temp_obj$UpperPreds[[PCn]]
  temp_lowerpreds <- temp_obj$LowerPreds[[PCn]]
  temp_xpreds <- temp_obj$XPreds

  # temp_CI <- Slope_table[PCn,"95%CI"]
  # upper_CI <- temp_slope + temp_CI
  # lower_CI <- temp_slope - temp_CI
  #
  # temp_int <- group_mean[2] - c(temp_slope*group_mean[1])

  plot(0, 0, type = "n",
       xlim = Xlim,
       ylim = Ylim,
       xlab = paste("Principal Component ", PCs[[1]], " (", round(100*Axes_data$pc.summary$importance[2,PCs[[1]]], digits = 1), "%)", sep = ""),
       ylab = paste("Principal Component ", PCs[[2]], " (", round(100*Axes_data$pc.summary$importance[2,PCs[[2]]], digits = 1), "%)", sep = ""),
       axes = FALSE,
       frame.plot = FALSE,
       asp=F)

  axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
  axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
  clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
  abline(h=0, lty=3)
  abline(v=0, lty=3)

  mtext(temp_group, at = -0.25)

  polygon(c(rev(temp_xpreds), temp_xpreds),
          c(rev(temp_upperpreds), temp_lowerpreds),
          col = alpha('grey',0.6), border = NA)

  clip(
    min(temp_xpreds),
    max(temp_xpreds),
    min(temp_lowerpreds),
    max(temp_upperpreds)
  )

  lines(temp_xpreds, temp_upperpreds, lty = 'dashed')
  lines(temp_xpreds, temp_lowerpreds, lty = 'dashed')
  abline(temp_int,temp_slope, lwd=3, lty=1)
}
}
###

### Calculates single group major axis slope and confidence interval
major.axis.lm <- function(group, transformed.ma, resampled_transformed.ma, MA_number = 1, PC_comp = 1){
  # group = name of group
  # transformed.ma = transformed major axis matrix of group
  # resampled_transformed.ma = list of resampled transformed major axis matrices of group
  # MA_number = which axis to compare, default is the first majro axis
  # PC_comp = the reference PC to use for model caclulation, default is PC1

  ##Load in variables
  group <- group
  transformed.major.axis <- transformed.ma
  resampled_transformed.major.axis <- resampled_transformed.ma
  MA_number <- MA_number
  PC_comp <- PC_comp
  ##

  CI_upper <- 0.975
  CI_lower <- 0.025

  col_names <- c("Slope","95%CI") # column names for matrix
  PCn <- ncol(transformed.major.axis) # grab all PCs except for dependent variable (PC_comp)

  # ##Construct arrays to capture pairwise slope differences, p-values, and the merged results table for each PC
  # Slope_Table <- matrix(data = NA,
  #                       dim = c(length(PCn),length(col_names),),
  #                       dimnames = list(paste("PC",PCn,sep = ""),col_names)
  # )
  # ##

  # make a named variable for the major axis of interest, if the value entered is a simple numberic
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }

  ma_rows <-   grep(MA_name, rownames(transformed.major.axis)) # grab the rows associated with the selected major axis

  ## calculate group slope and confidence interval
  temp_ma <- transformed.major.axis[ma_rows,] # select correct rows from matrix
  temp_lm <- lm(temp_ma[,-PC_comp] ~ temp_ma[,PC_comp]) # create linear model of major axis, relative to PC_comp
  temp_slope <- temp_lm$coefficients[2,] # measured group slope
  temp_int <- temp_lm$coefficients[1,]

  temp_x <- seq(min(temp_ma[,PC_comp]),max(temp_ma[,PC_comp]), length.out = 50)
  temp_preds <- (temp_slope*temp_x) + temp_int


  temp_resampled_ma <- lapply(resampled_transformed.major.axis, function(x) x[ma_rows,]) # select correct rows from all resampled matrices
  temp_resampled_lm <- resample.lm(temp_resampled_ma,PC_comp)
  temp_resampled_slopes <- temp_resampled_lm$resampled_lm_slopes # perform lm on each resampled major axis
  temp_resampled_ints <- temp_resampled_lm$resampled_lm_ints # calculate the 95% upper and lower quantile values
  temp_resampled_pred <- list()
  for (i in 1:ncol(temp_resampled_slopes)){
    temp_resampled_pred[[colnames(temp_resampled_slopes)[[i]]]] <- (temp_resampled_slopes[,i]%*%t(temp_x)) + temp_resampled_ints[,i]
  }

  temp_resampled_pred_CIupper <- lapply(temp_resampled_pred, function(x) apply(x,2,quantile, probs = CI_upper))
  temp_resampled_pred_CIlower <- lapply(temp_resampled_pred, function(x) apply(x,2,quantile, probs = CI_lower))
  temp_resampled_quant <- apply(temp_resampled_slopes, 2, quantile, probs = c(CI_lower,CI_upper), na.rm = TRUE) # calculate the 95% upper and lower quantile values
  temp_resampled_CIs <- apply(apply(temp_resampled_quant, 1, function(x) abs(x - temp_slope)), 1, "mean") # calculate the mean 95% CI value around the measured slope
  ##

  Slope_Table <- cbind(temp_slope,temp_resampled_CIs)
  colnames(Slope_Table) <- col_names

  out <- list(Slope_Table = Slope_Table, Intercepts = temp_int, XPreds = temp_x,
              UpperPreds = temp_resampled_pred_CIupper, LowerPreds = temp_resampled_pred_CIlower)
}
###

### Function to statistically test pairwise angular comparisons between sets of group major axes##
major.axis.comparison <- function(groups, Transformed.MA_list, resampled_transformed.ma_list, MA_number = 1,  PCs = c(1:4), PC_comp = 1){
  # groups = a list of groups to be compared
  # transformed.major.axis =
  # resampled_transformed.major.axis =
  # MA_number =

  ##Load in variables
  groups <- groups
  comp_groups <- unique.pairs(groups)
  transformed.major.axis <- Transformed.MA_list
  resampled_transformed.major.axis <- resampled_transformed.ma_list
  n_comp <- length(unlist(comp_groups))

  PCn <- grep(PC_comp,PCs,invert=TRUE) # grab all PCs except for dependent variable (PC_comp)
  col_n <- gsub("PC", "", colnames(transformed.major.axis[[1]])) # get column names and extract PC numbers for matching

  MA_number <- MA_number
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }

  MA_rows.list <- lapply(lapply(transformed.major.axis, FUN = "rownames"), FUN = "grep", pattern = MA_name)

  ##Construct arrays to capture pairwise slope differences, p-values, and the merged results table for each PC
  Slope_Diff_Tables <- array(data = NA,
                             dim = c(length(groups),length(groups),length(PCn)),
                             dimnames = list(groups,groups,paste("PC",PCn,sep = ""))
  )
  Slope_P_Tables <- array(data = NA,
                          dim = c(length(groups),length(groups),length(PCn)),
                          dimnames = list(groups,groups,paste("PC",PCn,sep = ""))
  )
  Adj_Slope_P_Tables <- array(data = NA,
                          dim = c(length(groups),length(groups),length(PCn)),
                          dimnames = list(groups,groups,paste("PC",PCn,sep = ""))
  )
  Results_Tables <- array(data = NA,
                          dim = c(length(groups),length(groups),length(PCn)),
                          dimnames = list(groups,groups,paste("PC",PCn,sep = ""))
  )
  ##

  ##Create indeces for lower triagonal, upper triagonal, and diagonals
  results_lower <- lower.tri(Results_Tables[,,1])
  results_upper <- upper.tri(Results_Tables[,,1])
  results_diagonal <- as.logical(diag(nrow(Results_Tables[,,1])))
  #

  # ##Construct table to capture pairwise slope differences, p-values, and the merged results table summed across PCs
  # Results_Table <- matrix(data = NA , nrow = length(groups), ncol = length(groups), dimnames=(list(groups,groups)))
  # #

  # ##Save list of normality test of the resampled diffrences
  # normality_tests <- list()
  # ##

  ##Make progress bar object
  message("Calculating pairwise differences in slopes for MA ", MA_number)
  pb = txtProgressBar(min = 0, max = n_comp, style = 3)
  step_n = 1
  ##

  ##Look for making pairwise comparisons among the two sets of groups
  for (i in 1:length(comp_groups)){
    temp_group <- groups[[i]] # take the name of group
    temp_MA_rows <- MA_rows.list[[i]] #
    temp_MA <- transformed.major.axis[[temp_group]][temp_MA_rows,]
    temp_resampled_MA <- lapply(resampled_transformed.major.axis[[temp_group]], function(x) x[temp_MA_rows,])

    # my failure to not need a second for loop
    # temp_comp_groups <- groups[comp_groups[[i]]]
    # temp_comp_MA.list <- lapply(comp_groups[[i]], function(x) transformed.major.axis[[x]][MA_rows.list[[x]],])
    # temp_comp_resampled_MA.list <- lapply(comp_groups[[i]], function(x) resampled_transformed.major.axis[[x]][MA_rows.list[[x]],])

    # calculate group slope and within group permuted difference in slope
    temp_lm <- lm(temp_MA[,-PC_comp] ~ temp_MA[,PC_comp])
    temp_slope <- temp_lm$coefficients[2,] # measured group slope

    temp_resampled_slopes <- resample.lm(temp_resampled_MA,PC_comp)$resampled_lm_slopes
    temp_slope_dist <- temp_resampled_slopes[sample(nrow(temp_resampled_slopes)),] - temp_resampled_slopes # distribution of differences in slope derived from within group permutation
    #

    temp_comp_groups <- comp_groups[[i]]

    for (j in 1:length(temp_comp_groups)){
      setTxtProgressBar(pb,step_n)
      step_n <- step_n+1

      temp_comp_n <- temp_comp_groups[j]
      temp_comp_group <- groups[temp_comp_n]

      temp_comp_MA_rows <- MA_rows.list[[temp_comp_n]] #
      temp_comp_MA <- transformed.major.axis[[temp_comp_n]][temp_comp_MA_rows,]
      temp_comp_resampled_MA <- lapply(resampled_transformed.major.axis[[temp_comp_n]], function(x) x[temp_comp_MA_rows,])

      # calculate group slope and within group permuted difference in slope
      temp_comp_lm <- lm(temp_comp_MA[,-PC_comp] ~ temp_comp_MA[,PC_comp])
      temp_comp_slope <- temp_comp_lm$coefficients[2,] # group slope

      temp_comp_resampled_slopes <- resample.lm(temp_comp_resampled_MA,PC_comp)$resampled_lm_slopes
      temp_comp_slope_dist <- temp_comp_resampled_slopes[sample(nrow(temp_comp_resampled_slopes)),] - temp_comp_resampled_slopes # distribution of differences in slope derived from within group permutation
      #

      # generate pooled distribution of slope differences
      Pooled_slope_dist <- rbind(temp_slope_dist,temp_comp_slope_dist)

      # calculate difference in slopes between groups for each PC
      temp_slope_diff <- abs(temp_slope - temp_comp_slope)

      # determine the number of cases where the resambled distribution showed
      # a greater difference in slope than the measured difference, by PC
      Pooled_slope_sigs <- apply(Pooled_slope_dist, 1, function(x) abs(temp_slope_diff) < abs(x))

      # calculate the percentage of greater differences to estimate p-value
      temp_slope_P_val <- apply(Pooled_slope_sigs, 1, function(x) sum(x) / length(x))

      # # calculate the combined mean percentage of greater differences, to estimate p-value across all PCs
      # temp_combined_P_val <- mean(temp_slope_P_val)

      # add values to correct cells across arrays
      Slope_Diff_Tables[temp_comp_n,i,] <- temp_slope_diff
      Slope_P_Tables[i,temp_comp_n,] <- temp_slope_P_val
      Adj_Slope_P_Tables[i,temp_comp_n,] <- p.adjust(temp_slope_P_val, n=n_comp)
    }
  }

  Results_Tables[results_lower] <- Slope_Diff_Tables[results_lower]
  Results_Tables[results_upper] <- Adj_Slope_P_Tables[results_upper]
  Results_Tables[results_diagonal] <- 1
  # for(p in 1:dim(Results_Tables)[3]){diag(Results_Tables[,,p])<-1}

#
#   Results_Table[lower.tri(Results_Table)] <- Slope_Diff_Table[lower.tri(Results_Table)]
#   Results_Table[upper.tri(Results_Table)] <- p.adjust(Slope_P_Table[upper.tri(Results_Table)], n=n_comp)
#   diag(Results_Table) <- 1

  out <- list(Slopes = Slope_Diff_Tables,
              Pvalues = Slope_P_Tables,
              Results = Results_Tables)
  # normality_tests = normality_tests

  out

}
###

##Function to statistically test pairwise angular comparisons between sets of group major axes##
posthoc.major.axis.comparison <- function(MASA_result, Comp_MASA_result, group_list, MA_number = 1,  PCs = c(1:4), PC_comp = 1){
  # group_list = a list of groups to be compared, must match at least one group from both MASA result files
  # MASA_result =
  # Comp_MASA_result =
  # MA_number =
  # PCs =
  # PC_comp =

  ##Load in variables
  MASA_result <- MASA_result
  Comp_MASA_result <- Comp_MASA_result

  ##Set the two groups of MA to be compared
  groups <- MASA_result$groups
  Measured_MA <- MASA_result$Transformed.MA
  Resampled_MA <- MASA_result$resampled_transformed.MA
  Comp_groups <- Comp_MASA_result$groups
  Comp_Measured_MA <- Comp_MASA_result$Transformed.MA
  Comp_Resampled_MA <- Comp_MASA_result$resampled_transformed.MA
  ##

  ##subset groups to match group_list, if it is called
  if (!missing(group_list)){
    groups <- groups[groups %in% group_list]
    Comp_groups <- Comp_groups[Comp_groups %in% group_list]
  }
  ##


  # if (!(class(MASA_result) == "ConfIntList" & class(Comp_MASA_result) == "ConfIntList")){
  #   print("Object(s) must be confidence interval list(s).")
  # }

  PCn <- grep(PC_comp,PCs,invert=TRUE) # grab all PCs except for dependent variable (PC_comp)
  col_n <- gsub("PC", "", colnames(Measured_MA[[1]])) # get column names and extract PC numbers for matching
  n_comp <- length(groups)*length(Comp_groups)

  MA_number <- MA_number
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }

  MA_rows.list <- lapply(lapply(Measured_MA, FUN = "rownames"), FUN = "grep", pattern = MA_name)
  Comp_MA_rows.list <- lapply(lapply(Comp_Measured_MA, FUN = "rownames"), FUN = "grep", pattern = MA_name)

  # if(is.character(MA_number) == TRUE){
  # sub('.*-([0-9]+).*','\\1', MA_number)
  # }
  # #

  # ##Construct tables to capture angular differences and p-values for pairwise comparisons, and the merged results table
  # Angular_Diff_Table <- matrix(data = NA , nrow = length(groups), ncol = length(Comp_groups), dimnames=(list(groups,Comp_groups)))
  # Angular_P_Table <- matrix(data = NA , nrow = length(groups), ncol = length(Comp_groups), dimnames=(list(groups,Comp_groups)))
  # Results_Table <- matrix(data = NA , nrow = length(groups), ncol = length(Comp_groups), dimnames=(list(groups,Comp_groups)))
  # ##

  ##Construct arrays to capture pairwise slope differences, p-values, and the merged results table for each PC
  Slope_Diff_Tables <- array(data = NA,
                             dim = c(length(groups),length(Comp_groups),length(PCn)),
                             dimnames = list(groups,Comp_groups,paste("PC",PCn,sep = ""))
  )
  Slope_P_Tables <- array(data = NA,
                          dim = c(length(groups),length(Comp_groups),length(PCn)),
                          dimnames = list(groups,Comp_groups,paste("PC",PCn,sep = ""))
  )
  Comp_Slope_P_Tables <- array(data = NA,
                          dim = c(length(groups),length(Comp_groups),length(PCn)),
                          dimnames = list(groups,Comp_groups,paste("PC",PCn,sep = ""))
  )
  Adj_Slope_P_Tables <- array(data = NA,
                              dim = c(length(groups),length(Comp_groups),length(PCn)),
                              dimnames = list(groups,Comp_groups,paste("PC",PCn,sep = ""))
  )
  Results_Tables <- array(data = NA,
                          dim = c(length(groups),length(Comp_groups),length(PCn)),
                          dimnames = list(groups,Comp_groups,paste("PC",PCn,sep = ""))
  )
  ##

  ##Create indeces for lower triagonal, upper triagonal, and diagonals
  results_lower <- lower.tri(Results_Tables[,,1])
  results_upper <- upper.tri(Results_Tables[,,1])
  results_diagonal <- as.logical(diag(nrow(Results_Tables[,,1])))
  ##

  # ##Save list of normality test of the resampled diffrences
  # normality_tests <- list()
  # ##

  ##Make progress bar object
  message("Calculating pairwise differences in slopes for MA ", MA_number)
  pb = txtProgressBar(min = 0, max = n_comp, style = 3)
  step_n = 1
  ##

  ##Look for making pairwise comparisons among the two sets of groups
  for (i in 1:length(groups)){
    temp_group <- groups[[i]] # take the name of group
    temp_MA_rows <- MA_rows.list[[i]] #
    temp_MA <- Measured_MA[[temp_group]][temp_MA_rows,]
    temp_resampled_MA <- lapply(Resampled_MA[[temp_group]], function(x) x[temp_MA_rows,])

    # calculate group slope and within group permuted difference in slope
    temp_lm <- lm(temp_MA[,-PC_comp] ~ temp_MA[,PC_comp])
    temp_slope <- temp_lm$coefficients[2,] # measured group slope

    temp_resampled_slopes <- resample.lm(temp_resampled_MA,PC_comp)$resampled_lm_slopes
    temp_slope_dist <- temp_resampled_slopes[sample(nrow(temp_resampled_slopes)),] - temp_resampled_slopes # distribution of differences in slope derived from within group permutation
    #

    for (j in 1:length(Comp_groups)){
      setTxtProgressBar(pb,step_n)
      step_n <- step_n+1

      temp_comp_group <- Comp_groups[[j]]
      temp_comp_MA_rows <- Comp_MA_rows.list[[temp_comp_group]] #
      temp_comp_MA <- Comp_Measured_MA[[temp_comp_group]][temp_comp_MA_rows,]
      temp_comp_resampled_MA <- lapply(Comp_Resampled_MA[[temp_comp_group]], function(x) x[temp_comp_MA_rows,])

      # calculate group slope and within group permuted difference in slope
      temp_comp_lm <- lm(temp_comp_MA[,-PC_comp] ~ temp_comp_MA[,PC_comp])
      temp_comp_slope <- temp_comp_lm$coefficients[2,] # group slope

      temp_comp_resampled_slopes <- resample.lm(temp_comp_resampled_MA,PC_comp)$resampled_lm_slopes
      temp_comp_slope_dist <- temp_comp_resampled_slopes[sample(nrow(temp_comp_resampled_slopes)),] - temp_comp_resampled_slopes # distribution of differences in slope derived from within group permutation
      #

      # generate pooled distribution of slope differences
      Pooled_slope_dist <- rbind(temp_slope_dist,temp_comp_slope_dist)

      # calculate difference in slopes between groups for each PC
      temp_slope_diff <- abs(temp_slope - temp_comp_slope)

      # determine the number of cases where the resambled distribution showed
      # a greater difference in slope than the measured difference, by PC
      Pooled_slope_sigs <- apply(Pooled_slope_dist, 1, function(x) abs(temp_slope_diff) < abs(x))

      # calculate the percentage of greater differences to estimate p-value
      temp_slope_P_val <- apply(Pooled_slope_sigs, 1, function(x) sum(x) / length(x))

      # determine the number of cases where the resambled distribution showed
      # a greater difference in slope than the measured difference, by PC
      comp_slope_sigs <- apply(temp_comp_slope_dist, 1, function(x) abs(temp_slope_diff) < abs(x))

      # calculate the percentage of greater differences to estimate p-value
      temp_comp_slope_P_val <- apply(comp_slope_sigs, 1, function(x) sum(x) / length(x))

      # # calculate the combined mean percentage of greater differences, to estimate p-value across all PCs
      # temp_combined_P_val <- mean(temp_slope_P_val)

      # add values to correct cells across arrays
      Slope_Diff_Tables[i,j,] <- temp_slope_diff
      Slope_P_Tables[i,j,] <- temp_slope_P_val
      Comp_Slope_P_Tables[i,j,] <- temp_comp_slope_P_val
      Adj_Slope_P_Tables[i,j,] <- p.adjust(temp_slope_P_val, n=n_comp)
    }
  }

  # for(l in 1:dim(Results_Tables)[3]){
  #   Results_Tables[,,l] <- cbind(Slope_Diff_Tables[,,l],Slope_P_Tables[,,l])
  # }

  # Results_Table[lower.tri(Results_Table)] <- Angular_Diff_Table[lower.tri(Results_Table)]
  # Results_Table[upper.tri(Results_Table)] <- p.adjust(Angular_P_Table[upper.tri(Results_Table)], n=n_comp)
  # diag(Results_Table) <- 1

  out <- list(Slopes = Slope_Diff_Tables,
              Pvalues = Slope_P_Tables,
              Comp_Pvalues = Comp_Slope_P_Tables,
              Adj_Pvalues = Adj_Slope_P_Tables,
              Results = Results_Tables)
  # normality_tests = normality_tests

  out

}
##

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


Major.Axis <- function(X, PCs = c(1:4), PC_comp = 1, MA_number = 1, method = c("bootstrap","jack-knife"), iter = 999, alpha = 0.05){

  resampled_ma <- resample.major.axis(X, PCs, MA_number, method, iter, alpha)

  ma_slopes <- list()
  for (i in 1:length(resampled_ma$groups)){
    temp_group <- resampled_ma$groups[i]
    ma_slopes[[temp_group]] <- major.axis.lm(temp_group, resampled_ma$Transformed.MA[i], resampled_ma$resampled_transformed.ma[i], MA_number, PC_comp)
  }

  ma_results <- major.axis.comparison(resampled_ma$groups, resampled_ma$Transformed.MA, resampled_ma$resampled_transformed.ma, MA_number, PCs, PC_comp)

  out <- list(groups <- groups, PCs <- PCs, axis <- MA_number,
              Slopes <- ma_slopes,
              R.squared <- resampled_ma$R.squared,
              Loadings <- resampled_ma$Loadings,
              Comparisons <- ma_results$Results
              )

  out
}


## old angular MA stuff
###
# # Creat vector of the change in each PC score along MA in question.
# #
# Y_MA <- temp_comp_MA[2,] - temp_comp_MA[1,]
# resampled_Y_MA <- t(sapply(temp_comp_resampled_MA, function(x) x[2,]-x[1,]))
# ###
#
# temp_measured_angle <- angle(t(X_MA),Y_MA)
# temp_resampled_angles <- c()

# ##code for intergroup comparions, with subtracted mean##
# for (k in 1:nrow(resampled_Y_MA)){
#   temp_X_MA <- resampled_X_MA[k,]
#   temp_Y_MA <- resampled_Y_MA[k,]
#
#   temp_angle <- angle(temp_X_MA,temp_Y_MA)
#
#   temp_resampled_angles <- append(temp_resampled_angles,temp_angle)
# }
# temp_angular_var <- temp_resampled_angles - mean(temp_resampled_angles)
# ##
# #code for intragroup comparions without mean subtraction##
# shuffled_X_MA <- resampled_X_MA[sample(nrow(resampled_X_MA)),]
# shuffled_Y_MA <- resampled_Y_MA[sample(nrow(resampled_Y_MA)),]
# for (k in 1:nrow(resampled_Y_MA)){
#   temp_X_MA <- resampled_X_MA[k,]
#   temp_shuffled_X_MA <- shuffled_X_MA[k,]
#   temp_Y_MA <- resampled_Y_MA[k,]
#   temp_shuffled_Y_MA <- shuffled_Y_MA[k,]
#
#   temp_X_angle <- angle(t(temp_X_MA),temp_shuffled_X_MA)
#   temp_Y_angle <- angle(t(temp_Y_MA),temp_shuffled_Y_MA)
#
#   temp_resampled_angles <- append(temp_resampled_angles,temp_X_angle,temp_Y_angle)
# }
# temp_angular_var <- temp_resampled_angles
# ##

# normality_tests[paste(temp_group, "&", temp_comp_group, "variance of angular differences", sep = " ")] <- shapiro.test(temp_angular_var)[[2]]




