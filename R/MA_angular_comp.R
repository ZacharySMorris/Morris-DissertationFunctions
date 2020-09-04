###Code for functions to calculate angular differences between multidimentional axes###

##Basic angle calculation function##
angle <- function(x,y){
  # x = first vector
  # y = second vector
  # dot.prod <- x%*%y
  # norm.x <- norm(x,type="2")
  # norm.y <- norm(y,type="2")
  # theta <- acos(dot.prod / (norm.x * norm.y))
  # as.numeric(theta)
  dot.prod <- as.numeric(x%*%y )
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  atheta<-dot.prod / (norm.x * norm.y)
  theta <- acos(pmin(pmax(atheta,-1.0),1.0))
  as.numeric(theta)

}
##

unique.pairs <- function(v){
  # v = a vector of group names to be paired

  n <- length(v)
  n_v <- c(1:n)

  temp.v <- lapply(n_v, FUN = function(x) (x+1):n)
  unique.v <- temp.v[-n]
}


##Function to statistically test pairwise angular comparisons between sets of group major axes##
MA_angular_comp <- function(groups, Transformed.MA_list, resampled_transformed.ma_list, MA_number = 1){
  # groups =
  # transformed.major.axis =
  # resampled_transformed.major.axis =
  # MA_number =

  ##Load in variables
  groups <- groups
  comp_groups <- unique.pairs(groups)
  transformed.major.axis <- Transformed.MA_list
  resampled_transformed.major.axis <- resampled_transformed.ma_list
  n_comp <- length(unlist(comp_groups))

  MA_number <- MA_number
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }

  MA_rows.list <- lapply(lapply(transformed.major.axis, FUN = "rownames"), FUN = "grep", pattern = MA_name)

  ##Construct tables to capture angular differences and p-values for pairwise comparisons, and the merged results table
  Angular_Diff_Table <- matrix(data = NA , nrow = length(groups), ncol = length(groups), dimnames=(list(groups,groups)))
  Angular_P_Table <- matrix(data = NA , nrow = length(groups), ncol = length(groups), dimnames=(list(groups,groups)))
  Results_Table <- matrix(data = NA , nrow = length(groups), ncol = length(groups), dimnames=(list(groups,groups)))
  ##

  ##Create indeces for lower triagonal, upper triagonal, and diagonals
  results_lower <- lower.tri(Results_Table)
  results_upper <- upper.tri(Results_Table)
  results_diagonal <- as.logical(diag(nrow(Results_Table)))
  #

  ##Save list of normality test of the resampled diffrences
  normality_tests <- list()
  ##

  ##Make progress bar object
  message("Calculating pairwise angular differences among major axes", MA_number)
  pb = txtProgressBar(min = 0, max = n_comp, style = 3)
  step_n = 1
  ##

  ##Look for making pairwise comparisons among the two sets of groups
  for (i in 1:length(comp_groups)){
    temp_group <- groups[[i]]
    temp_MA_rows <- MA_rows.list[[i]]
    temp_MA <- transformed.major.axis[[temp_group]][temp_MA_rows,]
    temp_resampled_MA <- lapply(resampled_transformed.major.axis[[temp_group]], function(x) x[temp_MA_rows,])

    # temp_comp_groups <- comp_groups[[i]]
    # temp_comp_MA.list <- lapply(temp_comp_groups, function(x) transformed.major.axis[[x]][MA_rows.list[[x]],])
    # temp_comp_resampled_MA.list <- t(sapply(temp_resampled_MA, function(x) x[2,]-x[1,]))

    ###
    # To calculate the angular difference we want all MA centered.
    # Thus, we want to subtract the minimum from the maximum point in PC space.
    # This gives us a single vector of PC scores reflecting the change in each PC score
    # that happens along the MA in question. We will do this for each resampled MA too.
    #
    X_MA <- temp_MA[2,] - temp_MA[1,]
    resampled_X_MA <- t(sapply(temp_resampled_MA, function(x) x[2,]-x[1,]))
#
#     Y_MA.list <- t(sapply(temp_comp_MA.list, function(x) x[2,]-x[1,]))
#     resampled_Y_MA.list <- lapply(comp_groups[[i]], function(x) temp_comp_resampled_MA.list[[x]][MA_rows.list[[x]],])
#     ###

    temp_comp_groups <- comp_groups[[i]]

    for (j in 1:length(temp_comp_groups)){
      setTxtProgressBar(pb,step_n)
      step_n <- step_n+1

      temp_comp_n <- temp_comp_groups[j]
      temp_comp_group <- groups[temp_comp_n]

      temp_comp_MA_rows <- MA_rows.list[[temp_comp_n]] #
      temp_comp_MA <- transformed.major.axis[[temp_comp_n]][temp_comp_MA_rows,]
      temp_comp_resampled_MA <- lapply(resampled_transformed.major.axis[[temp_comp_n]], function(x) x[temp_comp_MA_rows,])
      ###
      # Creat vector of the change in each PC score along MA in question.
      #
      Y_MA <- temp_comp_MA[2,] - temp_comp_MA[1,]
      resampled_Y_MA <- t(sapply(temp_comp_resampled_MA, function(x) x[2,]-x[1,]))
      ###

      temp_measured_angle <- angle(t(X_MA),Y_MA)
      temp_resampled_angles <- c()

      ##code for intergroup comparions, with subtracted mean##
      for (k in 1:nrow(resampled_Y_MA)){
        temp_X_MA <- resampled_X_MA[k,]
        temp_Y_MA <- resampled_Y_MA[k,]

        temp_angle <- angle(temp_X_MA,temp_Y_MA)

        temp_resampled_angles <- append(temp_resampled_angles,temp_angle)
      }
      temp_angular_var <- temp_resampled_angles - mean(temp_resampled_angles)
      ##
      #code for intragroup comparions without mean subtraction##
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

      temp_angular_P_val <- sum(temp_measured_angle < abs(temp_angular_var)) / length(temp_angular_var)

      Angular_Diff_Table[temp_comp_n,i] <- temp_measured_angle
      Angular_P_Table[i,temp_comp_n] <- temp_angular_P_val

      # # add values to correct cells across arrays
      # Angular_Diff_Table[temp_comp_n,i,] <- temp_measured_angle
      # Angular_P_Table[i,temp_comp_n,] <- temp_angular_P_val
      # Adj_Slope_P_Tables[i,temp_comp_n,] <- p.adjust(temp_angular_P_val, n=n_comp)
    }
  }

  Results_Table[lower.tri(Results_Table)] <- Angular_Diff_Table[lower.tri(Results_Table)]
  Results_Table[upper.tri(Results_Table)] <- p.adjust(Angular_P_Table[upper.tri(Results_Table)], n=n_comp)
  diag(Results_Table) <- 1

  out <- list(Angles =Angular_Diff_Table,
              P_values = Angular_P_Table,
              Results = Results_Table)
              # normality_tests = normality_tests

  out

}
##

##Function to statistically test pairwise angular comparisons between sets of group major axes##
posthoc_MASA_comp <- function(MASA_result, Comp_MASA_result, group_list, MA_number = 1){
  # ConfIntList =
  ##Load in variabless
  MASA_result <- MASA_result
  Comp_MASA_result <- Comp_MASA_result


  # if (!(class(MASA_result) == "ConfIntList" & class(Comp_MASA_result) == "ConfIntList")){
  #   print("Object(s) must be confidence interval list(s).")
  # }

  MA_number <- MA_number
  if(is.numeric(MA_number) == TRUE){
    MA_name <- paste("MA", MA_number, "_", sep = "")
  }
  # if(is.character(MA_number) == TRUE){
  # sub('.*-([0-9]+).*','\\1', MA_number)
  # }
  # #


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

  n_comp <- length(groups)*length(Comp_groups)

  MA_rows.list <- lapply(lapply(Measured_MA, FUN = "rownames"), FUN = "grep", pattern = MA_name)
  Comp_MA_rows.list <- lapply(lapply(Comp_Measured_MA, FUN = "rownames"), FUN = "grep", pattern = MA_name)


  ##Construct tables to capture angular differences and p-values for pairwise comparisons, and the merged results table
  Angular_Diff_Table <- matrix(data = NA , nrow = length(groups), ncol = length(Comp_groups), dimnames=(list(groups,Comp_groups)))
  Angular_P_Table <- matrix(data = NA , nrow = length(groups), ncol = length(Comp_groups), dimnames=(list(groups,Comp_groups)))
  # Results_Table <- matrix(data = NA , nrow = length(groups), ncol = 2*length(Comp_groups),
  #                         dimnames=(list(groups,c(paste("Angle with ", Comp_groups, sep = ""),paste("p-value with", Comp_groups, sep = "")))))
  ##

  ##Save list of normality test of the resampled diffrences
  normality_tests <- list()
  ##

  ##Make progress bar object
  message("Calculating pairwise angular differences among major axes", MA_number)
  pb = txtProgressBar(min = 0, max = n_comp, initial = 0, style = 3)
  step_n = 1
  ##Look for making pairwise comparisons among the two sets of groups
  for (i in 1:length(groups)){
    temp_group <- groups[[i]]
    temp_MA_rows <- MA_rows.list[[i]] #
    temp_MA <- Measured_MA[[temp_group]][temp_MA_rows,]
    temp_resampled_MA <- lapply(Resampled_MA[[temp_group]], function(x) x[temp_MA_rows,])

    ###
    # To calculate the angular difference we want all MA centered.
    # Thus, we want to subtract the minimum from the maximum point in PC space.
    # This gives us a single vector of PC scores reflecting the change in each PC score
    # that happens along the MA in question. We will do this for each resampled MA too.
    #
    X_MA <- temp_MA[2,] - temp_MA[1,]
    resampled_X_MA <- t(sapply(temp_resampled_MA, function(x) x[2,]-x[1,]))
    ###

    for (j in 1:length(Comp_groups)){
      setTxtProgressBar(pb,step_n)
      step_n <- step_n+1

      temp_comp_group <- Comp_groups[[j]]
      temp_comp_MA_rows <- Comp_MA_rows.list[[temp_comp_group]] #
      temp_comp_MA <- Comp_Measured_MA[[temp_comp_group]][temp_comp_MA_rows,]
      temp_comp_resampled_MA <- lapply(Comp_Resampled_MA[[temp_comp_group]], function(x) x[temp_comp_MA_rows,])

      ###
      # Creat vector of the change in each PC score along MA in question.
      #
      Y_MA <- temp_comp_MA[2,] - temp_comp_MA[1,]
      resampled_Y_MA <- t(sapply(temp_comp_resampled_MA, function(x) x[2,]-x[1,]))
      ###

      comp_PC_n <- min(length(X_MA), length(Y_MA))

      temp_measured_angle <- angle(t(X_MA[1:comp_PC_n]),Y_MA[1:comp_PC_n])
      temp_resampled_angles <- c()

      ##code for intergroup comparions, with subtracted mean##
      # for (k in 1:nrow(resampled_Y_MA)){
      #   temp_X_MA <- resampled_X_MA[k,]
      #   temp_Y_MA <- resampled_Y_MA[k,]
      #
      #   temp_angle <- angle(temp_X_MA,temp_Y_MA)
      #
      #   temp_resampled_angles <- append(temp_resampled_angles,temp_angle)
      # }
      # temp_angular_var <- temp_resampled_angles - mean(temp_resampled_angles)
      ##
      #code for intragroup comparions without mean subtraction##
      shuffled_X_MA <- resampled_X_MA[sample(nrow(resampled_X_MA)),]
      shuffled_Y_MA <- resampled_Y_MA[sample(nrow(resampled_Y_MA)),]
      for (k in 1:nrow(resampled_Y_MA)){
        temp_X_MA <- resampled_X_MA[k,]
        temp_shuffled_X_MA <- shuffled_X_MA[k,]
        temp_Y_MA <- resampled_Y_MA[k,]
        temp_shuffled_Y_MA <- shuffled_Y_MA[k,]

        comp_X_PC_n <- min(length(temp_X_MA), length(temp_shuffled_X_MA))
        comp_Y_PC_n <- min(length(temp_Y_MA), length(temp_shuffled_Y_MA))

        temp_X_angle <- angle(t(temp_X_MA[1:comp_X_PC_n]),temp_shuffled_X_MA[1:comp_X_PC_n])
        temp_Y_angle <- angle(t(temp_Y_MA[1:comp_Y_PC_n]),temp_shuffled_Y_MA[1:comp_Y_PC_n])

        temp_resampled_angles <- append(temp_resampled_angles,temp_X_angle,temp_Y_angle)
      }
      temp_angular_var <- temp_resampled_angles
      # ##

      # normality_tests[paste(temp_group, "&", temp_comp_group, "variance of angular differences", sep = " ")] <- shapiro.test(temp_angular_var)[[2]]

      temp_angular_P_val <- sum(temp_measured_angle < abs(temp_angular_var)) / length(temp_angular_var)

      Angular_Diff_Table[i,j] <- temp_measured_angle
      Angular_P_Table[i,j] <- temp_angular_P_val
    }
  }
  Adjusted_P_Table <- Angular_P_Table
  Adjusted_P_Table[1:length(groups),] <- p.adjust(Angular_P_Table)

  # Results_Table[lower.tri(Results_Table)] <- Angular_Diff_Table[lower.tri(Results_Table)]
  # Results_Table[upper.tri(Results_Table)] <- p.adjust(Angular_P_Table[upper.tri(Results_Table)], n=n_comp)
  # diag(Results_Table) <- 1

  out <- list(Angular_Differences = Angular_Diff_Table,
              P_values = Angular_P_Table,
              Adj_P_values = Adjusted_P_Table)
  # normality_tests = normality_tests

  out

}
##
