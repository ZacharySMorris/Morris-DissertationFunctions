###Script for using covariates to create values for plotting PCAs###

PlottingValues <- function(X,ColorGroup,ShapeGroup)
{
  set.seed(001)
  #Creates a vector with the species colors for each specimens for plotting
  gp <- as.factor(X[,ColorGroup])
  col.gp <- NULL
  col.gp <- colorRampPalette(brewer.pal(length(levels(gp)),"Paired"))(length(levels(X[,ColorGroup])))
  names(col.gp) <- sample(levels(gp)) ## assign group labels to colours
  col.gp <- col.gp[match(gp, names(col.gp))]

  #Creating vector of colors for the legend
  legend.col <- NULL
  legend.col <- unique(col.gp)
  names(legend.col) <- unique(gp)
  legend.col <- legend.col[match(unique(names(col.gp)), names(legend.col))] #Not needed?

  #Creates a vector with the size of specimens for plot
  age.gp <- NULL
  age.gp<-X$Size
  names(age.gp) <- X$Species

  #Creates a vector with the genus shape for each specimens for plotting
  gp <- as.factor(X[,ShapeGroup])
  shape.gp <- NULL
  shape.gp<-rep(25,length(X[[ShapeGroup]])) #again make vector length number specimens
  shape.gp[which(X[[ShapeGroup]] == levels(X[[ShapeGroup]])[[1]])] <- 21
  shape.gp[which(X[[ShapeGroup]] == levels(X[[ShapeGroup]])[2])] <- 22
  shape.gp[which(X[[ShapeGroup]] == levels(X[[ShapeGroup]])[3])] <- 23
  # shape.gp[which(X$Genus=='Tomistoma')]<-24

  names(shape.gp) <- gp ## assign group labels to colors
  shape.gp <- shape.gp[match(gp, names(shape.gp))]

  #Creating vector of shapes for the legend
  #Almost works just right, a few things matched up wrong
  legend.shape <- NULL
  legend.shape <- shape.gp[match(unique(names(shape.gp)), names(shape.gp))]
  legend.shape <- as.numeric(gsub(21, 19, legend.shape))
  legend.shape <- as.numeric(gsub(22, 15, legend.shape))
  legend.shape <- as.numeric(gsub(23, 18, legend.shape))
  legend.shape <- as.numeric(gsub(24, 17, legend.shape))
  names(legend.shape) <- unique(gp)

  #Assign values#
  assign("col.gp", col.gp, envir = .GlobalEnv)
  assign("legend.col", legend.col, envir = .GlobalEnv)
  assign("age.gp", age.gp, envir = .GlobalEnv)
  assign("shape.gp", shape.gp, envir = .GlobalEnv)
  assign("legend.shape", legend.shape, envir = .GlobalEnv)

  out <- list(color = col.gp, legendcolor = legend.col, size = age.gp, shape = shape.gp, legendshape = legend.shape)

  out

}

Simple_PVs <- function(X,covariate,values,name,Xappend=FALSE){
  #X is a covariate dataset, should be a data.frame
  #covariate is the name of a specific covariate from the data frame that should be used to create a plotting value
  #values is a list of values to replace the factors found in the data frame
  #name is the name of the new plotting value to be saved
  #Xappend is a logical value (T/F) of whether the output should be output as a single vector or appended to the original covariate dataset

  tmp_X <- X
  tmp_cov <- as.factor(tmp_X[[covariate]])
  tmp_values <- values
  values_class <- class(tmp_values)
  tmp_names <- append(names(tmp_X),name)

  if (!length(levels(tmp_cov)) == length(tmp_values)){
    break
    warning("The number of values supplied does not match the number of factors to be replaced.")
  }else{
    tmp_ch <- as.character(tmp_cov)
    for (i in 1:length(levels(tmp_cov))){
      tmp_level <- as.character(levels(tmp_cov)[[i]])
      tmp_ch[tmp_ch == tmp_level] <- tmp_values[[i]]
      out_values <- tmp_ch
    }
  }

  if (Xappend == TRUE){
    X_out <- cbind(tmp_X,tmp_ch)
    names(X_out) <- tmp_names
  }else{
    # X_out <- tmp_ch
    # names(X_out) <- name
    X_out <- data.frame(tmp_ch)
    colnames(X_out) <- name
  }

  out <- X_out
  out

}


