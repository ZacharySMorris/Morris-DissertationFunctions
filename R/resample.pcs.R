###Single function to resample PCA dataset###
resample.pcs <- function(pc.scores, method = c("bootstrap","jack-knife"), iter=999){
  # pc.scores = a matrix of PC scores to be resampled (ouput from princomp)
  # method = the method of resampling to be used (either bootstrap or jack-knife replacement)
  # iter = the number of resampled datasets to create

  group_to_resample <- pc.scores
  specimen_length <- c(1:dim(group_to_resample)[1])

  if (method == "bootstrap"){ToReplace = TRUE}
  if (method == "jack-knife"){ToReplace = FALSE}

  resamples <- list()
  resampled_PCs <- list()
  for (i in 1:iter){
    repeat{
      specimen_numbers <- c(sample(specimen_length, replace = ToReplace))
      if (length(unique(specimen_numbers)) > 3) break
    }
    resamples[[i]] <- specimen_numbers
    resampled_PCs[[i]] <- group_to_resample[resamples[[i]],]
  }

  out <- list(resampled_specimen = resamples, resampled_PCs = resampled_PCs,
              iter = iter, original_PCs = group_to_resample)
  class(out) <- "ResampledPCAList"
  message("Generated ", iter, " resampled versions of dataset via ", method, " replacement.", appendLF=TRUE)
  out
}
###

