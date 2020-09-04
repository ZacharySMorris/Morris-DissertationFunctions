## Functions to easily transform between GPA coordinates and PC scores given a rotation matrix ##

## Coords2PC takes a list of procrustes aligned coordinates and projects them into a target PC space ##
Coords2PC<-function(coords,pca,gpa){
  #coords = a set of GPA coordinates (2D or 3D) to be translated into principal components
  #pca = the center and rotation matrix for the target principal components (output from princomp)
  #gpa = the dimentions of the procrustes aligned dataset (output from gpagen)

  k <- dim(gpa$coords)[2]
  p <- dim(gpa$coords)[1]
  n <- dim(gpa$coords)[3]

  Mshape <- pca$pc.summary$center
  Rotation <- pca$rotation

  x <- array(data= NA,dim=c(p,k,1))

  if (dim(coords)[2]==k){
    if (length(dim(coords))==2) {

      dim(coords)[3]<-1
      x <- two.d.array(coords)

    }

    if (length(dim(coords))==3) {

      x <- two.d.array(coords)

    }
  }
  if (dim(coords)[2]>k) x <- coords

  ##Remove the PCA center, but it doesn't really change anything##
  z <- sweep(x, 2, Mshape)
  w <- z %*% Rotation
  return(w)
}
##
## PC2Coords takes a set of principal compoenents scores and calculates procrustes aligned coordinates ##
PC2Coords<-function(pc,pca,gpa){
  #coords = a set of PC scores to be translated into procrustes aligned coordinates
  #pca = the center and rotation matrix for the principal components of interest (output from princomp)
  #gpa = the dimentions of the procrustes aligned dataset (output from gpagen)

  Mshape <- mshape(gpa$coords)
  Rotation <- pca$rotation
  k <- dim(gpa$coords)[2]
  p <- dim(gpa$coords)[1]
  n <- dim(gpa$coords)[3]

  w <- as.matrix( (pc %*%  t(Rotation)) + pca$center)
  w <- arrayspecs(w,p,k)
  w[1:3,1:3,1]
  coords[1:3,1:3,1]

  return(w)
}
##
