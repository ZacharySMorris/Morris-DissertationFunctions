#Rotate matrix

rotate.LM <- function(M,LM1,LM2,plotLM){

  #Specify matrix of GPA coordinates
  M <- M
  #shift points, so that turning point is (0,0)
  M2.1 <- t(t(M) - c(M[LM1,1],M[LM1,2]))
  #calculate rotation angle
  alpha <- -atan(M2.1[LM2,2]/M2.1[LM2,1])
  #rotation matrix
  rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
  #rotate
  M2.2 <- t(rotm %*% (t(M2.1)))
  #shift back
  M2.3 <- t(t(M2.2)+c(M[LM1,1],M[LM1,2]))


  if (plotLM == TRUE){
    #plot data
    plot(M,xlim=c(-1,1),ylim=c(-1,1))
    points(M2.1,col="blue")
    points(M2.3,col="red")
    points(M2.2,col="green")
  }

  out = M2.3
}
