## All the helper functions for the longitudinal PCA plots

## Three Mode plots ##

## PCA accounting for correlation structures

# this function will calculate the scores and loadings for Joint biplots; based off the source code from the ThreeWay package on Cran
joint<-function (K, A, B, C, fixmode, fixunit){
  r1 = ncol(A)
  n = nrow(A)
  r2 = ncol(B)
  m = nrow(B)
  r3 = ncol(C)
  p = nrow(C)
  if (fixmode == 1) {
    Gmat = matrix(t(K[fixunit, ]), r2, r3)
    Smat = B %*% Gmat %*% t(C)
    SVD = svd(Smat)
    Bmat = (m/p)^0.25 * SVD$u[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    Cmat = (p/m)^0.25 * SVD$v[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    xmax = max(Bmat, Cmat) + 0.1
    xmin = min(Bmat, Cmat) - 0.1
    list=list(Bmat,Cmat,SVD$d)
    names(list)=c("Bmat","Cmat","eigen")
    list=lapply(list,cbind.data.frame)
  }
  if (fixmode == 2) {
    K = ThreeWay::permnew(K, r1, r2, r3)
    Gmat = matrix(t(K[fixunit, ]), r3, r1)
    Smat = C %*% Gmat %*% t(A)
    SVD = svd(Smat)
    Cmat = (p/n)^0.25 * SVD$u[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    Amat = (n/p)^0.25 * SVD$v[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    xmax = max(Cmat, Amat) + 0.1
    xmin = min(Cmat, Amat) - 0.1
    list=list(Cmat,Amat,SVD$d)
    list=lapply(list,cbind.data.frame)
    names(list)=c("Cmat","Amat","eigen")
  }
  if (fixmode == 3) {
    K = ThreeWay::permnew(K, r2, r3, r1)
    Gmat = matrix(t(K[fixunit, ]), r1, r2)
    Smat = A %*% Gmat %*% t(B)
    SVD = svd(Smat)
    Amat = (p/n)^0.25 * SVD$u[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    Bmat = (n/p)^0.25 * SVD$v[, c(1,2)] %*% diag(SVD$d[c(1,2)])^0.5
    list=list(Amat,Bmat,SVD$d)
    names(list)=c("Amat","Bmat","eigen")
    list=lapply(list,cbind.data.frame)
  }
  list
}

# this will create a joint biplot in ggplot
jointplot<-function(data,x,y,data_seg,x_seg,y_seg,ylab,xlab,point_lab,vec_lab){

  x_min=min(c(sapply(x,min),sapply(x_seg,min)))
  x_max=max(c(sapply(x,max),sapply(x_seg,max)))
  y_min=min(c(sapply(y,min),sapply(y_seg,min)))
  y_max=max(c(sapply(y,max),sapply(y_seg,max)))

  ggplot2::ggplot(data, ggplot2::aes(x=x,y=y))+
    ggplot2::labs(x=xlab, y=ylab)+
    ggplot2::theme(text = ggplot2::element_text(size=10))+
    ggplot2::geom_text(label=point_lab,hjust = "inward")+
    ggplot2::geom_segment(data=as.data.frame(data_seg),
                          ggplot2::aes(x=0,xend=x_seg,y=0,yend=y_seg),
                 colour="#990000",
                 stat="identity",inherit.aes = FALSE) +
    ggplot2::geom_text(data=as.data.frame(data_seg),
                       ggplot2::aes(x=x_seg,y=y_seg,label=vec_lab,
                                    size=3,col="#990000"),hjust = "inward")+
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size=10))+
    ggplot2::xlim(x_min-15,x_max+15)+
    ggplot2::ylim(y_min-15,y_max+15)}

# this will save the last name in the taxonomy line
name.split<-function(names){
  names=names
  save=strsplit(names,"_")

  h=0
  for(i in seq(1,length(names))){
    h[i]=length(save[[i]])}

  i=0
  name.list=NULL
  for(i in seq(1,length(save))){
    name.list[i]=save[[i]][h[i]]}
  name.list
}

# this will create data for plots found in the T3 function
# in r but not in the T3func in R
T3_plots<-function(r1,r2,r3,T3,mode,plot_scores){
  A=T3$A
  B=T3$B
  C=T3$C
  Ha = T3$H
  Hb = ThreeWay::permnew(Ha, r1, r2, r3)
  Hc = ThreeWay::permnew(Hb, r2, r3, r1)
  CBplot = Matrix::kronecker(C,B) %*% t(Ha)
  ACplot = Matrix::kronecker(A, C) %*% t(Hb)
  BAplot = Matrix::kronecker(B, A) %*% t(Hc)
  Aplot=T3$A%*%diag(diag(T3$La)^.5,nrow=r1)
  Bplot=T3$B%*%diag(diag(T3$Lb)^.5,nrow=r2)
  Cplot=T3$C%*%diag(diag(T3$Lc)^.5,nrow=r3)

  A_eig=diag(T3$La)^.5
  B_eig=diag(T3$Lb)^.5
  C_eig=diag(T3$Lc)^.5

  if(plot_scores){
    if(mode == "ACplot"){
      ACplot <- ACplot %*% B_eig
    } else if(mode == "BAplot"){
      BAplot <- BAplot %*% C_eig
    } else if(mode == "CBplot"){
      CBplot <- CBplot %*% A_eig
    }
  }

  list=list(A,B,C,Ha,Hb,Hc,CBplot,ACplot,BAplot,Aplot,Bplot,Cplot,A_eig,B_eig,C_eig)
  names(list)=c("A","B","C","Ha","Hb","Hc","CBplot","ACplot","BAplot","Aplot","Bplot","Cplot","A_eig","B_eig","C_eig")
  list
}

U_matrices=function(X, method, Acomp, Bcomp, Ccomp, j, m, k){

  method=as.character(method)
  Acomp=as.numeric(Acomp)
  Bcomp=as.numeric(Bcomp)
  Ccomp=as.numeric(Ccomp)

  double_center<-function(mat,n){
    mat_sqr=as.matrix(mat)^2
    mat_center=diag(1,n)-1/n * matrix(1,ncol=n,nrow=n)
    mat_double_center= -(1/2)*(mat_center%*%mat_sqr%*%mat_center)
    mat_double_center
  }

  #### A mode
  tuck_A=as.matrix(vegan::vegdist(X,method = method))
  h=double_center(tuck_A,j)
  A=svd(h)$u[,seq(1,Acomp)]

  #### B mode
  tuck_A1=ThreeWay::permnew(X, j,m,k)
  tuck_B=as.matrix(vegan::vegdist(tuck_A1,method = method))
  n=dim(tuck_A1)[1]
  h=double_center(tuck_B,n)
  B=svd(h)$u[,seq(1,Bcomp)]

  ###C mode
  tuck_A2=ThreeWay::permnew(X, m, k,j)
  tuck_C=as.matrix(vegan::vegdist(tuck_A2,method = method))
  n=dim(tuck_A2)[1]
  h=double_center(tuck_C,n)
  C=svd(h)$u[,seq(1,Ccomp)]

  G = t(A) %*% as.matrix(X) %*% (kronecker(B,C))

  list=list(A,B,C,G)
  names(list)=c("A","B","C","G")
  list
}
