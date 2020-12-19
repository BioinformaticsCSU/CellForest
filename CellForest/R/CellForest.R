#' CellForest
#'
#' Cluster samples based on functional similarity using gene-set based local similarity
#' @param data  A expression data matrix,gene in rows and samples in columns.
#' @param k  The number of clusters to output.
#' @param ncores The number of cores to be used when the program running in parallel.
#' @return samplegeneweight The weights of genes calculated by function randomForest based on random permuted predict label.
#' @return geneImportance The weights of genes calculated by function randomForest.
#' @return selgenes The genes Cell Forest final select.
#' @return selgeneratio The proportion of the number of selected genes in the total genes.
#' @return runningtime The running time of the program.
#' @return cutoff Threshold for selecting genes.
#' @return thresholdMethod Method for selecting genes.
#' @return plabel Predicated clusters.
#' @keywords clustering, unsupervised feature selection
#' @export
#' @author Hongdong Li, lhdcsu@gmail.com, Central South University.
#' @examples
#' data(cfdemo)
#' result = CellForest(data,k = 3)



CellForest<-function(data,k,ncores=-1){
  timestart <- proc.time()
  
  nsample=ncol(data)
  ngene=nrow(data)
  geneweight0=rep(0,length=ngene)
  names(geneweight0)=rownames(data)
  thresholdMethod="mean3sd"
  # Main
  # step 1: generate labels by clustering anlaysis
  plabel = preClustering(data,k)
  
  # run Feature selection
  cat("Run randomForest for feature selection.\n")
  
  library(randomForest)
  randomLabel = sample(as.numeric(plabel))
  # rfRandom=randomForest(t(data),factor(randomLabel),importance = T)
  # grank=rfRandom$importance[,k+1]
  rfRandom = ranger(y = factor(randomLabel),x = t(data),importance = 'impurity')
  grank=rfRandom$variable.importance
  
  
  samplegeneweight=geneweight0
  samplegeneweight[names(grank)]=grank
  
  max=samplegeneweight[which.max(samplegeneweight)]
  mean=mean(samplegeneweight)
  sd=sd(samplegeneweight)
  msd=mean+sd
  m3sd=mean+3*sd
  Uquartile=quantile(samplegeneweight, probs = 0.75)
  
  # rf=randomForest(t(data),factor(plabel),importance = T)
  # grank=rf$importance[,k+1]
  rfRandom = ranger(y = factor(plabel),x = t(data),importance = 'impurity')
  grank=rfRandom$variable.importance
  
  geneweight=geneweight0
  geneweight[names(grank)]=grank
  
  if(thresholdMethod == "zero"){
    goodgenes=names(geneweight[geneweight>0])
    cutoff = 0
  }else if(thresholdMethod == "max"){
    goodgenes=names(geneweight[geneweight>max])
    cutoff = max
  }else if(thresholdMethod == "meansd"){
    goodgenes=names(geneweight[geneweight>msd])
    cutoff = msd
  }else if(thresholdMethod == "mean3sd"){
    goodgenes=names(geneweight[geneweight>m3sd])
    cutoff = m3sd
  }else if(thresholdMethod == "UpperQuartile"){
    goodgenes=names(geneweight[geneweight>Uquartile])
    cutoff = Uquartile
  }
  
  nselgenes = length(goodgenes)
  selgeneratio = length(nselgenes)/ngene
  cat(nselgenes,' genes with imp >= ',cutoff,'\n')
  
  runningtime <- proc.time() - timestart
  
  return(list(samplegeneweight = samplegeneweight,
              geneImportance = geneweight,
              selgenes = goodgenes,
              selgeneratio = selgeneratio,
              runningtime = runningtime,
              cutoff = cutoff,
              plabel = plabel))
}


preClustering<-function(data,k){
  # cluster label prediction
  # X: A expression data matrix,gene in rows and samples in columns.
  nsamples = ncol(data)
  dist1 = as.matrix(compute_dist(data, distmethod = "pearson"))
  p1 <- prcomp(dist1,center = TRUE, scale. = TRUE)
  sdev1 <- p1$sdev[1:20]
  
  meanresiduals = c()
  for (dim in 2:15) {
    lm1 = lm(sdev1[1:dim]~c(1:dim))
    lm2 = lm(sdev1[(dim+1):20]~c((dim+1):20))
    meanresiduals = c(meanresiduals,(mean(lm1$residuals)+mean(lm2$residuals)))
  }
  estimate.K1 = c(2:15)[which.min(meanresiduals)]
  
  plabel = kmeans(p1$rotation[,1:estimate.K1],centers = k,iter.max = 1e+09,nstart = 1000)$cluster
  return(plabel)
}


compute_dist<-function(X,distmethod=c("pearson","euclidean","spearman")){
  # calculate distance between columns
  # X: A expression data matrix,gene in rows and samples in columns.
  
  if (distmethod %in% c("pearson","spearman")){
    corr = cor(X,method=distmethod)
    D=as.dist(1-corr)
  }else if(distmethod == "euclidean"){
    D=dist(t(X),method=distmethod)
  }
  return(D)
}




