
#' Title:basic functions
#'
#' @param data  A data.frame or matrix of expression counts
#' @return A data.frame or matrix.
simu.scale=function(data){
  new.data=apply(data, 2, function(y){y/max(y)})
  new.data
}


#' Title:scale.count  (original name: simu.count)
#' Transformation of expression counts.
#' To adjust the effects from library size between sample spots and outliers of gene expression.
#' @param data A data.frame or matrix of expression counts (spots x genes)
#' @param qh  A numeric value should be between 0 and 1.
#' @return A data.frame or matrix of adjusted expression counts (spots x genes)
scale.count=function(data,qh=0.985){
  counts_mat=apply(data, 1, function(y){z=y/max(y)})  ## 以后改成稀疏矩阵的计算
  ## 类似于调节测序深度
  counts_mat=apply(counts_mat, 2, function(y)
  {y[which(y>quantile(y,qh))]<-quantile(y,qh);
  y[which(y<quantile(y,0.25))]<-0;y})
  # counts_frame=data.frame(counts_mat)
  # counts_frame=data.frame(cell=rownames(counts_mat),
  #                         counts_mat)
  # rm(counts_mat)
  # colnames()
  # counts_mat=apply(counts_mat, 2, function(y){z=y/max(y)})  ## 以后改成稀疏矩阵的计算

  t(counts_mat)
}


#' Title:simu_zi
#' The basic function of simulation generation functions.
#'
#' @param family 'ZINB' or 'ZIP',zero-inflated negative binomial or zero-inflated Poisson distribution.
#'
#' @param subject.n the number of spots,because each gene follow a zero-infalted P/NB distribution
#' @param zi.p The zero proportion of the zero generation process.
#' if  zi.p=0, it generates data which follows Poisson or NB distribution.
#' @param mu  the lambda para in the poisson distribution or the mu para in the NB distribution
#' @param size the size para in the NB distribution
#'
#' @return  A numeric vector.
simu_zi=function(family,subject.n,zi.p=0.5,
                 mu=0.5,size=0.25){
  ##' family: ZIP or ZINB
  ##' subject.n: the number of spots,because each gene follow a zero-infalted P/NB distribution
  ##' zi.p: zero proportion of the zero inflated model
  ##' if  zi.p=0, it generates data which follows Poisson or NB distribution.
  ##' mu: the lambda para in the poisson distribution or the mu para in the NB distribution
  ##' size: the size para in the NB distribution

  Y = rep(NA,subject.n)
  # set.seed(sim.seed)
  ind.mix <- rbinom(length(Y), 1, zi.p)
  ##' p的比例为成功（不为0）的比例，然后通过以下code,将这一步骤中的非0值改为0
  ##' 0值改为非0值

  if (family=='ZIP'){
    Y[which(ind.mix!=0)]<-0
    Y[which(ind.mix==0)]= rpois(n=sum(ind.mix==0), lambda=mu)
  }

  if (family=='ZINB'){
    Y[which(ind.mix!=0)]<-0
    Y[which(ind.mix==0)]= rnbinom(n=sum(ind.mix==0),size=size,mu=mu)
  }

  return(Y)
}





##  other  func ----------------
#' Title   filter.gene
#' This function filters the non-expressed genes and low-proportion expressed genes.
#' @param object : A seurat object.
#' @param thr A number.
#' If thr is between 0 and 1, it indicates the least ratio of each gene expressed in the spots.
#' If thr is greater than 1, it indicates the minimum number of each gene expressed in the spots.
#' @return A seurat object.
filter.gene=function(object,thr){
  a=ifelse((thr>=0&thr<1),thr*nrow(object),thr)
  object@assays$RNA@data=object@assays$RNA@counts[rowSums(object@assays$RNA@counts!=0,na.rm=T)>=a,]
  object
}


#### seurat scale---------------
## 从seurat.object转换为data.frame

##  seurat.trans  :(原来的scale.seu.spk)
#' Title seurat.trans
#'
#' @param object :A seurat object.
#' @param counts A logic value.T/F
#' If counts==T, the expression matrix is 'counts' in the 'assays' slots.
#' If counts==F, the expression matrix is 'data' in the 'assays' slots.
#' @param sps A logic value.T/F
#' If sps==T, the outpot is a sparse matrix.
#' @return A data.frame of matrix of counts and coordinates.
seurat.trans=function(object,counts = F,sps=T){
  if (counts==T){
    c2=t(object@assays$RNA@counts)
  } else{
    c2=t(object@assays$RNA@data)
  }

  counts_frame=data.frame(cell=rownames(c2),
                          c2)
  coor=data.frame(object@images[[Images(object)]]@coordinates,
                  cell= rownames(object@images[[Images(object)]]@coordinates))
  coor=coor[c('cell','row','col')]
  data=merge(coor,counts_frame,by='cell')
  rownames(data)<-data$cell
  data=data[,-which(colnames(data)=='cell')]

  if (sps==T){
    dt=as(as.matrix(data), "sparseMatrix")
    rm(data)
  }else{
    dt=data
    rm(data)
  }
  dt
}





#' Title index2 Calculating performances indexes
#' @param data A data.frame has at least two columns.
#' (gene name and p-values (or adjusted p-values))
#' @param var  The colname of p-values or adjusted p-values.
#' @param alpha Threshold.
#' @return A numeric vector.
cal.index=function(data,var,alpha){
  data=as.data.frame(data)
  data$g.label=do.call(rbind,str_split(data$gene,pattern = '\\.',n=2))[,1]
  data$is.svg=ifelse(data$g.label=='se',1,0)

  new=data[data[var]<=alpha,]  ## 符合要求的GENE

  other=data[data[var]>alpha,]

  all=nrow(data)
  t=sum(data$is.svg==1,na.rm=T)
  f=sum(data$is.svg==0,na.rm=T)
  tp=sum(new$is.svg==1,na.rm=T)
  fp=sum(new$is.svg==0,na.rm=T)
  tn=sum(other$is.svg==0,na.rm=T)
  fn=sum(other$is.svg==1,na.rm=T)

  # t=sum(data$g.label=='se',na.rm=T)
  # f=sum(data$g.label=='ns',na.rm=T)
  # tp=sum(new$g.label=='se',na.rm=T)
  # fp=sum(new$g.label=='ns',na.rm=T)
  # tn=sum(other$g.label=='ns',na.rm=T)
  # fn=sum(other$g.label=='se',na.rm=T)

  recall=tp/t
  prec=tp/(tp+fp)
  spec=tn/(tn+fp)

  fdr=fp/(tp+fp)

  f1_score=2*recall*prec/(prec+recall)

  out=data.frame(recall,precision=prec,specificity=spec,
                 fdr,f1_score,tp,fp,fn,all)

  out



}






















