

## 1. ts.test原始版本    ts.0603.c4
#' Title  ts.svg1  Finding SVG.
#'
#' @param data A data.frame with counts and coordinates.(spots x genes)
#' The first two columns are coordinates. The first two colnames must be 'row' and 'col'.
#' @param scale A logic value.T/F
#' If scale==T, the expression counts matrix runs function scale.count.
#' @param noise  A logic value.T/F.
#' If noise==F, the SVG-rank considers only the p-value and adjusted p-value.
#' If noise==T, the SVG-rank takes into account the p-value,adjusted p-valus
#'  and the distribution characteristics of the genes.
#' @param padj.m A character string.
#' The adjusted method of p-values. The default value is 'holm'.
#' @return A data.frame. All genes' p values, adjusted p values and SVG-rank.
ts.svg1=function(data,scale=T,noise=T,padj.m='holm'){

  if (scale==T){new=cbind(data[,1:2],scale.count(data[,-c(1:2)]))}
  if (scale==F){new=data}
  new=as.data.frame(new)
  # new=data
  locus_in=new[,c(1,2)]
  counts=new[,-c(1,2)]

  l1='row'
  l2='col'
  z1.group=ceiling(log(diff(range(locus_in[l1])))) # 向上取整
  z2.group=ceiling(log(diff(range(locus_in[l2])))) # 向上取整

  locus_in1=data.frame(locus_in,
                       n.row=cut(x=as.matrix(locus_in[l1]),
                                 breaks = seq(from=min(locus_in[l1]),to=max(locus_in[l1]),
                                              length=z1.group+1),
                                 labels =1:z1.group,
                                 include.lowest = T,right = T))

  locus_in2=data.frame(locus_in,
                       n.col=cut(x=as.matrix(locus_in[l2]),
                                 breaks = seq(from=min(locus_in[l2]),to=max(locus_in[l2]),
                                              length=z2.group+1),
                                 labels =1:z2.group,
                                 include.lowest = T,right = T))

  new1=cbind(locus_in1[colnames(locus_in1)!=l1],counts)  # n.row
  new2=cbind(locus_in2[colnames(locus_in2)!=l2],counts)  #n.col
  colnames(new1)[1:2]<-c('x','coor.z')
  colnames(new2)[1:2]<-c('x','coor.z')

  new.row=aggregate(new1[,-c(1:2)],list(new1$x,new1$coor.z),mean)

  new.col=aggregate(new2[,-c(1:2)],list(new2$x,new2$coor.z),mean)

  ####   不匹配的问题出现在了这里  --0709
  zero.p=apply(new[,-c(1:2)], 2,function(y){sum(y!=0,na.rm = T)/length(y)})
  mean=apply(new[,-c(1:2)], 2,
             function(y){mean(y[y!=0],na.rm=T)})
  sum=data.frame(gene=names(mean),
                 zero.p,
                 mean)

  #### 0630test版本: apply(new.row[,-c(1:3)]改为[,-c(1:2)]
  row.t=apply(new.row[,-c(1:2)], 2,
              function(y){z=ifelse(sum(y!=0)==0,1,Box.test(y,lag=z1.group)$p.value);z})
  col.t=apply(new.col[,-c(1:2)], 2,
              function(y){z=ifelse(sum(y!=0)==0,1,Box.test(y,lag=z2.group)$p.value);z})


  new.x=aggregate(new[,-c(1:2)],list(new$row),mean)
  new.y=aggregate(new[,-c(1:2)],list(new$col),mean)

  x.t=apply(new.x[,-1], 2,
            function(y){Box.test(y,lag=z1.group)$p.value})
  y.t=apply(new.y[,-1], 2,
            function(y){Box.test(y,lag=z2.group)$p.value})

  test=data.frame(row.t,col.t,
                  x.t,y.t,
                  gene = names(row.t)
  )


  #within
  # ### 处理掉0和1
  # test[which(test==1,arr.ind = T)]=1-10^-6
  # test[which(test==0,arr.ind = T)]=10^-6

  test[,1:(ncol(test)-1)]=sapply(test[,1:(ncol(test)-1)],
                                 function(y){z=ifelse(is.na(y)==T,0.999,y)})

  test$min=apply(test[,1:(ncol(test)-1)],1,
                 function(y){poolr::stouffer(y)$p})

  test$adj_min=p.adjust(test$min,method=padj.m)   #

  #BY or BH，BH更加宽松
  test=merge(test,sum,by='gene')
  test=subset(test,zero.p>0)
  test$zero.p=test$zero.p/max(test$zero.p)
  test$mean=test$mean/max(test$mean)
  test$c2=(test$zero.p+test$mean)/2

  if(noise==T){
    # setorder(test,adj_min,min,comb,-zero.p,-mean)
    data.table::setorder(test,adj_min,min,-zero.p,-mean)
    # a1=c('gene','x.p_val','min','adj_min')
    a1=c('gene','min','adj_min')
    test2=test[a1]
    # b1=c('gene','x.p_val','min','p_adj')
    # b1=c('gene','min','p_adj')
    b1=c('gene','pval','p_adj')

    colnames(test2)<-b1
    test2$rank=1:nrow(test2)

  }

  if(noise==F){
    data.table::setorder(test,adj_min,min)
    # a1=c('gene','x.p_val','min','adj_min')
    a1=c('gene','min','adj_min')
    test2=test[a1]
    # b1=c('gene','x.p_val','min','p_adj')
    # b1=c('gene','min','p_adj')
    b1=c('gene','pval','p_adj')

    colnames(test2)<-b1
    test2$rank=1:nrow(test2)

  }

  test2

}




## 2. ts.test 并行版本------------------

#' Title  ts.svg2  Finding SVG (parallel).
#'
#' @param data A data.frame with counts and coordinates.(spots x genes)
#' The first two columns are coordinates. The first two colnames must be 'row' and 'col'.
#' @param scale A logic value.T/F
#' If scale==T, the expression counts matrix runs function scale.count.
#' @param noise  A logic value.T/F.
#' If noise==F, the SVG-rank considers only the p-value and adjusted p-value.
#' If noise==T, the SVG-rank takes into account the p-value,adjusted p-valus
#'  and the distribution characteristics of the genes.
#' @param padj.m A character string.
#' The adjusted method of p-values. The default value is 'holm'.
#'
#' @return A data.frame. All genes' p values, adjusted p values and SVG-rank.
ts.svg2=function(data,scale=T,noise=T,padj.m='holm'){

  if (scale==T){new=cbind(data[,1:2],scale.count(data[,-c(1:2)]))}
  if (scale==F){new=data}
  new=as.data.frame(new)

  # new=data
  locus_in=new[,c(1,2)]
  counts=new[,-c(1,2)]


  l1='row'
  l2='col'
  z1.group=ceiling(log(diff(range(locus_in[l1])))) # 向上取整
  z2.group=ceiling(log(diff(range(locus_in[l2])))) # 向上取整

  locus_in1=data.frame(locus_in,
                       n.row=cut(x=as.matrix(locus_in[l1]),
                                 breaks = seq(from=min(locus_in[l1]),to=max(locus_in[l1]),
                                              length=z1.group+1),
                                 labels =1:z1.group,
                                 include.lowest = T,right = T))

  locus_in2=data.frame(locus_in,
                       n.col=cut(x=as.matrix(locus_in[l2]),
                                 breaks = seq(from=min(locus_in[l2]),to=max(locus_in[l2]),
                                              length=z2.group+1),
                                 labels =1:z2.group,
                                 include.lowest = T,right = T))

  new1=cbind(locus_in1[colnames(locus_in1)!=l1],counts)  # n.row
  new2=cbind(locus_in2[colnames(locus_in2)!=l2],counts)  #n.col
  colnames(new1)[1:2]<-c('x','coor.z')
  colnames(new2)[1:2]<-c('x','coor.z')

  ##counts1 and counts2, for parallel
  counts1=data.frame(x=new$row,coor.z=1,
                     new[,-c(1:2)])

  counts2=data.frame(x=new$col,coor.z=1,
                     new[,-c(1:2)])


  ####   不匹配的问题出现在了这里  --0709
  zero.p=apply(new[,-c(1:2)], 2,function(y){sum(y!=0,na.rm = T)/length(y)})
  mean=apply(new[,-c(1:2)], 2,
             function(y){mean(y[y!=0],na.rm=T)})
  sum=data.frame(gene=names(mean),
                 zero.p,
                 mean)

  ## all.list and lag.ver for do parallel
  all.list=list(counts1,counts2,new1,new2)
  lag.vec=c(z1.group,z2.group,z1.group,z2.group)
  ###### 0716 parallel -------------
  # a2=Sys.time()
  clnum=parallel::detectCores(logical = F)
  cl.num=ifelse(clnum>1,clnum-1,clnum)
  cl <- parallel::makeCluster(cl.num)
  doParallel::registerDoParallel(cl)

  pooling=foreach(i=1:length(all.list),.verbose = F,.inorder = T,.combine = cbind) %dopar%{
    m1=all.list[[i]]
    sta=aggregate(m1[,-c(1:2)],list(m1$x,m1$coor.z),mean)

    pval=apply(sta[,-c(1:2)], 2,
               function(y){z=ifelse(sum(y!=0)==0,1,Box.test(y,lag=lag.vec[i])$p.value);z})


  }
  parallel::stopCluster(cl)   #别忘了结束并行
  # b2=difftime(Sys.time(),a2,units = 'secs')

  colnames(pooling)<-c('x.t','y.t','semi.xt','semi.yt')  ## 'x.t','y.t','row.t','col.t'

  test=data.frame(pooling,
                  gene=rownames(pooling))

  ##  parallel end --------------

  #within

  # ### 处理掉0和1
  # test[which(test==1,arr.ind = T)]=1-10^-6
  # test[which(test==0,arr.ind = T)]=10^-6

  test[,1:(ncol(test)-1)]=sapply(test[,1:(ncol(test)-1)],
                                 function(y){z=ifelse(is.na(y)==T,0.999,y)})

  test$min=apply(test[,1:(ncol(test)-1)],1,
                 function(y){poolr::stouffer(y)$p})
  test$adj_min=p.adjust(test$min,method=padj.m)   #

  #BY or BH，BH更加宽松
  test=merge(test,sum,by='gene')
  # test=test[order(test$adj_min),]
  test=subset(test,zero.p>0)
  test$zero.p=test$zero.p/max(test$zero.p)
  test$mean=test$mean/max(test$mean)
  test$c2=(test$zero.p+test$mean)/2

  if(noise==T){
    # setorder(test,adj_min,min,comb,-zero.p,-mean)
    data.table::setorder(test,adj_min,min,-zero.p,-mean)
    # a1=c('gene','x.p_val','min','adj_min')
    a1=c('gene','min','adj_min')
    test2=test[a1]
    # b1=c('gene','x.p_val','min','p_adj')
    # b1=c('gene','min','p_adj')
    b1=c('gene','pval','p_adj')
    colnames(test2)<-b1
    test2$rank=1:nrow(test2)

  }

  if(noise==F){
    data.table::setorder(test,adj_min,min)
    # a1=c('gene','x.p_val','min','adj_min')
    a1=c('gene','min','adj_min')
    test2=test[a1]
    # b1=c('gene','x.p_val','min','p_adj')
    # b1=c('gene','min','p_adj')
    b1=c('gene','pval','p_adj')

    colnames(test2)<-b1
    test2$rank=1:nrow(test2)

  }

  test2

}






