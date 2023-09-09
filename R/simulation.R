


## 1 random

#' Title   random.zi: Simulation data. Genes with random pattern.
#'
#' @param lambda The spots detection efficiency.It should be between 0 and 1.
#' The default value is 0.5.
#' @param spots The number of spots.
#' @param g  The number of genes.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param g.p A number. The zero proportion of the zero generation process.
#' @param g.size The size para in the NB distribution.
#' @param g.mu the lambda para in the poisson distribution or the mu para in the NB distribution
#' @param outlier  A numeric value should be between 0 and 1.
#' @return A data.frame of gene expression counts and coordinates.
random.zi=function(lambda=0.5,spots,g,type='ZINB',  ## type对应simu_zi中的 family 参数
                   g.p,  ## 对应simu_zi中的zi.p
                   g.size=0.25,g.mu,outlier=F){

  ##' spots: the sample size
  ##' g : the number of genes
  ##' type: corresponding to the 'family' parameter of the simu_zi function
  ##' g.p: corresponding to the 'zi.p' parameter of the simu_zi function

  ## spot generation
  win=ceiling(sqrt(spots/lambda))
  coor=spatstat.random::rpoispp(lambda = lambda,
                                win=spatstat.geom::owin(c(0,win),
                                                        c(0,win)))
  coor.s1=data.frame(row=round(coor$x,0),
                     col=round(coor$y,0),
                     row.names = paste0('c-',1:length(coor$x)))


  # exp=matrix(rnbinom(g*nrow(coor.s1),size=g.size,mu=g.mu),
  #            ncol=g,dimnames = list(NULL,paste0('ns.',1:g)))## nb(r=size,p=prob)

  ##  新的  simu_zi
  # 用nrow(coor.s1)而非spots是因为随机生成位置时，会存在波动
  exp=matrix(NA,nrow=nrow(coor.s1),
             ncol=g,dimnames = list(NULL,paste0('ns.',1:g)))## nb(r=size,p=prob)

  exp=apply(exp, 2,
            function(y){z=simu_zi(family = type,subject.n = nrow(exp),
                                  zi.p=g.p,mu=g.mu,size=g.size);z})





  all=cbind(coor.s1,exp)



  # rownames(all)<-all$num
  # all=all[,-1]
  if(outlier==FALSE){
    end=all
  }
  if(outlier>=1 | outlier<=0){
    cat("# outlier parameter is wrong! \n")
    end=all
  }
  if(outlier>=0 & outlier <1){
    ind=sample(nrow(all),round(nrow(all)*outlier));
    out.para=10;
    # all[ind,-c(1:2)]<-matrix(simu_zi(family = type,subject.n = (length(ind)*(ncol(all)-2)),
    #                                  zi.p=g.p/2,
    #                                  size=g.size*out.para,
    #                                  mu=g.mu*out.para),
    #                          nrow = length(ind))

    out=sample(c(max(exp),max(exp)*2,max(exp)*1.5,g.mu*out.para),
               length(ind)*(ncol(all)-2),replace = T)

    all[ind,-c(1:2)]<-out

    # for (i in ind) {
    #   # j=sample(3:ncol(all),round(ncol(all)*outlier))
    #   # out=sample(c(max(exp),max(exp)*2,max(exp)*1.5,g.mu*out.para),
    #   #            length(j),replace = T)
    #   # all[i,j]<-out
    #
    #   # j=sample(3:ncol(all),round(ncol(all)*outlier))
    #   out=sample(c(max(exp),max(exp)*2,max(exp)*1.5,g.mu*out.para),
    #              (ncol(all)-2),replace = T)
    #
    #   all[i,3:ncol(all)]<-out
    #
    #
    #
    #
    # }



    end=all
  }
  end

}



## 2 streak

#' Title Simulation data. Genes with streak expression pattern.
#'
#' @param lambda The spots detection efficiency.It should be between 0 and 1.
#' The default value is 0.5.
#' @param spots The number of spots.
#' @param prop  The proportion of streak area.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param se.size  The size para of SVGs in the NB distribution.
#' @param se.mu  For SVGs, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return A data.frame of gene expression counts and coordinates.
streak.zi=function(lambda=0.5,spots,prop,se, ns, ## prop,pattern占比，se,ns是SVG和non-SVG数量
                   type, ## family
                   se.p,se.size=0.2,se.mu,  ##SVG high 参数
                   ns.p,ns.size=0.2,ns.mu,   ## non-SVG 以及SVG low参数
                   outlier=F){
  ## spot generation
  win=ceiling(sqrt(spots/lambda))
  coor=spatstat.random::rpoispp(lambda = lambda,
                                win=spatstat.geom::owin(c(0,win),
                                                        c(0,win)))
  coor.s1=data.frame(row=round(coor$x,0),
                     col=round(coor$y,0),
                     row.names = paste0('c-',1:length(coor$x)))

  ###  spots  marks
  ## 以row为标准 取样
  max=max(coor.s1$row)
  cen=sample(0:win,1)
  if ((cen/max-prop/2)>=0 & (cen/max+prop/2)<=1){
    # if ((cen-max*prop/2)>=0 & (cen+max*prop/2)<=max){
    num.streak=data.frame(row=(floor(cen-max*prop/2)):(floor(cen+max*prop/2)))
    num.random=data.frame(row=setdiff(unique(coor.s1$row),num.streak$row))
  }
  if((cen/max-prop/2)<=0){
    num.streak=data.frame(row=0:ceiling(max*prop))
    num.random=data.frame(row=setdiff(unique(coor.s1$row),num.streak$row))
  }

  if((cen/max+prop/2)>=1){
    num.streak=data.frame(row=max:(max-floor(max*prop)))
    num.random=data.frame(row=setdiff(unique(coor.s1$row),num.streak$row))
  }

  ## coordinate
  coor.streak=data.frame(merge(coor.s1,num.streak,by='row'),
                         row.names = paste0('s-',1:nrow(merge(coor.s1,num.streak,by='row'))))

  coor.random=data.frame(merge(coor.s1,num.random,by='row'),
                         row.names = paste0('r-',1:nrow(merge(coor.s1,num.random,by='row'))))

  ## corresponding expression
  ## HIGH
  # exp.streak=matrix(rnbinom(se*nrow(coor.streak),size=se.size,mu=se.mu),
  #                   ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
  ## 新ZI版本

  exp.streak=matrix(NA,nrow=nrow(coor.streak),ncol=se,
                    dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  exp.streak=apply(exp.streak, 2,
                   function(y){z=simu_zi(family = type,subject.n = nrow(exp.streak),
                                         zi.p=se.p,mu=se.mu,size=se.size);z})

  ## LOW
  # exp.random=matrix(rnbinom(se*nrow(coor.random),size=ns.size,mu=ns.mu),
  #                   ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  ## 新ZI版本
  exp.random=matrix(NA,nrow=nrow(coor.random),ncol=se,
                    dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  exp.random=apply(exp.random, 2,
                   function(y){z=simu_zi(family = type,subject.n = nrow(exp.random),
                                         zi.p=ns.p,mu=ns.mu,size=ns.size);z})

  exp.svg=rbind(exp.streak,exp.random)
  ## NON svg
  non.coor=rbind(coor.streak,coor.random)
  # exp.non=matrix(rnbinom(ns*nrow(non.coor),size=ns.size,mu=ns.mu),
  #                ncol=ns,dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)

  ## 新ZI版本
  exp.non=matrix(NA,nrow=nrow(non.coor),ncol=ns,
                 dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)

  exp.non=apply(exp.non, 2,
                function(y){z=simu_zi(family = type,subject.n = nrow(exp.non),
                                      zi.p=ns.p,mu=ns.mu,size=ns.size);z})

  all=cbind(non.coor,exp.svg,exp.non)
  if(outlier==FALSE){
    end=all
  }
  if(outlier>=1 | outlier<=0){
    cat("# outlier parameter is wrong! \n")
    end=all
  }
  if(outlier>=0 & outlier <1){
    ind=sample(nrow(all),round(nrow(all)*outlier));
    out.para=10;

    # out=sample(c(max(exp.svg),max(exp.svg)*2,max(exp.svg)*1.5,se.mu*out.para),
    #            length(ind)*(ncol(all)-2),replace = T)
    #
    # all[ind,-c(1:2)]<-out

    # out=sample(c(max(exp.svg),max(exp.svg)*2,max(exp.svg)*1.5,se.mu*out.para),
    #            length(ind)*(ncol(all)-2),replace = T)
    # out=sample(c(max(exp.svg),max(exp.svg)*0.5,max(exp.svg)*1.5,max(exp.svg)*2),
    #            length(ind)*(ncol(all)-2),replace = T)
    out=sample((se.mu*c(5,10,15)),
               length(ind)*(ncol(all)-2),replace = T)

    all[ind,-c(1:2)]<-out




    end=all
  }
  end

}


## 3 hotspot

#' Title Simulation data. Genes with hotspot expression pattern.
#' @param lambda The spots detection efficiency.It should be between 0 and 1.
#' The default value is 0.5.
#' @param spots The number of spots.
#' @param prop  The proportion of streak area.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param se.size  The size para of SVGs in the NB distribution.
#' @param se.mu  For SVGs, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.

hotspot.zi=function(lambda=0.5,spots,prop,se, ns, ## prop,pattern占比，se,ns是SVG和non-SVG数量
                    type,
                    se.p,se.size=0.2,se.mu,  ##SVG high 参数
                    ns.p,ns.size=0.2,ns.mu, ## non-SVG 以及SVG low参数
                    outlier=FALSE){

  ## 1 spot generation
  win=ceiling(sqrt(spots/lambda))
  coor=spatstat.random::rpoispp(lambda = lambda,
                                win=spatstat.geom::owin(c(0,win),
                                                        c(0,win)))
  coor.s1=data.frame(num=paste0('c-',1:length(coor$x)),
                     row=round(coor$x),# row=round(coor$x,0),
                     col=round(coor$y),# col=round(coor$y,0),
                     row.names = paste0('c-',1:length(coor$x)))

  ## 2  spots  marks
  ##  取圆心
  max.x=max(coor.s1$row)
  max.y=max(coor.s1$col)
  ## radius of hotspot
  r=sqrt(prop/pi)  ## 单位圆上，用于确定比例
  r1=floor(r*max.x)   ## spots上的radius

  cen.x=floor(sample((r*max.x):((1-r)*max.x),1))
  cen.y=floor(sample((r*max.y):((1-r)*max.y),1))

  ##  range
  range.x=(cen.x-floor(r*max.x)):(cen.x+floor(r*max.x))
  y.l=ceiling(cen.y-sqrt(r1^2-range.x^2-cen.x^2+2*range.x*cen.x))
  y.u=floor(cen.y+sqrt(r1^2-range.x^2-cen.x^2+2*range.x*cen.x))

  coor.hot=foreach(i=1:length(y.l),.combine = rbind) %do%{
    if (y.l[i]<y.u[i]){
      z=y.l[i]:y.u[i]
    } else if (y.l[i]==y.u[i]){
      z=y.l[i]
    }
    end=data.frame(row=range.x[i],col=z)
  }
  coor.hot=data.frame(coor.hot,
                      row.names = paste0('h-',1:nrow(coor.hot)))

  # num.hotspot=data.frame(num=unique(num.hotspot$num)) ## 去除重复
  #
  coor.hotspot=merge(coor.s1,coor.hot,
                     by=c('row','col'),no.dups = T)

  num.random=data.frame(num=setdiff(coor.s1$num,coor.hotspot$num))
  coor.random=merge(coor.s1,num.random,
                    by='num')
  ## 3 corresponding expression
  #### HIGH
  # exp.hotspot=matrix(rnbinom(se*nrow(coor.hotspot),size=se.size,mu=se.mu),
  #                    ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  ##  新ZI 版本
  exp.hotspot=matrix(NA,nrow=nrow(coor.hotspot),ncol=se,
                     dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  exp.hotspot=apply(exp.hotspot, 2,
                    function(y){z=simu_zi(family = type,subject.n =nrow(exp.hotspot),
                                          zi.p=se.p, size=se.size,mu=se.mu);z})
  #### LOW
  # exp.random=matrix(rnbinom(se*nrow(coor.random),size=ns.size,mu=ns.mu),
  #                   ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

  ##  新ZI 版本
  exp.random=matrix(NA,nrow=nrow(coor.random), ncol=se,
                    dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
  exp.random=apply(exp.random, 2,
                   function(y){z=simu_zi(family=type,subject.n = nrow(exp.random),
                                         zi.p=ns.p,size=ns.size,mu=ns.mu);z})

  exp.svg=rbind(exp.hotspot,exp.random)
  ### NON svg
  non.coor=rbind(coor.hotspot,coor.random)
  # exp.non=matrix(rnbinom(ns*nrow(non.coor),size=ns.size,mu=ns.mu),
  #                ncol=ns,dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)

  ## 新ZI版本
  exp.non=matrix(NA,nrow=nrow(non.coor),ncol=ns,
                 dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)
  exp.non=apply(exp.non, 2,
                function(y){z=simu_zi(family = type,subject.n = nrow(exp.non),
                                      zi.p=ns.p, size=ns.size,mu=ns.mu);z})


  all=cbind(non.coor[c('row','col')],exp.svg,exp.non)
  # rownames(all)<-all$num

  if(outlier==FALSE){
    end=all
  }
  if(outlier>=1 | outlier<=0){
    cat("# outlier parameter is wrong! \n")
    end=all
  }
  if(outlier>=0 & outlier <1){
    ind=sample(nrow(all),round(nrow(all)*outlier));
    out.para=5;
    # all[ind,-c(1:2)]<-matrix(rnbinom(n=(length(ind)*(ncol(all)-2)),
    #                                  size=se.size*out.para,
    #                                  mu=se.mu*out.para),
    #                          nrow = length(ind))

    ## 新 ZI版本
    all[ind,-c(1:2)]<-matrix(simu_zi(family = type,subject.n = (length(ind)*(ncol(all)-2)),
                                     zi.p = se.p/2,
                                     size=se.size*out.para,
                                     mu=se.mu*out.para),
                             nrow = length(ind))
    end=all
  }
  end
}



## 4 Gradient
## 原来的Gradient.zi2




#' Title Simulation data. Genes with gradient expression pattern.
#' @param lambda The spots detection efficiency.It should be between 0 and 1.
#' The default value is 0.5.
#' @param spots The number of spots.
#' @param prop  The proportion of streak area.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' @param a   The strength parameter of gradient pattern.
#' @param b  The strength parameter of gradient pattern.
#' @param grid  The strength parameter of gradient pattern.
#' @param grid.part  The proportion of gradient area.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
gradient.zi=function(lambda=0.5,spots,se,ns, ##se,ns是SVG和non-SVG数量
                     type='ZINB',ns.p,
                     ns.size=0.2,ns.mu,
                     a=1.5,b=0,grid=10,grid.part=0.5, ##     grid.part between 0 and 1
                     outlier=FALSE){   ##non-SVG 参数,并据此生成其他


  ## 1 spot generation
  win=ceiling(sqrt(spots/lambda))
  coor=spatstat.random::rpoispp(lambda = lambda,
                                win=spatstat.geom::owin(c(0,win),
                                                        c(0,win)))
  coor.s1=data.frame(#num=paste0('c-',1:length(coor$x)),
    row=round(coor$x),# row=round(coor$x,0),
    col=round(coor$y),# col=round(coor$y,0),
    row.names = paste0('c-',1:length(coor$x)))
  coor.s1=coor.s1[order(coor.s1$row),]  ## 按row排序

  # ## 2  spots  marks+SVG
  # max.x=max(coor.s1$row)
  # part=round(seq(0,max.x,length.out=grid+1))
  # exp.gradient=c()
  # for (i in 1:grid) {
  #   if(i<grid){
  #     coor.i=subset(coor.s1,row>=part[i] & row<part[i+1])
  #   } else{
  #     coor.i=subset(coor.s1,row>=part[i] & row<=part[i+1])
  #   }
  #   # exp.i=matrix(rnbinom(se*nrow(coor.i),size=a*i*ns.size,mu=a*i*ns.mu),
  #   #              ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
  #   # exp.i=matrix(rnbinom(se*nrow(coor.i),size=ns.size,mu=a*i*ns.mu),
  #   #              ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
  #
  #   ## 新ZI 版本
  #   exp.i=matrix(NA,nrow=nrow(coor.i),ncol=se,
  #                dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
  #   exp.i=apply(exp.i, 2,
  #               function(y){z=simu_zi(family = type,
  #                                     subject.n = nrow(exp.i),
  #                                     zi.p=ns.p*(grid-i+1)/grid,
  #                                     size=ns.size,mu=a*i*ns.mu);z})
  #
  #
  #   data.i=cbind(coor.i,exp.i)
  #   exp.gradient=rbind(exp.gradient,data.i)
  #
  # }


  ## 2  spots  marks+SVG   0712 可以调节mark area版本


  if  (grid.part>1| grid.part<0) {cat('grid.part is wrong. grid.part should be between 0 and 1.')}

  max.x=max(coor.s1$row)
  part=round(seq(0,max.x,length.out=grid+1))

  exp.gradient=c()
  for (i in ceiling(grid*(1-grid.part)):grid) {
    if(i<grid){
      coor.i=subset(coor.s1,row>=part[i] & row<part[i+1])
    } else{
      coor.i=subset(coor.s1,row>=part[i] & row<=part[i+1])
    }
    # exp.i=matrix(rnbinom(se*nrow(coor.i),size=a*i*ns.size,mu=a*i*ns.mu),
    #              ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
    # exp.i=matrix(rnbinom(se*nrow(coor.i),size=ns.size,mu=a*i*ns.mu),
    #              ncol=se,dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)

    ## 新ZI 版本
    exp.i=matrix(NA,nrow=nrow(coor.i),ncol=se,
                 dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)
    exp.i=apply(exp.i, 2,
                function(y){z=simu_zi(family = type,
                                      subject.n = nrow(exp.i),
                                      zi.p=ns.p*(grid-i+1)/grid,
                                      size=ns.size,  ## mu=a*i*ns.mu,  老版本
                                      mu=ns.mu*a*(i/grid));z})

    data.i=cbind(coor.i,exp.i)
    exp.gradient=rbind(exp.gradient,data.i)

  }

  coor.other=subset(coor.s1,row<part[ceiling(grid*(1-grid.part))])
  exp.other=matrix(NA,nrow=nrow(coor.other),ncol=se,
                   dimnames = list(NULL,paste0('se.',1:se)))## nb(r=size,p=prob)


  exp.other=apply(exp.other, 2,
                  function(y){z=simu_zi(family = type,subject.n = nrow(exp.other),
                                        zi.p = ns.p,size=ns.size,mu=ns.mu);z})
  exp.other=cbind(coor.other,exp.other)
  exp.gradient=rbind(exp.other,exp.gradient)

  ##non-SVG
  # exp.non=matrix(rnbinom(ns*nrow(coor.s1),size=ns.size,mu=ns.mu),
  #                ncol=ns,dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)

  ## 新ZI 版本

  exp.non=matrix(NA,nrow=nrow(coor.s1),ncol=ns,
                 dimnames = list(NULL,paste0('ns.',1:ns)))## nb(r=size,p=prob)
  exp.non=apply(exp.non, 2,
                function(y){z=simu_zi(family = type,subject.n = nrow(exp.non),
                                      zi.p = ns.p,size=ns.size,mu=ns.mu);z})

  all=cbind(exp.gradient,exp.non)
  # rownames(all)<-all$num
  # all=all[,-1]
  if(outlier==FALSE){
    end=all
  }
  if(outlier>=1 | outlier<=0){
    cat("# outlier parameter is wrong! \n")
    end=all
  }
  if(outlier>=0 & outlier <1){
    ind=sample(nrow(all),round(nrow(all)*outlier));
    out.para=5;
    # all[ind,-c(1:2)]<-matrix(rnbinom(n=(length(ind)*(ncol(all)-2)),
    #                                  size=a*grid*ns.size,
    #                                  mu=a*grid*ns.mu),
    #                          nrow = length(ind))

    all[ind,-c(1:2)]<-matrix(simu_zi(family = type,subject.n = (length(ind)*(ncol(all)-2)),
                                     zi.p = ns.p/2,
                                     size=a*grid*ns.size,
                                     mu=a*grid*ns.mu),
                             nrow = length(ind))

    end=all
  }
  end



}


