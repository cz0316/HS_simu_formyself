## 1  pattern from png, spots from the Poisson process.

#' Title new.zi
#' @description
#' 20230907 new simulation, patterns from png. spots from the Poisson process.
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
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
#' @export
new.zi=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2, 
                se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                outlier = FALSE) {
  
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)), 
                       row = round(coor$x)+1, col = round(coor$y)+1, 
                       row.names = paste0("c_",1:length(coor$x)))
  
  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image) 
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1
  
  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1
  
  ## 
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)
  
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se, 
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark), 
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se, 
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  
  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, 
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.non), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  if (outlier == FALSE) {
    end = all
  } else  if (outlier >= 1 | outlier <= 0) {
    cat("# outlier parameter is wrong! \n")
    end = all
  } else  if (outlier >= 0 & outlier < 1) {
    ind = sample(nrow(all), round(nrow(all) * outlier))
    out.para = 5
    all[ind, -c(1:2)] <- matrix(simu_zi(family = type, subject.n = (length(ind) * 
                                                                      (ncol(all) - 2)), zi.p = se.p/2, size = se.size * 
                                          out.para, mu = se.mu * out.para), nrow = length(ind))
    end = all
  }
  
  
  
  end
}


## 2  pattern from png, spots from the sample() function process.

#' Title new.zi2
#' new simulation, patterns from png. spots from the sample() function process.
#' @param lambda The probability of the sample() function process.
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
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
#' @export
new.zi2=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2, 
                 se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                 outlier = FALSE) {
  
  # win = ceiling(sqrt(spots/lambda))
  win=sqrt(spots/lambda)
  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image) 
  re_img=imager::imresize(imager::grayscale(image) ,
                          scale =round(win/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1
  
  ## 合并，并区分marked area
  # coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  coor.s1=img_coor
  
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(row,col)
  coor.mark=coor.mark[sample(nrow(coor.mark),1.5*lambda*nrow(coor.mark)),]
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(row,col)
  coor.random=coor.random[sample(nrow(coor.random),lambda*nrow(coor.random)),]
  
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se, 
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark), 
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se, 
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  
  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, 
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.non), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  if (outlier == FALSE) {
    end = all
  } else  if (outlier >= 1 | outlier <= 0) {
    cat("# outlier parameter is wrong! \n")
    end = all
  } else  if (outlier >= 0 & outlier < 1) {
    ind = sample(nrow(all), round(nrow(all) * outlier))
    out.para = 5
    all[ind, -c(1:2)] <- matrix(simu_zi(family = type, subject.n = (length(ind) * 
                                                                      (ncol(all) - 2)), zi.p = se.p/2, size = se.size * 
                                          out.para, mu = se.mu * out.para), nrow = length(ind))
    end = all
  }
  
  
  
  end
}



## 3  pattern from png, spots from the sample() function process.non-svg is different
## non-svg generates from the randomization of the SVG.

#' Title new.zi3
#' new simulation, patterns from png. spots from the sample() function process.
#' @param lambda The probability of the sample() function process.
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
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
#' @export

new.zi3=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2, 
                 se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                 outlier = FALSE) {
  
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)), 
                       row = round(coor$x)+1, col = round(coor$y)+1, 
                       row.names = paste0("c_",1:length(coor$x)))
  
  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image) 
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1
  
  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1
  
  ## 
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)
  
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se, 
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark), 
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se, 
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  
  exp.svg = rbind(exp.mark, exp.random)
  
  mm=apply(exp.svg,2,function(z){return(sample(z,length(z)))})  ## randomization
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, 
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  
  exp.non=mm[,sample(ncol(mm),size =ns,replace = T)]
  colnames(exp.non)<-paste0("ns.", 1:ns)
  # exp.non=apply(exp.non,2,function(z){z=mm[,sample(length(mm),size =1)]})
  
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  if (outlier == FALSE) {
    end = all
  } else  if (outlier >= 1 | outlier <= 0) {
    cat("# outlier parameter is wrong! \n")
    end = all
  } else  if (outlier >= 0 & outlier < 1) {
    ind = sample(nrow(all), round(nrow(all) * outlier))
    out.para = 5
    all[ind, -c(1:2)] <- matrix(simu_zi(family = type, subject.n = (length(ind) * 
                                                                      (ncol(all) - 2)), zi.p = se.p/2, size = se.size * 
                                          out.para, mu = se.mu * out.para), nrow = length(ind))
    end = all
  }
  
  
  
  end
}





## 3  irregular pattern png, spots from the Poisson process.

#' Title irreg.zi
#' new simulation, patterns from irregular png. spots from the Poisson process.
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
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
#' @export

irreg.zi=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2, 
                  se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                  outlier = FALSE) {
  
  win1 = ceiling(sqrt(spots/lambda))
  win=win*sqrt(2)
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)), 
                       row = round(coor$x)+1, col = round(coor$y)+1, 
                       row.names = paste0("c_",1:length(coor$x)))
  
  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image) 
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win1)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1
  
  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1
  
  ## 
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)
  
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se, 
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark), 
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se, 
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  
  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, 
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.non), 
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  if (outlier == FALSE) {
    end = all
  } else  if (outlier >= 1 | outlier <= 0) {
    cat("# outlier parameter is wrong! \n")
    end = all
  } else  if (outlier >= 0 & outlier < 1) {
    ind = sample(nrow(all), round(nrow(all) * outlier))
    out.para = 5
    all[ind, -c(1:2)] <- matrix(simu_zi(family = type, subject.n = (length(ind) * 
                                                                      (ncol(all) - 2)), zi.p = se.p/2, size = se.size * 
                                          out.para, mu = se.mu * out.para), nrow = length(ind))
    end = all
  }
  
  
  
  end
}

## 4 pattern from png, zero-inflater+normal distribution, spots from the Poisson process.


#' Title new.norm
#' @description
#'  pattern from png, zero-inflater+normal distribution, 
#'  spots from the Poisson process.
#' @param lambda The spots detection efficiency.It should be between 0 and 1.
#' The default value is 0.5.
#' @param spots The number of spots.
#' @param prop  The proportion of streak area.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param se.sd  the sd of the SVG.
#' @param se.mu  the mean of the SVG.
#' @param ns.sd the sd of the non-SVG.
#' @param ns.mu the mean of the non-SVG.
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param outlier  A numeric value should be between 0 and 1.
#' @return  A data.frame of gene expression counts and coordinates.
#' @export
new.norm=function(lambda = 0.7, spots, se, ns, se.p, ns.p,
                  se.sd = 0.2, se.mu, ns.sd = 0.2, ns.mu, ptn,png_dir,
                  outlier = FALSE) {
  
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)), 
                       row = round(coor$x)+1, col = round(coor$y)+1, 
                       row.names = paste0("c_",1:length(coor$x)))
  
  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image) 
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1
  
  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1
  
  ## 
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)
  
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se, 
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_norm(subject.n = nrow(exp.mark), 
                  zi.p = se.p, sd = se.sd, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se, 
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_norm( subject.n = nrow(exp.random), 
                   zi.p = ns.p, sd = ns.sd, mu = ns.mu)
    z
  })
  
  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, 
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_norm(subject.n = nrow(exp.non), 
                  zi.p = ns.p, sd = ns.sd, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  if (outlier == FALSE) {
    end = all
  } else  if (outlier >= 1 | outlier <= 0) {
    cat("# outlier parameter is wrong! \n")
    end = all
  } else  if (outlier >= 0 & outlier < 1) {
    ind = sample(nrow(all), round(nrow(all) * outlier))
    out.para = 5
    all[ind, -c(1:2)] <- matrix(simu_zi(family = type, subject.n = (length(ind) * 
                                                                      (ncol(all) - 2)), zi.p = se.p/2, size = se.size * 
                                          out.para, mu = se.mu * out.para), nrow = length(ind))
    end = all
  }
  
  
  
  end
}





