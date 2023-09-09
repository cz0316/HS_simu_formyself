

#' Title toGCO_2
#' @description
#' Generating tsv for scGCO.
#' @param data input, a data.frame, top 2 columns are row and col.
#' @return a .tsv file for scGCO.
#' @export
toGCO_2=function(data){
  counts=data[,-c(1:2)]
  coor=paste0(round(data$row,1),'x',round(data$col,1))
  # rownames(counts)=coor
  counts=cbind(coor,counts)
  counts
  
}
