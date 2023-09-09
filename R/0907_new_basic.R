

#' Title to_da
#' @description
#' Transformation of p*n matrix to a data frame with 3 columns (row,col,exp).
#' 
#' @param matrix_data p*n matrix. From a .png file.
#'
#' @return A data frame with 3 columns (row,col,exp).
#' @export
to_dt=function(matrix_data){
  n_rows <- nrow(matrix_data)
  n_cols <- ncol(matrix_data)
  data_frame <- data.frame(
    row=rep(1:n_rows,times=n_cols),
    col=rep(1:n_cols,each=n_rows),
    value = as.vector(matrix_data)  # as.vector是按列来的拉长的
  )
  data_frame
}



#' Title simu_norm
#' @description
#' A intermediate step in the new.norm.
#' 
#' @param subject.n sample size
#' @param zi.p The zero proportion of the zero generation process.
#' if  zi.p=0, it generates data from normal distribution.
#' @param mu the mean of the normal distribution.
#' @param sd the sd of the naromal distribution.
#' @return simulation data.
#' @export
simu_norm=function (subject.n, zi.p = 0.5, mu = 0.5, sd = 0.25) {
  Y = rep(NA, subject.n)
  ind.mix <- rbinom(length(Y), 1, zi.p)
  Y[which(ind.mix != 0)] <- 0
  Y[which(ind.mix == 0)] = round(abs(rnorm(n = sum(ind.mix == 0), 
                                           mean=mu,sd=sd)))
  return(Y)
}


