## Function recod
##
##' @title Factor recoding
##'
##' @description
##' The functions of `ProbBreed` works with alphanumeric factors (for e.g., "G01").
##' `recod` was built to perform the recoding of a factor's levels in the dataset. It
##' adds a new column with the newly coded levels, keeping the original names for
##' reference.
##'
##' @param data  A data frame containing the observations
##' @param name.fact A string. The name of the column containing the factor to be recoded
##' @param cod A string. The letter that will appear before the numbers in the new code. For e.g.: "G01"
##'
##' @return The function returns the same data frame, but with a new column
##' with the prefix "cod." containing the recoded levels
##'
##' @export
##'
##' @examples
##' toy = data.frame(
##'             gen = rep(1:5, each =4),
##'             env = rep(1:4, times = 5),
##'             y = rnorm(20, mean = 20, sd = 5)
##'         )
##' new.toy = recod(toy, name.fact = 'gen', cod = 'G')
##' toy

recod = function(data, name.fact, cod){
  stopifnot("cod must be a string" = is.character(cod))
  stopifnot("name.fact must be a string" = is.character(name.fact))

  num = dim(as.matrix(unique(data[,name.fact])))[1]
  data$ord = 1:nrow(data)

  a = data.frame(
    "JJJJ" = unique(data[,name.fact]),
    'XYWFX' = paste0(cod, sprintf(paste0('%0', nchar(num),'d'), seq(1:num)))
  )

  colnames(a)[which(colnames(a) == 'JJJJ')] = name.fact

  data = merge(data, a, by = name.fact, sort=F)
  colnames(data)[which(colnames(data) == "XYWFX")] = paste("cod",name.fact,sep = '.')

  data = data[order(data$ord),]
  data = data[,-which(colnames(data)=='ord')]
  rownames(data) = NULL

  return(data)

}
