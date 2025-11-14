## Function prob_sup
##
##' @title Bayesian Probabilistic Selection Index (BPSI)
##'
##' @description
##' This function estimates the genotypic merit for multiple traits using the
##' probabilities of superior performance across environments.
##'
##' @param problist A list of object of class `probsup`, obtained from the [ProbBreed::prob_sup] function
##' @param increase Logical vector of amount of traits used in the same order of problist.`TRUE` (default) if the selection is for increasing the trait value, `FALSE` otherwise.
##' @param lambda A numeric representing the weight of each trait, the default is 1, the trait with more
##' economic interest should be greater.
##' @param int A numeric representing the selection intensity
##' (between 0 and 1)
##'##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##' completed steps. Defaults to `FALSE`.
##'
##' @return The function returns an object of class `BPSI`, which contains two lists,
##' one with the `BPSI`- Bayesian Probabilistic Selection Index, and another with the original `data`-
##' with across-environments probabilities of superior performance for each trait.
##'
##'
##'
##'
##' @details
##' Probabilities provide the risk of recommending a selection candidate for a
##'  multitrait ideotype in a target population of environments.
##'   `BPSI` computes the probabilities of superior performance for multitrai selection:
##'
##' \itemize{\item Bayesian Probabilistic Selection Index}
##'
##'
##' \deqn{BPSI_i = \sum_{m=1}^{t} \frac{\gamma_{pt} -\gamma_{it} }{(1/\lambda_t)}}
##'
##' where \eqn{\gamma_p} is the probability of superior performance of the worst genotype for the trait \eqn{t},
##' \eqn{\gamma} is the probability of superior performance of genotype  \eqn{i} for trait \eqn{t},
##'  \eqn{t} is the total number of traits evaluated,
##'   \eqn{\left(m = 1, 2, ..., t \right)},
##' and \eqn{\lambda} is the weight for each trait \eqn{t}.
##'
##'
##' More details about the usage of `BPSI` can be found at \url{https://tiagobchagas.github.io/BPSI/}.
##'
##' @references
##'
##'##' Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães, L. J. M.,
##' Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging probability concepts
##' for cultivar recommendation in multi-environment trials. \emph{Theoretical and
##' Applied Genetics}, 133(2):443-455. \doi{10.1007/s00122-022-04041-y}
##'
##' Chaves, S. F. S., Krause, M. D., Dias, L. A. S., Garcia, A. A. F., & Dias, K. O. G. (2024).
##' ProbBreed: a novel tool for calculating the risk of cultivar recommendation in multienvironment trials.
##' \emph{G3 Genes|Genomes|Genetics}, 14(3). \doi{10.1093/g3journal/jkae013}
##'
##' Chagas, J. T. B., Dias, K. O. das G., Quintão Carneiro, V., de Oliveira, L. M. C., Nunes, N. X., Júnior, J. D. P., Carneiro, P. C. S., & Carneiro, J. E. de S. (2025).
##' Bayesian probabilistic selection index in the selection of common bean families.
##'  \emph{Crop Science, 65(3). \doi{10.1002/CSC2.70072}
##'
##'
##' @import ggplot2
##' @importFrom utils write.csv combn
##' @importFrom stats reshape median quantile na.exclude model.matrix aggregate
##' @importFrom rlang .data
##'
##' @seealso [ProbBreed::plot.BPSI]
##'
##' @export
##'
##' @examples
#' \donttest{
##'
##'
##' met_df=read.csv("https://raw.githubusercontent.com/tiagobchagas/BPSI/refs/heads/main/Data/blues_long.csv",header=T)
##'
##' mod = bayes_met(data = met_df,
##'                 gen = "gen",
##'                 loc = "env",
##'                 repl = NULL,
##'                 trait = "PH",
##'                 reg = NULL,
##'                 year = NULL,
##'                 res.het = T,
##'                 iter =400, cores = 4, chain = 4) #recommended run at least 4k iterations
##'
##'
##' mod2 = bayes_met(data = met_df,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl = NULL,
##'                  trait = "GY",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = T,
##'                  iter = 400, cores = 4, chain = 4) #recommended run at least 4k iterations
##'
##' mod3 = bayes_met(data = met_df,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl =  NULL,
##'                  trait = "NDM",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = T,
##'                  iter = 400, cores = 4, chain = 4) #recommended run at least 4k iterations
##'
##'
##'
##' models=list(mod,mod2,mod3)
##' names(models) <- c("PH","GY","NDM")
##' inc=c(FALSE,TRUE,FALSE)
##' names(inc) <- names(models)
##'
##' results <- lapply(names(models), function(model_name) {
##'   x <- models[[model_name]]  # actual model object
##'
##'   outs <- extr_outs(model = x,
##'                     probs = c(0.05, 0.95),
##'                     verbose = TRUE)
##'
##'   a <- prob_sup(extr = outs,
##'                 int = .2,
##'                 increase = inc[[model_name]],  # ← now model_name is a character!
##'                 save.df = FALSE,
##'                 verbose = TRUE)
##'
##'   return(a)
##' })
##' names(results) <- names(models)
##'
##'
##'
##'
##'
##'
##' bpsi=BPSI(problist=results,
##'           increase = c(FALSE,TRUE,FALSE),
##'           int = 0.1,
##'           lambda=c(1,2,1),
##'           save.df = F)
##'
##'
##'
##' plot(bpsi,category = "Ranks")
##'
##' plot(bpsi,category = "BPSI")
##'
##' df=print(bpsi)
##'
##'
#' }
##'



BPSI = function(problist,increase, lambda, int, save.df = FALSE){

  # stopifnot("Please, provide a valid vector of selection ideotype" !=is.null(increase))
  stopifnot("Please, provide a list of objects from class probsup" = all(sapply(problist, inherits, "probsup")))
  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })

  geno=sapply(problist,function(x){   length(x$across$perfo$ID)}  )
  gen.min=which.min(geno)
  traits=sapply(problist, function(x)  attr(x,"control")$trait  )
  inc=sapply(problist, function(x) attr(x,"control")$increase)
  traits=t(traits)
  names(traits) <- traits
  prob_l <- lapply(seq_along(problist), function(x){
    problist[[x]]$across$perfo[which(problist[[x]]$across$perfo$ID %in% problist[[gen.min]]$across$perfo$ID),]
  })
  names(prob_l) <- names(traits)


  df <- data.frame(matrix(NA,nrow = length(prob_l[[1]]$ID),
                          ncol = length(prob_l)))

  rownames(df) <- prob_l[[1]]$ID
  colnames(df) <- names(traits)
  rownames(df)=sort(row.names(df),decreasing = F)

  df[names(traits)] <- lapply(seq_along(prob_l), function(i){
    prob_l[[i]]$prob[match(rownames(df), prob_l[[i]]$ID)]
  })

  bpsi <- as.data.frame(apply(df, 2, function(x) {
    rank(-x, ties.method = "min", na.last = "keep")
  }))

  if(is.null(lambda)) lambda=1 else lambda=lambda
  if(is.null(increase)) increase=inc else increase
  if(all(increase==inc)) {
    bpsi <- apply(bpsi, 2, function(x) {
      max(x) - x})
  }else{
    bpsi <- sapply(seq_along(bpsi), function(x) {
      if(increase[x]==TRUE){
        max(bpsi[,x])-bpsi[,x]
      }else{
        bpsi[,x]-min(bpsi[,x]) } })
    colnames(bpsi) <- colnames(df)
    rownames(bpsi) <- row.names(df)
  }




  bpsi <- bpsi / (1/lambda[col(bpsi)])

  bpsi <- cbind(bpsi, rowSums(bpsi))
  colnames(bpsi)[ncol(bpsi)] <- "bpsi"

  bpsi <- as.data.frame(bpsi)
  gen_sel = bpsi[order(bpsi$bpsi,decreasing = T),]
  gen_i = round((length(gen_sel$bpsi) * int),0)

  selected=gen_sel[1:gen_i,]

  bpsi$sel=ifelse(rownames(bpsi) %in% rownames(selected) , "Selected","Not_Selected")
  bpsi=bpsi[order(bpsi$bpsi,decreasing = T),]
  bpsi$gen=rownames(bpsi)
  # if(verbose) message('-> Genotypes selected using BPSI')
  message(paste0('-> BPSI selected ',int*100, "% of genotypes"))
  attr(bpsi, "control")=int
  output = list(BPSI=bpsi,df=df)
  class(output)="BPSI"



  ## Save data frames in the work directory -----------------
  if(save.df){
    dir.create(path = paste0(getwd(),'/BPSI'))
    utils::write.csv(output$BPSI,
                     file = paste0(getwd(),'/BPSI/bpsi.csv'),
                     row.names = F)
    utils::write.csv(output$df,
                     file = paste0(getwd(),'/BPSI/prob_df.csv'),
                     row.names = F)
  }
  # Final output -----------
  # if(verbose) message('Process completed!')
  message('Process completed!')
  return(output)


}



#' Plots for the `BPSI` object
#'
#' Build plots using the outputs stored in the `BPSI` object.
#'
#'
#' @param problist A list of object of class `probsup`, obtained from the [ProbBreed::prob_sup] function
##' @param increase Logical vector of amount of traits used in the same order of problist.`TRUE` (default) if the selection is for increasing the trait value, `FALSE` otherwise.
##' @param lambda A numeric representing the weigth of each trait, the default is 1, the trait with more
##' economic interest should be greater.
##' @param int A numeric representing the selection intensity
##' (between 0 and 1)
##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##' @param verbose A logical value. If `TRUE`, the function will indicate the
##' completed steps. Defaults to `FALSE`.
##'
##' @return The function returns an object of class `BPSI`, which contains two lists,
##' one with the `BPSI`- Bayesian Probabilistic Selection Index, and another with the original `data`-
##' with across-environments probabilities for each trait.
##'
##' \itemize{\item Probability of superior performance}
##'
##'
##' \deqn{BPSI_i = \sum_{m=1}^{t} \frac{\gamma_{pt} -\gamma_{it} }{(1/\lambda_t)}}
##'
##' where \eqn{\gamma_p} is the probability of superior performance of the worst genotype for the trait \eqn{m},
##' \eqn{\gamma} is the probability of superior performance of genotype  \eqn{i} for trait \eqn{t},
##'  \eqn{t} is the total number of traits evaluated,
##'   \eqn{\left(m = 1, 2, ..., t \right)},
##' and \eqn{\lambda} is the weight for each trait \eqn{t}.
##'
##'
##' More details about the usage of `BPSI`, as well as the other function of
##' the `ProbBreed` package can be found at \url{https://saulo-chaves.github.io/ProbBreed_site/}.
##'
##' @references
##'
##'
##' Chagas, J. T. B., Dias, K. O. das G., Quintão Carneiro, V., de Oliveira, L. M. C., Nunes, N. X., Júnior, J. D. P., Carneiro, P. C. S., & Carneiro, J. E. de S. (2025).
##' Bayesian probabilistic selection index in the selection of common bean families.
##'  \emph{Crop Science}, 65(3).\doi{https://doi.org/10.1002/CSC2.70072}
#'
#'
#' @seealso  [ProbBreed::prob_sup]
#'
##' @import ggplot2
##' @importFrom stats reshape na.exclude
##' @importFrom rlang .data
##'
#' @rdname plot.probsup
#' @export
#'
##' @examples
#' \donttest{
##'
##'
##' met_df=read.csv("https://raw.githubusercontent.com/tiagobchagas/BPSI/refs/heads/main/Data/blues_long.csv",header=T)
##'
##' mod = bayes_met(data = met_df,
##'                 gen = "gen",
##'                 loc = "env",
##'                 repl = NULL,
##'                 trait = "PH",
##'                 reg = NULL,
##'                 year = NULL,
##'                 res.het = T,
##'                 iter =400, cores = 2, chain = 4) #recommended run at least 4k iterations and 4 cores
##'
##'
##' mod2 = bayes_met(data = met_df,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl = NULL,
##'                  trait = "GY",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = T,
##'                  iter = 400, cores = 2, chain = 4) #recommended run at least 4k iterations and 4 cores
##'
##' mod3 = bayes_met(data = met_df,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl =  NULL,
##'                  trait = "NDM",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = T,
##'                  iter = 400, cores = 2, chain = 4) #recommended run at least 4k iterations and 4 cores
##'
##'
##'
##' models=list(mod,mod2,mod3)
##' names(models) <- c("PH","GY","NDM")
##' inc=c(FALSE,TRUE,FALSE)
##' names(inc) <- names(models)
##'
##' results <- lapply(names(models), function(model_name) {
##'   x <- models[[model_name]]  # actual model object
##'
##'   outs <- extr_outs(model = x,
##'                     probs = c(0.05, 0.95),
##'                     verbose = TRUE)
##'
##'   a <- prob_sup(extr = outs,
##'                 int = .2,
##'                 increase = inc[[model_name]],  # ← now model_name is a character!
##'                 save.df = FALSE,
##'                 verbose = TRUE)
##'
##'   return(a)
##' })
##' names(results) <- names(models)
##'
##'
##'
##'
##'
##'
##' bpsi=BPSI(problist=results,
##'           increase = c(FALSE,TRUE,FALSE),
##'           int = 0.1,
##'           lambda=c(1,2,1),
##'           save.df = F)
##'
##'
##'
##' plot(bpsi,category = "Ranks")
##'
##' plot(bpsi,category = "BPSI")
##'
##' df=print(bpsi)
##'
##'
#' }
#'


plot.BPSI = function(BPSI_result, ..., category = "BPSI"){

  obj = BPSI_result


  # Namespaces
  requireNamespace('ggplot2')

  stopifnot("Object is not of class 'BPSI'" = class(obj) == "BPSI")

  obj=BPSI_result[[1]]
  #control = attr(obj, "control")

  #ord_gen = factor(obj$across$perfo$ID, levels = obj$across$perfo$ID)
  #retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]
  obj$gen=rownames(obj)

  obj$id <- 1:nrow(obj)
  obj$angle <- 90 - 360 * (obj$id - 0.5) / nrow(obj)
  obj$hjust <-  ifelse(obj$angle < -90, 1, 0)
  obj$angle = ifelse(obj$angle < -90, obj$angle + 180, obj$angle)


  traits <- colnames(obj)[!colnames(obj) %in% c("bpsi","sel","gen", "angle", "hjust","id")]
  obja=obj


  obja <- reshape(obj, direction = "long",varying =list(traits),
                  v.names = "rank",timevar = "trait",idvar="gen",times = traits )




  selected=obja[which(obja$sel %in% "Selected"),]





  # Rank plot --------------
  if(category == "Ranks"){



    ggplot() +
      geom_col( aes(x = .data[["gen"]],y =.data[["rank"]], fill=.data[["sel"]]),data=obja)+

      facet_wrap(~.data[["trait"]], scales = "free_x") +
      theme(
        axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 8,face = "bold")) +
      scale_fill_manual(
        values = c("Selected" = "blue3", "Not_Selected" = "grey"),
        breaks = c("Not_Selected", "Selected"),
        labels = c("Not selected", "Selected")
      ) +
      labs(
        x = "Genotypes",
        y = "Ranking of superior performance",
        fill = expression(bold(Pr(g %in% Omega)))
      )




  }


  # BPSI plot --------------
  else if(category == "BPSI"){



    ggplot(obj, aes(x = as.factor(.data[["id"]]), y =  .data[["bpsi"]], fill = .data[["sel"]])) +  # Reverse values
      geom_col() +
      scale_fill_manual(values = c("Selected" = "blue3",
                                   "Not_Selected" = "grey"),
                        breaks = c("Not_Selected","Selected"),
                        labels = c("Not selected","Selected")) +
      geom_point(data = data.frame(Sel = c("Not_Selected","Selected")),
                 aes(x = 0, y = 0, color = .data[["Sel"]]),
                 inherit.aes = FALSE, size = 4, alpha = 0) +
      scale_color_manual(values = c("Selected" = "blue3",
                                    "Not_Selected" = "grey"),
                         name = NULL,
                         guide = "none") +
      coord_polar(start = 0) +
      geom_text(aes(x = .data[["id"]], y = max(.data[["bpsi"]]) + 10,  # Keep original positioning
                    label = .data[["gen"]], angle = .data[["angle"]], hjust = .data[["hjust"]]),
                size = 2.8) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top") +
      labs(fill = NULL) +
      # Remove any axis labels that might show transformed values
      scale_y_continuous(labels = NULL)





  }

}


#' Print an object of class `BPSI`
#'
#' Print a `BPSI` object in R console
#'
#' @param x An object of class `BPSI`
#' @param ... currently not used
#' @method print BPSI
#'
#' @seealso [ProbBreed::BPSI]
#'
#' @export
#'

print.BPSI = function(BPSI_result, ...){
  obj = BPSI_result
  message("==> Considering an intensity of ", attr(obj[[1]], "control") *100,'%, here are the candidates selected:')
  obj=obj[[1]]
  obj$gen=rownames(obj)
  rownames(obj) = NULL
  print(obj)
}

