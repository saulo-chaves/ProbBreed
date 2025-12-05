##' @title Bayesian Probabilistic Selection Index (BPSI)
##'
##' @description
##' This function estimates the genotype's merit for multiple traits using the
##' probabilities of superior performance across environments.
##'
##' @param problist A list of object of class `probsup`, obtained from the [ProbBreed::prob_sup] function
##' @param increase Optional logical vector with size corresponding to the number of traits
##' of `problist`, in the same order.`TRUE` if the selection is for increasing the trait value, `FALSE` otherwise.
##' If not declared, `bpsi` will consider the information provided in [ProbBreed::prob_sup]
##' @param lambda A numeric representing the weight of each trait. Defaults to 1 (equal weights).
##' The trait with more economic interest should be greater.
##' @param int A numeric representing the selection intensity (between 0 and 1), considering the selection index.
##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##'
##' @return
##' The function returns an object of class `bpsi`, which contains two lists,
##' one with the BPSI- Bayesian Probabilistic Selection Index, and another with the original `data`-
##' with across-environments probabilities of superior performance for each trait.
##'
##' @details
##' \itemize{\item Bayesian Probabilistic Selection Index}
##'
##' \deqn{BPSI_i = \sum_{m=1}^{t} \frac{\gamma_{pt} -\gamma_{it} }{(1/\lambda_t)}}
##'
##' where \eqn{\gamma_p} is the probability of superior performance of the worst genotype for the trait \eqn{t},
##' \eqn{\gamma} is the probability of superior performance of genotype \eqn{i} for trait \eqn{t},
##'  \eqn{t} is the total number of traits evaluated,
##'   \eqn{\left(m = 1, 2, ..., t \right)},
##' and \eqn{\lambda} is the weight for each trait \eqn{t}.
##'
##' More details about the usage of `bpsi` can be found at \url{https://tiagobchagas.github.io/BPSI/}.
##'
##' @references
##' Chagas, J. T. B., Dias, K. O. G., Carneiro, V. Q., Oliveira, L. M. C., Nunes, N. X.,
##' Pereira Júnior, J. D., Carneiro, P. C. S., & Carneiro, J. E. S. (2025).
##' Bayesian probabilistic selection index in the selection of common bean families.
##' \emph{Crop Science}, 65(3). \doi{10.1002/CSC2.70072}
##'
##'
##'
##' @import ggplot2
##' @importFrom utils write.csv combn
##' @importFrom stats reshape median quantile na.exclude model.matrix aggregate
##' @importFrom rlang .data
##'
##' @seealso [ProbBreed::plot.bpsi()]
##'
##' @author José Tiago Barroso Chagas
##'
##' @export
##'
##' @examples
##' \donttest{
##'
##'
##'
##' mod = bayes_met(data = soy_pat,
##'                 gen = "gen",
##'                 loc = "env",
##'                 repl = NULL,
##'                 trait = "PH",
##'                 reg = NULL,
##'                 year = NULL,
##'                 res.het = TRUE,
##'                 iter = 2000, cores = 2, chain = 4)
##'
##'
##' mod2 = bayes_met(data = soy_pat,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl = NULL,
##'                  trait = "GY",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = TRUE,
##'                  iter = 2000, cores = 2, chain = 4)
##'
##' mod3 = bayes_met(data = soy_pat,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl =  NULL,
##'                  trait = "NDM",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = TRUE,
##'                  iter = 2000, cores = 2, chain = 4)
##'
##'
##'
##' models=list(mod,mod2,mod3)
##' names(models) <- c("PH","GY","NDM")
##' increase = c(FALSE,TRUE,FALSE)
##' names(increase) <- names(models)
##'
##' probs = list()
##' for (i in names(models)) {
##'   outs <- extr_outs(model = models[[i]],
##'                     probs = c(0.05, 0.95),
##'                     verbose = TRUE)
##'   probs[[i]] <- prob_sup(
##'     extr = outs,
##'     int = .2,
##'     increase = increase[[i]],
##'     save.df = FALSE,
##'     verbose = TRUE
##'   )
##'
##' }
##'
##' index = bpsi(
##'   problist = probs,
##'   increase = increase,
##'   int = 0.1,
##'   lambda = c(1, 2, 1),
##'   save.df = FALSE
##' )
##' }

bpsi = function(problist, increase = NULL, lambda = NULL, int, save.df = FALSE){


  # stopifnot("Please, provide a valid vector of selection ideotype" !=is.null(increase))
  stopifnot("Please, provide a list of objects from class probsup" = all(sapply(problist, inherits, "probsup")))
  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })

  geno = sapply(problist, function(x) {
    length(x$across$perfo$ID)
  })
  gen.min=which.min(geno)
  traits = sapply(problist, function(x)
    attr(x, "control")$trait)
  if(is.null(increase)){
    inc=sapply(problist, function(x) attr(x,"control")$increase)
  } else inc = increase

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

  if(is.null(lambda)) lambda = rep(1, times = length(traits)) else lambda=lambda

  if (all(inc)) {
    bpsi <- apply(bpsi, 2, function(x) {
      max(x) - x
    })
  } else{
    bpsi <- sapply(seq_along(bpsi), function(x) {
      if (inc[x] == TRUE) {
        max(bpsi[, x]) - bpsi[, x]
      } else{
        bpsi[, x] - min(bpsi[, x])
      }
    })
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
  # message(paste0('-> BPSI selected ',int*100, "% of genotypes"))
  attr(bpsi, "control")=int
  output = list(BPSI=bpsi,df=df)
  class(output)="bpsi"





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
  # message('Process completed!')

  return(output)


}


##' Plots for the `bpsi` object
##'
##' Build plots using the outputs stored in the `bpsi` object.
##'
##' @param x An object of class `bpsi`.
##' @param ... currently not used
##' @param category A string indicating which plot to build. There are currently two
##' types of visualizations. Set "Ranks" for bar plots along each trait and "BPSI" (default) for circular bar plots multitrait.
##' @method plot bpsi
##'
##' @references
##' Chagas, J. T. B., Dias, K. O. das G., Quintão Carneiro, V., de Oliveira, L. M. C., Nunes, N. X., Júnior, J. D. P., Carneiro, P. C. S., & Carneiro, J. E. de S. (2025).
##' Bayesian probabilistic selection index in the selection of common bean families.
##'  \emph{Crop Science}, 65(3).\doi{https://doi.org/10.1002/CSC2.70072}
##'
##' @author José Tiago Barroso Chagas
##'
##' @seealso [ProbBreed::bpsi]
##'
##' @import ggplot2
##' @importFrom stats reshape na.exclude
##' @importFrom rlang .data
##'
##' @rdname plot.bpsi
##'
##' @export
##'
##' @examples
##' \donttest{
##'
##'
##'
##' mod = bayes_met(data = soy_pat,
##'                 gen = "gen",
##'                 loc = "env",
##'                 repl = NULL,
##'                 trait = "PH",
##'                 reg = NULL,
##'                 year = NULL,
##'                 res.het = TRUE,
##'                 iter = 2000, cores = 2, chain = 4)
##'
##'
##' mod2 = bayes_met(data = soy_pat,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl = NULL,
##'                  trait = "GY",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = TRUE,
##'                  iter = 2000, cores = 2, chain = 4)
##'
##' mod3 = bayes_met(data = soy_pat,
##'                  gen = "gen",
##'                  loc = "env",
##'                  repl =  NULL,
##'                  trait = "NDM",
##'                  reg = NULL,
##'                  year = NULL,
##'                  res.het = TRUE,
##'                  iter = 2000, cores = 2, chain = 4)
##'
##'
##'
##' models=list(mod,mod2,mod3)
##' names(models) <- c("PH","GY","NDM")
##' increase = c(FALSE,TRUE,FALSE)
##' names(increase) <- names(models)
##'
##' probs = list()
##' for (i in names(models)) {
##'   outs <- extr_outs(model = models[[i]],
##'                     probs = c(0.05, 0.95),
##'                     verbose = TRUE)
##'   probs[[i]] <- prob_sup(
##'     extr = outs,
##'     int = .2,
##'     increase = increase[[i]],
##'     save.df = FALSE,
##'     verbose = TRUE
##'   )
##'
##' }
##'
##' index = bpsi(
##'   problist = probs,
##'   increase = increase,
##'   int = 0.1,
##'   lambda = c(1, 2, 1),
##'   save.df = FALSE
##' )
##'
##' plot(index, category = "BPSI")
##' plot(index, category = "Ranks")
##' }


plot.bpsi = function(x, ..., category = "BPSI"){
  # Namespaces
  requireNamespace('ggplot2')

  stopifnot("Object is not of class 'bpsi'" = class(x) == "bpsi")

  obj=x[["BPSI"]]

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
      geom_col(aes(x = .data[["gen"]], y = .data[["rank"]], fill = .data[["sel"]]), data =
                 obja) +
      facet_wrap( ~ .data[["trait"]], scales = "free_x") +
      theme(
        axis.text.x = element_text(
          size = 4,
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        panel.background = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 8, face = "bold")
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


##' Print an object of class `bpsi`
##'
##' Print a `bpsi` object in R console
##'
##' @param x An object of class `bpsi`
##' @param ... currently not used
##' @method print bpsi
##'
##' @seealso [ProbBreed::bpsi]
##'
##' @author José Tiago Barroso Chagas
##'
##'
##' @export
##'

  print.bpsi = function(x, ...){
    obj = x
    message("==> Considering an intensity of ", attr(obj[[1]], "control") *100,'%, here are the selected candidates:')
    obj=obj[[1]]
    obj$gen=rownames(obj)
    rownames(obj) = NULL
    print(obj)
  }
