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
##' @param omega A numeric representing the weight of each trait, the default is 1, the trait with more
##' economic interest should be greater.
##' @param int A numeric representing the selection intensity
##' (between 0 and 1)
##'##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##' @param verbose A logical value. If `TRUE`, the function will indicate the
##' completed steps. Defaults to `FALSE`.
##'
##' @return The function returns an object of class `BPSI`, which contains two lists,
##' one with the `BPSI`- Bayesian Probabilistic Selection Index, and another with the original `data`-
##' with across-environments probabilities for each trait.
##'
##'
##'
##'
##' @details
##' Probabilities provide the risk of recommending a selection candidate for a target
##' population of environments or for a specific environment. `prob_sup`
##' computes the probabilities of superior performance and the probabilities of superior stability:
##'
##' \itemize{\item Bayesian Probabilistic Selection Index}
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
##'
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
##' library(ProbBreed)
##' library(tidyverse)
##'
##' source("Data/bpsi_fun.R")
##' met_df=read.csv("https://raw.githubusercontent.com/tiagobchagas/BPSI/refs/heads/main/Data/blues_long.csv",header=T)
##'
##' head(met_df)
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
##' outs= vector("list",length(models))
##' names(outs) <- c("PH","GY","NDM")
##'
##' for (mods in 1:length(models)) {
##'   a= extr_outs(model = models[[mods]],
##'                probs = c(0.05, 0.95),
##'                verbose = TRUE)
##'   outs[[mods]]=a
##'   rm(a)
##' }
##'
##' results= vector("list",length(outs))
##' names(results) <- c("PH","GY","NDM")
##' inc=c(FALSE,TRUE,FALSE)
##'
##' for (mods in 1:length(outs)) {
##'   a = prob_sup(extr = outs[[mods]],
##'                int = .2,
##'                increase = inc[[mods]],
##'                save.df = FALSE,
##'                verbose = TRUE)
##'   results[[mods]]=a
##'   rm(a)
##' }
##' source("R/bpsi.R")
##'
##'
##' setdiff(results[[2]]$across$g_hpd$gen,results[[1]]$across$g_hpd$gen)
##'
##' length(results[[3]]$across$g_hpd$gen)
##'
##' results[[2]]$across$perfo <- results[[2]]$across$perfo[-which(results[[2]]$across$perfo$ID %in% c("G10", "G36", "G37")),]
##' results[[3]]$across$perfo <- results[[3]]$across$perfo[-which(results[[3]]$across$perfo$ID %in% c("G10", "G36", "G37")),]
##'
##'
##'
##' bpsi=BPSI(problist=results,
##'           increase = c(FALSE,TRUE,FALSE),
##'           int = 0.1,
##'           omega=c(1,2,1),
##'           save.df = F,
##'           verbose = T)
##'
##'
##'
##' plot(bpsi,category = "Ranks")
##'
##' plot(bpsi,category = "BPSI")
##'
##' df=print(bpsi)
##'
#' }
##'


BPSI = function(problist, increase = c(TRUE,FALSE,FALSE), omega, int, save.df = FALSE, verbose = FALSE){

  stopifnot("Please, provide a valid vector of selection direction" !=is.null(increase))
  stopifnot("Please, provide names to the list of models" !=is.null(names(problist)))
  stopifnot("Please, provide a list of objects from class probsup" =class(problist[[1]])=="probsup")
  stopifnot("Please provide a list with same number of genotypes for all traits" =
              all(dim(problist[[1]]$across$perfo) == dim(problist[[2]]$across$perfo),
                  dim(problist[[2]]$across$perfo) == dim(problist[[3]]$across$perfo)))
  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })


  df <- data.frame(matrix(nrow = length(problist[[1]]$across$perfo$ID),
                          ncol = length(problist)))
  rownames(df) <- problist[[1]]$across$perfo$ID
  colnames(df) <- names(problist)
  rownames(df)=sort(row.names(df),decreasing = F)
  # Preparation
  for (traits in 1:length(problist)) {
    tname=names(problist[traits])
    a=problist[[traits]]$across$perfo[order(problist[[traits]]$across$perfo$ID,decreasing =F),]
    df[paste0(tname)]=a$prob
    rm(a)
  }

  if(is.null(omega)) omega=1 else omega=omega

  bpsi <- apply(df, 2, function(x) {
    rank(-x, ties.method = "min", na.last = "keep")
  })


  bpsi <- apply(bpsi, 2, function(x) {
    max(x) - x})


  bpsi <- bpsi / (1/omega[col(bpsi)])


  bpsi <- cbind(bpsi, rowSums(bpsi))
  colnames(bpsi)[ncol(bpsi)] <- "bpsi"

  bpsi <- as.data.frame(bpsi)
  gen_sel = bpsi[order(bpsi$bpsi,decreasing = T),]

  gen_i = round((length(gen_sel$bpsi) * int),0)

  selected=gen_sel[1:gen_i,]

  bpsi$sel=ifelse(rownames(bpsi) %in% rownames(selected) , "Selected","Not_Selected")
  bpsi=bpsi[order(bpsi$bpsi,decreasing = T),]
  bpsi$gen=rownames(bpsi)
  if(verbose) message('-> Genotypes selected using BPSI')
  attr(bpsi, "control")=int
  output = list(BPSI=bpsi,df=df)
  class(output)="BPSI"

  library(tidyverse)





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
  if(verbose) message('Process completed!')
  return(output)


}



#' Plots for the `BPSI` object
#'
#' Build plots using the outputs stored in the `BPSI` object.
#'
#'
#' @param problist A list of object of class `probsup`, obtained from the [ProbBreed::prob_sup] function
##' @param increase Logical vector of amount of traits used in the same order of problist.`TRUE` (default) if the selection is for increasing the trait value, `FALSE` otherwise.
##' @param omega A numeric representing the weigth of each trait, the default is 1, the trait with more
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
#' @examples
#' \donttest{
##' source("Data/bpsi_fun.R")
##'
##'
##' load("res_PL_year.rda")
##' load("res_PH_year.rda")
##' load("res_GY_year.rda")
##'
##' source("data/bpsi_fun.R")
##'
##' models= vector("list",length(3))
##' models[[1]] = res_GY;
##' models[[2]] = res_PH
##' models[[3]] = res_PL
##'
##' bpsi=BPSI(problist=results,
##'  increase = c(TRUE,FALSE,FALSE),
##'  int = 0.1,
##'  omega=c(2,1,1),
##'  save.df = F,
##'  verbose = T)
##'
##'
##'
##' plot(bpsi,category = "Ranks")
##'
##' plot(bpsi,category = "BPSI")
##'
##' df=print(BPSI_soy)
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


  obj <- obj |>
    mutate(id = row_number(),
           angle = 90 - 360 * (id - 0.5) / n(),
           hjust = ifelse(angle < -90, 1, 0),
           angle = ifelse(angle < -90, angle + 180, angle))

  traits <- colnames(obj)[!colnames(obj) %in% c("bpsi","sel","gen", "angle", "hjust","id")]
  obja=obj

  library(tidyverse)
  library(ggrepel)
  library(gghighlight)
  obja=obj |> pivot_longer(paste(traits,sep = ":"),names_to = "trait",values_to = "rank")

  selected=obja[which(obja$sel %in% "Selected"),]





  # Rank plot --------------
  if(category == "Ranks"){


    # geom_bar(aes(x = as.factor(gen), y = max(rank) - rank, fill = "Selected"),
    #          data = subset(obja, gen %in% selected$gen),
    #          stat = "identity") +
    # geom_bar(aes(x = as.factor(gen), y = max(rank) - rank, fill = "Not_Selected"),
    #          data = subset(obja, !gen %in% selected$gen),
    #          stat = "identity") +
    ggplot() +
      geom_col( aes(x = gen,y =rank, fill=sel),data=obja)+

      # geom_hline(data = subset(obja, gen %in% selected$gen) %>% group_by(trait) %>% summarize(max_rank = max(rank)),
      #            aes(yintercept = max_rank), linetype = "dashed", color = "red", size = 0.5) +

      facet_wrap(~trait, scales = "free_x") +
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



    ggplot(obj, aes(x = as.factor(id), y =  bpsi, fill = sel)) +  # Reverse values
      geom_col() +
      scale_fill_manual(values = c("Selected" = "blue3",
                                   "Not_Selected" = "grey"),
                        breaks = c("Not_Selected","Selected"),
                        labels = c("Not selected","Selected")) +
      geom_point(data = data.frame(Sel = c("Not_Selected","Selected")),
                 aes(x = 0, y = 0, color = Sel),
                 inherit.aes = FALSE, size = 4, alpha = 0) +
      scale_color_manual(values = c("Selected" = "blue3",
                                    "Not_Selected" = "grey"),
                         name = NULL,
                         guide = "none") +
      coord_polar(start = 0) +
      geom_segment(aes(x = id, xend = id, y = max(bpsi) - bpsi, yend = max(bpsi) + 5),  # Adjusted
                   color = "grey100", linewidth = 0.3) +
      geom_text(aes(x = id, y = max(bpsi) + 10,  # Keep original positioning
                    label = gen, angle = angle, hjust = hjust),
                size = 2.8) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top") +
      labs(fill = NULL) +
      # Remove any axis labels that might show transformed values
      scale_y_continuous(labels = NULL)


    # ggplot() +
    #   geom_col(
    #     aes(y = bpsi, x = reorder(gen,bpsi), fill = gen),
    #     data = subset(obj, gen %in% selected$gen)
    #   )+geom_col(
    #     aes(y = bpsi, x = reorder(gen,bpsi)),
    #     fill = "grey",
    #     data = subset(obj, !gen %in% selected$gen)
    #   ) +geom_hline(yintercept = max(selected$bpsi))+
    #   scale_fill_manual(values = viridis::turbo(length(unique(selected$gen)))) +
    #   theme_minimal() +
    #   labs(x = "Generation", y = "BPSI", fill = "Generation") +
    #   geom_label_repel(
    #     data = subset(obj, gen %in% selected$gen),
    #     aes(x = gen, y = bpsi, label = gen),  # Specify x and y aesthetics
    #     fill = "lightblue", size = 2.25,
    #     box.padding = 0.1, max.overlaps = 20,label.padding = 0.1,
    #     label.size = 0.1,max.time = 1.0,
    #     nudge_x = 7.0,nudge_y = 7.0,force=10,force_pull = 10,max.iter = 20000
    #   ) +
    #   gghighlight(gen %in% selected$gen) +
    #   geom_hline(
    #     aes(yintercept = max(bpsi, na.rm = TRUE), linetype = "BPSI"),  # Handle NA values
    #     color = "black"
    #   ) +
    #   labs(
    #     x = "Genotypes", y = "BPSI", , linetype = "Selected"
    #   ) +guides(fill = "none")+
    #   theme(
    #     axis.text.x = element_blank(),
    #     legend.position = "top",
    #     axis.title.y = element_text(size = 12),
    #     axis.title.x = element_text(size = 12),
    #     axis.text.y = element_text(size = 12)
    #   ) +
    #   scale_y_reverse(n.breaks=10)





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
  message("==> Considering an intensity of ", attr(obj[[1]], "control") *100,'%, here are the selected candidates:')
  obj=obj[[1]]
  obj$gen=rownames(obj)
  rownames(obj) = NULL
  print(obj)
}

