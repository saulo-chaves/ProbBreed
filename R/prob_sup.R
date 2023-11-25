## Function prob_sup
##
##' @title Probabilities of superior performance and stability
##'
##' @description
##' This function estimates the probabilities of superior performance and stability
##' across environments (`marginal` output). It also computes the probabilities
##' of superior performance within environments (`conditional` output).
##'
##' @param data A data frame containing the phenotypic data
##' @param trait,gen,loc A string. The name of the columns that correspond to
##' the trait, genotype and location information, respectively. If
##' the environment is a combination of other factors (for instance, location-year),
##' the name of the column that contains this information must be attributed to `loc`.
##' @param reg A string or NULL. If the dataset has information about regions,
##' `reg` will be a string with the name of the column that corresponds to the
##' region information. Otherwise, `reg = NULL` (default).
##' @param year A string or NULL. If the data set has information about time-related
##' environmental factors (years, seasons...), `year` will be a string with the
##' name of the column that corresponds to the time information. Otherwise, `year = NULL` (default).
##' @param mod.output An object from the [extr_outs()] function
##' @param int A number representing the selection intensity
##' (between 0 and 1)
##' @param increase Logical. Indicates the direction of the selection.
##' `TRUE` (default) for increasing the trait value, `FALSE` otherwise.
##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##' @param interactive Logical. Should ggplots be converted into interactive plots?
##' If `TRUE`, the function loads the `plotly` package and uses the [plotly::ggplotly()]
##' command.
##' @param verbose A logical value. If `TRUE`, the function will indicate the
##' completed steps. Defaults to `FALSE`.
##'
##' @return The function returns two lists, one with the `marginal` probabilities, and
##' another with the `conditional` probabilities.
##'
##' The `marginal` list has:
##' \itemize{
##' \item \code{df} : A list of data frames containing the calculated probabilities:
##' \itemize{
##' \item \code{perfo}: the probabilities of superior performance.
##' \item \code{pair_perfo}: the pairwise probabilities of superior performance.
##' \item \code{stabi}: the probabilities of superior stability. Can be `stabi_gl`,
##' `stabi_gm` (when `reg` is not `NULL`) or `stabi_gt` (when `year` is not `NULL`).
##' \item \code{pair_stabi}: the pairwise probabilities of superior stability.
##' Can be `pair_stabi_gl`, `pair_stabi_gm` (when `reg` is not `NULL`) or
##' `pair_stabi_gt` (when `year` is not `NULL`).
##' \item \code{joint_prob}: the joint probabilities of superior performance and stability.
##' }
##' \item \code{plot} : A list of ggplots illustrating the outputs:
##' \itemize{
##' \item \code{g_hpd}: a caterpillar plot representing the marginal genotypic value of
##' each genotype, and their respective highest posterior density interval (95% represented by the
##' thick line, and 97.5% represented by the thin line).
##' \item \code{perfo}: a bar plot illustrating the probabilities of superior performance
##' \item \code{pair_perfo}: a heatmap representing the pairwise probability of superior
##' performance (the probability of genotypes at the \emph{x}-axis being superior
##' to those on the \emph{y}-axis).
##' \item \code{stabi}: a bar plot with the probabilities of superior stability.
##' Different plots are generated for `stabi_gl`, `stabi_gm` and `stabi_gt` if `reg`
##' or/and `year` are not `NULL`.
##' \item \code{pair_stabi}: a heatmap with the pairwise probabilities of superior stability.
##' Different plots are generated for `stabi_gl`, `stabi_gm` and `stabi_gt` if `reg`
##' or/and `year` are not `NULL`. This plot represents
##' the probability of genotypes at the \emph{x}-axis being superior
##' to those on \emph{y}-axis.
##' \item \code{joint_prob}: a plot with the probabilities of superior performance,
##' probabilities of superior stability and the joint probabilities of superior
##' performance and stability.
##' }
##' }
##'
##' The `conditional` list has:
##' \itemize{
##' \item \code{df} : A list with:
##' \itemize{
##' \item \code{prob}: data frames containing the probabilities of superior performance
##' within environments. Can be `prob_loc`, `prob_reg` (if `reg` is not `NULL`), and
##' `prob_year` (if `year` is not `NULL`).
##' \item \code{pwprob}: lists with the pairwise probabilities of superior performance
##' within environments. Can be `pwprob_loc`, `pwprob_reg` (if `reg` is not `NULL`), and
##' `pwprob_year` (if `year` is not `NULL`).
##' }
##' \item \code{plot} : A list with:
##' \itemize{
##' \item \code{prob}: heatmaps with the probabilities of superior performance within
##' environments. Can be `prob_loc`, `prob_reg` (if `reg` is not `NULL`), and
##' `prob_year` (if `year` is not `NULL`).
##' \item \code{pwprob}: a list of heatmaps representing the pairwise probability of superior
##' performance within environments. Can be `pwprob_loc`, `pwprob_reg` (if `reg` is not `NULL`), and
##' `pwprob_year` (if `year` is not `NULL`). The interpretation is the same as in the
##' `pair_perfo` in the `marginal` list: the probability of genotypes at the \emph{x}-axis being superior
##' to those on \emph{y}-axis.
##' }
##' }
##'
##' @details
##' Probabilities provide the risk of recommending a selection candidate for a target
##' population of environments or for a specific environment. The function `prob_sup()`
##' computes the probabilities of superior performance and the probabilities of superior stability:
##'
##' \itemize{\item Probability of superior performance}
##'
##' Let \eqn{\Omega} represent the subset of selected genotypes based on their
##' performance across environments. A given genotype \eqn{j} will belong to \eqn{\Omega}
##' if its genotypic marginal value (\eqn{\hat{g}_j}) is high or low enough compared to
##' its peers. `prob_sup()` leverages the Monte Carlo discretized sampling
##' from the posterior distribution to emulate the occurrence of \eqn{S} trials. Then,
##' the probability of the \eqn{j^{th}} genotype belonging to \eqn{\Omega} is the
##' ratio of success (\eqn{\hat{g}_j \in \Omega}) events and the total number of sampled events,
##' as follows:
##'
##' \deqn{Pr(\hat{g}_j \in \Omega \vert y) = \frac{1}{S}\sum_{s=1}^S{I(\hat{g}_j^{(s)} \in \Omega \vert y)}}
##'
##' where \eqn{S} is the total number of samples (\eqn{s = 1, 2, ..., S}),
##' and \eqn{I(g_j^{(s)} \in \Omega \vert y)} is an indicator variable that can assume
##' two values: (1) if \eqn{\hat{g}_j^{(s)} \in \Omega} in the \eqn{s^{th}} sample,
##' and (0) otherwise. \eqn{S} is conditioned to the number of iterations and chains
##' previously set at [ProbBreed::bayes_met()].
##'
##' Similarly, the conditional probability of superior performance can be applied to
##' individual environments. Let \eqn{\Omega_k} represent the subset of superior
##' genotypes in the \eqn{k^{th}} environment, so that the probability of the
##' \eqn{j^{th} \in \Omega_k} can calculated as follows:
##'
##' \deqn{Pr(\hat{g}_{jk} \in \Omega_k \vert y) = \frac{1}{S} \sum_{s=1}^S I(\hat{g}_{jk}^{(s)} \in \Omega_k \vert y)}
##'
##' where \eqn{I(\hat{g}_{jk}^{(s)} \in \Omega_k \vert y)} is an indicator variable
##' mapping success (1) if \eqn{\hat{g}_{jk}^{(s)}} exists in \eqn{\Omega_k}, and
##' failure (0) otherwise, and \eqn{\hat{g}_{jk}^{(s)} = \hat{g}_j^{(s)} + \widehat{ge}_{jk}^{(s)}}.
##' Note that when computing conditional probabilities (i.e., conditional to the
##' \eqn{k^{th}} environment or mega-environment), we are accounting for
##' the interaction of the \eqn{j^{th}} genotype with the \eqn{k^{th}}
##' environment.
##'
##' The pairwise probabilities of superior performance can also be calculated across
##' or within environments. This metric assesses the probability of the \eqn{j^{th}}
##' genotype being superior to another experimental genotype or a commercial check.
##' The calculations are as follows, across and within environments, respectively:
##'
##' \deqn{Pr(\hat{g}_{j} > \hat{g}_{j^\prime} \vert y) = \frac{1}{S} \sum_{s=1}^S I(\hat{g}_{j}^{(s)} > \hat{g}_{j^\prime}^{(s)} \vert y)}
##'
##' or
##'
##' \deqn{Pr(\hat{g}_{jk} > \hat{g}_{j^\prime k} \vert y) = \frac{1}{S} \sum_{s=1}^S I(\hat{g}_{jk}^{(s)} > \hat{g}_{j^\prime k}^{(s)} \vert y)}
##'
##' These equations are set for when the selection direction is positive. If
##' `increase = FALSE`, \eqn{>} is simply switched by \eqn{<}.
##'
##'
##' \itemize{\item Probability of superior stability}
##'
##' Probabilities of superior performance highlight experimental genotypes with
##' high agronomic stability. For ecological stability (invariance), the probability
##' of superior stability is the more adequate. Making a direct analogy with the
##' method of Shukla (1972), a stable genotype is the one that has a low variance
##' of the GEI (genotype-by-environment interaction) effects \eqn{[var(\widehat{ge})]}.
##' Using the same probability principles previously described, the probability
##' of superior stability is given as follows:
##'
##' \deqn{Pr[var(\widehat{ge}_{jk}) \in \Omega \vert y] = \frac{1}{S} \sum_{s=1}^S I[var(\widehat{ge}_{jk}^{(s)}) \in \Omega \vert y]}
##'
##' where \eqn{I[var(\widehat{ge}_{jk}^{(s)}) \in \Omega \vert y]} indicates if
##' \eqn{var(\widehat{ge}_{jk}^{(s)})} exists in \eqn{\Omega} (1) or not (0).
##' Pairwise probabilities of superior stability are also possible in this context:
##'
##' \deqn{Pr[var(\widehat{ge}_{jk}) < var(\widehat{ge}_{j^\prime k}) \vert y] = \frac{1}{S} \sum_{s=1}^S I[var(\widehat{ge}_{jk})^{(s)} < var(\widehat{ge}_{j^\prime k})^{(s)} \vert y]}
##'
##' Note that \eqn{j} will be superior to \eqn{j^\prime} if it has a \strong{lower}
##' variance of the genotype-by-environment interaction effect. This is true regardless
##' if `increase` is set to `TRUE` or `FALSE`.
##'
##' The joint probability independent events is the product of the individual probabilities.
##' The estimated genotypic main effects and the variances of GEI effects are independent
##' by design, thus the joint probability of superior performance and stability as follows:
##'
##' \deqn{Pr[\hat{g}_j \in \Omega \cap var(\widehat{ge}_{jk}) \in \Omega] = Pr(\hat{g}_j \in \Omega) \times Pr[var(\widehat{ge}_{jk}) \in \Omega]}
##'
##' The estimation of these probabilities are strictly related to some key questions that
##' constantly arises in plant breeding:
##' \itemize{
##' \item \strong{What is the risk of recommending a selection candidate for a target population of environments?}
##' \item \strong{What is the probability of a given selection candidate having good performance if
##' recommended to a target population of environments? And for a specific environment?}
##' \item \strong{What is the probability of a given selection candidate having better performance
##' than a cultivar check in the target population of environments? And in specific environments?}
##' \item \strong{How probable is it that a given selection candidate performs similarly across environments?}
##' \item \strong{What are the chances that a given selection candidate is more stable
##' than a cultivar check in the target population of environments?}
##' \item \strong{What is the probability that a given selection candidate having a
##' superior and invariable performance across environments?}
##' }
##'
##' More details about the usage of `prob_sup`, as well as the other function of
##' the `ProbBreed` package can be found at \url{https://saulo-chaves.github.io/ProbBreed_site/}.
##'
##' @references
##'
##' Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães, L. J. M.,
##' Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging probability concepts
##' for cultivar recommendation in multi-environment trials. \emph{Theoretical and
##' Applied Genetics}, 133(2):443-455. \doi{10.1007/s00122-022-04041-y}
##'
##' Shukla, G. K. (1972) Some statistical aspects of partioning genotype environmental
##' componentes of variability. \emph{Heredity}, 29:237-245. \doi{10.1038/hdy.1972.87}
##'
##'
##' @import ggplot2
##' @importFrom utils write.csv combn
##' @importFrom stats reshape median quantile na.exclude model.matrix aggregate
##' @importFrom rlang .data
##'
##' @export
##'
##' @examples
##' \donttest{
##' mod = bayes_met(data = maize,
##'                 gen = "Hybrid",
##'                 loc = "Location",
##'                 repl = c("Rep", "Block"),
##'                 year = NULL,
##'                 reg = 'Region',
##'                 res.het = FALSE,
##'                 trait = 'GY',
##'                 iter = 6000, cores = 4, chains = 4)
##'
##' outs = extr_outs(data = maize, trait = "GY", model = mod,
##'                  probs = c(0.05, 0.95),
##'                  check.stan.diag = TRUE,
##'                  verbose = TRUE)
##'
##' results = prob_sup(data = maize,
##'                    trait = "GY",
##'                    gen = "Hybrid",
##'                    loc = "Location",
##'                    reg = 'Region',
##'                    year = NULL,
##'                    mod.output = outs,
##'                    int = .2,
##'                    increase = TRUE,
##'                    save.df = FALSE,
##'                    interactive = FALSE,
##'                    verbose = FALSE)
##' }
##'

prob_sup = function(data, trait, gen, loc, reg = NULL, year = NULL, mod.output, int,
                    increase = TRUE, save.df = FALSE, interactive = FALSE,
                    verbose = FALSE){

  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })

  # Namespaces
  requireNamespace('ggplot2')

  # Preparation
  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = mod.output

  # Transforming again
  names(mod$post)[which(names(mod$post) == 'location')] = "l"
  names(mod$post)[which(names(mod$post) == 'genotype')] = "g"
  names(mod$post)[which(names(mod$post) == 'gen.loc')] = "gl"
  if("replicate" %in% names(mod$post)){names(mod$post)[which(names(mod$post) == 'replicate')] = "r"}
  if("block" %in% names(mod$post)){names(mod$post)[which(names(mod$post) == 'block')] = "b"}
  if("region" %in% names(mod$post)){
    names(mod$post)[which(names(mod$post) == 'region')] = "m"
    names(mod$post)[which(names(mod$post) == 'gen.reg')] = "gm"
  }
  if("year" %in% names(mod$post)){
    names(mod$post)[which(names(mod$post) == 'year')] = "t"
    names(mod$post)[which(names(mod$post) == 'gen.year')] = "gt"
  }

  # Conditions
  if(!all(grepl('[A-Za-z]', data[, gen]))){data[,gen] = paste("G", data[,gen], sep = "_")}
  if(!all(grepl('[A-Za-z]', data[, loc]))){data[,loc] = paste("E", data[,loc], sep = "_")}

  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.loc = levels(factor(data[,loc]))
  num.loc = nlevels(factor(data[,loc]))
  num.sim = nrow(mod$post$g)

  if(!is.null(year)) # With year info ------------------
  {
    if(!all(grepl('[A-Za-z]', data[, year]))){data[,year] = paste("Y", data[, year], sep = "_")}

    name.year = levels(factor(data[,year]))
    num.year = nlevels(factor(data[,year]))

    colnames(mod$post$g) = paste0(name.gen, '_@#')
    colnames(mod$post$gt) = paste('Gen',rep(name.gen,  times = num.year),
                                  'year',rep(name.year,  each = num.gen), sep = '_@#')

    if(!is.null(reg))  # With region info --------------------
    {
      if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

      name.reg = levels(factor(data[,reg]))
      num.reg = nlevels(factor(data[,reg]))
      colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                    'Reg',rep(name.reg,  each = num.gen), sep = '_@#')
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.loc),
                                    "loc", rep(name.loc,  each = num.gen), sep = '_@#')

      aux = unique(data[,c(gen,loc,year,reg)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])
      Z3 = stats::model.matrix(~-1 + aux[,gen]:aux[,year])
      Z4 = stats::model.matrix(~-1 + aux[,gen]:aux[,reg])

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        UP = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        up = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        DOWN = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        down = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g, y = reorder(.data$gen, .data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#33a02c')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------

      if(increase){
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      } else {
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      }

      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g) +
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90)) # + ylim(0,1)

      if(verbose) message('1. Probability of superior performance estimated')

      ord_gen = factor(prob_g$ID, levels = prob_g$ID)

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))

      if(increase){
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
          }
        }
      } else {
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
          }
        }
      }

      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_@#'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      if(increase){
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      } else {
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      }

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gl[,grep(paste0(x,'$'), colnames(staprob_gl))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gl = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gl = prob_gl[order(prob_gl$prob, decreasing = T),]
      prob_gl.plot = ggplot(prob_gl)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('3. Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gl1 = pwsprob_gl[match(prob_gl$ID, rownames(pwsprob_gl)),
                               match(prob_gl$ID, rownames(pwsprob_gl))]
      pwsprob_gl1[upper.tri(pwsprob_gl1, diag = T)] = NA
      pwsprob_gl1 = stats::reshape(
        data.frame(pwsprob_gl1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl1), times = colnames(pwsprob_gl1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gl.plot = ggplot(pwsprob_gl1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gl[x]) < var(gl[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gm),'_@#Reg'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gm[,grep(paste0(x,'$'), colnames(staprob_gm))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gm = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gm = prob_gm[order(prob_gm$prob, decreasing = T),]
      prob_gm.plot = ggplot(prob_gm)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('5. Probability of superior stability (GM) estimated')


      ## Pairwise probability of superior stability - Region -----------------
      pwsprob_gm = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gm)){
        for (j in colnames(pwsprob_gm)) {
          pwsprob_gm[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gm1 = pwsprob_gm[match(prob_gm$ID, rownames(pwsprob_gm)),
                               match(prob_gm$ID, rownames(pwsprob_gm))]
      pwsprob_gm1[upper.tri(pwsprob_gm1, diag = T)] = NA
      pwsprob_gm1 = stats::reshape(
        data.frame(pwsprob_gm1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gm1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm1), times = colnames(pwsprob_gm1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gm.plot = ggplot(pwsprob_gm1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gm[x]) < var(gm[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gm[upper.tri(pwsprob_gm, diag = T)] = NA
      pwsprob_gm = stats::reshape(
        data.frame(pwsprob_gm),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gm))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm), times = colnames(pwsprob_gm),
        new.row.names = 1:length(c(pwsprob_gm)), v.names = 'prob'
      )
      pwsprob_gm = stats::na.exclude(pwsprob_gm[order(pwsprob_gm$x),])

      if(verbose) message('6. Pairwise probability of superior stability (GM) estimated')

      ## Probability of superior stability - year --------------
      staprob_gt = mod$post$gt
      colnames(staprob_gt) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gt),'_@#year'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gt[,grep(paste0(x,'$'), colnames(staprob_gt))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gt = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gt = prob_gt[order(prob_gt$prob, decreasing = T),]
      prob_gt.plot = ggplot(prob_gt)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('7. Probability of superior stability (GT) estimated')


      ## Pairwise probability of superior stability - year -----------------
      pwsprob_gt = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gt)){
        for (j in colnames(pwsprob_gt)) {
          pwsprob_gt[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gt1 = pwsprob_gt[match(prob_gt$ID, rownames(pwsprob_gt)),
                               match(prob_gt$ID, rownames(pwsprob_gt))]
      pwsprob_gt1[upper.tri(pwsprob_gt1, diag = T)] = NA
      pwsprob_gt1 = stats::reshape(
        data.frame(pwsprob_gt1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gt1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gt1), times = colnames(pwsprob_gt1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gt.plot = ggplot(pwsprob_gt1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gt[x]) < var(gt[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gt[upper.tri(pwsprob_gt, diag = T)] = NA
      pwsprob_gt = stats::reshape(
        data.frame(pwsprob_gt),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gt))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gt), times = colnames(pwsprob_gt),
        new.row.names = 1:length(c(pwsprob_gt)), v.names = 'prob'
      )
      pwsprob_gt = stats::na.exclude(pwsprob_gt[order(pwsprob_gt$x),])

      if(verbose) message('8. Pairwise probability of superior stability (GT) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gm, by = 'ID'),
                     merge(prob_g, prob_gt, by = 'ID'))

      j_prob$joint = j_prob[,2] * j_prob[,3]
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Location', 'Region', 'Year'), each = num.gen)
      j_prob = stats::reshape(j_prob, direction = 'long', varying = list(2:4),
                              times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]

      j_prob = j_prob[j_prob$category == "Joint",]

      retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]

      j_prob.plot = ggplot(cbind(j_prob, V4 = paste(j_prob$ID, j_prob$level, sep = '@#_')),
                           aes(x = stats::reorder(.data$V4, -.data$prob), y = .data$prob)) +
        facet_wrap(.~.data$level, ncol = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_x_discrete(labels = retrieve) +
        geom_segment(aes(x = reorder(.data$V4, -.data$prob), xend = .data$V4,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = "Genotype", y = "Joint probability")

      if(verbose) message('9. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------------

      if(interactive){
        g_hpd = plotly::ggplotly(g_hpd)
        prob_g.plot = plotly::ggplotly(prob_g.plot)
        prob_gl.plot = plotly::ggplotly(prob_gl.plot)
        prob_gm.plot = plotly::ggplotly(prob_gm.plot)
        prob_gt.plot = plotly::ggplotly(prob_gt.plot)
        j_prob.plot = plotly::ggplotly(j_prob.plot)
        pwsprob_g.plot = plotly::ggplotly(
          ggplot(data = pwsprob_g1, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gm.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gm, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gm(x)) < var(gm(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gt.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gt, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gt(x)) < var(gt(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------------
      marg_prob = list(
        df = list(
          perfo = prob_g,
          pair_perfo = pwsprob_g,
          stabi_gl = prob_gl,
          pair_stabi_gl = pwsprob_gl,
          stabi_gm = prob_gm,
          pair_stabi_gm = pwsprob_gm,
          stabi_gt = prob_gt,
          pair_stabi_gt = pwsprob_gt,
          joint_prob = j_prob
        ),
        plots = list(
          g_hpd = g_hpd,
          perfo = prob_g.plot,
          pair_perfo = pwsprob_g.plot,
          stabi_gl = prob_gl.plot,
          pair_stabi_gl = pwsprob_gl.plot,
          stabi_gm = prob_gm.plot,
          pair_stabi_gm = pwsprob_gm.plot,
          stabi_gt = prob_gt.plot,
          pair_stabi_gt = pwsprob_gt.plot,
          joint_prob = j_prob.plot
        )
      )

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          utils::write.csv(marg_prob$df[[i]],
                           file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                           row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = apply(mod$post$g, 1, function(x){
        Z1 %*% as.matrix(x)
      }) +
        apply(mod$post$gl, 1, function(x){
          Z2 %*% as.matrix(x)
        }) +
        apply(mod$post$gt, 1, function(x){
          Z3 %*% as.matrix(x)
        }) +
        apply(mod$post$gm, 1, function(x){
          Z4 %*% as.matrix(x)
        })

      if(increase){
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
        }
      } else {
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
        }
      }

      ## Probability of superior performance by location ----------------
      cond_ggl = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,loc]),
        mean
      )
      cond_ggl = lapply(split(cond_ggl, f = cond_ggl[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggl, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggl = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggl)[-1] = name.loc

      prob_ggl.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggl, direction = 'long',
          varying = list(colnames(prob_ggl)[-1]),
          ids = name.gen, times = name.loc,
          v.names = 'prob',
          idvar = 'gen', timevar = 'loc'
        ),
        aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Locations", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5)) +
        scale_y_discrete(limits = rev)

      if(verbose) message('10. Probability of superior performance within locations estimated')

      ## Pairwise probability of superior performance per location ----------------

      if(increase){
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }


      if(verbose) message('11. Pairwise probability of superior performance within locations estimated')

      ## Probability of superior performance by Region ----------------
      cond_ggm = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,reg]),
        mean
      )
      cond_ggm = lapply(split(cond_ggm, f = cond_ggm[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggm, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggm = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggm)[-1] = name.reg

      prob_ggm.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggm, direction = 'long',
          varying = list(colnames(prob_ggm)[-1]),
          ids = name.gen, times = name.reg,
          v.names = 'prob',
          idvar = 'gen', timevar = 'reg'
        ),
        aes(x = .data$reg, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Regions", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))+
        scale_y_discrete(limits = rev)

      if(verbose) message('12. Probability of superior performance within regions estimated')

      ## Pairwise probability of superior performance per Region ----------------

      if(increase){
        pwprobs.reg = lapply(cond_ggm, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.reg = lapply(cond_ggm, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }

      if(verbose) message('13. Pairwise probability of superior performance within regions estimated')

      ## Probability of superior performance by year ----------------
      cond_ggt = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,year]),
        mean
      )
      cond_ggt = lapply(split(cond_ggt, f = cond_ggt[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggt, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggt = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggt)[-1] = name.year

      prob_ggt.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggt, direction = 'long',
          varying = list(colnames(prob_ggt)[-1]),
          ids = name.gen, times = name.year,
          v.names = 'prob',
          idvar = 'gen', timevar = 'year'
        ),
        aes(x = .data$year, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Year", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))+
        scale_y_discrete(limits = rev)

      if(verbose) message('14. Probability of superior performance within year estimated')

      ## Pairwise probability of superior performance per year ----------------

      if(increase){
        pwprobs.year = lapply(cond_ggt, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.year.plots = lapply(pwprobs.year, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.year = lapply(cond_ggt, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.year.plots = lapply(pwprobs.year, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }

      if(verbose) message('15. Pairwise probability of superior performance within year estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggl, direction = 'long',
              varying = list(colnames(prob_ggl)[-1]),
              ids = name.gen, times = name.loc,
              v.names = 'prob',
              idvar = 'gen', timevar = 'loc'
            ),
            aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Locations", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        prob_ggm.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggm, direction = 'long',
              varying = list(colnames(prob_ggm)[-1]),
              ids = name.gen, times = name.reg,
              v.names = 'prob',
              idvar = 'gen', timevar = 'reg'
            ),
            aes(x = .data$reg, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Region", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        prob_ggt.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggt, direction = 'long',
              varying = list(colnames(prob_ggt)[-1]),
              ids = name.gen, times = name.year,
              v.names = 'prob',
              idvar = 'gen', timevar = 'year'
            ),
            aes(x = .data$year, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Year", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        if(increase){
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

        if(increase){
          pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

        if(increase){
          pwprobs.year.plots = lapply(pwprobs.year, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.year.plots = lapply(pwprobs.year, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

      }

      ## Save the outputs in a list -----------------
      cond_prob = list(
        df = list(
          prob_loc = prob_ggl,
          prob_reg = prob_ggm,
          prob_year = prob_ggt,
          pwprob_loc = pwprobs.loc,
          pwprob_reg = pwprobs.reg,
          pwprob_year = pwprobs.year
        ),
        plots = list(
          prob_loc = prob_ggl.plot,
          prob_reg = prob_ggm.plot,
          prob_year = prob_ggt.plot,
          pwprob_loc = pwprobs.loc.plots,
          pwprob_reg = pwprobs.reg.plots,
          pwprob_year = pwprobs.year.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/cond_prob/prob_loc.csv'),
                         row.names = F)
        utils::write.csv(prob_ggm,
                         file = paste0(getwd(),'/cond_prob/prob_reg.csv'),
                         row.names = F)
        utils::write.csv(prob_ggt,
                         file = paste0(getwd(),'/cond_prob/prob_year.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_loc'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_loc/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_reg'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_reg/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_year'))
        for (i in names(pwprobs.year)){
          utils::write.csv(pwprobs.year[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_year/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)

    }
    else # Without region info --------------------
    {
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.loc),
                                    "loc", rep(name.loc,  each = num.gen), sep = '_@#')

      aux = unique(data[,c(gen,loc,year)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])
      Z3 = stats::model.matrix(~-1 + aux[,gen]:aux[,year])

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        UP = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        up = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        DOWN = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        down = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g, y = reorder(.data$gen, .data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#33a02c')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------

      if(increase){
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      } else {
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      }

      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g) +
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90)) # + ylim(0,1)

      if(verbose) message('1. Probability of superior performance estimated')

      ord_gen = factor(prob_g$ID, levels = prob_g$ID)

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))

      if(increase){
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
          }
        }
      } else {
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
          }
        }
      }

      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_@#'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      if(increase){
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      } else {
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      }

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gl[,grep(paste0(x,'$'), colnames(staprob_gl))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gl = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gl = prob_gl[order(prob_gl$prob, decreasing = T),]
      prob_gl.plot = ggplot(prob_gl)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('3. Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gl1 = pwsprob_gl[match(prob_gl$ID, rownames(pwsprob_gl)),
                               match(prob_gl$ID, rownames(pwsprob_gl))]
      pwsprob_gl1[upper.tri(pwsprob_gl1, diag = T)] = NA
      pwsprob_gl1 = stats::reshape(
        data.frame(pwsprob_gl1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl1), times = colnames(pwsprob_gl1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gl.plot = ggplot(pwsprob_gl1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gl[x]) < var(gl[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - year --------------
      staprob_gt = mod$post$gt
      colnames(staprob_gt) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gt),'_@#year'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gt[,grep(paste0(x,'$'), colnames(staprob_gt))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gt = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gt = prob_gt[order(prob_gt$prob, decreasing = T),]
      prob_gt.plot = ggplot(prob_gt)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('5. Probability of superior stability (GT) estimated')


      ## Pairwise probability of superior stability - year -----------------
      pwsprob_gt = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gt)){
        for (j in colnames(pwsprob_gt)) {
          pwsprob_gt[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gt1 = pwsprob_gt[match(prob_gt$ID, rownames(pwsprob_gt)),
                               match(prob_gt$ID, rownames(pwsprob_gt))]
      pwsprob_gt1[upper.tri(pwsprob_gt1, diag = T)] = NA
      pwsprob_gt1 = stats::reshape(
        data.frame(pwsprob_gt1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gt1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gt1), times = colnames(pwsprob_gt1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gt.plot = ggplot(pwsprob_gt1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gt[x]) < var(gt[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gt[upper.tri(pwsprob_gt, diag = T)] = NA
      pwsprob_gt = stats::reshape(
        data.frame(pwsprob_gt),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gt))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gt), times = colnames(pwsprob_gt),
        new.row.names = 1:length(c(pwsprob_gt)), v.names = 'prob'
      )
      pwsprob_gt = stats::na.exclude(pwsprob_gt[order(pwsprob_gt$x),])

      if(verbose) message('6. Pairwise probability of superior stability (GT) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gt, by = 'ID'))

      j_prob$joint = j_prob[,2] * j_prob[,3]
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Location', 'Year'), each = num.gen)
      j_prob = stats::reshape(j_prob, direction = 'long', varying = list(2:4),
                              times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]

      j_prob = j_prob[j_prob$category == "Joint",]

      retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]

      j_prob.plot = ggplot(cbind(j_prob, V4 = paste(j_prob$ID, j_prob$level, sep = '@#_')),
                           aes(x = stats::reorder(.data$V4, -.data$prob), y = .data$prob)) +
        facet_wrap(.~.data$level, ncol = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_x_discrete(labels = retrieve) +
        geom_segment(aes(x = reorder(.data$V4, -.data$prob), xend = .data$V4,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = "Genotype", y = "Joint probability")

      if(verbose) message('7. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------------

      if(interactive){
        g_hpd = plotly::ggplotly(g_hpd)
        prob_g.plot = plotly::ggplotly(prob_g.plot)
        prob_gl.plot = plotly::ggplotly(prob_gl.plot)
        prob_gt.plot = plotly::ggplotly(prob_gt.plot)
        j_prob.plot = plotly::ggplotly(j_prob.plot)
        pwsprob_g.plot = plotly::ggplotly(
          ggplot(data = pwsprob_g1, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )

        pwsprob_gt.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gt, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gt(x)) < var(gt(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------------
      marg_prob = list(
        df = list(
          perfo = prob_g,
          pair_perfo = pwsprob_g,
          stabi_gl = prob_gl,
          pair_stabi_gl = pwsprob_gl,
          stabi_gt = prob_gt,
          pair_stabi_gt = pwsprob_gt,
          joint_prob = j_prob
        ),
        plots = list(
          g_hpd = g_hpd,
          perfo = prob_g.plot,
          pair_perfo = pwsprob_g.plot,
          stabi_gl = prob_gl.plot,
          pair_stabi_gl = pwsprob_gl.plot,
          stabi_gt = prob_gt.plot,
          pair_stabi_gt = pwsprob_gt.plot,
          joint_prob = j_prob.plot
        )
      )

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          utils::write.csv(marg_prob$df[[i]],
                           file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                           row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = apply(mod$post$g, 1, function(x){
        Z1 %*% as.matrix(x)
      }) +
        apply(mod$post$gl, 1, function(x){
          Z2 %*% as.matrix(x)
        }) +
        apply(mod$post$gt, 1, function(x){
          Z3 %*% as.matrix(x)
        })

      if(increase){
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
        }
      } else {
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
        }
      }

      ## Probability of superior performance by location ----------------
      cond_ggl = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,loc]),
        mean
      )
      cond_ggl = lapply(split(cond_ggl, f = cond_ggl[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggl, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggl = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggl)[-1] = name.loc

      prob_ggl.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggl, direction = 'long',
          varying = list(colnames(prob_ggl)[-1]),
          ids = name.gen, times = name.loc,
          v.names = 'prob',
          idvar = 'gen', timevar = 'loc'
        ),
        aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Location", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5)) +
        scale_y_discrete(limits = rev)

      if(verbose) message('8. Probability of superior performance within locations estimated')

      ## Pairwise probability of superior performance per location ----------------

      if(increase){
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }


      if(verbose) message('9. Pairwise probability of superior performance within locations estimated')

      ## Probability of superior performance by year ----------------
      cond_ggt = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,year]),
        mean
      )
      cond_ggt = lapply(split(cond_ggt, f = cond_ggt[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggt, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggt = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggt)[-1] = name.year

      prob_ggt.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggt, direction = 'long',
          varying = list(colnames(prob_ggt)[-1]),
          ids = name.gen, times = name.year,
          v.names = 'prob',
          idvar = 'gen', timevar = 'year'
        ),
        aes(x = .data$year, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Year", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))+
        scale_y_discrete(limits = rev)

      if(verbose) message('10. Probability of superior performance within year estimated')

      ## Pairwise probability of superior performance per year ----------------

      if(increase){
        pwprobs.year = lapply(cond_ggt, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.year.plots = lapply(pwprobs.year, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.year = lapply(cond_ggt, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.year.plots = lapply(pwprobs.year, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }

      if(verbose) message('11. Pairwise probability of superior performance within year estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggl, direction = 'long',
              varying = list(colnames(prob_ggl)[-1]),
              ids = name.gen, times = name.loc,
              v.names = 'prob',
              idvar = 'gen', timevar = 'loc'
            ),
            aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Locations", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )


        prob_ggt.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggt, direction = 'long',
              varying = list(colnames(prob_ggt)[-1]),
              ids = name.gen, times = name.year,
              v.names = 'prob',
              idvar = 'gen', timevar = 'year'
            ),
            aes(x = .data$year, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Year", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        if(increase){
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }


        if(increase){
          pwprobs.year.plots = lapply(pwprobs.year, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.year.plots = lapply(pwprobs.year, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

      }

      ## Save the outputs in a list -----------------
      cond_prob = list(
        df = list(
          prob_loc = prob_ggl,
          prob_year = prob_ggt,
          pwprob_loc = pwprobs.loc,
          pwprob_year = pwprobs.year
        ),
        plots = list(
          prob_loc = prob_ggl.plot,
          prob_year = prob_ggt.plot,
          pwprob_loc = pwprobs.loc.plots,
          pwprob_year = pwprobs.year.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/cond_prob/prob_loc.csv'),
                         row.names = F)

        utils::write.csv(prob_ggt,
                         file = paste0(getwd(),'/cond_prob/prob_year.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_loc'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_loc/',i,'.csv'),
                           row.names = F)
        }

        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_year'))
        for (i in names(pwprobs.year)){
          utils::write.csv(pwprobs.year[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_year/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)

    }
  }
  else  # Without year info ----------------
  {
    colnames(mod$post$g) = paste0(name.gen, '_@#')

    if(!is.null(reg))  # With region info --------------------
    {
      if(!all(grepl('[A-Za-z]', data[, reg]))){data[,reg] = paste("R", data[,reg], sep = "_")}

      name.reg = levels(factor(data[,reg]))
      num.reg = nlevels(factor(data[,reg]))
      colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                    'Reg',rep(name.reg,  each = num.gen), sep = '_@#')
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.loc),
                                    "loc", rep(name.loc,  each = num.gen), sep = '_@#')

      aux = unique(data[,c(gen,loc,reg)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])
      Z4 = stats::model.matrix(~-1 + aux[,gen]:aux[,reg])

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        UP = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        up = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        DOWN = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        down = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g, y = reorder(.data$gen, .data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#33a02c')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------

      if(increase){
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      } else {
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      }

      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g) +
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90)) # + ylim(0,1)

      if(verbose) message('1. Probability of superior performance estimated')

      ord_gen = factor(prob_g$ID, levels = prob_g$ID)

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))

      if(increase){
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
          }
        }
      } else {
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
          }
        }
      }

      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_@#'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      if(increase){
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      } else {
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      }

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gl[,grep(paste0(x,'$'), colnames(staprob_gl))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gl = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gl = prob_gl[order(prob_gl$prob, decreasing = T),]
      prob_gl.plot = ggplot(prob_gl)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('3. Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gl1 = pwsprob_gl[match(prob_gl$ID, rownames(pwsprob_gl)),
                               match(prob_gl$ID, rownames(pwsprob_gl))]
      pwsprob_gl1[upper.tri(pwsprob_gl1, diag = T)] = NA
      pwsprob_gl1 = stats::reshape(
        data.frame(pwsprob_gl1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl1), times = colnames(pwsprob_gl1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gl.plot = ggplot(pwsprob_gl1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gl[x]) < var(gl[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gm),'_@#Reg'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gm[,grep(paste0(x,'$'), colnames(staprob_gm))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gm = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gm = prob_gm[order(prob_gm$prob, decreasing = T),]
      prob_gm.plot = ggplot(prob_gm)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('5. Probability of superior stability (GM) estimated')


      ## Pairwise probability of superior stability - Region -----------------
      pwsprob_gm = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gm)){
        for (j in colnames(pwsprob_gm)) {
          pwsprob_gm[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gm1 = pwsprob_gm[match(prob_gm$ID, rownames(pwsprob_gm)),
                               match(prob_gm$ID, rownames(pwsprob_gm))]
      pwsprob_gm1[upper.tri(pwsprob_gm1, diag = T)] = NA
      pwsprob_gm1 = stats::reshape(
        data.frame(pwsprob_gm1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gm1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm1), times = colnames(pwsprob_gm1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gm.plot = ggplot(pwsprob_gm1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gm[x]) < var(gm[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gm[upper.tri(pwsprob_gm, diag = T)] = NA
      pwsprob_gm = stats::reshape(
        data.frame(pwsprob_gm),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gm))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm), times = colnames(pwsprob_gm),
        new.row.names = 1:length(c(pwsprob_gm)), v.names = 'prob'
      )
      pwsprob_gm = stats::na.exclude(pwsprob_gm[order(pwsprob_gm$x),])

      if(verbose) message('6. Pairwise probability of superior stability (GM) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gm, by = 'ID'))

      j_prob$joint = j_prob[,2] * j_prob[,3]
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Location', 'Region'), each = num.gen)
      j_prob = stats::reshape(j_prob, direction = 'long', varying = list(2:4),
                              times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]

      j_prob = j_prob[j_prob$category == "Joint",]

      retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]

      j_prob.plot = ggplot(cbind(j_prob, V4 = paste(j_prob$ID, j_prob$level, sep = '@#_')),
                           aes(x = stats::reorder(.data$V4, -.data$prob), y = .data$prob)) +
        facet_wrap(.~.data$level, ncol = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_x_discrete(labels = retrieve) +
        geom_segment(aes(x = reorder(.data$V4, -.data$prob), xend = .data$V4,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = "Genotype", y = "Joint probability")

      if(verbose) message('7. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------------

      if(interactive){
        g_hpd = plotly::ggplotly(g_hpd)
        prob_g.plot = plotly::ggplotly(prob_g.plot)
        prob_gl.plot = plotly::ggplotly(prob_gl.plot)
        prob_gm.plot = plotly::ggplotly(prob_gm.plot)
        j_prob.plot = plotly::ggplotly(j_prob.plot)
        pwsprob_g.plot = plotly::ggplotly(
          ggplot(data = pwsprob_g1, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gm.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gm, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gm(x)) < var(gm(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------------
      marg_prob = list(
        df = list(
          perfo = prob_g,
          pair_perfo = pwsprob_g,
          stabi_gl = prob_gl,
          pair_stabi_gl = pwsprob_gl,
          stabi_gm = prob_gm,
          pair_stabi_gm = pwsprob_gm,
          joint_prob = j_prob
        ),
        plots = list(
          g_hpd = g_hpd,
          perfo = prob_g.plot,
          pair_perfo = pwsprob_g.plot,
          stabi_gl = prob_gl.plot,
          pair_stabi_gl = pwsprob_gl.plot,
          stabi_gm = prob_gm.plot,
          pair_stabi_gm = pwsprob_gm.plot,
          joint_prob = j_prob.plot
        )
      )

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          utils::write.csv(marg_prob$df[[i]],
                           file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                           row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = apply(mod$post$g, 1, function(x){
        Z1 %*% as.matrix(x)
      }) +
        apply(mod$post$gl, 1, function(x){
          Z2 %*% as.matrix(x)
        }) +
        apply(mod$post$gm, 1, function(x){
          Z4 %*% as.matrix(x)
        })

      if(increase){
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
        }
      } else {
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
        }
      }

      ## Probability of superior performance by location ----------------
      cond_ggl = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,loc]),
        mean
      )
      cond_ggl = lapply(split(cond_ggl, f = cond_ggl[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggl, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggl = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggl)[-1] = name.loc

      prob_ggl.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggl, direction = 'long',
          varying = list(colnames(prob_ggl)[-1]),
          ids = name.gen, times = name.loc,
          v.names = 'prob',
          idvar = 'gen', timevar = 'loc'
        ),
        aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Locations", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5)) +
        scale_y_discrete(limits = rev)

      if(verbose) message('8. Probability of superior performance within locations estimated')

      ## Pairwise probability of superior performance per location ----------------

      if(increase){
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }


      if(verbose) message('9. Pairwise probability of superior performance within locations estimated')

      ## Probability of superior performance by Region ----------------
      cond_ggm = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,reg]),
        mean
      )
      cond_ggm = lapply(split(cond_ggm, f = cond_ggm[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggm, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggm = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggm)[-1] = name.reg

      prob_ggm.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggm, direction = 'long',
          varying = list(colnames(prob_ggm)[-1]),
          ids = name.gen, times = name.reg,
          v.names = 'prob',
          idvar = 'gen', timevar = 'reg'
        ),
        aes(x = .data$reg, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Regions", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))+
        scale_y_discrete(limits = rev)

      if(verbose) message('10. Probability of superior performance within regions estimated')

      ## Pairwise probability of superior performance per Region ----------------

      if(increase){
        pwprobs.reg = lapply(cond_ggm, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.reg = lapply(cond_ggm, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }

      if(verbose) message('11. Pairwise probability of superior performance within regions estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggl, direction = 'long',
              varying = list(colnames(prob_ggl)[-1]),
              ids = name.gen, times = name.loc,
              v.names = 'prob',
              idvar = 'gen', timevar = 'loc'
            ),
            aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Locations", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        prob_ggm.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggm, direction = 'long',
              varying = list(colnames(prob_ggm)[-1]),
              ids = name.gen, times = name.reg,
              v.names = 'prob',
              idvar = 'gen', timevar = 'reg'
            ),
            aes(x = .data$reg, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Region", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        if(increase){
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

        if(increase){
          pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

      }

      ## Save the outputs in a list -----------------
      cond_prob = list(
        df = list(
          prob_loc = prob_ggl,
          prob_reg = prob_ggm,
          pwprob_loc = pwprobs.loc,
          pwprob_reg = pwprobs.reg
        ),
        plots = list(
          prob_loc = prob_ggl.plot,
          prob_reg = prob_ggm.plot,
          pwprob_loc = pwprobs.loc.plots,
          pwprob_reg = pwprobs.reg.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/cond_prob/prob_loc.csv'),
                         row.names = F)
        utils::write.csv(prob_ggm,
                         file = paste0(getwd(),'/cond_prob/prob_reg.csv'),
                         row.names = F)

        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_loc'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_loc/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_reg'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_reg/',i,'.csv'),
                           row.names = F)
        }

      }

      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)

    }
    else # Without region info --------------------
    {
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.loc),
                                    "loc", rep(name.loc,  each = num.gen), sep = '_@#')

      aux = unique(data[,c(gen,loc)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        UP = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        up = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        DOWN = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        down = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g, y = reorder(.data$gen, .data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#33a02c')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------

      if(increase){
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      } else {
        ind_post = apply(mod$post$g, 1, function(x){
          ifelse(name.gen %in%
                   unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                   split = '_@#')), 1, 0)
        })
      }

      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g) +
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90)) # + ylim(0,1)

      if(verbose) message('1. Probability of superior performance estimated')

      ord_gen = factor(prob_g$ID, levels = prob_g$ID)

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))

      if(increase){
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
          }
        }
      } else {
        for(i in colnames(mod$post$g)){
          for (j in colnames(mod$post$g)) {
            pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
          }
        }
      }

      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_@#'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      if(increase){
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      } else {
        pwsprob_g.plot = ggplot(pwsprob_g1,
                                aes(x = factor(.data$x,
                                               levels = unique(.data$x)),
                                    y = factor(.data$y,
                                               levels = unique(.data$y))))+
          geom_tile(aes(fill = .data$prob), colour = 'white') +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')+
          scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                               option = 'turbo')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5))
      }

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_g))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gl[,grep(paste0(x,'$'), colnames(staprob_gl))]
        ),
        function(x) apply(x, 1, var)
      ))
      colnames(probsta) = name.gen

      ind_post = apply(probsta, 1, function(x){
        ifelse(name.gen %in% names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_gl = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_gl = prob_gl[order(prob_gl$prob, decreasing = T),]
      prob_gl.plot = ggplot(prob_gl)+
        geom_segment(aes(x = factor(.data$ID, levels = .data$ID), xend = .data$ID,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(aes(x = factor(.data$ID, levels = .data$ID), y = .data$prob),
                   size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90)) #+ ylim(0,1)

      if(verbose) message('3. Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      pwsprob_gl1 = pwsprob_gl[match(prob_gl$ID, rownames(pwsprob_gl)),
                               match(prob_gl$ID, rownames(pwsprob_gl))]
      pwsprob_gl1[upper.tri(pwsprob_gl1, diag = T)] = NA
      pwsprob_gl1 = stats::reshape(
        data.frame(pwsprob_gl1),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl1))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl1), times = colnames(pwsprob_gl1),
        new.row.names = NULL, v.names = 'prob'
      )
      pwsprob_gl.plot = ggplot(pwsprob_gl1,
                               aes(x = factor(.data$x, levels = unique(.data$x)),
                                   y = factor(.data$y, levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white')+
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold("Pr["~var(gl[x]) < var(gl[y])~"]")))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                             option = 'turbo')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(data.frame(pwsprob_gl))),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'))

      j_prob$joint = j_prob[,2] * j_prob[,3]
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Location'), each = num.gen)
      j_prob = stats::reshape(j_prob, direction = 'long', varying = list(2:4),
                              times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]

      j_prob = j_prob[j_prob$category == "Joint",]

      retrieve = function(x) do.call(rbind, strsplit(x, '@#_'))[,1]

      j_prob.plot = ggplot(cbind(j_prob, V4 = paste(j_prob$ID, j_prob$level, sep = '@#_')),
                           aes(x = stats::reorder(.data$V4, -.data$prob), y = .data$prob)) +
        #facet_wrap(.~.data$level, ncol = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_x_discrete(labels = retrieve) +
        geom_segment(aes(x = reorder(.data$V4, -.data$prob), xend = .data$V4,
                         y = 0, yend = .data$prob), linewidth = 1) +
        geom_point(size = 2, color = 'black', fill = alpha('#33a02c', 0.3),
                   alpha = 0.7, shape = 21, stroke = .4) +
        labs(x = "Genotype", y = "Joint probability")

      if(verbose) message('5. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------------

      if(interactive){
        g_hpd = plotly::ggplotly(g_hpd)
        prob_g.plot = plotly::ggplotly(prob_g.plot)
        prob_gl.plot = plotly::ggplotly(prob_gl.plot)
        j_prob.plot = plotly::ggplotly(j_prob.plot)
        pwsprob_g.plot = plotly::ggplotly(
          ggplot(data = pwsprob_g1, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'turbo', limits = c(0,1),
                                 na.value = 'white', direction = 1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------------
      marg_prob = list(
        df = list(
          perfo = prob_g,
          pair_perfo = pwsprob_g,
          stabi_gl = prob_gl,
          pair_stabi_gl = pwsprob_gl,
          joint_prob = j_prob
        ),
        plots = list(
          g_hpd = g_hpd,
          perfo = prob_g.plot,
          pair_perfo = pwsprob_g.plot,
          stabi_gl = prob_gl.plot,
          pair_stabi_gl = pwsprob_gl.plot,
          joint_prob = j_prob.plot
        )
      )

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          utils::write.csv(marg_prob$df[[i]],
                           file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                           row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = apply(mod$post$g, 1, function(x){
        Z1 %*% as.matrix(x)
      }) +
        apply(mod$post$gl, 1, function(x){
          Z2 %*% as.matrix(x)
        })

      if(increase){
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
        }
      } else {
        supprob = function(vector, num.gen, int){
          ifelse(names(vector) %in%
                   names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
        }
      }

      ## Probability of superior performance by location ----------------
      cond_ggl = stats::aggregate(
        x = posgge,
        by = list(aux[,gen], aux[,loc]),
        mean
      )
      cond_ggl = lapply(split(cond_ggl, f = cond_ggl[,2]), function(x){
        rownames(x) = x[, 1]
        x = x[,-c(1,2)]
        x = t(x)
        x
      })

      probs = lapply(cond_ggl, function(x){
        y = apply(
          x, MARGIN = 1, FUN = supprob, num.gen = num.gen, int = int
        )
        rownames(y) = colnames(x)
        y = rowMeans(y)
        data.frame(gen = names(y),
                   probs = y)
      })

      prob_ggl = suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by = 'gen', all = T), probs))
      colnames(prob_ggl)[-1] = name.loc

      prob_ggl.plot = ggplot(
        data =  stats::reshape(
          data = prob_ggl, direction = 'long',
          varying = list(colnames(prob_ggl)[-1]),
          ids = name.gen, times = name.loc,
          v.names = 'prob',
          idvar = 'gen', timevar = 'loc'
        ),
        aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen), fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'turbo') +
        labs(x = "Locations", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5)) +
        scale_y_discrete(limits = rev)

      if(verbose) message('6. Probability of superior performance within locations estimated')

      ## Pairwise probability of superior performance per location ----------------

      if(increase){
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] >
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      } else {
        pwprobs.loc = lapply(cond_ggl, function(x){
          combs = data.frame(t(utils::combn(colnames(x), 2)))
          colnames(combs) = c('x', 'y')

          a = cbind(combs,
                    pwprob = apply(combs, 1, function(y){
                      mean(x[,grep(paste0(y[1], "$"), colnames(x))] <
                             x[,grep(paste0(y[2], "$"), colnames(x))])
                    }))
          a
        })

        pwprobs.loc.plots = lapply(pwprobs.loc, function(x){
          ggplot(data = x, aes(x = factor(.data$x, levels = unique(.data$x)),
                               y = factor(.data$y, levels = unique(.data$y)),
                               fill = .data$pwprob)) +
            geom_tile() +
            labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
            scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                 option = 'turbo')+
            guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                         title.hjust = .5)) +
            theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                  legend.position = c(.8,.15), legend.direction = 'horizontal')
        })
      }


      if(verbose) message('7. Pairwise probability of superior performance within locations estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  stats::reshape(
              data = prob_ggl, direction = 'long',
              varying = list(colnames(prob_ggl)[-1]),
              ids = name.gen, times = name.loc,
              v.names = 'prob',
              idvar = 'gen', timevar = 'loc'
            ),
            aes(x = .data$loc, y = factor(.data$gen, levels = ord_gen),
                fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = 1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'turbo') +
            labs(x = "Locations", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )


        if(increase){
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        } else {
          pwprobs.loc_plots = lapply(pwprobs.loc, function(x){
            plotly::ggplotly(
              ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
                geom_tile() +
                labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
                scale_fill_viridis_c(direction = 1, na.value = 'white',limits = c(0,1),
                                     option = 'turbo') +
                theme(axis.text.x = element_text(angle = 90),
                      panel.background = element_blank())
            )
          })
        }

      }

      ## Save the outputs in a list -----------------
      cond_prob = list(
        df = list(
          prob_loc = prob_ggl,
          pwprob_loc = pwprobs.loc
        ),
        plots = list(
          prob_loc = prob_ggl.plot,
          pwprob_loc = pwprobs.loc.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/cond_prob/prob_loc.csv'),
                         row.names = F)

        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_loc'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_loc/',i,'.csv'),
                           row.names = F)
        }

      }

      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)

    }
  }
}
