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
##' @param trait,gen,env A string. The name of the columns that corresponds to
##' the variable, genotype and environment information, respectively
##' @param reg A string or NULL. If the dataset has information about regions,
##' `reg` will be a string with the name of the column that corresponds to the
##' region information. Otherwise, `reg = NULL` (default).
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
##' \item \code{perfo}: the probabilities of superior performance
##' \item \code{pair_perfo}: the pairwise probabilities of superior performance
##' \item \code{stabi}: the probabilities of superior stability. When `reg` is not
##' `NULL`, `stabi` is divided into `stabi_gl` for the stability across environments, and
##' `stabi_gm` for the stability across regions
##' \item \code{pair_stabi}: the pairwise probabilities of superior stability, which
##' is also divided into `pair_stabi_gl` and `pair_stabi_gm` when `reg` is not `NULL`
##' \item \code{joint_prob}: the joint probabilities of superior performance and stability
##' }
##' \item \code{plot} : A list of ggplots illustrating the outputs:
##' \itemize{
##' \item \code{g_hpd}: a caterpillar plot representing the marginal genotypic value of
##' each genotype, and their respective highest posterior density interval (95% represented by the
##' thick line, and 97.5% represented by the thin line)
##' \item \code{perfo}: a bar plot illustrating the probabilities of superior performance
##' \item \code{pair_perfo}: a heatmap representing the pairwise probability of superior
##' performance (the probability of genotypes at the \emph{x}-axis being superior
##' to those on the \emph{y}-axis)
##' \item \code{stabi}: a bar plot with the probabilities of superior stability. Like the data frames,
##' when `reg` is not `NULL`, two different plots are generated, one for the stability across
##' environments (`stabi_gl`), and another for the stability across regions (`stabi_gm`)
##' \item \code{pair_stabi}: a heatmap with the pairwise probabilities of superior stability
##' (also divided into two different plots when `reg` is not `NULL`). This plot represents
##' the probability of genotypes at the \emph{x}-axis being superior
##' to those on \emph{y}-axis
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
##' \item \code{perfo}: a data frame containing the probabilities of superior performance
##' within environments. It also has the probabilities of superior performance within regions
##' if `reg` is not `NULL`.
##' \item \code{pair_perfo}: a list with the pairwise probabilities of superior performance
##' within environments. If `reg` is not `NULL`, two lists are generated.
##' }
##' \item \code{plot} : A list with:
##' \itemize{
##' \item \code{perfo}: a heatmap with the probabilities of superior performance within
##' environments. If `reg` is not `NULL`, there will be two heatmaps: `perfo_env`
##' for the probabilities of superior performance within environments, and `perfo_reg`
##' for the same probabilities within regions.
##' \item \code{pair_perfo}: a list of heatmaps representing the pairwise probability of superior
##' performance within environments. If `reg` is not `NULL`, there will be two lists: `pair_perfo_env`
##' for the probabilities of superior performance within environments, and `pair_perfo_reg`
##' for the same probabilities within regions. The interpretation is the same as in the
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
##' \eqn{j^{\text{th}} \in \Omega_k} can calculated as follows:
##'
##' \deqn{Pr(\hat{g}_{jk} \in \Omega_k \vert y) = \frac{1}{S} \sum_{s=1}^S I(\hat{g}_{jk}^{(s)} \in \Omega_k \vert y)}
##'
##' where \eqn{I(\hat{g}_{jk}^{(s)} \in \Omega_k \vert y)} is an indicator variable
##' mapping success (1) if \eqn{\hat{g}_{jk}^{(s)}} exists in \eqn{\Omega_k}, and
##' failure (0) otherwise, and \eqn{\hat{g}_{jk}^{(s)} = \hat{g}_j^{(s)} + \widehat{ge}_{jk}^{(s)}}.
##' Note that when computing conditional probabilities (i.e., conditional to the
##' \eqn{k^{\text{th}}} environment or mega-environment), we are accounting for
##' the interaction of the \eqn{j^{\text{th}}} genotype with the \eqn{k^{\text{th}}}
##' environment.
##'
##' The pairwise probabilities of superior performance can also be calculated across
##' or within environments. This metric assesses the probability of the \eqn{j^{\text{th}}}
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
##' @importFrom stats reshape median quantile na.exclude
##' @importFrom rlang .data
##'
##' @export
##'
##' @examples
##' \donttest{
##' mod = bayes_met(data = soy,
##'                 gen = "Gen",
##'                 env = "Env",
##'                 repl = NULL,
##'                 reg = "Reg",
##'                 res.het = FALSE,
##'                 trait = "Y",
##'                 iter = 2000, cores = 1, chains = 4)
##'
##' outs = extr_outs(data = soy, trait = "Y", gen = "Gen", model = mod,
##'                  effects = c('l','g','gl','m','gm'),
##'                  nenv = length(unique(soy$Env)),
##'                  probs = c(0.05, 0.95),
##'                  check.stan.diag = TRUE,
##'                  verbose = FALSE)
##'
##' results = prob_sup(data = soy,
##'                    trait = "Y",
##'                    gen = "Gen",
##'                    env = "Env",
##'                    mod.output = outs,
##'                    reg = 'Reg',
##'                    int = .2,
##'                    increase = TRUE,
##'                    save.df = FALSE,
##'                    interactive = FALSE,
##'                    verbose = FALSE)
##' }
##'

prob_sup = function(data, trait, gen, env, reg = NULL, mod.output, int,
                    increase = TRUE, save.df = FALSE, interactive = FALSE,
                    verbose = FALSE){

  # Conditions
  stopifnot("Each 'gen' and 'env' must be represented by a string (e.g., 'G01' or 'L25')" = {
    is.character(data[,gen])
    is.character(data[,env])
  })

  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })

  # Namespaces
  requireNamespace('ggplot2')

  # Preparation
  df = data
  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = mod.output
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))
  num.sim = nrow(mod$post$g)


  if(increase){ # Selection for increasing
    if(!is.null(reg)){# If there are breeding regions

      # Preparations

      stopifnot("Each 'reg' must be represented by a string (e.g., 'R08')" = {
        is.character(data[,reg])
      })

      name.reg = levels(factor(data[,reg]))
      num.reg = nlevels(factor(data[,reg]))
      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                    'Reg',rep(name.reg,  each = num.gen), sep = '_')
      name.env.reg = sort(paste('Env',unique(data[,c(env,reg)])[,1],
                                'Reg',unique(data[,c(env,reg)])[,2], sep = '_'))
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.env),
                                    rep(name.env.reg,  each = num.gen), sep = '_')

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
        geom_point(size = 4, color = '#781c1e')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------
      ind_post = apply(mod$post$g, 1, function(x){
        ifelse(name.gen %in%
                 unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                 split = '_')), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g, aes(x = factor(.data$ID, levels = .data$ID),
                                       y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90))

      if(verbose) message('1. Probability of superior performance estimated')

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                       dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))
      for(i in colnames(mod$post$g)){
        for (j in colnames(mod$post$g)) {
          pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
        }
      }
      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                         match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(pwsprob_g1)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      pwsprob_g.plot = ggplot(pwsprob_g1,
                              aes(x = factor(.data$x,
                                             levels = unique(.data$x)),
                                  y = factor(.data$y,
                                             levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white') +
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'plasma')+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(pwsprob_g)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gl),'_Env'))[,1]
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
      prob_gl.plot = ggplot(prob_gl,
                           aes(x = factor(.data$ID, levels = .data$ID),
                               y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

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
        varying = list(colnames(pwsprob_gl1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(pwsprob_gl)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gm),'_Reg'))[,1]
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
      prob_gm.plot = ggplot(prob_gm,
                           aes(x = factor(.data$ID, levels = .data$ID),
                               y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

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
        varying = list(colnames(pwsprob_gm1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gm[upper.tri(pwsprob_gm, diag = T)] = NA
      pwsprob_gm = stats::reshape(
        data.frame(pwsprob_gm),
        direction = 'long',
        varying = list(colnames(pwsprob_gm)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm), times = colnames(pwsprob_gm),
        new.row.names = 1:length(c(pwsprob_gm)), v.names = 'prob'
      )
      pwsprob_gm = stats::na.exclude(pwsprob_gm[order(pwsprob_gm$x),])

      if(verbose) message('6. Pairwise probability of superior stability (GM) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gm, by = 'ID'))
      j_prob$joint = j_prob$prob.x * j_prob$prob.y
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Environment', 'Region'), each = num.gen)
      j_prob = stats::reshape(j_prob, direction = 'long', varying = list(2:4),
              times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]
      j_prob.plot = ggplot(j_prob, aes(x = .data$ID, y = .data$prob)) +
        facet_wrap(.~.data$level, ncol = 1) +
        geom_segment(data = data.frame(
                   ID = name.gen,
                   prob = c(apply(merge(prob_g, prob_gl, by = 'ID')[,-1], 1,
                                 function(x) x[which.max(x)]),
                           apply(merge(prob_g, prob_gm, by = 'ID')[,-1], 1,
                                 function(x) x[which.max(x)])),
                   level = rep(c('Environment','Region'), each = num.gen)
                   ),
                 aes(x = .data$ID, y = 0, yend = .data$prob, xend = .data$ID),
                 linewidth = 1.2) +
        geom_point(aes(fill = .data$category, shape = .data$category), size = 2,
                   color = 'black') +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = 'top') +
        scale_fill_manual(
          label = c(
            'Joint' = 'Joint probability',
            'Performance' = 'Superior stability',
            'Stability' = 'Superior performance'
          ),
          values = c(
            'Joint' = '#1b9e77',
            'Performance' = '#d95f02',
            'Stability' = '#7570b3'
          )
        ) +
        scale_shape_manual(label = c(
          'Joint' = 'Joint probability',
          'Performance' = 'Superior stability',
          'Stability' = 'Superior performance'
        ),
        values = c(
          'Joint' = 21,
          'Performance' = 24,
          'Stability' = 25
        ))+
        ylim(0, 1) +
        labs(x = 'Genotype', y = 'Probabilities', fill = 'Probabilities',
             shape = 'Probabilities')

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
            scale_fill_viridis_c(option = 'plasma', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
          )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gm.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gm, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
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
      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl
      for (i in name.reg) {
        posgge[,grep(paste0(i,'$'), do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] =
          posgge[,grep(paste0(i,'$'), do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] +
          matrix(mod$post$gm[,grep(paste0(i,'$'), do.call(rbind, strsplit(colnames(mod$post$gm),'Reg'))[,2])],
                 nrow = num.sim, ncol = num.gen *
                   length(name.env.reg[grep(paste0(i,'$'), do.call(rbind,strsplit(name.env.reg,'Reg'))[,2])]))
      }

      ## Probability of superior performance ----------------
      supprob = function(vector, num.gen, int){
        ifelse(names(vector) %in%
                 names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
      }

      probs = apply(do.call(rbind, lapply(
        lapply(
          apply(
            posgge, 1, function(x){
              list(matrix(x, nrow = num.gen, ncol = num.env,
                          dimnames = list(name.gen, name.env.reg)))}
          ),
          Reduce, f = '+'
        ),
        function(x){
          apply(
            x, MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
          )}
      )), 2, function(x){
        tapply(x, rep(name.gen, num.sim), mean)
      })

      probs.df = stats::reshape(
        data = data.frame(probs), direction = 'long',
        varying = list(colnames(probs)),
        ids = name.gen, times = colnames(probs),
        new.row.names = 1:length(c(probs)), v.names = 'prob',
        idvar = 'gen', timevar = 'env'
      )
      probs.df = cbind(probs.df, do.call(rbind, strsplit(probs.df$env, '_Reg_')))
      probs.df = probs.df[,-1]; colnames(probs.df) = c('prob', 'gen', 'loc', 'reg')
      probs.df = probs.df[,c('reg', 'loc', 'gen', 'prob')]
      probs.df$loc = sub('Env_','', probs.df$loc)

      ### Per Location ----------------
      con_gl = ifelse(table(data[,gen], data[,env]) != 0, 1, NA)

      prob_ggl.plot = ggplot(
        data =  merge(
          x = probs.df,
          y = stats::reshape(
            data = data.frame(con_gl), direction = 'long',
            varying = list(colnames(con_gl)),
            ids = name.gen, times = name.env,
            v.names = 'freq',
            idvar = 'gen', timevar = 'loc'
          ), by= c('loc', 'gen')),
        aes(x = .data$loc, y = .data$gen, fill = .data$prob * .data$freq)
        ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Environment", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('8. Probability of superior performance within environments estimated')

      ### Per Region ----------------
      con_gm = ifelse(table(data[,gen], data[,reg]) != 0, 1, NA)

      prob_ggm.plot = ggplot(
        data =  merge(
          x = probs.df,
          y = stats::reshape(
            data = data.frame(con_gm), direction = 'long',
            varying = list(colnames(con_gm)),
            ids = name.gen, times = name.reg,
            v.names = 'freq',
            idvar = 'gen', timevar = 'reg'
          ), by= c('reg', 'gen')),
        aes(x = .data$reg, y = .data$gen, fill = .data$prob * .data$freq)
      )  +
        geom_tile(colour = 'white', na.rm = T) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Region", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('9. Probability of superior performance within regions estimated')

      ### Adjusting the data frame ----------------
      condprobs = merge(
        x = probs.df,
        y = stats::reshape(
          data = data.frame(con_gl), direction = 'long',
          varying = list(colnames(con_gl)),
          ids = name.gen, times = name.env,
          v.names = 'freq',
          idvar = 'gen', timevar = 'loc'
        ), by= c('loc', 'gen'))
      condprobs$prob = condprobs$prob * condprobs$freq
      condprobs = condprobs[,-5]
      colnames(condprobs) = c('env', 'gen', 'reg', 'prob')
      condprobs = condprobs[,c('env','reg','gen','prob')]

      ## Pairwise probability of superior performance ----------------
      ### Per Location -------------
      combs = data.frame(t(utils::combn(paste('Gen', name.gen, sep = '_'), 2)))
      colnames(combs) = c('x', 'y')
      pwprobs.env = lapply(
        sapply(paste('Env', name.env, sep = '_'),
               function(x) posgge[,grep(paste0(x,'_'), colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] >
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )

          a[,1] = sub('Gen_', '', a[,1])
          a[,2] = sub('Gen_', '', a[,2])

          a
        }
      )
      names(pwprobs.env) = sub('Env_', '', names(pwprobs.env))
      for (i in names(pwprobs.env)) {
        pwprobs.env[[i]] = merge(
          merge(pwprobs.env[[i]],
                data.frame(index = table(data[,gen],data[,env])[,i],
                           x = rownames(table(data[,gen],data[,env])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,env])[,i],
                     y = rownames(table(data[,gen],data[,env])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs.env[[i]] = pwprobs.env[[i]][which(pwprobs.env[[i]]$index.x != 0 &
                                                    pwprobs.env[[i]]$index.y != 0), ]
        pwprobs.env[[i]] = pwprobs.env[[i]][,-c(4,5)]
      }

      pwprobs.env.plots = lapply(pwprobs.env, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('10. Pairwise probability of superior performance within environments estimated')

      ### Per Region --------------
      pwprobs.reg = lapply(
        sapply(paste('Reg', name.reg, sep = '_'),
               function(x) posgge[,grep(paste0(x,'$'), colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] >
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )

          a[,1] = sub('Gen_', '', a[,1])
          a[,2] = sub('Gen_', '', a[,2])

          a
        }
      )
      names(pwprobs.reg) = sub('Reg_', '', names(pwprobs.reg))
      for (i in names(pwprobs.reg)) {
        pwprobs.reg[[i]] = merge(
          merge(pwprobs.reg[[i]],
                data.frame(index = table(data[,gen],data[,reg])[,i],
                           x = rownames(table(data[,gen],data[,reg])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,reg])[,i],
                     y = rownames(table(data[,gen],data[,reg])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs.reg[[i]] = pwprobs.reg[[i]][which(pwprobs.reg[[i]]$index.x != 0 &
                                                    pwprobs.reg[[i]]$index.y != 0), ]
        pwprobs.reg[[i]] = pwprobs.reg[[i]][,-c(4,5)]
      }

      pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('11. Pairwise probability of superior performance within regions estimated')


      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  merge(
              x = probs.df,
              y = stats::reshape(
                data = data.frame(con_gl), direction = 'long',
                varying = list(colnames(con_gl)),
                ids = name.gen, times = name.env,
                v.names = 'freq',
                idvar = 'gen', timevar = 'loc'
              ), by= c('loc', 'gen')),
            aes(x = .data$loc, y = .data$gen, fill = .data$prob * .data$freq)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Environment", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        prob_ggm.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  merge(
              x = probs.df,
              y = stats::reshape(
                data = data.frame(con_gm), direction = 'long',
                varying = list(colnames(con_gm)),
                ids = name.gen, times = name.reg,
                v.names = 'freq',
                idvar = 'gen', timevar = 'reg'
              ), by= c('reg', 'gen')),
            aes(x = .data$reg, y = .data$gen, fill = .data$prob * .data$freq)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Region", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        pwprobs.env_plots = lapply(pwprobs.env, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })
      }

      ## Save the outputs in a list -----------------

      cond_prob = list(
        df = list(
          perfo = condprobs,
          pair_perfo_env = pwprobs.env,
          pair_perfo_reg = pwprobs.reg
        ),
        plots = list(
          perfo_env = prob_ggl.plot,
          perfo_reg = prob_ggm.plot,
          pair_perfo_env = pwprobs.env.plots,
          pair_perfo_reg = pwprobs.reg.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(condprobs,
                         file = paste0(getwd(),'/cond_prob/perfo.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_env'))
        for (i in names(pwprobs.env)){
          utils::write.csv(pwprobs.env[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_env/perfo_',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_reg'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_reg/perfo_',i,'.csv'),
                           row.names = F)
        }
      }



      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)


    }else{ #If there is no breeding region

      # Preparation
      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                    rep(name.env,  each = num.gen), sep = '_')

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
        geom_point(size = 4, color = '#781c1e')

      # Marginal probabilities ----------------

      ## Probability of superior performance --------------
      ind_post = apply(mod$post$g, 1, function(x){
        ifelse(name.gen %in%
                 unlist(strsplit(names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]),
                                 split = '_')), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g, aes(x = factor(.data$ID, levels = .data$ID),
                                       y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90))

      cat('1. Probability of superior performance estimated')

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))
      for(i in colnames(mod$post$g)){
        for (j in colnames(mod$post$g)) {
          pwsprob_g[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
        }
      }
      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(pwsprob_g1)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      pwsprob_g.plot = ggplot(pwsprob_g1,
                              aes(x = factor(.data$x,
                                             levels = unique(.data$x)),
                                  y = factor(.data$y,
                                             levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white') +
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'plasma')+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(pwsprob_g)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability  -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gl),'_'))[,1]
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
      prob_gl.plot = ggplot(prob_gl,
                            aes(x = factor(.data$ID, levels = .data$ID),
                                y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

      if(verbose) message('3. Probability of superior stability estimated')

      ## Pairwise probability of superior stability  -------------
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
        varying = list(colnames(pwsprob_gl1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(pwsprob_gl)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability estimated')


      ## Joint probability of superior performance and stability -----------------
      j_prob = merge(prob_g, prob_gl, by = 'ID')
      j_prob$joint = j_prob$prob.x * j_prob$prob.y
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob = reshape(j_prob, direction = 'long', varying = list(2:4),
                       times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-4]
      colnames(j_prob) = c('ID', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$category,  j_prob$ID),]
      j_prob.plot = ggplot(j_prob, aes(x = .data$ID, y = .data$prob)) +
        geom_segment(data = data.frame(
          ID = name.gen,
          prob = apply(merge(prob_g, prob_gl, by = 'ID')[,-1], 1,
                         function(x) x[which.max(x)])
        ),
        aes(x = .data$ID, xend = .data$ID, y = 0, yend = .data$prob),
        linewidth = 1.2) +
        geom_point(aes(fill = .data$category, shape = .data$category), size = 2,
                   color = 'black') +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = 'top') +
        scale_fill_manual(
          label = c(
            'Joint' = 'Joint probability',
            'Performance' = 'Superior stability',
            'Stability' = 'Superior performance'
          ),
          values = c(
            'Joint' = '#1b9e77',
            'Performance' = '#d95f02',
            'Stability' = '#7570b3'
          )
        ) +
        scale_shape_manual(label = c(
          'Joint' = 'Joint probability',
          'Performance' = 'Superior stability',
          'Stability' = 'Superior performance'
        ),
        values = c(
          'Joint' = 21,
          'Performance' = 24,
          'Stability' = 25
        ))+
        ylim(0, 1) +
        labs(x = 'Genotype', y = 'Probabilities', fill = 'Probabilities',
             shape = 'Probabilities')

      if(verbose) message('5. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------

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
            scale_fill_viridis_c(option = 'plasma', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) > g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------

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

      ## Save data frames in the work directory -----------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          write.csv(marg_prob$df[[i]],
                    file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                    row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl

      ## Probability of superior performance ----------------
      supprob = function(vector, num.gen, int){
        ifelse(names(vector) %in%
                 names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
      }

      probs = apply(do.call(rbind, lapply(
        lapply(
          apply(
            posgge, 1, function(x){
              list(matrix(x, nrow = num.gen, ncol = num.env,
                          dimnames = list(name.gen, name.env)))}
          ),
          Reduce, f = '+'
        ),
        function(x){
          apply(
            x,MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
          )}
      )), 2, function(x){
        tapply(x, rep(name.gen, num.sim), mean)
      })
      probs = probs * ifelse(table(data[,gen],data[,env]) != 0, 1, NA)

      probs.df = stats::reshape(
        data = data.frame(probs), direction = 'long',
        varying = list(colnames(probs)),
        ids = name.gen, times = colnames(probs),
        new.row.names = 1:length(c(probs)), v.names = 'prob',
        idvar = 'gen', timevar = 'env'
      )
      probs.df = probs.df[,c('env', 'gen', 'prob')]

      prob_ggl.plot = ggplot(
        data =  probs.df,
        aes(x = .data$env, y = .data$gen, fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Environment", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('6. Probability of superior performance within environments estimated')

      ## Pairwise probability of superior performance ----------------
      ### Per Location -------------
      combs = data.frame(t(utils::combn(name.gen, 2)))
      colnames(combs) = c('x', 'y')
      pwprobs = lapply(
        sapply(name.env,
               function(x) posgge[,grep(paste0(x,'$'), colnames(posgge))],
               simplify = F),
        function(y){
          cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] >
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )
        }
      )

      for (i in names(pwprobs)) {
        pwprobs[[i]] = merge(
          merge(pwprobs[[i]],
                data.frame(index = table(data[,gen],data[,env])[,i],
                           x = rownames(table(data[,gen],data[,env])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,env])[,i],
                     y = rownames(table(data[,gen],data[,env])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs[[i]] = pwprobs[[i]][which(pwprobs[[i]]$index.x != 0 &
                                            pwprobs[[i]]$index.y != 0), ]
        pwprobs[[i]] = pwprobs[[i]][,-c(4,5)]
      }

      pwprobs.plots = lapply(pwprobs, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('7. Pairwise probability of superior performance within environments estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  probs.df,
            aes(x = .data$env, y = .data$gen, fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Environment", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        pwprobs.plots = lapply(pwprobs, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] > g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })
      }

      ## Save the outputs in a list -----------------

      cond_prob = list(
        df = list(
          perfo = probs.df,
          pair_perfo_env = pwprobs
        ),
        plots = list(
          perfo_env = prob_ggl.plot,
          pair_perfo_env = pwprobs.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(condprobs,
                         file = paste0(getwd(),'/cond_prob/perfo.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_env'))
        for (i in names(pwprobs.env)){
          utils::write.csv(pwprobs.env[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_env/perfo_',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)

    }
  }else{ # Selection for decreasing
    if(!is.null(reg)){# If there are breeding regions

      # Preparations

      stopifnot("Each 'reg' must be represented by a string (e.g., 'R08')" = {
        is.character(data[,reg])
      })

      name.reg = levels(factor(data[,reg]))
      num.reg = nlevels(factor(data[,reg]))
      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                    'Reg',rep(name.reg,  each = num.gen), sep = '_')
      name.env.reg = sort(paste('Env',unique(data[,c(env,reg)])[,1],
                                'Reg',unique(data[,c(env,reg)])[,2], sep = '_'))
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.env),
                                    rep(name.env.reg,  each = num.gen), sep = '_')

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

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g,
                                       y = reorder(.data$gen, -.data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#781c1e')


      # Marginal probabilities ----------------

      ## Probability of superior performance --------------
      ind_post = apply(mod$post$g, 1, function(x){
        ifelse(name.gen %in%
                 unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                 split = '_')), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g, aes(x = factor(.data$ID, levels = .data$ID),
                                       y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90))

      if(verbose) message('1. Probability of superior performance estimated')

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))
      for(i in colnames(mod$post$g)){
        for (j in colnames(mod$post$g)) {
          pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
        }
      }
      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(pwsprob_g1)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      pwsprob_g.plot = ggplot(pwsprob_g1,
                              aes(x = factor(.data$x,
                                             levels = unique(.data$x)),
                                  y = factor(.data$y,
                                             levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white') +
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'plasma')+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(pwsprob_g)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gl),'_Env'))[,1]
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
      prob_gl.plot = ggplot(prob_gl,
                            aes(x = factor(.data$ID, levels = .data$ID),
                                y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

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
        varying = list(colnames(pwsprob_gl1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(pwsprob_gl)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gm),'_Reg'))[,1]
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
      prob_gm.plot = ggplot(prob_gm,
                            aes(x = factor(.data$ID, levels = .data$ID),
                                y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

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
        varying = list(colnames(pwsprob_gm1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gm[upper.tri(pwsprob_gm, diag = T)] = NA
      pwsprob_gm = stats::reshape(
        data.frame(pwsprob_gm),
        direction = 'long',
        varying = list(colnames(pwsprob_gm)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gm), times = colnames(pwsprob_gm),
        new.row.names = 1:length(c(pwsprob_gm)), v.names = 'prob'
      )
      pwsprob_gm = stats::na.exclude(pwsprob_gm[order(pwsprob_gm$x),])

      if(verbose) message('6. Pairwise probability of superior stability (GM) estimated')

      ## Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gm, by = 'ID'))
      j_prob$joint = j_prob$prob.x * j_prob$prob.y
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Environment', 'Region'), each = num.gen)
      j_prob = reshape(j_prob, direction = 'long', varying = list(2:4),
                       times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-5]
      colnames(j_prob) = c('ID', 'level', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$level, j_prob$category,  j_prob$ID),]
      j_prob.plot = ggplot(j_prob, aes(x = .data$ID, y = .data$prob)) +
        facet_wrap(.~.data$level, ncol = 1) +
        geom_segment(data = data.frame(
          ID = name.gen,
          prob = c(apply(merge(prob_g, prob_gl, by = 'ID')[,-1], 1,
                         function(x) x[which.max(x)]),
                   apply(merge(prob_g, prob_gm, by = 'ID')[,-1], 1,
                         function(x) x[which.max(x)])),
          level = rep(c('Environment','Region'), each = num.gen)
        ),
        aes(x = .data$ID, y = 0, yend = .data$prob, xend = .data$ID),
        linewidth = 1.2) +
        geom_point(aes(fill = .data$category, shape = .data$category), size = 2,
                   color = 'black') +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = 'top') +
        scale_fill_manual(
          label = c(
            'Joint' = 'Joint probability',
            'Performance' = 'Superior stability',
            'Stability' = 'Superior performance'
          ),
          values = c(
            'Joint' = '#1b9e77',
            'Performance' = '#d95f02',
            'Stability' = '#7570b3'
          )
        ) +
        scale_shape_manual(label = c(
          'Joint' = 'Joint probability',
          'Performance' = 'Superior stability',
          'Stability' = 'Superior performance'
        ),
        values = c(
          'Joint' = 21,
          'Performance' = 24,
          'Stability' = 25
        ))+
        ylim(0, 1) +
        labs(x = 'Genotype', y = 'Probabilities', fill = 'Probabilities',
             shape = 'Probabilities')

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
            scale_fill_viridis_c(option = 'plasma', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) < g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gm.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gm, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
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
      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl
      for (i in name.reg) {
        posgge[,grep(paste0(i,'$'), do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] =
          posgge[,grep(paste0(i,'$'), do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] +
          matrix(mod$post$gm[,grep(paste0(i,'$'), do.call(rbind, strsplit(colnames(mod$post$gm),'Reg'))[,2])],
                 nrow = num.sim, ncol = num.gen *
                   length(name.env.reg[grep(paste0(i,'$'), do.call(rbind,strsplit(name.env.reg,'Reg'))[,2])]))
      }

      ## Probability of superior performance ----------------
      supprob = function(vector, num.gen, int){
        ifelse(names(vector) %in%
                 names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      }

      probs = apply(do.call(rbind, lapply(
        lapply(
          apply(
            posgge, 1, function(x){
              list(matrix(x, nrow = num.gen, ncol = num.env,
                          dimnames = list(name.gen, name.env.reg)))}
          ),
          Reduce, f = '+'
        ),
        function(x){
          apply(
            x, MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
          )}
      )), 2, function(x){
        tapply(x, rep(name.gen, num.sim), mean)
      })

      probs.df = stats::reshape(
        data = data.frame(probs), direction = 'long',
        varying = list(colnames(probs)),
        ids = name.gen, times = colnames(probs),
        new.row.names = 1:length(c(probs)), v.names = 'prob',
        idvar = 'gen', timevar = 'env'
      )
      probs.df = cbind(probs.df, do.call(rbind, strsplit(probs.df$env, '_Reg_')))
      probs.df = probs.df[,-1]; colnames(probs.df) = c('prob', 'gen', 'loc', 'reg')
      probs.df = probs.df[,c('reg', 'loc', 'gen', 'prob')]
      probs.df$loc = sub('Env_','', probs.df$loc)

      ### Per Location ----------------
      con_gl = ifelse(table(data[,gen], data[,env]) != 0, 1, NA)

      prob_ggl.plot = ggplot(
        data =  merge(
          x = probs.df,
          y = stats::reshape(
            data = data.frame(con_gl), direction = 'long',
            varying = list(colnames(con_gl)),
            ids = name.gen, times = name.env,
            v.names = 'freq',
            idvar = 'gen', timevar = 'loc'
          ), by= c('loc', 'gen')),
        aes(x = .data$loc, y = .data$gen, fill = .data$prob * .data$freq)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Environment", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('8. Probability of superior performance within environments estimated')

      ### Per Region ----------------
      con_gm = ifelse(table(data[,gen], data[,reg]) != 0, 1, NA)

      prob_ggm.plot = ggplot(
        data =  merge(
          x = probs.df,
          y = stats::reshape(
            data = data.frame(con_gm), direction = 'long',
            varying = list(colnames(con_gm)),
            ids = name.gen, times = name.reg,
            v.names = 'freq',
            idvar = 'gen', timevar = 'reg'
          ), by= c('reg', 'gen')),
        aes(x = .data$reg, y = .data$gen, fill = .data$prob * .data$freq)
      )  +
        geom_tile(colour = 'white', na.rm = T) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Region", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('9. Probability of superior performance within regions estimated')

      ### Adjusting the data frame ----------------
      condprobs = merge(
        x = probs.df,
        y = stats::reshape(
          data = data.frame(con_gl), direction = 'long',
          varying = list(colnames(con_gl)),
          ids = name.gen, times = name.env,
          v.names = 'freq',
          idvar = 'gen', timevar = 'loc'
        ), by= c('loc', 'gen'))
      condprobs$prob = condprobs$prob * condprobs$freq
      condprobs = condprobs[,-5]
      colnames(condprobs) = c('env', 'gen', 'reg', 'prob')
      condprobs = condprobs[,c('env','reg','gen','prob')]

      ## Pairwise probability of superior performance ----------------
      ### Per Location -------------
      combs = data.frame(t(utils::combn(paste('Gen', name.gen, sep = '_'), 2)))
      colnames(combs) = c('x', 'y')
      pwprobs.env = lapply(
        sapply(paste('Env', name.env, sep = '_'),
               function(x) posgge[,grep(paste0(x, '_'), colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] <
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )

          a[,1] = sub('Gen_', '', a[,1])
          a[,2] = sub('Gen_', '', a[,2])

          a
        }
      )
      names(pwprobs.env) = sub('Env_', '', names(pwprobs.env))
      for (i in names(pwprobs.env)) {
        pwprobs.env[[i]] = merge(
          merge(pwprobs.env[[i]],
                data.frame(index = table(data[,gen],data[,env])[,i],
                           x = rownames(table(data[,gen],data[,env])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,env])[,i],
                     y = rownames(table(data[,gen],data[,env])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs.env[[i]] = pwprobs.env[[i]][which(pwprobs.env[[i]]$index.x != 0 &
                                                    pwprobs.env[[i]]$index.y != 0), ]
        pwprobs.env[[i]] = pwprobs.env[[i]][,-c(4,5)]
      }

      pwprobs.env.plots = lapply(pwprobs.env, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('10. Pairwise probability of superior performance within environments estimated')

      ### Per Region --------------
      pwprobs.reg = lapply(
        sapply(paste('Reg', name.reg, sep = '_'),
               function(x) posgge[,grep(paste0(x, '$'), colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] <
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )

          a[,1] = sub('Gen_', '', a[,1])
          a[,2] = sub('Gen_', '', a[,2])

          a
        }
      )
      names(pwprobs.reg) = sub('Reg_', '', names(pwprobs.reg))
      for (i in names(pwprobs.reg)) {
        pwprobs.reg[[i]] = merge(
          merge(pwprobs.reg[[i]],
                data.frame(index = table(data[,gen],data[,reg])[,i],
                           x = rownames(table(data[,gen],data[,reg])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,reg])[,i],
                     y = rownames(table(data[,gen],data[,reg])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs.reg[[i]] = pwprobs.reg[[i]][which(pwprobs.reg[[i]]$index.x != 0 &
                                                    pwprobs.reg[[i]]$index.y != 0), ]
        pwprobs.reg[[i]] = pwprobs.reg[[i]][,-c(4,5)]
      }

      pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('11. Pairwise probability of superior performance within regions estimated')


      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  merge(
              x = probs.df,
              y = stats::reshape(
                data = data.frame(con_gl), direction = 'long',
                varying = list(colnames(con_gl)),
                ids = name.gen, times = name.env,
                v.names = 'freq',
                idvar = 'gen', timevar = 'loc'
              ), by= c('loc', 'gen')),
            aes(x = .data$loc, y = .data$gen, fill = .data$prob * .data$freq)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Environment", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        prob_ggm.plot = suppressWarnings(plotly::ggplotly(
          ggplot(
            data =  merge(
              x = probs.df,
              y = stats::reshape(
                data = data.frame(con_gm), direction = 'long',
                varying = list(colnames(con_gm)),
                ids = name.gen, times = name.reg,
                v.names = 'freq',
                idvar = 'gen', timevar = 'reg'
              ), by= c('reg', 'gen')),
            aes(x = .data$reg, y = .data$gen, fill = .data$prob * .data$freq)
          )  +
            geom_tile(colour = 'white', na.rm = T) +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Region", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        ))

        pwprobs.env_plots = lapply(pwprobs.env, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })

        pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })
      }

      ## Save the outputs in a list -----------------

      cond_prob = list(
        df = list(
          perfo = condprobs,
          pair_perfo_env = pwprobs.env,
          pair_perfo_reg = pwprobs.reg
        ),
        plots = list(
          perfo_env = prob_ggl.plot,
          perfo_reg = prob_ggm.plot,
          pair_perfo_env = pwprobs.env.plots,
          pair_perfo_reg = pwprobs.reg.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(condprobs,
                         file = paste0(getwd(),'/cond_prob/perfo.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_env'))
        for (i in names(pwprobs.env)){
          utils::write.csv(pwprobs.env[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_env/perfo_',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_reg'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_reg/perfo_',i,'.csv'),
                           row.names = F)
        }
      }



      # Final output -----------
      if(verbose) message('Process completed!')
      output = list(marginal = marg_prob, conditional = cond_prob)
      return(output)


    }else{ #If there is no breeding region

      # Preparation
      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                    rep(name.env,  each = num.gen), sep = '_')

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

      g_hpd = ggplot(data = g_hpd, aes(x = .data$g, y = reorder(.data$gen, -.data$g))) +
        geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
        geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
        labs(x = 'Genotypic main effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#781c1e')

      # Marginal probabilities ----------------

      ## Probability of superior performance --------------
      ind_post = apply(mod$post$g, 1, function(x){
        ifelse(name.gen %in%
                 unlist(strsplit(names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]),
                                 split = '_')), 1, 0)
      })
      rownames(ind_post) = name.gen
      prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
      prob_g = prob_g[order(prob_g$prob, decreasing = T),]
      prob_g.plot = ggplot(prob_g, aes(x = factor(.data$ID, levels = .data$ID),
                                       y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior performance') +
        theme(axis.text.x = element_text(angle = 90))

      if(verbose) message('1. Probability of superior performance estimated')

      ## Pairwise probability of superior performance ----------------
      pwsprob_g = matrix(NA, num.gen, num.gen,
                         dimnames = list(colnames(mod$post$g), colnames(mod$post$g)))
      for(i in colnames(mod$post$g)){
        for (j in colnames(mod$post$g)) {
          pwsprob_g[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
        }
      }
      colnames(pwsprob_g) = rownames(pwsprob_g) = unlist(strsplit(rownames(pwsprob_g), split = '_'))
      pwsprob_g1 = pwsprob_g[match(prob_g$ID, rownames(pwsprob_g)),
                             match(prob_g$ID, rownames(pwsprob_g))]
      pwsprob_g1[upper.tri(pwsprob_g1, diag = T)] = NA
      pwsprob_g1 = stats::reshape(
        data.frame(pwsprob_g1),
        direction = 'long',
        varying = list(colnames(pwsprob_g1)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g1), times = colnames(pwsprob_g1),
        new.row.names = NULL, v.names = 'prob'
      )

      pwsprob_g.plot = ggplot(pwsprob_g1,
                              aes(x = factor(.data$x,
                                             levels = unique(.data$x)),
                                  y = factor(.data$y,
                                             levels = unique(.data$y))))+
        geom_tile(aes(fill = .data$prob), colour = 'white') +
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')+
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'plasma')+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_g[upper.tri(pwsprob_g, diag = T)] = NA
      pwsprob_g = stats::reshape(
        data.frame(pwsprob_g),
        direction = 'long',
        varying = list(colnames(pwsprob_g)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_g), times = colnames(pwsprob_g),
        new.row.names = 1:length(c(pwsprob_g)), v.names = 'prob'
      )
      pwsprob_g = stats::na.exclude(pwsprob_g[order(pwsprob_g$x),])

      if(verbose) message('2. Pairwise probability of superior performance estimated')

      ## Probability of superior stability  -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gl),'_'))[,1]
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
      prob_gl.plot = ggplot(prob_gl,
                            aes(x = factor(.data$ID, levels = .data$ID),
                                y = .data$prob))+
        geom_bar(stat = 'identity', fill = '#781c1e', color = 'black')+
        labs(x = 'Genotypes', y = 'Probability of superior stability') +
        theme(axis.text.x = element_text(angle = 90))

      if(verbose) message('3. Probability of superior stability estimated')

      ## Pairwise probability of superior stability  -------------
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
        varying = list(colnames(pwsprob_gl1)),
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
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                             option = 'viridis')+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5))

      pwsprob_gl[upper.tri(pwsprob_gl, diag = T)] = NA
      pwsprob_gl = stats::reshape(
        data.frame(pwsprob_gl),
        direction = 'long',
        varying = list(colnames(pwsprob_gl)),
        idvar = 'y', timevar = 'x',
        ids = rownames(pwsprob_gl), times = colnames(pwsprob_gl),
        new.row.names = 1:length(c(pwsprob_gl)), v.names = 'prob'
      )
      pwsprob_gl = stats::na.exclude(pwsprob_gl[order(pwsprob_gl$x),])

      if(verbose) message('4. Pairwise probability of superior stability estimated')


      ## Joint probability of superior performance and stability -----------------
      j_prob = merge(prob_g, prob_gl, by = 'ID')
      j_prob$joint = j_prob$prob.x * j_prob$prob.y
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob = reshape(j_prob, direction = 'long', varying = list(2:4),
                       times = colnames(j_prob)[2:4], v.names = 'value')
      j_prob = j_prob[,-4]
      colnames(j_prob) = c('ID', 'category', 'prob')
      rownames(j_prob) = NULL
      j_prob = j_prob[order(j_prob$category,  j_prob$ID),]
      j_prob.plot = ggplot(j_prob, aes(x = .data$ID, y = .data$prob)) +
        geom_segment(data = data.frame(
          ID = name.gen,
          prob = apply(merge(prob_g, prob_gl, by = 'ID')[,-1], 1,
                       function(x) x[which.max(x)])
        ),
        aes(x = .data$ID, xend = .data$ID, y = 0, yend = .data$prob),
        linewidth = 1.2) +
        geom_point(aes(fill = .data$category, shape = .data$category), size = 2,
                   color = 'black') +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = 'top') +
        scale_fill_manual(
          label = c(
            'Joint' = 'Joint probability',
            'Performance' = 'Superior stability',
            'Stability' = 'Superior performance'
          ),
          values = c(
            'Joint' = '#1b9e77',
            'Performance' = '#d95f02',
            'Stability' = '#7570b3'
          )
        ) +
        scale_shape_manual(label = c(
          'Joint' = 'Joint probability',
          'Performance' = 'Superior stability',
          'Stability' = 'Superior performance'
        ),
        values = c(
          'Joint' = 21,
          'Performance' = 24,
          'Stability' = 25
        ))+
        ylim(0, 1) +
        labs(x = 'Genotype', y = 'Probabilities', fill = 'Probabilities',
             shape = 'Probabilities')

      if(verbose) message('5. Joint probability of superior performance and stability estimated')

      ## Transform into interactive plots -----------

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
            scale_fill_viridis_c(option = 'plasma', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[g(x) < g(y)]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
        pwsprob_gl.plot = plotly::ggplotly(
          ggplot(data = pwsprob_gl, aes(x = factor(.data$x, levels = unique(.data$x)),
                                        y = factor(.data$y, levels = unique(.data$y)),
                                        fill = .data$prob)) +
            geom_tile() +
            scale_fill_viridis_c(option = 'viridis', limits = c(0,1),
                                 na.value = 'white', direction = -1) +
            labs(x = 'Genotypes', y = 'Genotypes', fill = 'Pr[var(gl(x)) < var(gl(y))]') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank())
        )
      }

      ## Save the outputs in a list -----------

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

      ## Save data frames in the work directory -----------
      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          write.csv(marg_prob$df[[i]],
                    file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                    row.names = F)
        }
      }

      # Conditional probabilities ----------------
      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl

      ## Probability of superior performance ----------------
      supprob = function(vector, num.gen, int){
        ifelse(names(vector) %in%
                 names(vector[order(vector, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
      }

      probs = apply(do.call(rbind, lapply(
        lapply(
          apply(
            posgge, 1, function(x){
              list(matrix(x, nrow = num.gen, ncol = num.env,
                          dimnames = list(name.gen, name.env)))}
          ),
          Reduce, f = '+'
        ),
        function(x){
          apply(
            x,MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
          )}
      )), 2, function(x){
        tapply(x, rep(name.gen, num.sim), mean)
      })
      probs = probs * ifelse(table(data[,gen], data[,env]) != 0, 1, NA)

      probs.df = stats::reshape(
        data = data.frame(probs), direction = 'long',
        varying = list(colnames(probs)),
        ids = name.gen, times = colnames(probs),
        new.row.names = 1:length(c(probs)), v.names = 'prob',
        idvar = 'gen', timevar = 'env'
      )
      probs.df = probs.df[,c('env', 'gen', 'prob')]

      prob_ggl.plot = ggplot(
        data =  probs.df,
        aes(x = .data$env, y = .data$gen, fill = .data$prob)
      ) +
        geom_tile(colour = 'white') +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_blank(),
              legend.position = 'top') +
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                             limits = c(0,1), option = 'plasma') +
        labs(x = "Environment", y = 'Genotypes',
             fill = expression(bold(Pr(g %in% Omega)))) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                                     title.position = 'top',
                                     title.hjust = .5))

      if(verbose) message('6. Probability of superior performance within environments estimated')

      ## Pairwise probability of superior performance ----------------
      ### Per Location -------------
      combs = data.frame(t(utils::combn(name.gen, 2)))
      colnames(combs) = c('x', 'y')
      pwprobs = lapply(
        sapply(name.env,
               function(x) posgge[,grep(paste0(x, '$'), colnames(posgge))],
               simplify = F),
        function(y){
          cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(paste0(z[1],'_'), colnames(y))] <
                     y[,grep(paste0(z[2],'_'), colnames(y))])
            })
          )
        }
      )

      for (i in names(pwprobs)) {
        pwprobs[[i]] = merge(
          merge(pwprobs[[i]],
                data.frame(index = table(data[,gen],data[,env])[,i],
                           x = rownames(table(data[,gen],data[,env])),
                           row.names = NULL),
                by = 'x'),
          data.frame(index = table(data[,gen],data[,env])[,i],
                     y = rownames(table(data[,gen],data[,env])),
                     row.names = NULL),
          by = 'y'
        )
        pwprobs[[i]] = pwprobs[[i]][which(pwprobs[[i]]$index.x != 0 &
                                            pwprobs[[i]]$index.y != 0), ]
        pwprobs[[i]] = pwprobs[[i]][,-c(4,5)]
      }

      pwprobs.plots = lapply(pwprobs, function(x){
        ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                               option = 'plasma')+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(verbose) message('7. Pairwise probability of superior performance within environments estimated')

      ## Transform into interactive plots -----------------
      if(interactive){
        prob_ggl.plot = plotly::ggplotly(
          ggplot(
            data =  probs.df,
            aes(x = .data$env, y = .data$gen, fill = .data$prob)
          ) +
            geom_tile(colour = 'white') +
            theme(axis.text.x = element_text(angle = 90),
                  panel.background = element_blank(),
                  legend.position = 'top') +
            scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',
                                 limits = c(0,1), option = 'plasma') +
            labs(x = "Environment", y = 'Genotypes',
                 fill = expression(bold(Pr(g %in% Omega))))
        )

        pwprobs.plots = lapply(pwprobs, function(x){
          plotly::ggplotly(
            ggplot(data = x, aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
              geom_tile() +
              labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[x] < g[y]))))+
              scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1),
                                   option = 'plasma') +
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_blank())
          )
        })
      }

      ## Save the outputs in a list -----------------

      cond_prob = list(
        df = list(
          perfo = probs.df,
          pair_perfo_env = pwprobs
        ),
        plots = list(
          perfo_env = prob_ggl.plot,
          pair_perfo_env = pwprobs.plots
        )
      )

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/cond_prob'))
        utils::write.csv(condprobs,
                         file = paste0(getwd(),'/cond_prob/perfo.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/cond_prob/pairwise_env'))
        for (i in names(pwprobs.env)){
          utils::write.csv(pwprobs.env[[i]],
                           file = paste0(getwd(),'/cond_prob/pairwise_env/perfo_',i,'.csv'),
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
