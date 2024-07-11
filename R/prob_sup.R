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
##' @param mod.output An object of class `extr`, obtained from the [extr_outs()] function
##' @param int A number representing the selection intensity
##' (between 0 and 1)
##' @param increase Logical. Indicates the direction of the selection.
##' `TRUE` (default) for increasing the trait value, `FALSE` otherwise.
##' @param save.df Logical. Should the data frames be saved in the work directory?
##' `TRUE` for saving, `FALSE` (default) otherwise.
##' @param verbose A logical value. If `TRUE`, the function will indicate the
##' completed steps. Defaults to `FALSE`.
##'
##' @return The function returns an object of class `probsup`, which contains two lists,
##' one with the `across`-environment probabilities, and another with the `within`-environments probabilities.
##'
##' The `across` list has the following elements:
##' \itemize{
##' \item \code{g_hpd}: Highest posterior density (HPD) of the posterior genotypic main effects.
##' \item \code{perfo}: the probabilities of superior performance.
##' \item \code{pair_perfo}: the pairwise probabilities of superior performance.
##' \item \code{stabi}: a list with the probabilities of superior stability. It can contain the data frames `gl`,
##' `gm` (when `reg` is not `NULL`) and `gt` (when `year` is not `NULL`).
##' \item \code{pair_stabi}: a list with the pairwise probabilities of superior stability. It can contain the data frames `gl`,
##' `gm` (when `reg` is not `NULL`) and `gt` (when `year` is not `NULL`).
##' \item \code{joint_prob}: the joint probabilities of superior performance and stability.
##' }
##'
##' The `within` list has the following elements:
##' \itemize{
##' \item \code{perfo}: a list of data frames containing the probabilities of superior performance
##' within locations (`gl`), regions (`gm`) and years (`gt`).
##' \item \code{pair_perfo}: lists with the pairwise probabilities of superior performance
##' within locations (`gl`), regions (`gm`) and years (`gt`).
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
##' mod = bayes_met(data = soy,
##'                 gen = "Gen",
##'                 loc = "Loc",
##'                 repl = NULL,
##'                 year = NULL,
##'                 reg = NULL,
##'                 res.het = FALSE,
##'                 trait = 'Y',
##'                 iter = 6000, cores = 4, chains = 4)
##'
##' outs = extr_outs(data = soy, trait = "Y", model = mod,
##'                  probs = c(0.05, 0.95), plots = TRUE,
##'                  verbose = TRUE)
##'
##' results = prob_sup(data = soy,
##'                    trait = "Y",
##'                    gen = "Gen",
##'                    loc = "Loc",
##'                    reg = NULL,
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
                    increase = TRUE, save.df = FALSE, verbose = FALSE){

  stopifnot("Please, provide a valid selection intensity (number between 0 and 1)" = {
    is.numeric(int)
    int >= 0 & int <=1
  })

  # Namespaces
  requireNamespace('ggplot2')

  # Preparation
  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = mod.output
  output = list(across = list(), within = list())
  class(output) = "probsup"
  attr(output, "control") = data.frame(increase = increase, trait = trait, gen = gen,
                                       loc = loc, reg = ifelse(is.null(reg), 0, reg),
                                       year = ifelse(is.null(year), 0, year))

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
      if(any(colSums(Z2) == 0)){
        mod$post$gl = mod$post$gl[,-which(colSums(Z2) == 0)]
        Z2 = Z2[,-which(colSums(Z2) == 0)]
      }
      Z3 = stats::model.matrix(~-1 + aux[,gen]:aux[,year])
      if(any(colSums(Z3) == 0)){
        mod$post$gt = mod$post$gt[,-which(colSums(Z3) == 0)]
        Z3 = Z3[,-which(colSums(Z3) == 0)]
      }
      Z4 = stats::model.matrix(~-1 + aux[,gen]:aux[,reg])
      if(any(colSums(Z4) == 0)){
        mod$post$gm = mod$post$gm[,-which(colSums(Z4) == 0)]
        Z4 = Z4[,-which(colSums(Z4) == 0)]
      }

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        HPD95 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        HPD97.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        HPD5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        HPD7.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )
      output$across$g_hpd = g_hpd


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

      if(verbose) message('-> Probability of superior performance estimated')

      output$across$perfo = prob_g

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
      output$across$pair_perfo = pwsprob_g

      if(verbose) message('-> Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      if(any(table(colnames(staprob_gl)) == 1)){
        oncegeno =  names(table(colnames(staprob_gl))[which(table(colnames(staprob_gl)) == 1)])
        staprob_gl = staprob_gl[,-which(colnames(staprob_gl) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one location (environment), so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gl = prob_gl

      if(verbose) message('-> Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }

      output$across$pair_stabi$gl = pwsprob_gl

      if(verbose) message('-> Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gm),'_@#Reg'))[,1]
      )
      if(any(table(colnames(staprob_gm)) == 1)){
        oncegeno =  names(table(colnames(staprob_gm))[which(table(colnames(staprob_gm)) == 1)])
        staprob_gm = staprob_gm[,-which(colnames(staprob_gm) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one region, so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gm = prob_gm

      if(verbose) message('-> Probability of superior stability (GM) estimated')


      ## Pairwise probability of superior stability - Region -----------------
      pwsprob_gm = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gm)){
        for (j in colnames(pwsprob_gm)) {
          pwsprob_gm[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      output$across$pair_stabi$gm = pwsprob_gm

      if(verbose) message('-> Pairwise probability of superior stability (GM) estimated')

      ## Probability of superior stability - year --------------
      staprob_gt = mod$post$gt
      colnames(staprob_gt) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gt),'_@#year'))[,1]
      )
      if(any(table(colnames(staprob_gt)) == 1)){
        oncegeno =  names(table(colnames(staprob_gt))[which(table(colnames(staprob_gt)) == 1)])
        staprob_gt = staprob_gt[,-which(colnames(staprob_gt) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one year, so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gt = prob_gt

      if(verbose) message('-> Probability of superior stability (GT) estimated')


      ## Pairwise probability of superior stability - year -----------------
      pwsprob_gt = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gt)){
        for (j in colnames(pwsprob_gt)) {
          pwsprob_gt[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      output$across$pair_stabi$gt = pwsprob_gt

      if(verbose) message('-> Pairwise probability of superior stability (GT) estimated')

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

      output$across$joint = j_prob

      if(verbose) message('-> Joint probability of superior performance and stability estimated')

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/across_probs'))
        utils::write.csv(output$across$pair_perfo,
                         file = paste0(getwd(),'/across_probs/pair_perfo.csv'),
                         row.names = T)
        for (i in names(output$across)[-grep("stabi|pair", names(output$across))]){
          utils::write.csv(output$across[[i]],
                           file = paste0(getwd(),'/across_probs/',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$stabi)){
          utils::write.csv(output$across$stabi[[i]],
                           file = paste0(getwd(),'/across_probs/stabi_',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$pair_stabi)){
          utils::write.csv(output$across$pair_stabi[[i]],
                           file = paste0(getwd(),'/across_probs/pair_stabi_',i,'.csv'),
                           row.names = T)
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

      output$within$perfo$gl = prob_ggl

      if(verbose) message('-> Probability of superior performance within locations estimated')

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
      }

      output$within$pair_perfo$gl = pwprobs.loc

      if(verbose) message('-> Pairwise probability of superior performance within locations estimated')

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

      output$within$perfo$gm = prob_ggm

      if(verbose) message('-> Probability of superior performance within regions estimated')

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

      }

      output$within$pair_perfo$gm = pwprobs.reg

      if(verbose) message('-> Pairwise probability of superior performance within regions estimated')

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

      output$within$perfo$gt = prob_ggt

      if(verbose) message('-> Probability of superior performance within year estimated')

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

      }

      output$within$pair_perfo$gt = pwprobs.year

      if(verbose) message('-> Pairwise probability of superior performance within year estimated')


      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/within'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/within/prob_ggl.csv'),
                         row.names = F)
        utils::write.csv(prob_ggm,
                         file = paste0(getwd(),'/within/prob_ggm.csv'),
                         row.names = F)
        utils::write.csv(prob_ggt,
                         file = paste0(getwd(),'/within/prob_ggt.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/within/pairwise_ggl'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggl/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/within/pairwise_ggm'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggm/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/within/pairwise_ggt'))
        for (i in names(pwprobs.year)){
          utils::write.csv(pwprobs.year[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggt/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      return(output)

    }
    else # Without region info --------------------
    {
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.loc),
                                    "loc", rep(name.loc,  each = num.gen), sep = '_@#')

      aux = unique(data[,c(gen,loc,year)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])
      if(any(colSums(Z2) == 0)){
        mod$post$gl = mod$post$gl[,-which(colSums(Z2) == 0)]
        Z2 = Z2[,-which(colSums(Z2) == 0)]
      }
      Z3 = stats::model.matrix(~-1 + aux[,gen]:aux[,year])
      if(any(colSums(Z3) == 0)){
        mod$post$gt = mod$post$gt[,-which(colSums(Z3) == 0)]
        Z3 = Z3[,-which(colSums(Z3) == 0)]
      }

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        HPD95 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        HPD97.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        HPD5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        HPD7.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )
      output$across$g_hpd = g_hpd


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

      if(verbose) message('-> Probability of superior performance estimated')

      output$across$perfo = prob_g

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
      output$across$pair_perfo = pwsprob_g

      if(verbose) message('-> Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      if(any(table(colnames(staprob_gl)) == 1)){
        oncegeno =  names(table(colnames(staprob_gl))[which(table(colnames(staprob_gl)) == 1)])
        staprob_gl = staprob_gl[,-which(colnames(staprob_gl) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one location (environment), so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gl = prob_gl

      if(verbose) message('-> Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }

      output$across$pair_stabi$gl = pwsprob_gl

      if(verbose) message('-> Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - year --------------
      staprob_gt = mod$post$gt
      colnames(staprob_gt) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gt),'_@#year'))[,1]
      )
      if(any(table(colnames(staprob_gt)) == 1)){
        oncegeno =  names(table(colnames(staprob_gt))[which(table(colnames(staprob_gt)) == 1)])
        staprob_gt = staprob_gt[,-which(colnames(staprob_gt) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one year, so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gt = prob_gt

      if(verbose) message('-> Probability of superior stability (GT) estimated')


      ## Pairwise probability of superior stability - year -----------------
      pwsprob_gt = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gt)){
        for (j in colnames(pwsprob_gt)) {
          pwsprob_gt[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      output$across$pair_stabi$gt = pwsprob_gt

      if(verbose) message('-> Pairwise probability of superior stability (GT) estimated')

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

      output$across$joint = j_prob

      if(verbose) message('-> Joint probability of superior performance and stability estimated')

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/across_probs'))
        utils::write.csv(output$across$pair_perfo,
                         file = paste0(getwd(),'/across_probs/pair_perfo.csv'),
                         row.names = T)
        for (i in names(output$across)[-grep("stabi|pair", names(output$across))]){
          utils::write.csv(output$across[[i]],
                           file = paste0(getwd(),'/across_probs/',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$stabi)){
          utils::write.csv(output$across$stabi[[i]],
                           file = paste0(getwd(),'/across_probs/stabi_',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$pair_stabi)){
          utils::write.csv(output$across$pair_stabi[[i]],
                           file = paste0(getwd(),'/across_probs/pair_stabi_',i,'.csv'),
                           row.names = T)
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

      output$within$perfo$gl = prob_ggl

      if(verbose) message('-> Probability of superior performance within locations estimated')

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
      }

      output$within$pair_perfo$gl = pwprobs.loc

      if(verbose) message('-> Pairwise probability of superior performance within locations estimated')

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

      output$within$perfo$gt = prob_ggt

      if(verbose) message('-> Probability of superior performance within year estimated')

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

      }

      output$within$pair_perfo$gt = pwprobs.year

      if(verbose) message('-> Pairwise probability of superior performance within year estimated')


      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/within'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/within/prob_ggl.csv'),
                         row.names = F)
        utils::write.csv(prob_ggt,
                         file = paste0(getwd(),'/within/prob_ggt.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/within/pairwise_ggl'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggl/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/within/pairwise_ggt'))
        for (i in names(pwprobs.year)){
          utils::write.csv(pwprobs.year[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggt/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
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
      if(any(colSums(Z2) == 0)){
        mod$post$gl = mod$post$gl[,-which(colSums(Z2) == 0)]
        Z2 = Z2[,-which(colSums(Z2) == 0)]
      }
      Z4 = stats::model.matrix(~-1 + aux[,gen]:aux[,reg])
      if(any(colSums(Z4) == 0)){
        mod$post$gm = mod$post$gm[,-which(colSums(Z4) == 0)]
        Z4 = Z4[,-which(colSums(Z4) == 0)]
      }

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        HPD95 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        HPD97.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        HPD5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        HPD7.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )
      output$across$g_hpd = g_hpd


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

      if(verbose) message('-> Probability of superior performance estimated')

      output$across$perfo = prob_g

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
      output$across$pair_perfo = pwsprob_g

      if(verbose) message('-> Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      if(any(table(colnames(staprob_gl)) == 1)){
        oncegeno =  names(table(colnames(staprob_gl))[which(table(colnames(staprob_gl)) == 1)])
        staprob_gl = staprob_gl[,-which(colnames(staprob_gl) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one location (environment), so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gl = prob_gl

      if(verbose) message('-> Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }

      output$across$pair_stabi$gl = pwsprob_gl

      if(verbose) message('-> Pairwise probability of superior stability (GL) estimated')

      ## Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gm),'_@#Reg'))[,1]
      )
      if(any(table(colnames(staprob_gm)) == 1)){
        oncegeno =  names(table(colnames(staprob_gm))[which(table(colnames(staprob_gm)) == 1)])
        staprob_gm = staprob_gm[,-which(colnames(staprob_gm) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one region, so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gm = prob_gm

      if(verbose) message('-> Probability of superior stability (GM) estimated')


      ## Pairwise probability of superior stability - Region -----------------
      pwsprob_gm = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gm)){
        for (j in colnames(pwsprob_gm)) {
          pwsprob_gm[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }
      output$across$pair_stabi$gm = pwsprob_gm

      if(verbose) message('-> Pairwise probability of superior stability (GM) estimated')

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

      output$across$joint = j_prob

      if(verbose) message('-> Joint probability of superior performance and stability estimated')

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/across_probs'))
        utils::write.csv(output$across$pair_perfo,
                         file = paste0(getwd(),'/across_probs/pair_perfo.csv'),
                         row.names = T)
        for (i in names(output$across)[-grep("stabi|pair", names(output$across))]){
          utils::write.csv(output$across[[i]],
                           file = paste0(getwd(),'/across_probs/',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$stabi)){
          utils::write.csv(output$across$stabi[[i]],
                           file = paste0(getwd(),'/across_probs/stabi_',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$pair_stabi)){
          utils::write.csv(output$across$pair_stabi[[i]],
                           file = paste0(getwd(),'/across_probs/pair_stabi_',i,'.csv'),
                           row.names = T)
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

      output$within$perfo$gl = prob_ggl

      if(verbose) message('-> Probability of superior performance within locations estimated')

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
      }

      output$within$pair_perfo$gl = pwprobs.loc

      if(verbose) message('-> Pairwise probability of superior performance within locations estimated')

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

      output$within$perfo$gm = prob_ggm

      if(verbose) message('-> Probability of superior performance within regions estimated')

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

      }

      output$within$pair_perfo$gm = pwprobs.reg

      if(verbose) message('-> Pairwise probability of superior performance within regions estimated')

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/within'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/within/prob_ggl.csv'),
                         row.names = F)
        utils::write.csv(prob_ggm,
                         file = paste0(getwd(),'/within/prob_ggm.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/within/pairwise_ggl'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggl/',i,'.csv'),
                           row.names = F)
        }
        dir.create(path = paste0(getwd(),'/within/pairwise_ggm'))
        for (i in names(pwprobs.reg)){
          utils::write.csv(pwprobs.reg[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggm/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      return(output)

    }
    else # Without region info --------------------
    {
      aux = unique(data[,c(gen,loc)])
      Z1 = stats::model.matrix(~-1 + aux[,gen])
      colnames(Z1) = gsub("aux\\[, gen\\]",'', colnames(Z1))
      Z2 = stats::model.matrix(~-1 + aux[,gen]:aux[,loc])
      if(any(colSums(Z2) == 0 & ncol(Z2) == ncol(mod$post$gl))){
        mod$post$gl = mod$post$gl[,-which(colSums(Z2) == 0)]
        Z2 = Z2[,-which(colSums(Z2) == 0)]
      }else if(any(colSums(Z2) == 0)){
        Z2 = Z2[,-which(colSums(Z2) == 0)]
      }
      colnames(Z2) = gsub("aux\\[, gen\\]",'', gsub("aux\\[, loc\\]",'_@#', colnames(Z2)))
      colnames(mod$post$gl) = paste("Gen",aux[,1], "loc", aux[,2], sep = "_@#")

      # Genotypic effects and their HPD ------------
      g_hpd = data.frame(
        gen = name.gen,
        g = apply(mod$post$g, 2, stats::median),
        HPD95 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.95)),
        HPD97.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.975)),
        HPD5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.05)),
        HPD7.5 = apply(mod$post$g, 2, function(x) stats::quantile(x, probs = 0.025)),
        row.names = NULL
      )
      output$across$g_hpd = g_hpd


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

      if(verbose) message('-> Probability of superior performance estimated')

      output$across$perfo = prob_g

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
      output$across$pair_perfo = pwsprob_g

      if(verbose) message('-> Pairwise probability of superior performance estimated')

      ## Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_@#','',do.call(rbind,strsplit(colnames(staprob_gl),'_@#loc'))[,1]
      )
      if(any(table(colnames(staprob_gl)) == 1)){
        oncegeno =  names(table(colnames(staprob_gl))[which(table(colnames(staprob_gl)) == 1)])
        staprob_gl = staprob_gl[,-which(colnames(staprob_gl) %in% oncegeno)]
        warning("Some genotypes were evaluated in only one location (environment), so we could not compute their stability: ",
                paste(oncegeno, collapse = ", "))
      }
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

      output$across$stabi$gl = prob_gl

      if(verbose) message('-> Probability of superior stability (GL) estimated')

      ## Pairwise probability of superior stability - Location -------------
      pwsprob_gl = matrix(NA, num.gen, num.gen,
                          dimnames = list(colnames(probsta), colnames(probsta)))
      for(i in colnames(pwsprob_gl)){
        for (j in colnames(pwsprob_gl)) {
          pwsprob_gl[i,j] = mean(probsta[,j] < probsta[,i])
        }
      }

      output$across$pair_stabi$gl = pwsprob_gl

      if(verbose) message('-> Pairwise probability of superior stability (GL) estimated')

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

      output$across$joint = j_prob

      if(verbose) message('-> Joint probability of superior performance and stability estimated')

      ## Save data frames in the work directory -----------------
      if(save.df){
        dir.create(path = paste0(getwd(),'/across_probs'))
        utils::write.csv(output$across$pair_perfo,
                         file = paste0(getwd(),'/across_probs/pair_perfo.csv'),
                         row.names = T)
        for (i in names(output$across)[-grep("stabi|pair", names(output$across))]){
          utils::write.csv(output$across[[i]],
                           file = paste0(getwd(),'/across_probs/',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$stabi)){
          utils::write.csv(output$across$stabi[[i]],
                           file = paste0(getwd(),'/across_probs/stabi_',i,'.csv'),
                           row.names = F)
        }
        for (i in names(output$across$pair_stabi)){
          utils::write.csv(output$across$pair_stabi[[i]],
                           file = paste0(getwd(),'/across_probs/pair_stabi_',i,'.csv'),
                           row.names = T)
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

      output$within$perfo$gl = prob_ggl

      if(verbose) message('-> Probability of superior performance within locations estimated')

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
      }

      output$within$pair_perfo$gl = pwprobs.loc

      if(verbose) message('-> Pairwise probability of superior performance within locations estimated')

      ## Save data frames in the work directory -----------------

      if(save.df){
        dir.create(path = paste0(getwd(),'/within'))
        utils::write.csv(prob_ggl,
                         file = paste0(getwd(),'/within/prob_ggl.csv'),
                         row.names = F)
        dir.create(path = paste0(getwd(),'/within/pairwise_ggl'))
        for (i in names(pwprobs.loc)){
          utils::write.csv(pwprobs.loc[[i]],
                           file = paste0(getwd(),'/within/pairwise_ggl/',i,'.csv'),
                           row.names = F)
        }
      }

      # Final output -----------
      if(verbose) message('Process completed!')
      return(output)

    }
  }
}


#' Plots for the `probsup` object
#'
#' Build plots using the outputs stored in the `probsup` object.
#'
#'
#' @param object An object of class `probsup`.
#' @param category A string indicating which plot to build. See options in the Details section.
#' @param level A string indicating the information level to be used for building
#' the plots. Options are `"across"` for focusing on the probabilities across environments,
#' or `"within"` to focus on the within-environment effects. Defaults `"across"`.
#' @param ... Currently not used.
#'
#' @method plot probsup
#'
#'
#' @details The available options are:
##'   \itemize{
##'     \item \code{hpd} : a caterpillar plot representing the marginal genotypic value of
##'       each genotype, and their respective highest posterior density interval (95% represented by the
##'       thick line, and 97.5% represented by the thin line). Available only if `level = "across"`.
##'     \item \code{perfo} : if `level = "across"`, a bar plot illustrating the probabilities of superior performance.
##'       If `level = "within"`, heatmaps with the probabilities of superior performance within
##'       environments. Can be `prob_loc`, `prob_reg` (if `reg` is not `NULL`), and
##'       `prob_year` (if `year` is not `NULL`).
##'     \item \code{stabi}: a bar plot with the probabilities of superior stability.
##'       Different plots are generated (`stabi_gl`, `stabi_gm` and `stabi_gt`) if `reg`
##'       or/and `year` are not `NULL`. Available of if `level = "across"`.
##'     \item \code{pair_perfo} : if `level = "across"`, a heatmap representing the pairwise probability of superior
##'       performance (the probability of genotypes at the \emph{x}-axis being superior
##'       to those on the \emph{y}-axis). If `level = "within"`, a list of heatmaps representing the pairwise probability of superior
##'       performance within environments. Can be `pwprob_loc`, `pwprob_reg` (if `reg` is not `NULL`), and
##'       `pwprob_year` (if `year` is not `NULL`).
##'     \item \code{pair_stabi}: a heatmap with the pairwise probabilities of superior stability.
##'       Different plots are generated (`stabi_gl`, `stabi_gm` and `stabi_gt`) if `reg`
##'       or/and `year` are not `NULL`. This plot represents
##'       the probability of genotypes at the \emph{x}-axis being superior
##'       to those on \emph{y}-axis. Available of if `level = "across"`.
##'     \item \code{joint_prob}: a plot with the joint probabilities of superior
##'       performance and stability.
##'   }
#'
#'
#' @seealso  [ggplot2], [ProbBreed::prob_sup]
#'
##' @import ggplot2
##' @importFrom utils write.csv combn
##' @importFrom stats reshape median quantile na.exclude model.matrix aggregate
##' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#'
#'  }





#'
#' #' Print an object of class `probsup`
#' #'
#' #' Print a `probsup` object in R console
#' #'
#' #' @param obj An object of class `probsup`
#' #' @method print probsup
#' #'
#' #' @seealso [ProbBreed::prob_sup]
#' #'
#' #' @export
#' #'
#'
#' print.extr = function(obj){
#'
#'   o
#'
#' }
