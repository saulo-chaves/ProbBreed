## Function marg_prob
##
##' @title
##' Marginal probabilities of superior performance
##' @description
##' This function estimates the marginal probabilities of superior performance and
##' yields plots for better visualization. All the plots are customizable, using `ggplot2`
##'
##' @details
##' Using `marg_prob()`, you can estimate the probability of a genotype being amongst the
##' selected, based on a given selection intensity. If we let \eqn{\Omega} represent the
##' subset of superior genotypes, then the probability of the <i>jth</i> genotype belonging to
##' \eqn{\Omega} is:
##'
##' \deqn{Pr(g_j \in \Omega \vert y) = \frac{1}{S}\sum_{s=1}^S{I(g_j^{(s)} \in \Omega \vert y)}}
##'
##' where \eqn{S} is the total number of samples, and \eqn{I(g_j^{(s)} \in \Omega \vert y)}
##' is an indicator variable mapping success (1) if \eqn{g_j^{(s)}} exists in \eqn{\Omega},
##' and failure (0) otherwise.
##'
##' This function also estimates the pairwise probability of superior performance, i.e.
##' the probability of a given genotype <i>j</i> being superior to genotype <i>i</i>:
##'
##' \deqn{Pr(g_j > g_i \vert y) = \frac{1}{S}\sum_{s=1}^S{I(g_j^{(s)} > g_i^{(s)} \vert y)}}
##'
##' @references
##'
##' Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães, L. J. M.,
##' Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging probability concepts
##' for cultivar recommendation in multi-environment trials. <i>Theoretical and
##' Applied Genetics</i>, 133(2):443-455. https://doi.org/10.1007/s00122-022-04041-y
##'
##'
##' @param data A dataframe containing the observations
##' @param trait A character representing the name of the column that
##' corresponds to the analysed variable
##' @param gen A character representing the name of the column that
##' corresponds to the evaluated genotypes
##' @param env A character representing the name of the column that
##' corresponds to the environments
##' @param extr_outs An object from the `extr_outs` function
##' @param int An integer representing the selection intensity
##' (superior limit = 1)
##' @param increase Logical: TRUE if genotypes with higher trait values are better.
##' FALSE otherwise.
##' @param save.df A logical value indicating if the data frames with the marginal
##' probability of each genotype and the pairwise probabilities should be saved at
##' the work directory.
##' @param interactive A logical value indicating if the plots should be interactive.
##' If TRUE, the function loads the `plotly` package and uses the [plotly::ggplotly]
##' command.
##' @return The function returns a list with:
##' \itemize{
##' \item \code{g_hpd} : a caterpillar plot representing the performance of
##' each genotype, and their respective confidence interval (95% (thick) and 97.5% (thin))
##' \item \code{marg_prob.df} : a dataframe with the marginal probability of superior
##' performance of each genotype
##' \item \code{marg_prob.plot} : a bar plot with the marginal probability of superior
##' performance of each genotype
##' \item \code{pair_prob.df} : a dataframe with the pairwise marginal probability
##' of superior performance between genotypes
##' \item \code{pair_prob.plot} a heatmap representing the pairwise probability of
##' superior performance
##' }
##'
##' @seealso [ggplot2::ggplot()]
##'
##' @import ggplot2 dplyr
##' @importFrom utils write.csv
##' @importFrom tidyr pivot_longer
##' @importFrom tidyr separate
##' @importFrom tibble rownames_to_column
##' @importFrom rlang .data
##'
##' @export
##'
##' @examples
##' \dontrun{
##' mod = bayes_met(data = maize, gen = c("Hybrid", "normal", "cauchy"),
##'                 env = c("Location", "normal", "cauchy"),
##'                 rept = list(c("Rep", "normal", "cauchy"), c("Block", "normal", "cauchy")),
##'                 trait = "GY", hyperparam = "default", sigma.dist = c("cauchy", "cauchy"),
##'                 mu.dist = c("normal", "normal"), gei.dist = c("normal", "cauchy"),
##'                 reg = list(c("Region", "normal", "cauchy"), c("normal", "cauchy")),
##'                 res.het = T,
##'                 iter = 100, cores = 2, chain = 2)
##'                 #Remember, increase the number of iterations, cores and chains
##'
##' outs = extr_outs(data = maize, trait = "GY", gen = "Hybrid", model = mod,
##'                  effects = c("r", "b", "l", "m", "g", "gl", "gm"),
##'                  nenv = 16, res.het = T, check.stan.diag = TRUE)
##'
##' margs = marg_prob(data = maize, trait = "GY", gen = "Hybrid",env = "Location",
##'                   increase = TRUE,extr_outs = outs, int = .2,
##'                   save.df = TRUE, interactive = TRUE)
##'                   }
##'

marg_prob = function(data, trait, gen, env, extr_outs, int = .2, increase = T,
                     save.df = FALSE, interactive = FALSE){

  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = extr_outs
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))

  g_hpd = as.data.frame(matrix(mod$post$g,
                               dimnames=list(t(outer(colnames(mod$post$g),
                                                     rownames(mod$post$g),
                                                     FUN=paste)), NULL))) %>%
    dplyr::mutate(gen = rep(name.gen, each = nrow(mod$post$g))) %>%
    dplyr::group_by(.data$gen) %>%
    summarise(
      g = stats::median(.data$V1),
      UP = stats::quantile(.data$V1, probs = 0.95),
      up = stats::quantile(.data$V1, probs = 0.975),
      DOWN = stats::quantile(.data$V1, probs = 0.05),
      down = stats::quantile(.data$V1, probs = 0.025)
    ) %>%
    ggplot(aes(x = .data$g, y = reorder(.data$gen, .data$g))) +
    geom_errorbar(aes(xmin = .data$down, xmax = .data$up), width = 0)+
    geom_errorbar(aes(xmin = .data$DOWN, xmax = .data$UP), width = 0, linewidth = 2, alpha = .8) +
    labs(x = 'Values', y = 'Genotypes') +
    geom_point(size = 4, color = '#781c1e')

  if(increase){
  colnames(mod$post$g) = name.gen
  ind_post = apply(mod$post$g, 1, function(x){
    ifelse(name.gen %in%
             names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
  })
  rownames(ind_post) = name.gen

  prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
  prob_g = prob_g[order(prob_g$prob, decreasing = T),]
  prob_g$ID = factor(prob_g$ID, levels = prob_g$ID)

  psps.bar = ggplot(prob_g, aes(x = .data$ID, y = .data$prob))+
    geom_bar(stat = 'identity', fill = '#781c1e')+
    labs(x = 'Genotypes', y = 'Probability') +
    theme(axis.text.x = element_text(angle = 90))


  pwsprob = matrix(NA, num.gen, num.gen)

  for(i in 1:ncol(mod$post$g)){
    for (j in 1:ncol(mod$post$g)) {
      pwsprob[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
    }
  }
  colnames(pwsprob) = rownames(pwsprob) = name.gen
  pwsprob1 = pwsprob[match(prob_g$ID, rownames(pwsprob)),
                    match(prob_g$ID, rownames(pwsprob))]
  diag(pwsprob1) = NA
  pwsprob1[upper.tri(pwsprob1)] = NA
  pwsprob1 = as.data.frame(pwsprob1) %>%
    tibble::rownames_to_column(var = 'gen') %>%
    tidyr::pivot_longer(cols = c(colnames(pwsprob1)[1]:colnames(pwsprob1)[length(colnames(pwsprob1))])) %>%
    dplyr::rename('gen2' = .data$name, 'prob' = .data$value)

  pwsprob1$gen = factor(pwsprob1$gen, levels = unique(pwsprob1$gen))
  pwsprob1$gen2 = factor(pwsprob1$gen2, levels = unique(pwsprob1$gen))

  pwsp.heat = ggplot(pwsprob1, aes(x = .data$gen, y = .data$gen2))+
    geom_tile(aes(fill = .data$prob), colour = 'white')+
    labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
    theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
          legend.position = c(.8,.15), legend.direction = 'horizontal')+
    coord_flip()+
    scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                 title.hjust = .5))

  pwsprob = pwsprob[order(rownames(pwsprob)), order(colnames(pwsprob))]
  diag(pwsprob) = NA
  pwsprob[upper.tri(pwsprob)] = NA
  pwsprob = as.data.frame(pwsprob) %>%
    tibble::rownames_to_column(var = 'gen') %>%
    tidyr::pivot_longer(cols = c(colnames(pwsprob)[1]:colnames(pwsprob)[length(colnames(pwsprob))])) %>%
    dplyr::rename('gen2' = .data$name, 'prob' = .data$value)


  if(save.df){
    utils::write.csv(prob_g, file = paste0(getwd(),'/probsp.csv'), row.names = F)
    utils::write.csv(stats::na.exclude(pwsprob), file = paste0(getwd(),'/pwsprob.csv'), row.names = F)
  }

  if(interactive){

    g_hpd = plotly::ggplotly(g_hpd)
    psps.bar = plotly::ggplotly(ggplot(prob_g, aes(x = .data$ID, y = .data$prob))+
                          geom_bar(stat = 'identity', fill = '#781c1e')+
                          labs(x = 'Genotypes', y = 'Probability') +
                          theme(axis.text.x = element_text(angle = 90)))
    pwsp.heat = plotly::ggplotly(ggplot(pwsprob, aes(x = .data$gen, y = .data$gen2))+
                           geom_tile(aes(fill = .data$prob), colour = 'white')+
                           labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
                           theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                                 legend.position = 'right', legend.direction = 'vertical')+
                           coord_flip()+
                           scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1)))

  }

    out = list(g_hpd, prob_g, psps.bar, na.exclude(pwsprob), pwsp.heat)
    names(out) = c('g_hpd', 'marg_prob.df', 'marg_prob.plot', 'pair_prob.df', 'pair_prob.plot')

  return(out)

  }else{
    colnames(mod$post$g) = name.gen
    ind_post = apply(mod$post$g, 1, function(x){
      ifelse(name.gen %in%
               names(x[order(x, decreasing = F)][1:ceiling(int * num.gen)]), 1, 0)
    })
    rownames(ind_post) = name.gen

    prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
    prob_g = prob_g[order(prob_g$prob, decreasing = T),]
    prob_g$ID = factor(prob_g$ID, levels = prob_g$ID)

    psps.bar = ggplot(prob_g, aes(x = .data$ID, y = .data$prob))+
      geom_bar(stat = 'identity', fill = '#781c1e')+
      labs(x = 'Genotypes', y = 'Probability') +
      theme(axis.text.x = element_text(angle = 90))


    pwsprob = matrix(NA, num.gen, num.gen)

    for(i in 1:ncol(mod$post$g)){
      for (j in 1:ncol(mod$post$g)) {
        pwsprob[i,j] = mean(mod$post$g[,j] < mod$post$g[,i])
      }
    }
    colnames(pwsprob) = rownames(pwsprob) = name.gen
    pwsprob1 = pwsprob[match(prob_g$ID, rownames(pwsprob)),
                      match(prob_g$ID, rownames(pwsprob))]
    diag(pwsprob1) = NA
    pwsprob1[upper.tri(pwsprob1)] = NA
    pwsprob1 = as.data.frame(pwsprob1) %>%
      tibble::rownames_to_column(var = 'gen') %>%
      tidyr::pivot_longer(cols = c(colnames(pwsprob1)[1]:colnames(pwsprob1)[length(colnames(pwsprob1))])) %>%
      dplyr::rename('gen2' = .data$name, 'prob' = .data$value)

    pwsprob1$gen = factor(pwsprob1$gen, levels = unique(pwsprob1$gen))
    pwsprob1$gen2 = factor(pwsprob1$gen2, levels = unique(pwsprob1$gen))

    pwsp.heat = ggplot(pwsprob1, aes(x = .data$gen, y = .data$gen2))+
      geom_tile(aes(fill = .data$prob), colour = 'white')+
      labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
      theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
            legend.position = c(.8,.15), legend.direction = 'horizontal')+
      coord_flip()+
      scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                   title.hjust = .5))

    pwsprob = pwsprob[order(rownames(pwsprob)), order(colnames(pwsprob))]
    diag(pwsprob) = NA
    pwsprob[upper.tri(pwsprob)] = NA
    pwsprob = as.data.frame(pwsprob) %>%
      tibble::rownames_to_column(var = 'gen') %>%
      tidyr::pivot_longer(cols = c(colnames(pwsprob)[1]:colnames(pwsprob)[length(colnames(pwsprob))])) %>%
      dplyr::rename('gen2' = .data$name, 'prob' = .data$value)


    if(save.df){
      utils::write.csv(prob_g, file = paste0(getwd(),'/probsp.csv'), row.names = F)
      utils::write.csv(stats::na.exclude(pwsprob), file = paste0(getwd(),'/pwsprob.csv'), row.names = F)
    }

    if(interactive){

      requireNamespace('plotly')

      g_hpd = plotly::ggplotly(g_hpd)
      psps.bar = plotly::ggplotly(ggplot(prob_g, aes(x = .data$ID, y = .data$prob))+
                                    geom_bar(stat = 'identity', fill = '#781c1e')+
                                    labs(x = 'Genotypes', y = 'Probability') +
                                    theme(axis.text.x = element_text(angle = 90)))
      pwsp.heat = plotly::ggplotly(ggplot(pwsprob, aes(x = .data$gen, y = .data$gen2))+
                                     geom_tile(aes(fill = .data$prob), colour = 'white')+
                                     labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
                                     theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                                           legend.position = 'right', legend.direction = 'vertical')+
                                     coord_flip()+
                                     scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1)))

    }

    out = list(g_hpd, prob_g, psps.bar, na.exclude(pwsprob), pwsp.heat)
    names(out) = c('g_hpd', 'marg_prob.df', 'marg_prob.plot', 'pair_prob.df', 'pair_prob.plot')

    return(out)
  }

}

