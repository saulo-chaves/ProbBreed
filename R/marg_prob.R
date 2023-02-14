## Function marg_prob
##
##' This function computates the marginal probabilities of superior performance and
##' yields plots for better visualization. All the plots are customizable, using ggplot2
##'
##' @title Marginal probabilities of superior performance
##' @param data A dataframe containing the observations
##' @param trait A character representing the name of the column that
##' corresponds to the analysed variable
##' @param gen A character representing the name of the column that
##' corresponds to the evaluated genotypes
##' @param env A character representing the name of the column that
##' corresponds to the environments
##' @param extr_outs An object from the 'extr_outs' function
##' @param int An integer representing the selection intensity
##' (superior limit = 1)
##' @param save.df A logical value indicating if the data frames with the marginal
##' probability of each genotype and the pairwise probabilities should be saved at
##' the work directory.
##' @param interactive A logical value indicating if the plots should be interactive.
##' If TRUE, the function loads the "plotly" package and uses the "ggplotly" command.
##' @return The function returns a list with:
##' \itemize{
##' \item \code{g_hpd} : a catterpillar plot representing the performance of
##' each genotype, and their respective confidence interval (95% and 97.5%)
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
##' @export
##'
##' @import ggplot2 dplyr
##' @importFrom utils write.csv
##' @importFrom tidyr pivot_longer
##' @importFrom tidyr separate
##' @importFrom tibble rownames_to_column
##'
##' @examples
##' \dontrun{
##' margs = marg_prob(data = df, trait = "GY", gen = "H",env = "L", extr_outs = outs,
##'                   int = .2, save.df = TRUE, interactive = TRUE)
##'                   }
##'

marg_prob = function(data, trait, gen, env, extr_outs, int = .2,
                     save.df = FALSE, interactive = FALSE){

  requireNamespace('ggplot2')
  requireNamespace('tidyr')
  requireNamespace('dplyr')

  mod = extr_outs
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))

  g_hpd = as.data.frame(matrix(mod$post$g,
                               dimnames=list(t(outer(colnames(mod$post$g),
                                                     rownames(mod$post$g),
                                                     FUN=paste)), NULL))) %>%
    mutate(gen = rep(name.gen, each = nrow(mod$post$g))) %>%  group_by(gen) %>%
    summarise(
      g = median(V1),
      UP = quantile(V1, probs = 0.95),
      up = quantile(V1, probs = 0.975),
      DOWN = quantile(V1, probs = 0.05),
      down = quantile(V1, probs = 0.025)
    ) %>%
    ggplot(aes(x = g, y = reorder(gen,g))) +
    geom_errorbar(aes(xmin = down, xmax = up), width = 0)+
    geom_errorbar(aes(xmin = DOWN, xmax = UP), width = 0, linewidth = 2, alpha = .8) +
    labs(x = 'Values', y = 'Genotypes') +
    geom_point(size = 4, color = '#781c1e')


  colnames(mod$post$g) = name.gen
  ind_post = apply(mod$post$g, 1, function(x){
    ifelse(name.gen %in%
             names(x[order(x, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
  })
  rownames(ind_post) = name.gen

  prob_g = data.frame(ID = rownames(ind_post), prob = rowMeans(ind_post), row.names = NULL)
  prob_g = prob_g[order(prob_g$prob, decreasing = T),]
  prob_g$ID = factor(prob_g$ID, levels = prob_g$ID)

  psps.bar = ggplot(prob_g, aes(x = ID, y = prob))+
    geom_bar(stat = 'identity', fill = '#781c1e')+
    geom_text(aes(label = round(prob,2)), fontface = 'bold', vjust = -.4,
              angle = 45, hjust = -.3)+
    labs(x = 'Genotypes', y = 'Probability') +
    theme(axis.text.x = element_text(angle = 90))


  pwsprob = matrix(NA, num.gen, num.gen)

  for(i in 1:ncol(mod$post$g)){
    for (j in 1:ncol(mod$post$g)) {
      pwsprob[i,j] = mean(mod$post$g[,j] > mod$post$g[,i])
    }
  }
  colnames(pwsprob) = rownames(pwsprob) = name.gen
  pwsprob = pwsprob[match(prob_g$ID, rownames(pwsprob)),
                    match(prob_g$ID, rownames(pwsprob))]
  diag(pwsprob) = NA
  pwsprob[upper.tri(pwsprob)] = NA
  pwsprob = as.data.frame(pwsprob) %>%
    tibble::rownames_to_column(var = 'gen') %>%
    tidyr::pivot_longer(cols = c(colnames(pwsprob)[1]:colnames(pwsprob)[length(colnames(pwsprob))])) %>%
    rename('gen2' = name, 'prob' = value)

  pwsprob$gen = factor(pwsprob$gen, levels = unique(pwsprob$gen))
  pwsprob$gen2 = factor(pwsprob$gen2, levels = unique(pwsprob$gen))

  pwsp.heat = ggplot(pwsprob, aes(x = gen, y = gen2))+
    geom_tile(aes(fill = prob), colour = 'white')+
    labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
    theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
          legend.position = c(.8,.15), legend.direction = 'horizontal')+
    coord_flip()+
    scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                 title.hjust = .5))

  if(save.df){
    write.csv(prob_g, file = paste0(getwd(),'/probsp.csv'), row.names = F)
    write.csv(pwsprob, file = paste0(getwd(),'/pwsprob.csv'), row.names = F)
  }

  if(interactive){

    requireNamespace('plotly')

    g_hpd = plotly::ggplotly(g_hpd)
    psps.bar = plotly::ggplotly(ggplot(prob_g, aes(x = ID, y = prob))+
                          geom_bar(stat = 'identity', fill = '#781c1e')+
                          labs(x = 'Genotypes', y = 'Probability') +
                          theme(axis.text.x = element_text(angle = 90)))
    pwsp.heat = plotly::ggplotly(ggplot(pwsprob, aes(x = gen, y = gen2))+
                           geom_tile(aes(fill = prob), colour = 'white')+
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

