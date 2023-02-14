## Function cond_prob
##
##' This function computates the  probabilities of superior performance within
##' environments and yields plots for better visualization.
##' All the plots are customizable, using ggplot2
##'
##' @title Conditional probabilities of superior performance
##' @param data A dataframe containing the observations
##' @param trait A character representing the name of the column that
##' corresponds to the analysed variable
##' @param gen A character representing the name of the column that
##' corresponds to the evaluated genotypes
##' @param env A character representing the name of the column that
##' corresponds to the environments
##' @param reg A character representing the name of the column that
##' corresponds to the regions. NULL otherwise (default)
##' @param extr_outs An object from the 'extr_outs' function
##' @param int numeric A number representing the selection intensity
##' (superior limit = 1)
##' @param save.df A logical value indicating if the data frames with the marginal
##' probability of each genotype and the pairwise probabilities should be saved at
##' the work directory.
##' @param interactive A logical value indicating if the plots should be interactive.
##' If TRUE, the function loads the "plotly" package and uses the "ggplotly" command.
##' @return The function returns a list with:
##' \itemize{
##' \item \code{risk.plot} a plot representing the Bayesian Eskridge's Risk
##' index (Adapted from Eskridge (1990))
##' \item \code{psp_env_reg.mat} : a matrix containing the probability of superior
##' performance within environments and regions. Exists only if a string is provided
##' for "reg"
##' \item \code{psp_env.mat} : a matrix containing the probability of superior
##' performance within environments. Exists only "reg" is NULL
##' \item \code{psp_env.plot} : a heatmap representing the probability of superior
##' performance of each genotype, within each environment
##' \item \code{psp_reg.plot} : a heatmap representing the probability of superior
##' performance of each genotype, within each region. Exists only if a string is
##' provided for "reg"
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
##'
##' @examples
##' \dontrun{conds = marg_prob(data = df, trait = "GY", gen = "H",env = "L", extr_outs = outs,
##'                   reg = "Region", int = .2, save.df = TRUE, interactive = TRUE)}


cond_prob = function(data, trait, gen, env, reg = NULL, extr_outs, int = .2,
                     save.df = FALSE, interactive = FALSE){

  requireNamespace('ggplot2')
  requireNamespace('dplyr')
  mod = extr_outs
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))
  num.sim = nrow(mod$post$g)

  if(!is.null(reg)){

    # Eskridge

    name.reg = levels(factor(data[,reg]))
    num.reg = nlevels(factor(data[,reg]))

    V1 = apply(matrix(mod$map$gl, num.gen, num.env,
                      dimnames = list(name.gen, name.env)), 1, sd)
    V2 = apply(matrix(mod$map$gm, num.gen, num.reg,
                      dimnames = list(name.gen, name.reg)), 1, sd)
    Zi = quantile(mod$post$g, probs = .95)
    Risk = mod$map$g - (Zi * V1+V2)
    Risk = as.data.frame(Risk) %>% rownames_to_column(var = 'gen') %>%
      arrange(desc(Risk))
    Risk$gen = factor(Risk$gen, levels = Risk$gen)

    sfi = ggplot(Risk, aes(x = gen, y = Risk))+
      geom_segment(aes(x = gen, xend = gen, y = 0, yend = Risk), linewidth = 1)+
      geom_point(color = '#781c1e', size = 2)+
      labs(x = 'Genotypes', y = 'Safety-first Index')+
      theme(axis.text.x = element_text(angle = 90))

    # Probabiliies of superior performance within environments

    colnames(mod$post$g) = paste0(name.gen, '_')
    colnames(mod$post$gm) = paste(rep(name.gen,  times = num.reg),
                                  rep(name.reg,  each = num.gen), sep = '_')
    name.env.reg = sort(paste(unique(data[,c(env,reg)])[,1],
                              unique(data[,c(env,reg)])[,2], sep = '_'))
    colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                  rep(name.env.reg,  each = num.gen), sep = '_')


    posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl
    for (i in name.reg) {
      posgge[,grep(i, colnames(posgge))] = posgge[,grep(i, colnames(posgge))] +
        matrix(mod$post$gm[,grep(i, colnames(mod$post$gm))],
               nrow = num.sim, ncol = num.gen * length(name.env.reg[grep(i, name.env.reg)]))
    }

   supprob = function(vector, num.gen, int){
     ifelse(names(vector) %in%
            names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
   }

   probs = lapply(
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
            x,MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
            )}
          )

   probs = apply(do.call(rbind, probs), 2, function(x){
                 tapply(x, rep(name.gen, num.sim), mean)
                 })

   env.heat = as.data.frame(probs) %>% rownames_to_column(var = 'gen') %>%
     pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
     separate(name, into = c('envir','region'), sep = '_') %>%
     ggplot(aes(x = envir, y = reorder(gen, value), fill = value))+
     geom_tile(colour = 'white')+
     labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
     theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
           legend.position = 'right', legend.direction = 'vertical')+
     scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))

   reg.heat = as.data.frame(probs) %>% rownames_to_column(var = 'gen') %>%
     pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
     separate(name, into = c('envir','region'), sep = '_') %>%
     ggplot(aes(x = region, y = reorder(gen, value), fill = value))+
     geom_tile(colour = 'white')+
     labs(x = 'Regions', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
     theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
           legend.position = 'right', legend.direction = 'vertical')+
     scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))

   if(save.df){
     write.csv(probs, file = paste0(getwd(),'/psp_env_reg.mat.csv'), row.names = F)
   }

   if(interactive){
     env.heat = suppressWarnings(plotly::ggplotly(env.heat))
     reg.heat = suppressWarnings(plotly::ggplotly(reg.heat))
     sfi = suppressWarnings(plotly::ggplotly(sfi))
   }

   outs = list(sfi, probs, env.heat, reg.heat)
   names(outs) = c("risk.plot", 'psp_env_reg.mat', 'psp_env.plot','psp_reg.plot')

   return(outs)


  }else{

    # Eskridge

    Vi = apply(matrix(mod$map$gl, num.gen, num.env,
                      dimnames = list(name.gen, name.env)), 1, sd)
    Zi = quantile(mod$post$g, probs = .95)
    Risk = mod$map$g - (Zi * Vi)
    Risk = as.data.frame(Risk) %>% rownames_to_column(var = 'gen') %>%
      arrange(desc(Risk))
    Risk$gen = factor(Risk$gen, levels = Risk$gen)

    sfi = ggplot(Risk, aes(x = gen, y = Risk))+
      geom_segment(aes(x = gen, xend = gen, y = 0, yend = Risk), linewidth = 1)+
      geom_point(color = '#781c1e', size = 2)+
      labs(x = 'Genotypes', y = 'Safety-first Index')+
      theme(axis.text.x = element_text(angle = 90))


    # Probability of superior performance

    colnames(mod$post$g) = paste0(name.gen, '_')
    colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                  rep(name.env,  each = num.gen), sep = '_')


    posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl

    supprob = function(vector, num.gen, int){
      ifelse(names(vector) %in%
               names(vector[order(vector, decreasing = T)][1:ceiling(int * num.gen)]), 1, 0)
    }

    probs = lapply(
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
    )

    probs = apply(do.call(rbind, probs), 2, function(x){
      tapply(x, rep(name.gen, num.sim), mean)
    })

    env.heat = as.data.frame(probs) %>% rownames_to_column(var = 'gen') %>%
      pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))]),
                  names_to = 'envir') %>%
      ggplot(aes(x = envir, y = reorder(gen, value), fill = value))+
      geom_tile(colour = 'white')+
      labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
      theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
            legend.position = 'right', legend.direction = 'vertical')+
      scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))


    if(save.df){
      write.csv(probs, file = paste0(getwd(),'/psp_env.mat.csv'), row.names = F)
    }

    if(interactive){
      requireNamespace('plotly')
      env.heat = suppressWarnings(plotly::ggplotly(env.heat))
      sfi = suppressWarnings(plotly::ggplotly(sfi))
    }

    outs = list(sfi, probs, env.heat)
    names(outs) = c("risk.plot", 'psp_env.mat', 'psp_env.plot')

    return(outs)
  }

}

