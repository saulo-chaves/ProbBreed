## Function cond_prob
##
##' @title
##' Conditional probabilities of superior performance
##'
##' @description
##' This function estimates the  probabilities of superior performance within
##' environments and yields plots for better visualization.
##' All the plots are customizable, using `ggplot2`
##'
##' @details
##' Using `cond_prob()`, you can estimate the probability of a genotype being amongst the
##' selected, based on a given selection intensity in a specific environment or
##' breeding region. If we let \eqn{\Omega_k} represent the subset of superior
##' genotypes in the <i>kth</i> environment, then the probability of the <i>jth</i>
##' genotype belonging to \eqn{\Omega_k} is:
##'
##' \deqn{Pr(g_{jk} \in \Omega_k \vert y) = \frac{1}{S}\sum_{s=1}^S{I(g_{jk}^{(s)} \in \Omega \vert y)}}
##'
##' where \eqn{S} is the total number of samples, and \eqn{I(g_{jk}^{(s)} \in \Omega \vert y)}
##' is an indicator variable mapping success (1) if \eqn{g_{jk}^{(s)}} exists in \eqn{\Omega},
##' and failure (0) otherwise.
##'
##' @param data A data frame containing the observations
##' @param trait A character representing the name of the column that
##' corresponds to the analysed variable
##' @param gen A character representing the name of the column that
##' corresponds to the evaluated genotypes
##' @param env A character representing the name of the column that
##' corresponds to the environments
##' @param reg A character representing the name of the column that
##' corresponds to the regions. `NULL` otherwise (default)
##' @param extr_outs An object from the `extr_outs` function
##' @param int A number representing the selection intensity
##' (superior limit = 1)
##' @param increase Logical: `TRUE` if genotypes with higher trait values are better.
##' FALSE otherwise.
##' @param save.df A logical value indicating if the data frames with the marginal
##' probability of each genotype and the pairwise probabilities should be saved at
##' the work directory.
##' @param interactive A logical value indicating if the plots should be interactive.
##' If `TRUE`, the function loads the `plotly` package and uses the [plotly::ggplotly()]
##' command.
##' @return The function returns a list with:
##' \itemize{
##' \item \code{conds_prob} : a matrix containing the probability of superior
##' performance within environments (and regions, if available).
##' \item \code{conds_pwprob} : a list with the pairwise probabilities of superior
##' performance within environments.
##' \item \code{conds_pwprob.reg} : a list with the pairwise probabilities of superior
##' performance within regions (if available).
##' \item \code{psp_env.plot} : a heatmap representing the probability of superior
##' performance of each genotype, within each environment
##' \item \code{psp_reg.plot} : a heatmap representing the probability of superior
##' performance of each genotype, within each region. Exists only if a string is
##' \item \code{pwprobs.plots} : a list containing a heatmap representing the pairwise
##'  probabilities of superior performance for each environment.
##'  \item \code{pwprobs.reg.plots} : a list containing a heatmap representing the pairwise
##'  probabilities of superior performance for each region (if available).
##' }
##'
##' @seealso [ggplot2::ggplot()]
##'
##' @references
##'
##' Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães, L. J. M.,
##' Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging probability concepts
##' for cultivar recommendation in multi-environment trials. <i>Theoretical and
##' Applied Genetics</i>, 133(2):443-455. https://doi.org/10.1007/s00122-022-04041-y
##'
##' @import ggplot2 dplyr
##' @importFrom utils write.csv combn
##' @importFrom tidyr pivot_longer
##' @importFrom tidyr separate
##' @importFrom tibble rownames_to_column
##' @importFrom rlang .data
##'
##' @export
##'
##' @examples
##' \dontrun{
##' mod = bayes_met(data = soy,
##'                 gen = "Gen",
##'                 env = "Env",
##'                 rept = NULL,
##'                 reg = "Reg",
##'                 res.het = F,
##'                 trait = "Y",
##'                 iter = 100, cores = 4, chain = 4)
##'                 #Remember, increase the number of iterations, cores and chains
##'
##' outs = extr_outs(data = soy, trait = "Y", gen = "Gen", model = mod,
##'                  effects = c('l','g','gl','m','gm'),
##'                  nenv = length(unique(soy$Env)), res.het = FALSE,
##'                  probs = c(0.05, 0.95), check.stan.diag = TRUE)
##'
##' conds = cond_prob(data = soy, trait = "Y", gen = "Gen", env = "Env",
##'                   extr_outs = outs, reg = "Reg", int = .2,
##'                   increase = TRUE, save.df = TRUE, interactive = TRUE)
##'                   }


cond_prob = function(data, trait, gen, env, reg = NULL, extr_outs, int = .2,
                     increase = TRUE, save.df = FALSE, interactive = FALSE){

  requireNamespace('ggplot2')
  requireNamespace('dplyr')
  df = data
  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = extr_outs
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))
  num.sim = nrow(mod$post$g)

  if(increase){
    if(!is.null(reg)){

       name.reg = levels(factor(data[,reg]))
       num.reg = nlevels(factor(data[,reg]))

      # Probabilities of superior performance within environments

      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                    'Reg',rep(name.reg,  each = num.gen), sep = '_')
      name.env.reg = sort(paste('Env',unique(data[,c(env,reg)])[,1],
                                'Reg',unique(data[,c(env,reg)])[,2], sep = '_'))
      colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.env),
                                    rep(name.env.reg,  each = num.gen), sep = '_')


      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl
      for (i in name.reg) {
        posgge[,grep(i, do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] =
          posgge[,grep(i, do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] +
          matrix(mod$post$gm[,grep(i, do.call(rbind,strsplit(colnames(mod$post$gm),'Reg'))[,2])],
                 nrow = num.sim, ncol = num.gen * length(name.env.reg[grep(i, do.call(rbind,strsplit(name.env.reg,'Reg'))[,2])]))
      }

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
            x,MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
          )}
      )), 2, function(x){
        tapply(x, rep(name.gen, num.sim), mean)
      })
      probs = probs * ifelse(table(data[,gen],data[,env]) != 0, 1, NA)


      combs = data.frame(t(utils::combn(paste('Gen', name.gen, sep = '_'), 2)))
      colnames(combs) = c('x', 'y')
      pwprobs = lapply(
        sapply(paste('Env', name.env, sep = '_'),
               function(x) posgge[,grep(x, colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(z[1], colnames(y))] > y[,grep(z[2], colnames(y))])
            })
          )

          a[,1] = sub('Gen_', '', a[,1])
          a[,2] = sub('Gen_', '', a[,2])

          a
        }
      )
      names(pwprobs) = sub('Env_', '', names(pwprobs))
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



      pwprobs.reg = lapply(
        sapply(paste('Reg', name.reg, sep = '_'),
               function(x) posgge[,grep(x, colnames(posgge))],
               simplify = F),
        function(y){

          a = cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(z[1], colnames(y))] > y[,grep(z[2], colnames(y))])
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



      env.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
        tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
        tidyr::separate(.data$name, into = c('envir','region'), sep = '_Reg_') %>%
        dplyr::mutate(
          envir = sub('Env_',"",.data$envir)
        ) %>%
        ggplot(aes(x = .data$envir, y = reorder(.data$gen, .data$value), fill = .data$value))+
        geom_tile(colour = 'white')+
        labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[jk] %in% Omega[k]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = 'right', legend.direction = 'vertical')+
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))

      reg.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
        tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
        tidyr::separate(.data$name, into = c('envir','region'), sep = '_Reg_') %>%
        dplyr::group_by(.data$gen,.data$region) %>%
        dplyr::summarise(value = mean(.data$value, na.rm=T), .groups = 'drop') %>%
        ggplot(aes(x = .data$region, y = reorder(.data$gen, .data$value), fill = .data$value))+
        geom_tile(colour = 'white')+
        labs(x = 'Regions', y = 'Genotypes', fill = expression(bold(Pr(g[jm] %in% Omega[m]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = 'right', legend.direction = 'vertical')+
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))

      pwprobs.plots = lapply(pwprobs, function(x){
        x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
        x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })


      if(save.df){
        utils::write.csv(probs, file = paste0(getwd(),'/conds_prob.csv'), row.names = F)
        dir.create(path = paste0(getwd(),'/cond_pwprob_env'))
        for (i in names(pwprobs)){
          write.csv(pwprobs[[i]], file = paste0(getwd(),'/cond_pwprob_env/cond_pwprob_env_',i,'.csv'),
                    row.names = F)
        }
        dir.create(path = paste0(getwd(),'/cond_pwprob_reg'))
        for (i in names(pwprobs.reg)){
          write.csv(pwprobs.reg[[i]], file = paste0(getwd(),'/cond_pwprob_reg/cond_pwprob_reg_',i,'.csv'),
                    row.names = F)
        }
      }

      if(interactive){
        env.heat = suppressWarnings(plotly::ggplotly(env.heat))
        reg.heat = suppressWarnings(plotly::ggplotly(reg.heat))
      }

      outs = list(probs, pwprobs, pwprobs.reg, env.heat, reg.heat,
                  pwprobs.plots, pwprobs.reg.plots)
      names(outs) = c('conds_prob', 'conds_pwprobs.env', 'conds_pwprobs.reg',
                      'psp_env.plot','psp_reg.plot', 'pwprobs.plots',
                      'pwprobs.reg.plots')
      return(outs)


    }else{

      # Probability of superior performance

      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                    rep(name.env,  each = num.gen), sep = '_')


      posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl

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


      combs = data.frame(t(utils::combn(name.gen, 2)))
      colnames(combs) = c('x', 'y')

      pwprobs = lapply(
        sapply(name.env,
               function(x) posgge[,grep(x, colnames(posgge))],
               simplify = F),
        function(y){
          cbind(
            combs,
            pwprob = apply(combs, 1, function(z){
              mean(y[,grep(z[1], colnames(y))] > y[,grep(z[2], colnames(y))])
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


      env.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
        tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))]),
                            names_to = 'envir') %>%
        ggplot(aes(x = .data$envir, y = reorder(.data$gen, .data$value), fill = .data$value))+
        geom_tile(colour = 'white')+
        labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[jk] %in% Omega[k]))))+
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = 'right', legend.direction = 'vertical')+
        scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))

      pwprobs.plots = lapply(pwprobs, function(x){
        x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
          geom_tile() +
          labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] > g[j]))))+
          scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                       title.hjust = .5)) +
          theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
                legend.position = c(.8,.15), legend.direction = 'horizontal')
      })

      if(save.df){
        utils::write.csv(probs, file = paste0(getwd(),'/conds_prob.csv'), row.names = F)
        dir.create(path = paste0(getwd(),'/test'))
        for (i in names(pwprobs)){
          write.csv(pwprobs[[i]], file = paste0(getwd(),'/test/cond_pwprob_',i,'.csv'),
                    row.names = F)
        }
      }

      if(interactive){
        requireNamespace('plotly')
        env.heat = suppressWarnings(plotly::ggplotly(env.heat))
      }

      outs = list(probs, pwprobs, env.heat, pwprobs.plots)
      names(outs) = c('conds_prob', 'conds_pwprobs', 'psp_env.plot', 'pwprobs.plots')

      return(outs)
    }
  }else{

  if(!is.null(reg)){

     name.reg = levels(factor(data[,reg]))
     num.reg = nlevels(factor(data[,reg]))

    # Probabilities of superior performance within environments

    colnames(mod$post$g) = paste0(name.gen, '_')
    colnames(mod$post$gm) = paste('Gen',rep(name.gen,  times = num.reg),
                                  'Reg',rep(name.reg,  each = num.gen), sep = '_')
    name.env.reg = sort(paste('Env',unique(data[,c(env,reg)])[,1],
                              'Reg',unique(data[,c(env,reg)])[,2], sep = '_'))
    colnames(mod$post$gl) = paste('Gen',rep(name.gen,  times = num.env),
                                  'Env',rep(name.env.reg,  each = num.gen), sep = '_')


    posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl
    for (i in name.reg) {
      posgge[,grep(i, do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] =
        posgge[,grep(i, do.call(rbind,strsplit(colnames(posgge),'Reg'))[,2])] +
        matrix(mod$post$gm[,grep(i, do.call(rbind,strsplit(colnames(mod$post$gm),'Reg'))[,2])],
               nrow = num.sim, ncol = num.gen * length(name.env.reg[grep(i, do.call(rbind,strsplit(name.env.reg,'Reg'))[,2])]))
    }

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
         x,MARGIN = 2, FUN = supprob, num.gen = num.gen, int = .2
       )}
   )), 2, function(x){
                 tapply(x, rep(name.gen, num.sim), mean)
                 })
   probs = probs * ifelse(table(data[,gen],data[,env]) != 0, 1, NA)


   combs = data.frame(t(utils::combn(paste('Gen', name.gen, sep = '_'), 2)))
   colnames(combs) = c('x', 'y')
   pwprobs = lapply(
     sapply(paste('Env', name.env, sep = '_'),
            function(x) posgge[,grep(x, colnames(posgge))],
            simplify = F),
     function(y){

       a = cbind(
         combs,
         pwprob = apply(combs, 1, function(z){
           mean(y[,grep(z[1], colnames(y))] < y[,grep(z[2], colnames(y))])
         })
       )

       a[,1] = sub('Gen_', '', a[,1])
       a[,2] = sub('Gen_', '', a[,2])

       a
     }
   )
   names(pwprobs) = sub('Env_', '', names(pwprobs))
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

   pwprobs.reg = lapply(
     sapply(paste('Reg', name.reg, sep = '_'),
            function(x) posgge[,grep(x, colnames(posgge))],
            simplify = F),
     function(y){

       a = cbind(
         combs,
         pwprob = apply(combs, 1, function(z){
           mean(y[,grep(z[1], colnames(y))] < y[,grep(z[2], colnames(y))])
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

   env.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
     tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
     tidyr::separate(.data$name, into = c('envir','region'), sep = '_Reg_') %>%
     dplyr::mutate(
       envir = sub('Env_',"",.data$envir)
     ) %>%
     ggplot(aes(x = .data$envir, y = reorder(.data$gen, .data$value), fill = .data$value))+
     geom_tile(colour = 'white')+
     labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[jk] %in% Omega[k]))))+
     theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
           legend.position = 'right', legend.direction = 'vertical')+
     scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))

   reg.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
     tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))])) %>%
     tidyr::separate(.data$name, into = c('envir','region'), sep = '_Reg_') %>%
     dplyr::group_by(.data$gen,.data$region) %>%
     dplyr::summarise(value = mean(.data$value, na.rm=T), .groups = 'drop') %>%
     ggplot(aes(x = .data$region, y = reorder(.data$gen, .data$value), fill = .data$value))+
     geom_tile(colour = 'white')+
     labs(x = 'Regions', y = 'Genotypes', fill = expression(bold(Pr(g[jm] %in% Omega[m]))))+
     theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
           legend.position = 'right', legend.direction = 'vertical')+
     scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))


   pwprobs.plots = lapply(pwprobs, function(x){
     x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
       geom_tile() +
       labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] < g[j]))))+
       scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
       guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                    title.hjust = .5)) +
       theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
             legend.position = c(.8,.15), legend.direction = 'horizontal')
   })

   pwprobs.reg.plots = lapply(pwprobs.reg, function(x){
     x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
       geom_tile() +
       labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] < g[j]))))+
       scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
       guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                    title.hjust = .5)) +
       theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
             legend.position = c(.8,.15), legend.direction = 'horizontal')
   })


   if(save.df){
     utils::write.csv(probs, file = paste0(getwd(),'/conds_prob.csv'), row.names = F)
     dir.create(path = paste0(getwd(),'/test'))
     for (i in names(pwprobs)){
       utils::write.csv(pwprobs[[i]], file = paste0(getwd(),'/test/cond_pwprob_',i,'.csv'),
                 row.names = F)
     }
     dir.create(path = paste0(getwd(),'/cond_pwprob_reg'))
     for (i in names(pwprobs.reg)){
       write.csv(pwprobs.reg[[i]], file = paste0(getwd(),'/cond_pwprob_reg/cond_pwprob_reg_',i,'.csv'),
                 row.names = F)
     }
   }

   if(interactive){
     env.heat = suppressWarnings(plotly::ggplotly(env.heat))
     reg.heat = suppressWarnings(plotly::ggplotly(reg.heat))
   }

   outs = list(probs, pwprobs, pwprobs.reg, env.heat, reg.heat,
               pwprobs.plots, pwprobs.reg.plots)
   names(outs) = c('conds_prob', 'conds_pwprobs.env', 'conds_pwprobs.reg',
                   'psp_env.plot','psp_reg.plot', 'pwprobs.plots',
                   'pwprobs.reg.plots')
   return(outs)


  }else{

    # Probability of superior performance

    colnames(mod$post$g) = paste0(name.gen, '_')
    colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                  rep(name.env,  each = num.gen), sep = '_')


    posgge = matrix(mod$post$g, nrow = num.sim, ncol = num.env * num.gen) + mod$post$gl

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
    probs = probs * ifelse(table(data[,gen],data[,env]) != 0, 1, NA)


    combs = data.frame(t(utils::combn(name.gen, 2)))
    colnames(combs) = c('x', 'y')
    pwprobs = lapply(
      sapply(name.env,
             function(x) posgge[,grep(x, colnames(posgge))],
             simplify = F),
      function(y){
        cbind(
          combs,
          pwprob = apply(combs, 1, function(z){
            mean(y[,grep(z[1], colnames(y))] < y[,grep(z[2], colnames(y))])
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


    env.heat = as.data.frame(probs) %>% tibble::rownames_to_column(var = 'gen') %>%
      tidyr::pivot_longer(cols = c(colnames(probs)[1]:colnames(probs)[length(colnames(probs))]),
                  names_to = 'envir') %>%
      ggplot(aes(x = .data$envir, y = reorder(.data$gen, .data$value), fill = .data$value))+
      geom_tile(colour = 'white')+
      labs(x = 'Environments', y = 'Genotypes', fill = expression(bold(Pr(g[jk] %in% Omega[k]))))+
      theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
            legend.position = 'right', legend.direction = 'vertical')+
      scale_fill_viridis_c(direction = -1, na.value = '#D3D7DC',limits = c(0,1))

    pwprobs.plots = lapply(pwprobs, function(x){
      x |> ggplot(aes(x = .data$x, y = .data$y, fill = .data$pwprob)) +
        geom_tile() +
        labs(x = 'Genotypes', y = 'Genotypes', fill = expression(bold(Pr(g[i] < g[j]))))+
        scale_fill_viridis_c(direction = -1, na.value = 'white',limits = c(0,1))+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5, title.position = 'top',
                                     title.hjust = .5)) +
        theme(axis.text.x = element_text(angle = 90),panel.background = element_blank(),
              legend.position = c(.8,.15), legend.direction = 'horizontal')
    })

    if(save.df){
      utils::write.csv(probs, file = paste0(getwd(),'/conds_prob.csv'), row.names = F)
      dir.create(path = paste0(getwd(),'/test'))
      for (i in names(pwprobs)){
        utils::write.csv(pwprobs[[i]], file = paste0(getwd(),'/test/cond_pwprob_',i,'.csv'),
                  row.names = F)
      }
    }

    if(interactive){
      requireNamespace('plotly')
      env.heat = suppressWarnings(plotly::ggplotly(env.heat))
    }

    outs = list(probs, pwprobs, env.heat, pwprobs.plots)
    names(outs) = c('conds_prob', 'conds_pwprobs', 'psp_env.plot', 'pwprobs.plots')

    return(outs)
  }
  }

}

