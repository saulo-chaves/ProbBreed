prob_sup = function(data, trait, gen, env, reg = NULL, extr_outs, int = .2,
                    increase = TRUE, save.df = FALSE, interactive = FALSE){

  # Namespaces
  requireNamespace('ggplot2')
  requireNamespace('dplyr')

  # Preparation
  df = data
  data = if(any(is.na(data[,trait]))) data[-which(is.na(data[,trait])),] else data
  mod = extr_outs
  name.gen = levels(factor(data[,gen]))
  num.gen = nlevels(factor(data[,gen]))
  name.env = levels(factor(data[,env]))
  num.env = nlevels(factor(data[,env]))
  num.sim = nrow(mod$post$g)

  # Selection for increasing
  if(increase){
    if(!is.null(reg)){# If there are breeding regions

      # Preparations
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
        labs(x = 'Posterior effects (HPD)', y = 'Genotypes') +
        geom_point(size = 4, color = '#781c1e')

      # Marginal probabilities ----------------

      # Probability of superior performance --------------
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

      cat('1. Probability of superior performance estimated \n')

      # Pairwise probability of superior performance ----------------
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

      cat('2. Pairwise probability of superior performance estimated \n')

      # Probability of superior stability - Location -----------------
      staprob_gl = mod$post$gl
      colnames(staprob_gl) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gl),'_Env'))[,1]
        )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gl[,grep(x, colnames(staprob_gl))]
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

      cat('3. Probability of superior stability (GL) estimated \n')

      # Pairwise probability of superior stability - Location -------------
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

      cat('4. Pairwise probability of superior stability estimated (GL) \n')

      # Probability of superior stability - Region --------------
      staprob_gm = mod$post$gm
      colnames(staprob_gm) = sub(
        'Gen_','',do.call(rbind,strsplit(colnames(staprob_gm),'_Env'))[,1]
      )
      probsta = do.call(cbind, lapply(
        lapply(
          name.gen, function(x) staprob_gm[,grep(x, colnames(staprob_gm))]
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

      cat('5. Probability of superior stability (GM) estimated \n')


      # Pairwise probability of superior stability - Region -----------------
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

      cat('6. Pairwise probability of superior stability estimated (GL) \n')

      # Joint probability of superior performance and stability -----------------
      j_prob = rbind(merge(prob_g, prob_gl, by = 'ID'),
                     merge(prob_g, prob_gm, by = 'ID'))
      j_prob$joint = j_prob$prob.x * j_prob$prob.y
      colnames(j_prob) = c('ID', 'Performance', 'Stability', 'Joint')
      j_prob$lev = rep(c('Location', 'Region'), each = num.gen)
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
                   level = rep(c('Location','Region'), each = num.gen)
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

      cat('7. Joint probability of superior performance and stability estimated \n')


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


      if(save.df){
        dir.create(path = paste0(getwd(),'/marg_prob'))
        for (i in names(marg_prob$df)){
          write.csv(marg_prob$df[[i]],
                    file = paste0(getwd(),'/marg_prob/',i,'.csv'),
                    row.names = F)
        }
      }



    }else{ #If there is no breeding region

      # Preparation
      colnames(mod$post$g) = paste0(name.gen, '_')
      colnames(mod$post$gl) = paste(rep(name.gen,  times = num.env),
                                    rep(name.env,  each = num.gen), sep = '_')






    }
  }else{ # Selection for decreasing
    if(!is.null(reg)){# If there are breeding regions






    }else{ #If there is no breeding region






    }
  }

}
