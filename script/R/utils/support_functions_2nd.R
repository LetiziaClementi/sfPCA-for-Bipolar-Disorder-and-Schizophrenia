#funzione che preso in input cartella ritorna lista di mappe complete di ID sbj


join_maps_name <- function(maps_directory, ...){
  
  all_maps <-list()
  
  maps_data <- list.files(maps_directory)
  

  maps_files <- list()
  rois <- list()
  participants_files <- list()
  k <-1
  l <-1
  
  # individuo file soggetti e file mappe
  for(i in 1:length(maps_data)){
    if (substr(maps_data[i], 1, 6) == 'z_maps'){
      maps_files[k] <- maps_data[i]
      rois[k] <- substr((maps_data[i]), 8, nchar(maps_data[i])-4)
      k <- k+1
      
    } else if (substr(maps_data[i], 1, 4) == 'part'){
      participants_files[l] <- maps_data[i]
      #rois[l] <- substr((maps_data[i]), 8, nchar(maps_data[i])-4)
      l <- l+1
    }
  }
  


  
  # bisogna parallelizzare
  # apro mappe e corrispondenti soggetti e unisco
  for(i in 1:length(rois)){

    z_map <- read.csv2(file.path(maps_dir, paste0('z_maps_',rois[i],'.csv')),
                       header = F,
                       sep = ',',
                       dec = '.',
                       na.strings = 'NaN')
    #print(paste('opened z map for ', rois[i], '. Dimensions are: ', dim(z_map)))
    part_ofinterest <- read.csv2(file.path(maps_dir, paste0('part_ofinterest_', rois[i], '.csv')),
                                 header = T,
                                 sep = ',',
                                 dec = '.',
                                 na.strings = 'NaN')

    rownames(z_map) <- part_ofinterest$participant_id
    #print(paste('opened participants for ', rois[i], '. Dimensions are: ', dim(part_ofinterest)))

    all_maps[[i]] <- z_map

    print(paste('***', rois[i],' completed***'))
  }
  
  names(all_maps) <- rois
  return(all_maps)
}

drop_empty_sbj <- function(all_maps,...){
  for(i in 1:length(all_maps)){
    l <- rowMeans(all_maps[[i]])
    all_maps[[i]] <- all_maps[[i]][-which(l==0),]
  }
  
  return(all_maps)
}

compute_mean_maps <- function(all_maps,...){
  
  mean_maps <- list()
  #rois <- names(all_maps)
  
  for(i in 1:length(all_maps)){
    mean_maps[[names(all_maps[i])]] <- colMeans(all_maps[[i]], na.rm = T)
  }
  
  #names(mean_maps)<-names(all_maps)
  
  return(mean_maps)
}

plot_all_paraview <- function(mean_maps,FEMbasis, paraview_dir){
  require(fdaPDE)
  
  for (i in 1:length(mean_maps)) {
    write.vtu(FEM(mean_maps[[i]], FEMbasis), file = file.path(paraview_dir,paste0('mean_',names(mean_maps[i]),'.vtu')))
  }
}

means_by_group <- function(all_maps, participants, ...){
  groups <- unique(participants$diagnosis)
  means_bygroup <- list()
  
  for(i in 1:length(all_maps)){
    map <- my_maps[[i]]
    part_ofinterest <- read.csv2(file.path(maps_dir, paste0('part_ofinterest_', names(my_maps)[i], '.csv')),
                                 header = T,
                                 sep = ',',
                                 dec = '.',
                                 na.strings = 'NaN')
    part_ofinterest <- part_ofinterest[part_ofinterest$participant_id %in% rownames(map),]
    
    for (k in 1:length(groups)) {
      media <- colMeans(map[which(part_ofinterest$diagnosis == groups[k]), ], na.rm = T)
      means_bygroup[[names(my_maps)[i]]][[groups[k]]] <- media
    }
  }
  return(means_bygroup)
}

plot_bygroup_paraview <- function(means_bygroup, FEMbasis, paraview_dir, participants, task = 'rest'){
  require(fdaPDE)
  
  groups <- unique(participants$diagnosis)
  
  for (i in 1:length(mean_maps)) {
    for (k in 1:length(groups)) {
      write.vtu(FEM(means_bygroup[[i]][[groups[k]]], FEMbasis), 
                file = file.path(paraview_dir,
                                 paste0('mean_',names(mean_maps[i]), '_' ,groups[k],'_',task ,'.vtu')))
    }
    
  }
}

compute_fpca <- function(all_maps, FEMbasis, lambda = 1, nPC = 5, validation = 'KFold', GCVmethod = 'Exact'){
  require(fdaPDE)
  
  fpcas <- list()
  
  for (i in 1:length(all_maps)) {
    fpca <- FPCA.FEM(locations = NULL,
                              all_maps[[i]],
                              FEMbasis,
                              lambda = lambda,
                              nPC = nPC,
                              validation = validation)
    
    fpcas[[i]] <- fpca
    print(paste('***', names(all_maps[i]), ' done***'))
    
  }
  
  names(fpcas)<-names(all_maps)
  return(fpcas)
}

compute_mvpca <- function(all_maps, ncp = 5){
  require(FactoMineR)
  mvpcas <- list()
  
  for (i in 1:length(all_maps)) {
    pca <- PCA(all_maps[[i]], ncp = 5, graph = F)
    
    mvpcas[[i]] <- pca
    print(paste('***', names(all_maps[i]), ' done***'))
    
  }
  
  names(mvpcas)<-names(all_maps)
  return(mvpcas)
}

join_src_data <- function(fpcas, mvpca, part_ofinterest, with_mv = T, ...){
  
  nomi <- colnames(part_ofinterest[,1:4])
  all_data <- list()
  
  for (i in 1:dim(fpcas[[1]][['scores']])[2]) {
    nomi <- c(nomi, paste0('sfPC_',i))
  }
  
  if(with_mv){
  for (i in 1:dim(fpcas[[1]][['scores']])[2]) {
    nomi <- c(nomi, paste0('mvPC_',i))
  }}
  
  for(i in 1:length(fpcas)){
    
    fpca <- fpcas[[i]]
    all_data[[names(my_maps)[i]]] <- cbind(part_ofinterest[,1:4], fpca$scores)
    
    if(with_mv){
      mvpca <- mv_pcas[[i]]
      all_data[[names(my_maps)[i]]] <- cbind(all_data[[names(my_maps)[i]]], mvpca$svd$U)
    }
    
    
    
    names(all_data[[names(my_maps)[i]]]) <- nomi
  }

  return(all_data)
}

create_boxplots <- function(all_data, participants, nPC =5, criterion = 'diagnosis', with_mv = T){
  
  if (with_mv){
    m=2
  } else {
    m=1
  }
  
  for (i in 1:length(all_data)) {
    for (k in 1:(nPC*m)) {
      
      titolo = paste((names(all_data)[i]))
      ylab = paste(names(all_data[[i]][4+k]))     
               
      boxplot(all_data[[i]][,4+k]~ all_data[[i]][,criterion], 
              data = all_data[[i]],
              col = c('cornflowerblue', 'coral', 'palegreen'),
              main = titolo,
              ylab = ylab,
              xlab =  paste(criterion)
              )
      
      
      readline(prompt = 'press enter to show next plot')
      
    }
    
  }
  

  
  
}

subtract_mean <- function(all_maps){
  
  mean_zero <- list()
  
  for (i in 1:length(all_maps)) {
    map <- all_maps[[i]]
    # m <- colMeans(map, na.rm = T)
    # 
    # 
    # t <- rep(m, rep = nrow(map))
    # mean_zero[[i]] <- map 
    
    ones = rep(1, nrow(map))
    m = ones %*% t(colMeans(map, na.rm = T))
    mean_zero[[i]] <- map -m
    
  }
  
  names(mean_zero) <- names(all_maps)
  return(mean_zero)
  
}

select_mean <- function(all_maps, part_ofinterest, group_one = 'SCHZ', group_two = 'BIPOLAR'){
  
  mean_selected <- list()
  sbj_ofinterest <- part_ofinterest[which(part_ofinterest$diagnosis == group_one | part_ofinterest$diagnosis == group_two), 'participant_id']
  
  
  for (i in 1:length(all_maps)) {
    map <- all_maps[[i]]
    
    
    map <- map[sbj_ofinterest, ]
    
    mean_selected[[i]] <- map
    
  }
  
  names(mean_selected) <- names(all_maps)
  return(mean_selected)
  
}

boxplots_p <- function(all_data, participants, nPC=5, first_group = 'BIPOLAR', second_group = 'CONTROL', criterion = 'diagnosis', with_mv = T, save_dir, name_fpca, to_save = T){
  
  
  if (with_mv){
    m=2
  } else {
    m=1
  }
  
  
  for (i in 1:length(all_data)) {
    data <- all_data[[i]]
    data <- data[which(data$diagnosis==first_group | data$diagnosis==second_group), ]
    titolo = paste((names(all_data)[i]))
    
    for (k in 1:(nPC*m)) {
      
      a <- data[which(data$diagnosis == first_group), (4+k)]
      b <- data[which(data$diagnosis == second_group), (4+k)]
      p <- wilcox.test(a,b)
      
      
      measure <- names(data)[4+k]
      
      # ylab = paste(names(data[,4+k]))
      # xlab = str(criterion)
      # 
      # c <- data[,4+k]
      # d <- data[which(data$diagnosis==first_group | data$diagnosis==second_group), 'diagnosis']
      
      if(p$p.value <0.51){
        boxplot(data[,4+k]~data[,'diagnosis'],
                data = data,
                col = c('cornflowerblue', 'coral'),
                main = paste(titolo, 'p:', round(p$p.value, digits = 3), '***'),
                ylab = paste(measure),
                xlab = paste(criterion)
        )
      } else {
        boxplot(data[,4+k]~data[,'diagnosis'],
                data = data,
                col = c('cornflowerblue', 'coral'),
                main = paste(titolo, 'p:', round(p$p.value, digits = 3)),
                ylab = paste(measure),
                xlab = paste(criterion))}
      if(to_save){
        pdf(file = paste0(save_dir, name_fpca,'_', measure, '_boxplot.pdf'))
        dev.off()
      }

      readline(prompt = 'press enter to show next plot')
    }
  }
}

plot_paraview_fpcas <- function(fpcas, paraview_dir, task = 'rest'){
  
  for(i in 1:length(fpcas)){
    
    fpca <- fpcas[[i]]
    write.vtu(FEM(fpca$loadings.FEM$coeff, 
                  fpca$loadings.FEM$FEMbasis), 
              file = file.path(paraview_dir,
                               paste0('loadings_', names(fpcas[i]), '_' ,task ,'.vtu')))
    
  }

  
}


plot_paraview_mvpcas <- function(mv_pcas, FEMbasis, paraview_dir, task = 'rest'){
  
  for(i in 1:length(mv_pcas)){
    
    mvpca <- mv_pcas[[i]]
    write.vtu(FEM(mvpca$svd$V, 
                  FEMbasis), 
              file = file.path(paraview_dir,
                               paste0('mv_loadings_', names(mv_pcas[i]), '_' ,task ,'.vtu')))
    
  }
  
  
}


print_plots_cumsum <- function(fpcas_extracted, save_dir, nPC = 5, names_fpcas){
  cl <- rainbow(10)
  number_fpcas=length(fpcas_extracted)
  rois <- names(fpcas_extracted[[1]])
  
  for (i in 1:length(rois)) {
    roi <- rois[i]
    dat <- matrix(data = NA, nrow = number_fpcas, ncol = nPC)
    
    for (k in 1:number_fpcas) {
      try(dat[k, ] <- fpcas_extracted[[k]][[roi]][['cumsum_percentage']])
    }
    pdf(file = paste0(save_dir, '/', roi, '_cumsum_percentage.pdf'))
    
    matplot(t(dat), 
            type = c("b"),
            pch=1,
            col = 1:length(fpcas_extracted), 
            ylab = 'cumsum percentage',
            main = paste(rois[i]))
    legend("bottomright", legend = names_fpcas, col=1:length(fpcas_extracted), pch=1) 
    dev.off()
  }

  
}




