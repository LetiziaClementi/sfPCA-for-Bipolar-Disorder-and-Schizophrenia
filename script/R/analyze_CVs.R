rm(list = ls())
library(stringr)
library(ggplot2)
library(RColorBrewer)

setwd("~/WORK/sfPCA-for-Bipolar-Disorder-and-Schizophrenia/script/R")

path_cvs <- file.path('..', '..', 'results')
path_results <- file.path('..', '..', 'results')
path_maps <- file.path('..', '..', 'maps')
path_participants <- file.path('..', '..', 'participants')
path_data <- file.path('..', '..', 'data')


participants <- read.csv2(file = file.path(path_data,'participants.csv'),
                          header = T,
                          sep = '\t',
                          dec = '.')

rois <- list.files(path_maps)
rois <- str_sub(rois, 8, -5)
rois


lambdas <- c('lambda_m_3', 'lambda_m_2', 'lambda_m_1', 'lambda_zero', 
             'lambda_1', 'lambda_2', 'lambda_3')

sfpca_names <- vector()
for(i in 1:5){
  sfpca_names[i] <- paste0('sfPC_', i)
}

# import CV data
my_cvs <- list()
for(i in 1:length(rois)){
  cv <- read.csv2(file = file.path(path_cvs, paste0(rois[i], '_CV.csv')),
                  header = F,
                  sep = ',',
                  dec = '.',
                  row.names = sfpca_names,
                  col.names = lambdas)
  my_cvs[[i]] <- (cv)
  rm(cv)
}
names(my_cvs) <- rois

# import explained variance data
my_variances <- list()
for(i in 1:length(rois)){
  cv <- read.csv2(file = file.path(path_cvs, paste0(rois[i], '_exp_var.csv')),
                  header = F,
                  sep = ',',
                  dec = '.',
                  row.names = sfpca_names,
                  #col.names = lambdas
                  )
  my_variances[[i]] <- t(cv)*100
  rm(cv)
}
names(my_variances) <- rois


# plot CVs with matplot
par(mfrow=c(1,length(my_cvs)))
for(i in 1:length(my_cvs)){
  matplot(t(my_cvs[[i]]), 
          type = 'l',
          ylim = c(0.20, 0.50) ,
          ylab = ' ',
          main = rois[i])
}


# plot CVs with classic plot (vertical)
par(mfrow=c(1, length(my_cvs)))
for (j in 1:length(my_cvs)) {
  plot(as.double(my_cvs[[j]][1,]), 
       type = 'l', 
       ylim = c(0.20, 0.50),
       ylab = 'CV',
       xlab = 'lambda',
       xaxt = 'n', 
       main = names(my_cvs)[j],
       col = brewer.pal(5, 'Set2')[1],
       lwd = 3)
  axis(1, at = c(1:7),
       labels = c(-3:3)) 
  
  if (j==1){
    legend('bottom', 
           legend=sfpca_names, 
           lty=1, 
           col = brewer.pal(dim(my_cvs[[1]])[1], 'Set2'), 
           cex=1,
           lwd = 3)
  }

  
  for (i in 2:dim(my_cvs[[1]])[1]) {
    lines(as.double(my_cvs[[j]][i,]),  
          ylim = c(0.20, 0.50),
          col = brewer.pal(dim(my_cvs[[1]])[1], 'Set2')[i],
          cex = 1,
          pch = 10,
          lwd = 3
          )
  }
  
}




# plot explained variance with classic plot (vertical)
par(mfrow=c(1, 1))
plot(as.double(my_variances[[1]]), 
     type = 'b', 
     ylim = c(10, 30),
     ylab = 'explained variance [%]',
     xlab = 'n sfPC',
     lty = 2,
     main = "Explained Variance",
     col = brewer.pal(length(my_variances), 'Paired')[1],
     lwd = 3)
legend('topleft', 
       legend=names(my_variances), 
       lty=1, 
       col = brewer.pal(length(my_variances), 'Paired'), 
       cex=1,
       lwd = 3)
for (j in 2:length(my_variances)) {  
  lines(as.double(my_variances[[j]]),  
        ylim = c(10, 30),
        col = brewer.pal(length(my_variances), 'Paired')[j],
        cex = 1,
        pch = 10,
        lty = 2,
        type = 'b',
        lwd = 3
  )}




# boxplot of scores -------------------------------------------------------
# import U_norm data (aka SCORES)
my_U_norm <- list()
for(i in 1:length(rois)){
  U <- read.csv2(file = file.path(path_results, paste0(rois[i], '_U_normalized.csv')),
                  header = F,
                  sep = ',',
                  dec = '.',
                  #row.names = sfpca_names,
                  col.names = sfpca_names
                 )
  part <- read.csv2(file = file.path(path_participants, 
                                     paste0('part_ofinterest_',
                                            rois[i], 
                                            '.csv')),
                    header = T,
                    sep = ',',
                    dec = '.')
  
  part_lacking <- read.csv2(file = file.path(path_participants, 
                                     paste0('part_lacking_',
                                            rois[i], 
                                            '.csv')),
                    header = T,
                    sep = ',',
                    dec = '.')
  
  part <- part[-which(part$participant_id %in% part_lacking$participant_id),]
  rownames(U) <- part$participant_id
  
  my_U_norm[[i]] <- cbind(part[1:4], U)
  rm(U, part, part_lacking)
}
names(my_U_norm) <- rois



# BOXPLOT ALL
for(i in 1:length(my_U_norm)){
  data <- my_U_norm[[i]]
  
  pdf(file=file.path('..', '..', 'plots', 'U_norm', 'all', paste0('U_norm_', rois[i], '_ALL.pdf')),
      width = 15, height = 8)
  par(mfrow=c(1,5))
  
  for(j in 1:5){
    measure <- j+4
    boxplot(data[which(data$diagnosis == 'CONTROL'), measure], 
            data[which(data$diagnosis == 'BIPOLAR'), measure], 
            data[which(data$diagnosis == 'SCHZ'), measure],
            col = brewer.pal(3, 'Set2'), 
            names = c('CTRL', 'BIPOLAR', 'SCHZ'),
            xlab = paste('sfPC', j)
    )
    
  }
  mtext(paste(rois[i]), side = 3, line = -3, outer = TRUE)
  
  dev.off() 
  
  
  #readline()
}

#relevant_U_norm <- list()
relevant_U_norm_DF <- list()
par(mfrow=c(1, 5))
# BIP versus CTRL
for(i in 1:length(my_U_norm)){
  data <- my_U_norm[[i]]

  pdf(file=file.path('..', '..', 'plots', 'U_norm', 'BIP_CTRL', paste0('U_norm_', rois[i], '_BIP_CTRL.pdf')),
      width = 15, height = 8)
  par(mfrow=c(1,5))
  
  for(j in 1:5){
    measure <- j+4
    
    p <- wilcox.test(data[which(data$diagnosis == 'CONTROL'), measure], 
                     data[which(data$diagnosis == 'BIPOLAR'), measure])
    if(p$p.value < 0.05){
      #relevant_U_norm <- append(relevant_U_norm, paste0('BIP_CTRL_', rois[i], '_sfPC_', j))
      relevant_U_norm_DF <- rbind(relevant_U_norm_DF, paste0('BIP_CTRL', ';', rois[i],';', 'sfPC_', j))
    }
    
    boxplot(data[which(data$diagnosis == 'CONTROL'), measure], 
            data[which(data$diagnosis == 'BIPOLAR'), measure], 
            col = brewer.pal(3, 'Set2'), 
            names = c('CTRL', 'BIPOLAR'),
            xlab = paste('sfPC', j)
    )
    title(paste('p = ', round(p$p.value, digits = 3)), line = -2)
    
  }
  mtext(paste(rois[i]), side = 3, line = -2, outer = TRUE)
  
  dev.off() 
  
  
  #readline()
}

par(mfrow=c(1, 5))
# CTRL_SCHZ
for(i in 1:length(my_U_norm)){
  data <- my_U_norm[[i]]
  
  pdf(file=file.path('..', '..', 'plots', 'U_norm', 'CTRL_SCHZ', paste0('U_norm_', rois[i], '_CTRL_SCHZ.pdf')),
      width = 15, height = 8)
  par(mfrow=c(1,5))
  
  for(j in 1:5){
    measure <- j+4
    
    p <- wilcox.test(data[which(data$diagnosis == 'CONTROL'), measure], 
                     data[which(data$diagnosis == 'SCHZ'), measure]) 
    
    if(p$p.value < 0.05){
      #relevant_U_norm <- append(relevant_U_norm, paste0('CTRL_SCHZ_', rois[i], '_sfPC_', j))
      relevant_U_norm_DF <- rbind(relevant_U_norm_DF, paste0('CTRL_SCHZ', ';', rois[i], ';','sfPC_', j))
    }
    
    boxplot(data[which(data$diagnosis == 'CONTROL'), measure], 
            data[which(data$diagnosis == 'SCHZ'), measure], 
            col = c(brewer.pal(3, 'Set2')[1], brewer.pal(3, 'Set2')[3]), 
            names = c('CTRL', 'SCHZ'),
            xlab = paste('sfPC', j)
    )
    title(paste('p = ', round(p$p.value, digits = 3)), line = -2)
    
  }
  mtext(paste(rois[i]), side = 3, line = -2, outer = TRUE)
  
  dev.off() 
  
  
  #readline()
}


par(mfrow=c(1, 5))
# BIP_SCHZ
for(i in 1:length(my_U_norm)){
  data <- my_U_norm[[i]]
  
  pdf(file=file.path('..', '..', 'plots', 'U_norm', 'BIP_SCHZ', paste0('U_norm_', rois[i], '_BIP_SCHZ.pdf')),
      width = 15, height = 8)
  par(mfrow=c(1,5))
  
  for(j in 1:5){
    measure <- j+4
    
    p <- wilcox.test(data[which(data$diagnosis == 'BIPOLAR'), measure], 
                     data[which(data$diagnosis == 'SCHZ'), measure]) 
    if(p$p.value < 0.05){
      #relevant_U_norm <- append(relevant_U_norm, paste0('BIP_SCHZ_', rois[i], '_sfPC_', j))
      relevant_U_norm_DF <- rbind(relevant_U_norm_DF, paste0('BIP_SCHZ', ';',rois[i],';', 'sfPC_', j))
    }
    
    boxplot(data[which(data$diagnosis == 'BIPOLAR'), measure], 
            data[which(data$diagnosis == 'SCHZ'), measure], 
            col = c(brewer.pal(3, 'Set2')[2], brewer.pal(3, 'Set2')[3]), 
            names = c('CTRL', 'SCHZ'),
            xlab = paste('sfPC', j)
    )
    title(paste('p = ', round(p$p.value, digits = 3)), line = -2)
    
  }
  mtext(paste(rois[i]), side = 3, line = -2, outer = TRUE)
  
  dev.off() 
  
  
  #readline()
}


write.csv2(relevant_U_norm_DF, 
           file = file.path(path_results, 'relevant_U_norm_DF.csv'),
           row.names = F, 
           col.names = F)
rel <- read.csv2(file = file.path(path_results, 'relevant_U_norm_DF.csv'),
                 sep = ';',
                 header = F)
rel <- rel[-1,]
colnames(rel) <- c('COMPARISON', 'ROI', 'sfPC')
barplot(rel)


ggplot(rel, aes(x = COMPARISON)) +
  geom_bar()

ggplot(rel, aes(x = ROI)) +
  geom_bar()


ggplot(rel, aes(x = ROI, fill = sfPC)) +
  geom_bar()

ggplot(rel, aes(x = COMPARISON, fill = ROI)) +
  geom_bar()



