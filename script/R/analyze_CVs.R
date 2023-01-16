rm(list = ls())
library(stringr)

setwd("~/WORK/sfPCA-for-Bipolar-Disorder-and-Schizophrenia/script/R")

path_cvs <- file.path('..', '..', 'results')
path_maps <- file.path('..', '..', 'maps')


rois <- list.files(path_maps)
rois <- str_sub(rois, 8, -5)

lambdas <- c('lambda_m_3', 'lambda_m_2', 'lambda_m_1', 'lambda_zero', 
             'lambda_1', 'lambda_2', 'lambda_3')

sfpca_names <- vector()
for(i in 1:5){
  sfpca_names[i] <- paste0('sfPC_', i)
}


my_cvs <- list()
for(i in 1:length(rois)){
  cv <- read.csv2(file = file.path(path_cvs, paste0(rois[i], '_CV.csv')),
                  header = F,
                  sep = ',',
                  dec = '.',
                  row.names = sfpca_names,
                  col.names = lambdas)
  my_cvs[[i]] <- cv
}
names(my_cvs) <- rois




matplot(t(my_cvs$Amygdala_l_sum))








par(mfrow=c(1,8))
for(i in 1:length(my_cvs)){
  matplot(t(my_cvs[[i]]), 
          type = 'l',
          ylim = c(0.30, 0.50) ,
          main = rois[i])
}





