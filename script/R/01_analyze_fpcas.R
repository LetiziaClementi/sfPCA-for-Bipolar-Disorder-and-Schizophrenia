rm(list=ls())


library(fdaPDE)
library(stringr)

# change according to your needs
setwd("~/OneDrive - Politecnico di Milano/PhD.LAVORO/3zo_anno/UCLA_WORK/script/R")
source('utils/fdaPDE.write.vtk.R')
source('utils/support_functions_2nd.R')

load('../../results/env_maps_afterfpca.Rdata')


## join data from fpcas and participants
participants <-  read.csv2(file.path(data_dir, 'participants.csv'),
                           header = T,
                           sep = '\t',
                           dec = '.',
                           na.strings = 'n/a')
part_ofinterest <- participants[which(participants$participant_id %in% rownames(my_maps$Amygdala_l_sum)), ]



create_boxplots(all_data, part_ofinterest, 5)
create_boxplots(all_data_schz_bip, part_ofinterest, 5)

create_boxplots(all_data_gcv, part_ofinterest, 5)



# merge data --------------------------------------------------------------

all_data <- join_src_data(fpcas, mv_pcas, part_ofinterest)
all_data_gcv <- join_src_data(fpca = fpcas_gcv, mv_pcas, part_ofinterest)
#remember that join_src_data merges all subjects so you have to specify in the call the important ones
all_data_schz_bip <- join_src_data(fpca = fpcas_schz_bip, 
                                   part_ofinterest = part_ofinterest[which(part_ofinterest$diagnosis=='SCHZ' | part_ofinterest$diagnosis=='BIPOLAR'), ],
                                   with_mv = F)

all_data_schz_bip_GCV <- join_src_data(fpca = fpcas_schz_bip_GCV, 
                                   part_ofinterest = part_ofinterest[which(part_ofinterest$diagnosis=='SCHZ' | part_ofinterest$diagnosis=='BIPOLAR'), ],
                                   with_mv = F)

all_data_schz_ctrl_GCV <- join_src_data(fpca = fpcas_schz_bip_GCV, 
                                       part_ofinterest = part_ofinterest[which(part_ofinterest$diagnosis=='SCHZ' | part_ofinterest$diagnosis=='BIPOLAR'), ],
                                       with_mv = F)

# save boxplots -----------------------------------------------------------

save_dir <- '../../plots/fpcas_performance/boxplots/'

boxplots_p(all_data, part_ofinterest, save_dir = save_dir, name_fpca = 'fpca_mvpca')
boxplots_p(all_data, part_ofinterest, save_dir = save_dir, second_group = 'SCHZ', name_fpca = 'fpca_mvpca_2groups')


boxplots_p(all_data_gcv, part_ofinterest, save_dir = save_dir, first_group = 'BIPOLAR', second_group = 'SCHZ', name_fpca = 'GCV_all')
boxplots_p(all_data_schz_bip, part_ofinterest, save_dir = save_dir, first_group = 'BIPOLAR', second_group = 'SCHZ', name_fpca = 'GCV_2', with_mv = F)



# paraview ----------------------------------------------------------------
plot_paraview_fpcas(fpcas,  paraview_dir = file.path(paraview_dir, 'loadings'))
plot_paraview_mvpcas(mv_pcas, FEMbasis, paraview_dir = file.path(paraview_dir, 'loadings'))
plot_paraview_fpcas(fpcas_gcv,  paraview_dir = file.path(paraview_dir, 'loadings', 'gcv'))


plot_paraview_fpcas(fpcas_schz_bip_GCV,  paraview_dir = file.path(paraview_dir, 'loadings', 'gcv'))
plot_paraview_fpcas(fpcas_schz_ctrl_GCV,  paraview_dir = file.path(paraview_dir, 'loadings', 'gcv'))


# plot fpca performance ---------------------------------------------------

fpcas_extracted <- list(fpcas, fpcas_gcv, fpcas_schz_bip, fpcas_schz_bip_GCV, fpcas_schz_ctrl_GCV)
names_fpcas <- c('fpcas', 'fpcas_gcv', 'fpcas_schz_bip', 'fpcas_schz_bip_GCV', 'fpcas_schz_ctrl_GCV')
names(fpcas_extracted) <- names_fpcas


print_plots_cumsum(fpcas_extracted, save_dir = '~/OneDrive - Politecnico di Milano/PhD.LAVORO/3zo_anno/UCLA_WORK/plots/fpcas_performance',
                   nPC = 5,
                   names_fpcas)





