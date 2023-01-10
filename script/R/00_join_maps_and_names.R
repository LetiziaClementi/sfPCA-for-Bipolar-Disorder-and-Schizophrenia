rm(list=ls())


library(fdaPDE)
library(stringr)

# change according to your needs
setwd("~/WORK/UCLA_WORK/script/R")
setwd("~/OneDrive - Politecnico di Milano/PhD.LAVORO/3zo_anno/UCLA_WORK/script/R")
source('utils/fdaPDE.write.vtk.R')
source('utils/support_functions_2nd.R')

maps_directory <- file.path('../../maps/')
data_dir <- file.path('../../data/')
paraview_dir <- file.path('../../plots/paraview/')
results_dir <- file.path('../../results/')


participants <-  read.csv2(file.path(data_dir, 'participants.csv'),
                           header = T,
                           sep = '\t',
                           dec = '.',
                           na.strings = 'n/a')
part_ofinterest <- participants[which(participants$participant_id %in% rownames(my_maps$Amygdala_l_sum)), ]
## this needs to be checked bc the null subjects are not dropped atm

## mesh
nodes <-read.csv('../../data/mesh/node_spaziodati_onlygm.csv', header = F)
elem <-read.csv('../../data/mesh/elem_onlygm.csv', header = F)
labels <- read.csv('../../data/mesh/labels_onlygm.csv', header = F)
mesh <- create.mesh.3D(nodes = nodes, tetrahedrons = elem[,1:4])
FEMbasis <- create.FEM.basis(mesh = mesh)

my_maps <- join_maps_name(maps_directory)

my_maps <- drop_empty_sbj(my_maps)
save(my_maps, file = file.path(data_directory, 'maps_list.RData'))

mean_maps <- list()
mean_maps <- compute_mean_maps(my_maps)
save(mean_maps, file = file.path(data_directory, 'mean_maps.RData'))

means_bygroup <- list()
means_bygroup <- means_by_group(my_maps, part_ofinterest)
#save(means_bygroup, file = file.path(data_directory, 'means_bygroup.RData'))

mean_zero <- list()
mean_zero <- subtract_mean(my_maps)


plot_bygroup_paraview(means_bygroup, FEMbasis, file.path(paraview_dir, 'mean'))

s <- seq(-3,3)
lambda <- 10^s

fpcas <- list()
fpcas <- compute_fpca(my_maps, FEMbasis, lambda = lambda)
save(fpcas, file = file.path(results_dir, 'fpcas.RData'))

fpcas_zero <- list()
fpcas_zero <- compute_fpca(mean_zero, FEMbasis, lambda = lambda)
save(fpcas, file = file.path(mean_zero, 'fpcas_zeromean.RData'))


mvpcas <- list()
mvpcas <- compute_mvpca(my_maps, FEMbasis, lambda = lambda)
save(mv_pcas, file = file.path(results_dir, 'mvpcas.RData'))


mvpcas_zero <- list()
mvpcas_zero <- compute_mvpca(mean_zero, FEMbasis, lambda = lambda)
save(mv_pcas_zero, file = file.path(results_dir, 'mvpcas_zeromean.RData'))



# 2 by 2 ------------------------------------------------------------------

maps_schz_bip <- select_mean(mean_zero, part_ofinterest)
maps_schz_ctrl <- select_mean(mean_zero, part_ofinterest, 'SCHZ', 'CONTROL')
maps_bip_ctrl <- select_mean(mean_zero, part_ofinterest, 'BIPOLAR', 'CONTROL')

fpcas_schz_bip <- list()
fpcas_schz_bip <- compute_fpca(maps_schz_bip[3:6], FEMbasis, lambda = lambda)
fpcas_schz_bip_GCV <- compute_fpca(maps_schz_bip[3:6], FEMbasis, lambda = lambda, validation = 'GCV')
save(fpcas_schz_bip, file = file.path(results_dir, 'fpcas_schz_bip.RData'))
save(fpcas_schz_bip_GCV, file = file.path('fpcas_schz_bip_GCV.RData'))

fpcas_schz_ctrl <- list()
fpcas_schz_ctrl <- compute_fpca(maps_schz_ctrl, FEMbasis, lambda = lambda)
fpcas_schz_ctrl_GCV <- compute_fpca(maps_schz_ctrl, FEMbasis, lambda = lambda, validation = 'GCV')
save(fpcas_schz_ctrl, file = file.path(results_dir, 'fpcas_schz_ctrl.RData'))
save(fpcas_schz_ctrl_GCV, file = file.path(results_dir, 'fpcas_schz_ctrl_GCV.RData'))


fpcas_bip_ctrl <- list()
fpcas_bip_ctrl <- compute_fpca(maps_bip_ctrl, FEMbasis, lambda = lambda)
fpcas_bip_ctrl_GCV <- compute_fpca(maps_bip_ctrl, FEMbasis, lambda = lambda, validation = 'GCV')
save(fpcas_bip_ctrl, file = file.path(results_dir, 'fpcas_bip_ctrl.RData'))
save(fpcas_bip_ctrl_GCV, file = file.path(results_dir, 'fpcas_bip_ctrl_GCV.RData'))



