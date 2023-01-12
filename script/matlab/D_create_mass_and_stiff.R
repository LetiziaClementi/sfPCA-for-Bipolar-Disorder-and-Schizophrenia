rm(list=ls())

setwd("~/OneDrive - Politecnico di Milano/PhD.LAVORO/2ndo_anno/PRJ_PRINCIPALE/03_MARZO/script/")

library(fdaPDE)

nodes <-read.csv('../data/mesh/node_spaziodati_onlygm.csv', header = F)
elem <-read.csv('../data/mesh/elem_onlygm.csv', header = F)
mesh <- create.mesh.3D(nodes = nodes, tetrahedrons = elem[,1:4])
FEMbasis <- create.FEM.basis(mesh = mesh)
plot(mesh)

mass <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
stiff <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEMbasis)

writeMM(mass, file = '../data/mass.txt')
writeMM(stiff, file = '../data/stiff.txt')
# NB: ricordare di rielaborare a mano