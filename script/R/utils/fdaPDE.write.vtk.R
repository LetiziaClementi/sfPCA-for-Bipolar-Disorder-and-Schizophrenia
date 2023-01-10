
#### A function to export in .vtu format ####
write.vtu.old <- function(FEMobj, file = ""){
  
  if(class(FEMobj) != "FEM"){
    stop("FEMobj is not of class FEM.")
  }
  
  if(!is.character(file)){
    stop("file is not of type character.")
  }
  
  ## Header
  cat('<?xml version="1.0"?>\n', file = file)
  cat('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n', file = file, append = TRUE)
  cat('<UnstructuredGrid>\n', file = file, append = TRUE)
  
  
  p   = matrix(FEMobj$FEMbasis$mesh$nodes, nrow = 3, ncol = FEMobj$FEMbasis$mesh$nnodes)
  dim = 3
  t   = matrix(FEMobj$FEMbasis$mesh$tetrahedrons, nrow = 4, ncol = FEMobj$FEMbasis$mesh$ntetrahedrons)
  
  t   = t - 1
  
  nnodes = FEMobj$FEMbasis$mesh$nnodes
  nelems = FEMobj$FEMbasis$mesh$ntetrahedrons
  
  ## Header for <Piece>
  cat(paste0('<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',nelems,'">\n'), file = file, append = TRUE)

  ## Print grid
  print_grid(file, dim, p, nnodes, t, nelems);
  
  ## Print Data
  print_data_points(file, FEMobj$coeff, nnodes)
  
  ## Footer for <Piece>
  cat('</Piece>\n', file = file, append = TRUE)
  
  ## Footer
  cat('</UnstructuredGrid>\n', file = file, append = TRUE)
  cat('</VTKFile>', file = file, append = TRUE)
}



print_grid <- function(file, dim, p, nnodes, t, nelems){
  eltype = 10
  
  ## VTK-Points (mesh nodes)
  cat('<Points>\n', file = file, append = TRUE)
  cat('<DataArray type="Float64" Name="Array" NumberOfComponents="3" format="ascii">\n', file = file, append = TRUE)
  write.table(t(p),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('</Points>\n', file = file, append = TRUE)
  
  ## VTK-Cells (mesh elements)
  cat('<Cells>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="connectivity" format="ascii">\n', file = file, append = TRUE)
  write.table(format(t(t), scientific = FALSE),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="offsets" format="ascii">\n', file = file, append = TRUE)
  write.table(format(matrix(seq(dim+1,(dim+1)*nelems,by=(dim+1)),nrow = 1, byrow = TRUE),  scientific = FALSE),  
              file = file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="types" format="ascii">\n', file = file, append = TRUE)
  write.table(format(matrix(rep(eltype,nelems),nrow = 1, byrow = TRUE), scientific = FALSE),  file = file, append = TRUE, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('</Cells>\n', file = file, append = TRUE)
}



# ## Print DataPoints
# print_data_points <- function(file, nodedata, nnodes){
# 
#   ## # of data to print in 
#   ## <PointData> field
#   nvdata = length(nodedata)  
#   
#   cat('<PointData>\n', file = file, append = TRUE)
#   cat('<DataArray type="Float64" Name="Data" NumberOfComponents="1" format="ascii">\n', file = file, append = TRUE)
#   write.table(matrix(nodedata,nrow = 1),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
#   cat('</DataArray>\n', file = file, append = TRUE)
#   cat('</PointData>\n', file = file, append = TRUE)
# 
# }

## Print DataPoints
print_data_points <- function(file, nodedata, nnodes){
  
  ## # of data to print in 
  ## <PointData> field
  nvdata = ncol(nodedata)  
  
  cat('<PointData>\n', file = file, append = TRUE)
  for (ind in 1:nvdata) {
    cat(paste0('<DataArray type="Float64" Name="Data',ind,'" NumberOfComponents="1" format="ascii">\n'), file = file, append = TRUE)
    write.table(matrix(nodedata[,ind],nrow = 1),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
    cat('</DataArray>\n', file = file, append = TRUE)
  }
  cat('</PointData>\n', file = file, append = TRUE)
  
}



############################################################
write.vtu <- function(FEMobj, file = ""){
  
  if(class(FEMobj) != "FEM"){
    stop("FEMobj is not of class FEM.")
  }
  
  if(!is.character(file)){
    stop("file is not of type character.")
  }
  
  ## Header
  cat('<?xml version="1.0"?>\n', file = file)
  cat('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n', file = file, append = TRUE)
  cat('<UnstructuredGrid>\n', file = file, append = TRUE)
  
  
  p   = t(FEMobj$FEMbasis$mesh$nodes)
  dim = 3
  t   = t(FEMobj$FEMbasis$mesh$tetrahedrons)
  
  t   = t - 1
  
  nnodes = nrow(FEMobj$FEMbasis$mesh$nodes)
  nelems = nrow(FEMobj$FEMbasis$mesh$tetrahedrons)
  
  ## Header for <Piece>
  cat(paste0('<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',nelems,'">\n'), file = file, append = TRUE)
  
  ## Print grid
  print_grid(file, dim, p, nnodes, t, nelems);
  
  ## Print Data
  print_data_points(file, FEMobj$coeff, nnodes)
  
  ## Footer for <Piece>
  cat('</Piece>\n', file = file, append = TRUE)
  
  ## Footer
  cat('</UnstructuredGrid>\n', file = file, append = TRUE)
  cat('</VTKFile>', file = file, append = TRUE)
}