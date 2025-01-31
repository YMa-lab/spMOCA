#' The function to create a spatial dependency matrix given spatial coordiante
#' @param spatialCoord a location-by-2 spatial coordinate matrix
#' @param normCount a gene-by-location normalized matrix count 
#' @param kernelType the type of Spatial Kernel to be constructed. The default is Gaussian
#' @param kernelBw for constructing a spatial kernel, which type of bandwidth estimation will be used. The default is Rule-Of-Thumb
#' @return Return a spatial kernel matrix

createSpatialKernel = function(spatialCoord,
                               normCount,
                               kernelType = "Gaussian",
                               kernelBw = "Rule-Of-Thumb"){
  
  colnames(spatialCoord) = c("x","y")
  spatialCoord = scale(spatialCoord)
  L_dist = as.matrix(dist(spatialCoord,upper = T)^2)
  
  if(kernelType == "Gaussian"){
    
    if(kernelBw == "Rule-Of-Thumb"){
      normCount_ScaleGene = normCount
      for(i in 1:nrow(normCount_ScaleGene)){
        rowexp = normCount_ScaleGene[i,]
        rowexp = as.numeric(scale(rowexp))
        normCount_ScaleGene[i,] = rowexp
      }
      geneBw = apply(normCount_ScaleGene,1,function(x){tryCatch({bw.nrd0(x)},error = function(e){NA})})
      estBw = median(geneBw)
    }
    
    V = exp(-L_dist/estBw) # Gaussian kernel matrix
    
  }
  
  rownames(V) = colnames(V) = rownames(spatialCoord)
  return(V)
}

#' The function to estimate the gene correlation by maximizing log-likelihood given a spatial dependency
#' @param normCount a gene-by-location normalized count matrix
#' @param spatialKernel a spatial kernel
#' @return Return a maximized likelihood gene coexpression correlation estimation.

spMOCAest = function(normCount,
                     spatialKernel){
  
  V_R = as.matrix(nearPD(spatialKernel)$mat) 
  V_Rinv = solve(V_R) 
  X = normCount
  G = nrow(X)
  N = ncol(X)
  Vec_1G = matrix(1,nrow = G)
  Vec_1N = matrix(1,nrow = N)
  ### Mu
  mu_hat_num = matmult(matmult(X,V_Rinv),Vec_1N)
  mu_hat_dom = matmult(matmult(t(Vec_1N),V_Rinv),Vec_1N)
  mu_hat = mu_hat_num/as.numeric(mu_hat_dom)
  
  ### X_Center
  X_Center = X - matmult(mu_hat,t(Vec_1N))
  ### Sigma G
  XV_R = matmult(X_Center,V_Rinv)
  SigmaG = (1/N)*matmult(XV_R,t(X_Center))
  
  ### Corr RG
  W_invsq = matrix(0,nrow = G,ncol = G)
  diag(W_invsq) = 1/sqrt(diag(SigmaG))
  WXV_R = matmult(matmult(W_invsq,X_Center),V_Rinv)
  XW = matmult(t(X_Center),W_invsq)
  RG = (1/N)*matmult(WXV_R,XW)
  
  return(RG)
}