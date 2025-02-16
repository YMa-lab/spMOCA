setClass("spMOCA",
         slots = list(
           spatial_countMat = "ANY",
           spatial_normCount = "ANY",
           feature_select = "vector",
           spatial_location = "data.frame",
           spatial_kernel = "matrix",
           corr_est = "matrix",
           network = "data.frame",
           project = "character"
         ))

#' The function to create a spMOCA object for inferring gene co-expression network for SRT
#' @param spatial_countMat a SRT gene-by-locatiion raw count matrix
#' @param spatial_location a location-by-2 spatial coordinate data fra,e
#' @param feature_list a list of selected features
#' @return Return a spMOCA object with a raw count data, a location data and initial selected features

createSPMOCAObject = function(spatial_countMat,
                              spatial_location,
                              feature_list){
  object = new(Class = "spMOCA",
               spatial_countMat = spatial_countMat,
               spatial_location = spatial_location,
               feature_select = feature_list,
               project = "coexpression")
  cat(paste0("## spMOCA object created! ...\n"))
  return(object)
}

#' The function to normalize the count matrix in a spMOCA object
#' @param object a spMOCA object
#' @param norm.method Normalization method. The default is a log-transformation with library size normalization
#' @return Return a spMOCA object with a normalized expression matrix

normalizeSPMOCA = function(object,
                                   norm.method = "logNorm"){
  if(norm.method == "logNorm"){
    cat(paste0("## Performing log-mormalization to count data ...\n"))
    count_data = object@spatial_countMat
    normCount = sweep(count_data,2,colSums(count_data),"/")
    normCount = log1p(normCount * 10000)
    normCount = as.matrix(normCount)
    object@spatial_normCount = normCount
  }
  return(object)
}

#' The function to perform a feature selection
#' @param object a spMOCA object
#' @param feature.selection Feature selection method. The default is SPARK for spatial variable genes selection. 
#' If feature.selection == Null, no selection will be performed. 
#' @param feature.list A list of pre-specified selected features. 
#' Users can pre-specify their features list and updated the spMOCA object directly.
#' @return Return a spMOCA object with updated selected features
#' 
featureSelectionSPMOCA = function(object,
                           feature.selection = "SPARK",
                           feature.list = NULL){
  if(!is.null(feature.selection)){
    if(feature.selection == "SPARK"){
      cat(paste0("## Performing feature selection ...\n"))
      location = object@spatial_location
      count_data = object@spatial_countMat
      
      location_coord = as.data.frame(location[,c("x","y")])
      location_coord[,1] = as.numeric(location_coord[,1])
      location_coord[,2] = as.numeric(location_coord[,2])
      
      spark <- CreateSPARKObject(counts=count_data,
                                 location=location_coord,
                                 percentage = 0.1,
                                 min_total_counts = 5)
      spark@lib_size <- apply(spark@counts, 2, sum)
      spark <- spark.vc(spark,
                        covariates = NULL,
                        lib_size = spark@lib_size,
                        num_core = parallel::detectCores(),
                        verbose = F)
      spark <- spark.test(spark,
                          check_positive = T,
                          verbose = F)
      
      sign.res = spark@res_mtest %>%
        dplyr::filter(adjusted_pvalue <= 0.05)
      gene_spark = rownames(sign.res)
      object@feature_select = gene_spark
      cat(paste0("## Feature selection done! ...\n"))
    }
  }
  
  
  if(is.null(feature.selection)){
    if(!is.null(feature.list)){
      cat(paste0("## Updated Selected Features! ...\n"))
      object@feature_select = feature.list
    }
  }
  return(object)
}

#' The function to create a spatial dependency matrix given spatial coordiante
#' @param object a spMOCA object
#' @param kernelType the type of Spatial Kernel to be constructed. The default is Gaussian
#' @param kernelBw for constructing a spatial kernel, which type of bandwidth estimation will be used. The default is Rule-Of-Thumb
#' @return Return a spMOCA object a spatial kernel matrix
#' 
createSpatialKernel = function(object,
                               kernelType = "Gaussian",
                               kernelBw = "Rule-Of-Thumb"){
  spatialCoord = object@spatial_location
  normCount = object@spatial_normCount
  feature = object@feature_select
  normCount = normCount[feature,]
  
  colnames(spatialCoord) = c("x","y")
  spatialCoord = scale(spatialCoord)
  L_dist = as.matrix(dist(spatialCoord,upper = T)^2)
  
  if(kernelType == "Gaussian"){
    cat(paste0("## Constructing Gaussian Spatial Kernel ...\n"))
    if(kernelBw == "Rule-Of-Thumb"){
      normCount_ScaleGene = normCount
      for(i in 1:nrow(normCount_ScaleGene)){
        rowexp = normCount_ScaleGene[i,]
        rowexp = as.numeric(scale(rowexp))
        normCount_ScaleGene[i,] = rowexp
      }
      geneBw = apply(normCount_ScaleGene,1,function(x){tryCatch({bw.nrd0(x)},error = function(e){NA})})
      estBw = median(geneBw,na.rm = T)
    }
    
    V = exp(-L_dist/estBw) # Gaussian kernel matrix
    rownames(V) = colnames(V) = rownames(spatialCoord)
    object@spatial_kernel = V
  }
  return(object)
}

#' The function to estimate the gene correlation by maximizing log-likelihood given a spatial dependency
#' @param object a spMOCA object
#' @return Return a spMOCA object with maximized likelihood estimated gene-gene correlation.

estimatespMOCA = function(object){
  spatialKernel = object@spatial_kernel
  normCount = object@spatial_normCount
  feature = object@feature_select
  normCount = normCount[feature,]
  
  cat(paste0("## Calculating Inverse of Spatial Kernel ...\n"))
  V_R = as.matrix(nearPD(spatialKernel)$mat) 
  V_Rinv = solve(V_R) 
  X = normCount
  G = nrow(X)
  N = ncol(X)
  Vec_1G = matrix(1,nrow = G)
  Vec_1N = matrix(1,nrow = N)
  ### Mu
  cat(paste0("## Estimating Gene Correlation ...\n"))
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
  
  rownames(RG) = colnames(RG) = feature
  object@corr_est = RG
  return(object)
}

#' The function to perform gene module detection
#' @param object a spMOCA object
#' @param module.method module detection method. The default is WGCNA
#' @param res_parameter resolution parameter for module detection method
#' @return Return a maximized likelihood gene coexpression correlation estimation.

moduleDetectionspMOCA = function(object,
                                 module.method = "WGCNA",
                                 res_parameter = 400){
  REst = object@corr_est
  if(module.method == "WGCNA"){
    cat(paste0("## Performing WGCNA ...\n"))
    # Rounding correlation matrix to make it runable for WGCNA
    corr.mat = round(REst,6)
    corr.mat[lower.tri(corr.mat)] = t(corr.mat)[lower.tri(corr.mat)]
    
    # Compute the adjacency matrix based on the co-expression matrix
    adj.mat = WGCNA::adjacency.fromSimilarity(corr.mat, power = 1,type = "signed")
    # Compute the topological overlap matrix
    TOM = WGCNA::TOMsimilarity(adj.mat)
    dissTOM = 1-TOM
    rownames(dissTOM) <- colnames(dissTOM) <- rownames(corr.mat)
    hclust_dist = hclust(as.dist(dissTOM)) 
    wcgna_label = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist,
                                                distM = dissTOM,
                                                deepSplit = T,
                                                pamRespectsDendro = FALSE,
                                                minClusterSize = res_parameter)
    
    ##### Calculate Gene's Global and Within-Module Connectivity
    gene.module.all = intramodularConnectivity( (corr.mat + 1)/2,
                                                colors = wcgna_label)
    gene.module.all$k = wcgna_label
    gene.module.all$gene = object@feature_select
    
    object@network = gene.module.all
  }
  return(object)
}
