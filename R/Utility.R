#' The function to calculate module score given gene modules
#' @param object a spMOCA object
#' @param gene.use the features used to calculate module scores. The default is hub genes.
#' @return Return a location-by-module module score matrix.

calModuleScore = function(object,
                                 gene.use = "hubgene"){
  cat(paste0("## Calculating Module Score ...\n"))
  gene.module.all = object@network
  normCount = object@spatial_normCount
  if(gene.use == "hubgene"){
    ##### Calculate Module Score 
    hubgene.module.all = gene.module.all %>%
      dplyr::group_by(k) %>%
      dplyr::top_frac(0.05,kWithin)
    num_cluster_wcgna = length(unique(gene.module.all$k)[which(unique(gene.module.all$k) != 0)])
    
    ms_prop = NULL
    for(kk in 1:num_cluster_wcgna){
      module_net = hubgene.module.all[hubgene.module.all$k == kk,]
      genesModule = module_net$gene
      sub_count = normCount[genesModule,]
      ms_prop = cbind(ms_prop,colMeans(sub_count))
    }
    colnames(ms_prop) = paste0("Mod",1:num_cluster_wcgna)
  }
  return(ms_prop)
}



#' The function to perform gene set analysis
#' @param object a spMOCA object
#' @param gene.use the features used to calculate module scores. The default is hub genes.
#' @param fgsea_sets a list of gene set 
#' @return Return a list of odds ratio and p-value tables for each module per each gene set

geneSetAnalysis = function(object,
                          gene.use = "hubgene",
                          fgsea_sets){
  gene.module.all = object@network
  bg_gene = object@feature_select
  cat(paste0("## Performing Gene Set Analysis ...\n"))
  if(gene.use == "hubgene"){
    hubgene.module.all = gene.module.all %>%
      dplyr::group_by(k) %>%
      dplyr::top_frac(0.05,kWithin)
    
    gene_cluster = hubgene.module.all
    gene_num_modular = length(unique(hubgene.module.all$k))
    
    binary_gene_set = list()
    count_gene_set = list()
    
    odds_gene_set = matrix(0,nrow = gene_num_modular,ncol = length(fgsea_sets))
    pval_gene_set = matrix(0,nrow = gene_num_modular,ncol = length(fgsea_sets))
    
    module = gtools::mixedsort(unique(gene.module.all$k))
    colnames(odds_gene_set) =
      colnames(pval_gene_set) = names(fgsea_sets)
    rownames(odds_gene_set) = rownames(odds_gene_set) = module
    
    for(k in 1:gene_num_modular){
      print(k)
      member = gene_cluster$gene[gene_cluster$k == module[k]]
      Geneset = matrix(0,nrow = length(member),ncol = length(fgsea_sets))
      rownames(Geneset) = member
      colnames(Geneset) = names(fgsea_sets)
      
      
      for(ig in 1:length(fgsea_sets)){
        if(sum(member %in% fgsea_sets[[ig]]) > 0){
          # entry of binary matrix
          Geneset[member %in% fgsea_sets[[ig]],ig] = 1
          # entry of odds ratio matrix
          binary_all_gene_in_module = ifelse(bg_gene %in% member,1,0) # in-module-indicators of selected genes 
          binary_all_gene_in_set = ifelse(bg_gene %in% fgsea_sets[[ig]],1,0) # in-gene-set-indicators of selected genes
          cont.table = table(binary_all_gene_in_module,binary_all_gene_in_set)
          cont.table = cont.table[c(2,1),c(2,1)]
          odds_gene_set[k,ig] = as.numeric(fisher.test(cont.table)$estimate)
          pval_gene_set[k,ig] = as.numeric(fisher.test(cont.table)$p.value)
        }
      }
      binary_gene_set[[module[k]]] = Geneset
    }
    
    # exclude gene set
    count_in_set = matrix(0,nrow = gene_num_modular,ncol = length(fgsea_sets))
    for(k in 1:gene_num_modular){
      count_in_set[k,] = colSums(binary_gene_set[[k]])
    }
    exclude_by_count = which(colSums(count_in_set) == 0) # gene sets have no gene in all modules
    exclude_by_odds = which(is.infinite(colSums(odds_gene_set))) # gene sets have infinite odd ratios
    exclude_col = union(exclude_by_count,exclude_by_odds) 
    
    odds_gene_set[which(odds_gene_set ==0,arr.ind = T)] =  NA
    pval_gene_set[which(pval_gene_set ==0,arr.ind = T)] =  NA
  }
  return(list(odds = odds_gene_set,
              pval = pval_gene_set))
}