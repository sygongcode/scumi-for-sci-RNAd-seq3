library(cluster)

# Given a set of marker genes (examples in the marker_gene folder), 
# and the raw UMI count matrix (or TPM matrix for Smart-seq2 data), 
# turning parameter to produce clustering 
# that maximize the avarge AUCs

SelectClusterModel = function(umi.count, marker, 
                              number.neighbor=c(15, 30), 
                              resolution = c(0.8, 1.0, 1.2, 1.5),
                              num.pc = c(20, 30, 50, 100),
                              variable.gene = c(TRUE, FALSE),
                              min.cluster=3, 
                              max.cluster=20, normalize=TRUE, 
                              min.cells=3, min.genes=1, ...) {
  
  seu.obj = CreateSeuratObject(raw.data = umi.count, min.cells = min.cells, 
                               min.genes = min.genes)
  
  if (normalize == TRUE) {
    seu.obj = NormalizeData(object = seu.obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)
  } else {
    seu.obj@data = as.matrix(log1p(seu.obj@raw.data))
    seu.obj@calc.params$NormalizeData$normalization.method = 'TPM'
  }

  gene = unique(gsub('.{1}$', '', unlist(marker)))
  gene = intersect(gene, rownames(seu.obj@data))
  dis = dist(t(seu.obj@data[gene, , drop=FALSE]))^2
  
  umi.count = umi.count[, names(seu.obj@ident)]
  
  out = vector('list', length(number.neighbor) * 
                length(resolution) * length(num.pc) * length(variable.gene))
  para = out
  
  iter = 0
  for (v in variable.gene) 
    for (k in number.neighbor) 
      for (r in resolution) 
        for (n in num.pc) { 
          iter = iter + 1
          cat(iter, '\n')
          
          para[[iter]] = list(number.neighbor=k, variable.gene=v, 
                              resolution=r, num.pc=n)
          
          res = ClusterSeurat(umi.count, 
                              number.neighbor = k, 
                              variable.gene = v, 
                              resolution = r, 
                              num.pc = n, do.plot = FALSE, 
                              normalize=normalize, 
                              min.cells=min.cells, min.genes=min.genes,
                              ...)
          
          cluster = as.numeric(res@ident)
          
          cluster.count = table(cluster)
          num.cluster = length(cluster.count)
          
          sil = -1
          auc = 0
          if (num.cluster >= min.cluster & num.cluster <= max.cluster) {
            sil = silhouette(cluster, dis)
            sil = mean(sil[, 3])
            
            min.cell.in.cluster = max(3, k / 5.0)
            cluster.keep = names(cluster.count[cluster.count >= min.cell.in.cluster])
            cell.keep = which(cluster %in% as.numeric(cluster.keep))
            
            auc = ComputeClusterAuc(umi.count[, cell.keep, drop=FALSE], 
                                    cluster[cell.keep], marker)[[1]]
          }
          
          out[[iter]] = list(cluster=cluster, auc=auc, sil=sil, 
                             num.cluster=num.cluster)
        }
  
  # Finall results, parameter sets, clustering, AUCs 
  res = ExtractResults(para, out, min.cluster=min.cluster, 
                       max.cluster = max.cluster)
  
  para.sil = para[[res$model.sil]]
  para.auc = para[[res$model.auc]]
  
  # 
  auc.sil = CalculateAuc(umi.count, 
                         cluster=out[[res$model.sil]][[1]], 
                         marker = marker)$auc
  
  auc.auc = CalculateAuc(umi.count, 
                         cluster=out[[res$model.auc]][[1]], 
                         marker = marker)$auc
  
  # 
  auc.final = auc.auc
  para.final = para.auc
  
  if ((sum(auc.sil) - 0.5 * length(auc.sil))  > 
      (sum(auc.auc) - 0.5 * length(auc.auc))) {
    auc.final = auc.sil
    para.final = para.sil
  }
  
  model = ClusterSeurat(umi.count, 
                        number.neighbor = para.final$number.neighbor, 
                        variable.gene = para.final$variable.gene, 
                        resolution = para.final$resolution, 
                        num.pc = para.final$num.pc, do.plot = FALSE, 
                        normalize=normalize, 
                        min.cells=min.cells, min.genes=min.genes,
                        ...)
  TSNEPlot(model, do.label = TRUE)
  
  return(list(auc.final=auc.final, para.final=para.final, model.final=model, 
              para=para, model.out=out))
}


ExtractResults = function(parameter, model.out, 
                          min.cluster=5, max.cluster=20) {
  num.cluster = unlist(sapply(model.out, function(z) z[[4]]))
  sil = unlist(sapply(model.out, function(z) z[[3]]))
  # auc = unlist(sapply(model.out, function(z) mean(z[[2]])))
  
  id = which(num.cluster >= min.cluster & 
               num.cluster <= max.cluster)
  
  if (length(id) == 0) {
    id = seq(length(num.cluster))
  }
  
  model.sil = which.max(sil[id])
  # model.auc = which.max(auc)
  
  model.auc.avg = sapply(model.out, 
             function(z) AverageAUC(z[[2]]))
  model.auc = which.max(model.auc.avg[id])
  
  return(list(model.sil=id[model.sil], model.auc=id[model.auc]))
}


AverageAUC = function(x) {
  cell.type = unique(names(x))
  num.cluster = length(cell.type) 
  
  if (length(x) >= length(cell.type) * 2) {
    return(0)
  }
  
  res = sapply(cell.type, function(z) {
    id = names(x) == z
    mean(x[id]) - 0.7 
  }, simplify = TRUE)

  if (is.list(res)) {
    return(0)
  }
  
  sum(res)
}


