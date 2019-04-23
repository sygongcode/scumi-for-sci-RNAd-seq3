# =================================================================
# 
# Jiarui Ding (jding@broadinstitute.org)
# 
# Gvien a gene by cell UMI raw count matrix 
# (or TPM matrix for Smart-seq2 data), the clustering of cells, 
# and a list of marker genes (examples are in the marker_gene folder), 
# automatically assigning cell types to clusters and 
# computing the AUCs for the assigned cell types
#

library(ROCR)
 
CalculateAuc = function(umi.count, cluster, marker, merge=TRUE, 
                        unassigned.threshold=0.65) {
  res = ComputeClusterAuc(umi.count, as.character(cluster), marker)
  auc = res[[1]]
  
  out = ComputeMergeClusterAuc(res[[3]], 
    as.character(cluster), auc, res[[2]], marker=marker, 
    merge=merge, unassigned.threshold=unassigned.threshold)
  
  out 
}

ComputeClusterAuc = function(x, cluster, marker) {
  marker = SplitMarker(marker)
  
  cell.type.score = sapply(marker, function(z) {
    res = ComputeMarkerScore(x=x, marker.list = z)
  })
  
  marker.name = names(marker)
  uniq.cluster = unique(cluster)
  
  auc = setNames(vector('list', length(marker.name)), marker.name)
  for (cell.type in marker.name) {
    res = sapply(uniq.cluster, function(z) {
      ScoreAuc(cell.type.score[, cell.type], cluster == z)
    })
    auc[[cell.type]] = unlist(res)
  }
  
  auc = do.call(rbind, auc)
  auc.assign = AssignCluster(auc)
  
  list(auc.assign, auc, cell.type.score)
}

AssignCluster = function(auc.table) {
  res = apply(auc.table, 2, function(z) {
    id = which.max(z)
    res = z[id]
    names(res) = rownames(auc.table)[id]
    
    list(res, names(res))
  })
  
  auc.assign = sapply(res, function(z) z[[1]])
  name = sapply(res, function(z) z[[2]])
  
  names(auc.assign) = name
  
  auc.assign
}

ReassignCluster = function(x, cluster, auc, auc.table, marker) {
  marker = SplitMarker(marker)
  
  cluster.name = colnames(auc.table)
  cell.type = names(auc)
  names(cell.type) = cluster.name
  
  cell.type.candidate = unique(cell.type)
  
  auc.table.update = auc.table[cell.type.candidate, , drop=FALSE]
  auc.table = auc.table.update
  
  for (z in cluster.name) {
    res = sapply(cell.type.candidate, function(z1) {
      if (max(auc.table.update[, z]) > 0.9) {
        return(NULL)
      }
      
      auc.one = MergeDetector(x, z, z1, cluster, cell.type)
      if (sum(cell.type==cell.type[z])==1 & auc.one > 0.55) {
        return(NULL)
      }
      
      id.rm = which(cell.type == z1)
      cluster.rm = cluster.name[id.rm]
      cluster.rm = setdiff(cluster.rm, z)
      
      id.keep = !(cluster %in% cluster.rm)
      
      cluster.keep = cluster[id.keep]
      cluster.score = x[id.keep, z1]
      
      auc.two = ScoreAuc(score = cluster.score, cluster.keep %in% z)
    }, simplify = TRUE)
    
    id.keep = sapply(res, length) > 0
    auc.table.update[id.keep, z] = unlist(res[id.keep])
  }
  
  auc.update = AssignCluster(auc.table.update)
  
  return(list(auc=auc.update, auc.table=auc.table.update))
}

MergeDetector = function(x, cluster, cell.type, cluster.all, cell.type.all) {
  if (cell.type.all[cluster] == cell.type) {
    return (0.5)
  }
  
  id.rm = which(cell.type.all == cell.type)
  
  cluster.name = names(cell.type.all)
  cluster.rm = cluster.name[id.rm]
  
  id.keep = cluster.all %in% c(cluster, cluster.rm)
  
  cluster.keep = cluster.all[id.keep]
  cluster.score = x[id.keep, cell.type.all[cluster]]
  
  auc = ScoreAuc(score = cluster.score, cluster.keep %in% cluster)
}

ScoreAuc = function(score, label) {
  if(length(table(label)) != 2) {
    return(0.5)
  }
  
  pred = ROCR::prediction(predictions=score, label)
  
  roc.perf = performance(pred, measure = 'tpr', 
                         x.measure = 'fpr')
  
  auc.perf = performance(pred, measure = 'auc')
  res = auc.perf@y.values
  
  res[[1]]
}

ComputeMergeClusterAuc = function(x, cluster, auc, auc.table, 
                                  marker, merge=TRUE, 
                                  unassigned.threshold=0.65) {
  
  if (merge == TRUE) {
    res = ReassignCluster(x, cluster, auc, auc.table, marker = marker)
    auc.update = res[[1]]
    auc.table.update = res[[2]]
    
    num.iter = 1
    while (sum(names(auc) != names(auc.update)) > 0 & num.iter <= 5) {
      num.iter = num.iter + 1
      
      auc = auc.update
      auc.table = auc.table.update
      
      res = ReassignCluster(x, cluster, auc, auc.table, marker = marker)
      auc.update = res[[1]]
      auc.table.update = res[[2]]
    }
    auc = auc.update
    auc.table = auc.table.update
  }
  
  cluster.name = colnames(auc.table)
  cell.type = names(auc)

  res = apply(auc.table, 2, function(z) max(z) < unassigned.threshold)
  id.unassigned = which(res)
  
  if (length(id.unassigned) > 0) {
    cell.type[id.unassigned] = 'Unassigned'
    auc.unassigned = mean(auc[id.unassigned])
    cluster.unassigned = names(id.unassigned)
  }
  
  auc = sapply(setdiff(unique(cell.type), 'Unassigned'), function(z) {
    id = which(cell.type == z)
    cluster.merge = cluster.name[id]
    
    cluster.score = x[, z]
    ScoreAuc(cluster.score, cluster %in% cluster.merge)
  })
  auc = auc[!duplicated(names(auc))]
  
  label = ConvertClusterLabel(cluster = cluster, 
                              auc = AssignCluster(auc.table), 
                              auc.table = auc.table)
  
  if (length(id.unassigned) > 0) {
    auc['Unassigned'] = auc.unassigned
    label[cluster %in% cluster.unassigned] = 'Unassigned'
  }
  
  return(list(auc=auc, label=label))
}

ComputeMarkerScore = function(x, marker.list) {
  total = Matrix::colSums(x)
  
  gene.pos = intersect(marker.list[[1]], rownames(x))
  gene.neg = intersect(marker.list[[2]], rownames(x))
  
  if (length(gene.pos) > 0 & length(gene.neg) > 0) {
    out = Matrix::colSums(x[gene.pos, , drop=FALSE]) -
      Matrix::colSums(x[gene.neg, , drop=FALSE])
    
  } else if (length(gene.pos) > 0) {
    out = Matrix::colSums(x[gene.pos, , drop=FALSE])
    
  } else if (length(gene.neg) > 0) {
    out = Matrix::colSums(x[gene.neg, , drop=FALSE])
    out = max(out) - out
  }
  
  res = out / total * 10000
  res[res < 0] = 0
  
  res = suppressWarnings(log1p(res))
  
  res
}

SplitMarker = function(marker) {
  res = sapply(marker, function(z) {
    gene = gsub('.{1}$', '', z)
    
    id = grep(z, pattern = '-$')
    negative.marker = gene[id]
    
    id = grep(z, pattern = '+$')
    positive.marker = gene[id]
    
    positive.marker = setdiff(positive.marker, negative.marker)
    
    res = list(positive.marker, negative.marker)
    res
  }, simplify = FALSE)
  
  res
}


ConvertClusterLabel = function(cluster, auc, auc.table) {
  cell.type = names(auc)
  cluster.name = colnames(auc.table)
  
  names(cell.type) = cluster.name
  label = cell.type[cluster]
  
  label
}

