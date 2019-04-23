# =====
# Some utility functions (using Seurat v2.3.4) to clustering scRNA-seq data
# and remove cluster than are likey consiting of low-quality cells
# 
# umi.count is a raw UMI count matrix (or TPM matrix for Smart-seq2 data)
# 

library(ROCR)
library(Seurat)
library(dplyr)


# =====
ClusterSeurat = function(umi.count, name='PBMC1_Drop-seq', number.neighbor=30, 
                         variable.gene=FALSE, resolution=1.0, num.pc=50, 
                         species='human', do.plot=FALSE, normalize=TRUE, 
                         min.cells=3, min.genes=1) {
  
  seu.obj = CreateSeuratObject(raw.data = umi.count, min.cells = min.cells, 
                               min.genes = min.genes, project = name)
  
  if (species == 'human') {
    mito.genes = grep('MT-', rownames(umi.count), value = TRUE)
  } else if (species == 'mouse') {
    mito.genes = grep('mt-', rownames(umi.count), value = TRUE)
  } else {
    stop('species must be human or mouse!')
  }
  
  mito.ratio = Matrix::colSums(umi.count[mito.genes, ]) / Matrix::colSums(umi.count)
  
  seu.obj = AddMetaData(object = seu.obj, metadata = mito.ratio, 
                        col.name = "percent.mito")
  
  # seu.obj@meta.data[, 'scale.ngene'] = scale(seu.obj@meta.data$nGene)
  
  if (TRUE == do.plot) {
    VlnPlot(object = seu.obj, features.plot = c("nGene", "nUMI", "percent.mito"), 
            nCol = 3)
    
    op = par(mfrow = c(1, 2))
    GenePlot(object = seu.obj, gene1 = "nUMI", gene2 = "percent.mito")
    GenePlot(object = seu.obj, gene1 = "nUMI", gene2 = "nGene")
    par(op)
  }
  
  if (TRUE == normalize) {
    seu.obj = NormalizeData(object = seu.obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
  } else {
    # Smart-seq2 data 
    seu.obj@data = as.matrix(log1p(seu.obj@raw.data))
    seu.obj@calc.params$NormalizeData$normalization.method = 'TPM'
  }
  
  seu.obj = FindVariableGenes(object = seu.obj, mean.function = ExpMean,
                              dispersion.function = LogVMR, do.plot = do.plot,
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
  # ===================
  if (species == 'human') {
    mt.rb.gene = grep('^MT-|^RPL|^RPS|^MRPS|^MRPL', rownames(seu.obj@raw.data), value = TRUE)
  } else if (species == 'mouse') {
    mt.rb.gene = grep('^mt-|^Rpl|^Rps|^Mrps|^Mrpl', rownames(seu.obj@raw.data), value = TRUE)
  } else {
    stop('species must be human or mouse!')
  }
  
  # seu.obj = ScaleData(object = seu.obj, vars.to.regress = c("nUMI", "percent.mito"))
  seu.obj = ScaleData(object = seu.obj)
  # seu.obj = ScaleData(object = seu.obj, vars.to.regress = c("nUMI"))
  
  if (variable.gene) {
    seu.obj = RunPCA(object = seu.obj,
                     do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute=num.pc)
  } else {
    seu.obj = RunPCA(object = seu.obj, pc.genes = setdiff(rownames(seu.obj@raw.data), c(mt.rb.gene)), 
                     do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute=num.pc)
  }
  
  
  PrintPCA(object = seu.obj, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  
  if (TRUE == do.plot) {
    VizPCA(object = seu.obj, pcs.use = 1:2)
    
    PCAPlot(object = seu.obj, dim.1 = 1, dim.2 = 2)
  }
  
  seu.obj = ProjectPCA(object = seu.obj, do.print = FALSE)
  
  # PCHeatmap(object = seu.obj, pc.use = 1, cells.use = 500, 
  #           do.balanced = TRUE, label.columns = FALSE)
  
  seu.obj = FindClusters(object = seu.obj, reduction.type = "pca", dims.use = 1:num.pc, 
                         resolution = resolution, print.output = 0, 
                         save.SNN = TRUE, k.param=number.neighbor)
  
  seu.obj = RunTSNE(object = seu.obj, dims.use = 1:num.pc, do.fast = TRUE)
  if (TRUE == do.plot) {
    TSNEPlot(object = seu.obj, do.label=TRUE)
  }
  
  seu.obj
} 

SelectCluster = function(x, gene.to.remove, gene.to.ignore, 
                         num.top.gene=15, ratio.th=0.7, 
                         min.percent=0.25, max.percent=1.0) {
  top.gene = x %>% filter(p_val_adj <= 0.01) %>%  
    filter(!(gene %in% gene.to.ignore)) %>% 
    filter(pct.2 <= max.percent) %>% 
    group_by(cluster) %>% dplyr::slice(1:num.top.gene) %>% 
    filter(!(gene %in% gene.to.remove)) %>% 
    filter(pct.1 > min.percent) %>% 
    as.data.frame()
  
  count.marker.gene = table(top.gene$cluster)
  ratio = count.marker.gene / num.top.gene
  
  cluster = levels(x$cluster)
  ratio = ratio[cluster]
  cluster.keep = names(ratio[which(ratio > ratio.th)])
  
  return(list(cluster=cluster.keep, ratio=ratio))
}

# =====
ResetGeneName = function(umi.count) {
  gene_id_name = rownames(umi.count)
  gene = strsplit(gene_id_name, '_')
  
  id = sapply(gene, function(z) length(z) > 2)
  gene[id] = sapply(gene[id], function(z) {
    # z = z[z != ""]
    num.field = length(z)
    zz = c(z[1], NULL)
    zz[2] = paste(z[2:num.field], collapse ='_')
  })
  
  gene = do.call(rbind, gene)[, 2]
  gene[duplicated(gene)] = gene_id_name[duplicated(gene)]
  
  rownames(umi.count) = gene
  
  umi.count
} 

# =====
FilterByMito = function(umi.count, return.ratio=FALSE, species='human') {
  num.read = Matrix::colSums(umi.count)
  
  if (species == 'human') {
    mito.genes = grep('MT-', rownames(umi.count), value = TRUE)
  } else if (species == 'mouse') {
    mito.genes = grep('mt-', rownames(umi.count), value = TRUE)
  } else {
    stop('species must be human or mouse!')
  }
  
  mito.ratio = Matrix::colSums(umi.count[mito.genes, ]) / num.read
  mito.th = quantile(mito.ratio, 0.75) + 3 *IQR(mito.ratio)
  
  if (mito.th <= 0.02) {
    print('Mito-th \n')
    mito.th = max(max(mito.ratio) / 2, mito.th)
  }
  
  plot(log10(num.read), mito.ratio)
  abline(h=mito.th)
  
  if (TRUE == return.ratio) {
    return(list(umi.count[, mito.ratio <= mito.th], mito.ratio))
  }
  
  umi.count[, mito.ratio <= mito.th]
}

SelectClusterObj = function(seu.obj, cluster.marker, test.use = "wilcox", 
                            latent.vars=NULL, num.top.gene=15, 
                            species='human', ...) {
  if (species == 'human') {
    rb.gene = grep('RPL|RPS|MRPS|MRPL', rownames(seu.obj@raw.data), value = TRUE)
    mt.gene = grep('MT-', rownames(seu.obj@raw.data), value = TRUE)
    
    prob.gene = grep('MALAT1|^MTRNR', rownames(seu.obj@raw.data), value = TRUE)
  } else if (species == 'mouse') {
    rb.gene = grep('Rpl|Rps|Mrps|Mrpl', rownames(seu.obj@raw.data), value = TRUE)
    mt.gene = grep('mt-', rownames(seu.obj@raw.data), value = TRUE)
    
    prob.gene = grep('Malat1', rownames(seu.obj@raw.data), value = TRUE)
  } else {
    stop('species must be human or mouse!')
  }
  
  if (missing(cluster.marker)) {
    cluster.marker = FindAllMarkers(object = seu.obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25, test.use = test.use, 
                                    latent.vars=latent.vars) 
  }
  
  cluster = SelectCluster(cluster.marker, gene.to.remove = mt.gene, 
                          gene.to.ignore = c(rb.gene, prob.gene), 
                          num.top.gene = num.top.gene, ...)
  
  selected.cell = names(seu.obj@ident[as.character(seu.obj@ident) %in% cluster$cluster])
  
  return(list(cell=selected.cell, marker=cluster.marker, cluster=cluster))
} 
