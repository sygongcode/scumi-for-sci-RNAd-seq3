FitLogCount = function(x, show.plot=TRUE, loglik=FALSE, filter.ratio=0.1, 
                       xlab="Expression", ylab="Density", take.log=TRUE, ...) {
  # Given an input sparse gene by cell matrix x (raw UMI count matrix or read matrix), 
  # use Student's t mixture model to remove likely low-quality cells 
  # that have small numbers of UMIs or reads 
  # (if x is the raw uncollapsed gene by cell read matrix)
  # 
  # First install the xseq R package (in this folder) from terminal:
  # R CMD build xseq
  # R CMD install xseq_0.2.2.tar.gz  
  #
 
  x.cell = Matrix::colSums(x)
  
  if (take.log==TRUE) {
    expr.quantile = log10(x.cell)
  } else {
    expr.quantile = x.cell
  }
  
  th = quantile(expr.quantile, filter.ratio)
  
  mu = c(mean(expr.quantile[expr.quantile<=th]), 
         mean(expr.quantile[expr.quantile>th]))
  
  sigma  = c(sd(expr.quantile[expr.quantile<=th]), 
             sd(expr.quantile[expr.quantile>th]))
  
  lambda = c(filter.ratio, 1 - filter.ratio)
  
  prior = list()
  prior$lambda = lambda
  prior$mu     = mu
  prior$sigma  = sigma
  
  prior$alpha = c(5, 5)
  
  prior$kapp = 0.2
  prior$dof  = 3 # max(length(x.cell) * 0.05, 2)
  
  ## 
  model = xseq:::MixStudentFitEM(expr.quantile, lambda=lambda, prior = prior, 
                                 mu=mu, sigma=sigma, K=2, nu.equal = TRUE)
  if (show.plot == TRUE) {
    xseq:::MixModelPlot(model,  xlab2=xlab, ylab2=ylab, 
                        breaks=40, loglik=loglik, ...)
  }
  
  range.x = range(expr.quantile)
  x.range = seq(range.x[1], range.x[2], diff(range.x) / 1000)
  
  expec = xseq:::MixStudentProb(x.range, K=2, model$lambda, 
                                model$mu, model$sigma, model$nu)
  
  post.prob = expec$tau
  id.cross = which(post.prob[, 1] > post.prob[, 2])
  
  id = which(diff(id.cross) > 1)
  if (length(id) >= 1) {
    id.cross = id.cross[id[1]]
  } else {
    id.cross = id.cross[length(id.cross)]
  }
  
  abline(v = x.range[id.cross], lwd=2.5, col='dodgerblue')
  
  id = which(expr.quantile >= x.range[id.cross])
  if (length(id) > 0) {
    x = x[, id]
  }
  
  x
}


