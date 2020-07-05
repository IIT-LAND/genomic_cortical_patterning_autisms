ac_sparsePCA_gene_finder <- function(X,Y,PCnum){
  result1cv <- acPCAtuneLambda(X=X, Y=Y, nPC=2, lambdas=seq(0, 20, 0.1),
                               anov=T, kernel = "linear", quiet=T)
  result1 <- acPCA(X=X, Y=Y, lambda=result1cv$best_lambda, kernel="linear", nPC=2)
  v_ini <- as.matrix(result1$v[,PCnum])
  # par(mfrow=c(1,2), pin=c(2.5,2.5), mar=c(4.1, 3.9, 3.2, 1.1))
  c2s <- seq(1, 0, -0.1)*sum(abs(v_ini))
  resultcv_spc1_coarse <- acSPCcv( X=X, Y=Y, c2s=c2s, v_ini=v_ini,
                                   kernel="linear", quiet=T, fold=10, plot=FALSE)
  c2s <- seq(0.9, 0.7, -0.02)*sum(abs(v_ini))
  resultcv_spc1_fine <- acSPCcv( X=X, Y=Y, c2s=c2s, v_ini=v_ini,
                                 kernel="linear", quiet=T, fold=10, plot=FALSE)
  result_spc1 <- acSPC( X=X, Y=Y, c2=resultcv_spc1_fine$best_c2,
                        v_ini=v_ini, kernel="linear")
  v1 <- result_spc1$v
  sum(v1!=0)
  
  gene_mask = v1!=0
  return(gene_mask)
}

