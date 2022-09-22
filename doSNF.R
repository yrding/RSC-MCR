#' doSNF
#'
#' @param data List of matrices.
#' @param K Number of clusters
#' @param sigma Variance for local model
#' @param K_n Number of nearest neighbors
#'
#' @return a list of \code{clust} the clustering of samples and
#' \code{fit} the results of the method SNF
#' @export
#'
#' @examples
#' set.seed(333)
#' c_1 <- simulateY(J=1000, prop=0.1, noise=1)
#' c_2 <- simulateY(J=2000, prop=0.1, noise=1)
#' c_3 <- simulateY(J=500, prop=0.1,  noise=0.5)
#' data <- list(c_1$data , c_2$data , c_3$data)
#' res <- doSNF(data,K=4, K_n=10, sigma=0.5)
#' @import SNFtool
#' @importFrom dplyr %>%
#' @export
library(SNFtool)
library(dplyr)
library(fpc)
lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc=t(lnc)
mi=t(mi)
m=t(m)
lnc1=read.table("p_lnc.txt",sep="\t",header=F);
mi1=read.table("p_mi.txt",sep="\t",header=F);
m1=read.table("p_m.txt",sep="\t",header=F);
data=list(lnc,mi,m)
data1=list(as.matrix(lnc1),as.matrix(mi1),as.matrix(m1))
doSNF <- function (data, K,  K_n=10, sigma=0.5){
  dat <- lapply(data, function (dd){
    dd <- dd %>% as.matrix
    W <- dd %>% dist2(dd) %>% affinityMatrix(K=K_n, sigma=sigma)
  })
  W <-  SNF(dat, K_n, K_n)
  clust.SNF = W %>% spectralClustering(K)
  res <- list(clust=clust.SNF, fit= W)
  return(res)
}

res <- doSNF(data,K=2, K_n=10, sigma=0.5)
res1 <- doSNF(data1,K=2, K_n=10, sigma=0.5)
label=res[["clust"]]
b=res[["fit"]]
stats <- cluster.stats(b,label)
stats$avg.silwidth
calDBI(b,label)

label1=res1[["clust"]]
b1=res1[["fit"]]
stats <- cluster.stats(b1,label1)
stats$avg.silwidth
calDBI(b1,label1)


Dist1 = dist2(as.matrix(dff),as.matrix(dff));
K=10
alpha=0.5
W1 = affinityMatrix(Dist1, K, alpha)
utrue_label=spectralClustering(W1,2)
stats <- cluster.stats(W1,utrue_label)
stats$avg.silwidth
calDBI(W1,utrue_label)

data=list(ICA$S[,1],ICA$S[,2],ICA$S[,3])
res <- doSNF(data,K=2, K_n=10, sigma=0.5)
label=res[["clust"]]
b=res[["fit"]]
stats <- cluster.stats(b,label)
stats$avg.silwidth
calDBI(b,label)