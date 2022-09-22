install.packages("fastICA")
library(fastICA)
lnc=read.table("lnc_trans.txt",sep = '\t',header = F)
mi=read.table("mi_trans.txt",sep = '\t',header = F)
m=read.table("m_trans.txt",sep = '\t',header = F)
lnc=t(lnc)
mi=t(mi)
m=t(m)
data=cbind(lnc,mi,m)
n.comp=3
X=data
ICA=fastICA(X, n.comp, alg.typ = c("parallel","deflation"),
        fun = c("logcosh","exp"), alpha = 1.0, method = c("R","C"),
        row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE,
        w.init = NULL)
b1=b$S[,1]
b2=b$S[,2]
b3=b$S[,3]
