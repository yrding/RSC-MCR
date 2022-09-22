library(fpc)
source("FUM_function.R")
source('DB.R')
lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc=t(lnc)
mi=t(mi)
m=t(m)
lnc1=read.table("result_lnc.txt",sep="\t",header=F);
mi1=read.table("result_mi.txt",sep="\t",header=F);
m1=read.table("result_m.txt",sep="\t",header=F);
data=cbind(lnc,mi,m)
data1=cbind(lnc1,mi1,m1)

result <- FCM(data,K=2)
labels=result[["cluster_id"]]
b=result[["u"]]
dis=dist(b)
stats <- cluster.stats(dis, labels)
stats$avg.silwidth
calDBI(b, labels)


result <- FCM(data1,K=2)
labels=result[["cluster_id"]]
b=result[["u"]]
dis=dist(b)
stats <- cluster.stats(dis, labels)
stats$avg.silwidth
calDBI(b, labels)

result1 = FCM(dff,K=2)
b1=result1[["u"]]
dis1=dist(b1)
labels1=result1[["cluster_id"]]
stats <- cluster.stats(dis1, labels1)
stats$avg.silwidth
calDBI(b1, labels1)

result1 = FCM(ICA$S,K=2)
b1=result1[["u"]]
dis1=dist(b1)
labels1=result1[["cluster_id"]]
stats <- cluster.stats(dis1, labels1)
stats$avg.silwidth
calDBI(b1, labels1)
