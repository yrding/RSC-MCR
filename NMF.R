lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
install.packages("NMF")
library(NMF)
lnc=t(lnc)
mi=t(mi)
m=t(m)
data=cbind(lnc,mi,m)
data=t(data)
res=nmf(data,2)


b=res@fit@H
b=t(b)
labels=t(group)
labels=as.numeric(labels)

library(fpc)
stats=cluster.stats(dist(b),labels)
stats[["avg.silwidth"]]
dbi <- calDBI(b, labels)

group <- predict(res)

group <- as.data.frame(group)

#group$group <- paste0('Cluster',group$group)

#group$sample <- rownames(group)

#group <- group[order(group$group),]

table(group$group) 
