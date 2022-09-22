

lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc=t(lnc)
mi=t(mi)
m=t(m)
lnc=scale(lnc)
mi=scale(mi)
m=scale(m)
data=cbind(lnc,mi,m)
df=as.data.frame(data)
df.pr<- prcomp(df, center = F,scale. = F)
d=summary(df.pr) 
dff=df.pr$x[,1:160]
#85
write.table(d[["importance"]],"ddd.txt",row.names = F, col.names = F, quote = F,sep="\t")
