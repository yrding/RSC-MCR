install.packages("Hmisc")
library(Hmisc)

lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc=t(lnc)
mi=t(mi)
m=t(m)
data=cbind(lnc,mi,m)
res2 <- rcorr(as.matrix(data))
#write.table(r,"r.txt",row.names = F, col.names = F, quote = F,sep="\t");
res2$P[is.na(res2$P)] <- 0
res2$r[is.na(res2$r)] <- 0
p=res2$P
r=res2[["r"]]
r[p>0.05]=0

A=r[1:1492,1493:3199] #lnc-mi+m
B1=r[1493:1684,1:1492]
B2=r[1493:1684,1685:3199]
B=cbind(B1,B2)  #mi-lnc+m
C=r[1685:3199,1:1684]#m-lnc+mi

a=t(A)
b=t(B)
c=t(C)

write.table(a,"a1.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(b,"b1.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(c,"c1.txt",row.names = F, col.names = F, quote = F,sep="\t")


lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc1=read.table("a1.txt",sep="\t",header=F);
mi1=read.table("b1.txt",sep="\t",header=F);
m1=read.table("c1.txt",sep="\t",header=F);
lnc1=t(lnc1)
mi1=t(mi1)
m1=t(m1)
lnc=t(lnc)
mi=t(mi)
m=t(m)
sss=function(dataset){
  marg.center <- apply(dataset, 2, median)
  marg.scale  <- apply(dataset, 2, mad)
  X <- scale(dataset, center = marg.center, scale = marg.scale)
  return(X)
}

lnc=sss(lnc) 
mi=sss(mi) 
m=sss(m) 
lnc1=sss(lnc1) 
mi1=sss(mi1) 
m1=sss(m1) 

projectRes1=lnc %*% lnc1
projectRes2=mi %*% mi1
projectRes3=m %*% m1

lncRNA=lnc-  projectRes1
miRNA=mi- projectRes2  
mRNA=m-projectRes3
write.table(lncRNA,"p_lnc.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(miRNA,"p_mi.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(mRNA,"p_m.txt",row.names = F, col.names = F, quote = F,sep="\t")

