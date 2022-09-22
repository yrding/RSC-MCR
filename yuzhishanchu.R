#library(latticeExtra)
#install.packages("latticeExtra")
lnc=read.table("lnc_trans.txt",sep="\t",header=F);
mi=read.table("mi_trans.txt",sep="\t",header=F);
m=read.table("m_trans.txt",sep="\t",header=F);
lnc=t(lnc)
mi=t(mi)
m=t(m)
data=cbind(lnc,mi,m)

i=as.list.numeric_version(i0)
res=cor(data,method = "pearson")
res[res>0.5]=0
res[res<-0.5]=0
n0 <- apply(res == 0, 1, sum)
# 统计0元素大于2个的行号
i0 <- which(n0 > 1)
# 删除0元素大于2个的行
dat=t(data)
#rr=dat[-i0, ]
i1=i0[1:1050]
i2=i0[1051:1159]-1050
i3=i0[1160:2574]-1160
lnc_delet=lnc[,-i1]
mi_delet=mi[,-i2]
m_delet=m[,-i3]

write.table(rr,"yuzhishanchu.txt",row.names = F, col.names = F, quote = F,sep="\t")
#write.table(i0,"ii.txt",row.names = F, col.names = F, quote = F,sep="\t")
