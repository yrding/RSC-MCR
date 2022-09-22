#这是一种鲁棒相关矩阵估计器的实现，在高维数据的情况下有效地消除了偏差。
#采用交叉验证技术实现了自适应阈值技术，对相关矩阵进行了正则化。
library(parallel) 
source("cormad.R")
source("thresholding.R")
source("bin_matrix.R")
source("RSC_thresholding.R")
source("robust_scaling.R")

#导入已经做好标准化、筛选、批次处理的数据

K = 17;		
alpha = 0.5;
T =20; 
library('readr')
library("SNFtool");
library(fpc)
data=read.table("yangjing_result.txt",header=T,row.names=1,sep="\t")#

data=t(data)
data=apply(data,2,as.numeric);
sample_size = 39
#多核运算，看电脑配置进行修改（只改变速度不改变结果）
n_cores = 8
nsplits=25
chips = sample(x = 1:nrow(data),size = sample_size,replace = FALSE)
red_data = data[chips,]
print("Evaluating RSC Correlation...\n")

RSC = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores,monitor = 0)
R_vec = cormad.vec(red_data,jobs = n_cores)
R = cormad.vec2mat(R_vec)
R_bin = bin_matrix(R,RSC$threshold)


#RSC矩阵：gene-gene
temp=RSC$avgloss         #RSC=R_bin
RSC1=RSC_thresholding(temp,R)         #RSC=RSC_thresholding(FLOSSES,R)
write.table(RSC1,"RSC.txt",row.names = 1, col.names = 1, quote = F,sep="\t")


RSC1=read.table("RSC.txt",sep="\t",header=F);
A=RSC1[1:127,128:382] #lnc-mi+m
B1=RSC1[128:254,1:127]
B2=RSC1[128:254,255:382]
B=cbind(B1,B2)  #mi-lnc+m
C=RSC1[255:382,1:254]#m-lnc+mi

a=t(A)
b=t(B)
c=t(C)

write.table(a,"a1.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(b,"b1.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(c,"c1.txt",row.names = F, col.names = F, quote = F,sep="\t")

asum=colSums(abs(A))
A1=rbind(asum,A)
A2=t(A1)
A2=A2[order(A2[,1]),]
A2=A2[129:255,2:128]

asum=colSums(abs(B))
B1=rbind(asum,B)
B2=t(B1)
B2=B2[order(B2[,1]),]
B2=B2[129:255,2:128]

asum=colSums(abs(C))
C1=rbind(asum,C)
C2=t(C1)
C2=C2[order(C2[,1]),]
C2=C2[127:254,2:129]

lnc2=t(lnc)
mi2=t(mi)
m2=t(m)

sss=function(dataset){
  marg.center <- apply(dataset, 2, median)
  marg.scale  <- apply(dataset, 2, mad)
  X <- scale(dataset, center = marg.center, scale = marg.scale)
  return(X)
}

lnc1=sss(lnc1) 
mi1=sss(mi1) 
m1=sss(m1) 

lnc1=apply(lnc1, 1, as.numeric)
mi1=apply(mi1, 1, as.numeric)
m1=apply(m1, 1, as.numeric)


  projectRes1=lnc1 %*% A2
  projectRes2=mi1 %*% B2
  projectRes3=m1 %*% C2
  
lncRNA=lnc-  projectRes1
miRNA=mi- projectRes2
mRNA=m-projectRes3


n0 <- apply(Data2 == 0, 2, sum)
i0 <- which(n0 > 142)
d2<- Data2[, -i0]

dim(d1)
dim(d2)
dim(d3)
d=cbind(Data1,Data2,Data3)
dim(d)
write.table(d,"input.txt",row.names = F, col.names = F, quote = F,sep="\t")


