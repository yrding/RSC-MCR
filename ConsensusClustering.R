library("ConsensusClusterPlus")
library(dplyr)
library(fpc)
library(survival)
library(survminer)
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

datalist=list(lnc,mi,m);
datalist=list(lnc1,mi1,m1);
data=datalist;
data=list(dff)
data=list(ICA$S)
K=2;
d <- do.call(cbind, data) %>% t
d = sweep(d,1, apply(d,1,median,na.rm=T))
if(K==2){
    Kmax <- 3
}else{
    Kmax <- K
}
fit <-  ConsensusClusterPlus(d,maxK=Kmax,plot = "png")
clust <- fit[[2]]$consensusClass
res5 <- list(clust=clust, fit=fit)
b=res5[["fit"]][[2]][["consensusMatrix"]]
label=res5[["fit"]][[2]][["consensusClass"]]
stats <- cluster.stats(b, label)
stats$avg.silwidth
calDBI(b,label)



#生存分析
survivalData = read.table("ConsensusClustering_survial.txt",sep="\t",header=TRUE)
ClustGroup =res5[["clust"]]
SurvDat = survivalData

colors <- c("red", "blue", "orange", "pink","green")
coxFit <- coxph(
  Surv(time = Survival, Death) ~ as.factor(ClustGroup),
  data = SurvDat,
  ties = "exact"
)

mfit <- survfit(Surv(Survival, Death) ~ as.factor(ClustGroup), data = SurvDat)
plot(
  mfit, col = colors , main = "Survival curves for KIRC, level 2",cex.main="0.7",
  xlab = "Days", ylab = "Survival",lwd = 2,pch=c(50,30),
)
legend("bottomright",
       text.font=0.1,
       text.width = 2650,
       legend = paste(
         "Cox p-value:", round(summary(coxFit)$sctest[3], digits = 5), sep = ""
       )
)
text.legend = paste("Group ", levels(factor(ClustGroup)), ": ",
                    table(ClustGroup)[levels(factor(ClustGroup))], sep ="")
legend(
  "topright",
  fill = colors,
  ncol=2,
  text.font=0.1,
  text.width =1600,
  legend =text.legend,
  )

