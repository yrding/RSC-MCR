library(LRAcluster)
library(dplyr)

library(survival)
library(survminer)

getwd
getwd()
setwd("/Users/mac/Downloads/text/膀胱癌")

lnc=read.table("lnc_trans.txt");
mi=read.table("mi_trans.txt");
m=read.table("m_trans.txt");

lnc1=read.table("result_lnc.txt");
mi1=read.table("result_mi.txt");
m1=read.table("result_m.txt");

nr1=nrow(Data1);
nc1=ncol(Data1);
Data1=Data1[2:nr1,2:nc1];

nr2=nrow(Data2);
nc2=ncol(Data2);
Data2=Data2[2:nr2,2:nc2];

nr3=nrow(Data3);
nc3=ncol(Data3);
Data3=Data3[2:nr3,2:nc3];

Data1=t(Data1);
Data2=t(Data2);
Data3=t(Data3);

lnc=apply(lnc,2,as.numeric);
mi=apply(mi,2,as.numeric);
m=apply(m,2,as.numeric);

data=read.table("yuzhishanchu.txt");

datalist=list(as.matrix(lnc),as.matrix(mi),as.matrix(m))
datalist1=list(as.matrix(lnc1),as.matrix(mi1),as.matrix(m1))

datalist=list(as.matrix(lnc_dd),as.matrix(mi_dd),as.matrix(m_dd))

data=datalist;

lnc_dd=t(lnc_delet)
mi_dd=t(lnc_delet)
m_dd=t(m_delet)


#dat <- lapply(data, t);
K=2;
type= rep("gaussian", length(data));
fit <- LRAcluster(data, dimension = K, type= type);
hcl <-fit$coordinate %>% t %>% dist %>% hclust(method="ward.D2")
clust <-cutree(hcl,k=K);
res6 <- list(clust=clust, fit=fit);

labels=res6[["clust"]]
b=t(res6[["fit"]][["coordinate"]])
library(fpc)
stats=cluster.stats(dist(b),labels)
stats[["avg.silwidth"]]
dbi <- calDBI(b, labels)

stats1=cluster.stats(dist(b1),labels1)
stats1[["avg.silwidth"]]
calDBI <- function(x=data,labels=labesls)
  
  ##data必须行为样本，列为特征
{
  clusters_n <- length(unique(labels))
  cluster_k <- list()
  for (i in c(1:clusters_n)) {
    cluster_k[[i]] <- x[which(labels==i),]
  }
  
  centroids <- list()
  for (i in c(1:clusters_n)) {
    centroids[[i]] <- apply(cluster_k[[i]],2,mean)
  }
  
  s <- list()
  for (i in c(1:clusters_n)) {
    a <- c()
    for (j in c(1:nrow(cluster_k[[i]]))) {
      b <- dist(rbind(cluster_k[[i]][j,],centroids[[i]]),method = "euclidean")
      a <- c(a,b)
    }
    s[[i]] <- mean(a)
  }
  
  Ri <- list()
  for (i in c(1:clusters_n)){
    r <- c()
    for (j in c(1:clusters_n)){
      if (j!=i){
        h <- (s[[i]]+s[[j]])/dist(rbind(centroids[[i]],centroids[[j]]),method = "euclidean")
        r <- c(r,h)
      }
    }
    Ri[[i]] <- max(r)
  }
  dbi <- mean(unlist(Ri))
  return(dbi)
}

dbi1 <- calDBI(b1, labels1)#x为样本——特征矩阵（行为样本，列为特征），labels为聚类结果
dbi <- calDBI(b, labels)

calCH <- function(X,labels){ 
  ##X必须行为样本，列为特征
  labels_n <- length(unique(labels))
  samples_n <- nrow(X)
  X_mean <- apply(X,2,mean)
  ex_disp <- c()
  in_disp <- c()
  for (i in c(1:labels_n)) {
    cluster_k <- X[which(labels==i),]
    mean_k <- apply(cluster_k,2,mean)
    a1 <- nrow(cluster_k)*sum((mean_k-X_mean)^2)
    ex_disp <- c(ex_disp,a1)
    a2 <- sum((t(t(cluster_k)-mean_k))^2)
    in_disp <- c(in_disp,a2)
  }
  k1<- sum(ex_disp)
  k2<- sum(in_disp)
  if(k2==0)
  {
    return(1)
  }
  else
  {
    return((k1*(samples_n-labels_n))/(k2*(labels_n-1)))
  }
}
#sample
ch<- calCH(b,labels)#X为样本——特征矩阵（行为样本，列为特征），labels为聚类结果
ch1<- calCH(b1,labels1)

#生存分析
survivalData = read.table("LRAcluster_survial.txt",sep="\t",header=TRUE)
ClustGroup =res3[["clust"]]
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


