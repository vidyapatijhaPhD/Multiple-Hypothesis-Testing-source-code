#install.packages("ADMM")
#install.packages("cauchypca")
#install.packages("ggplot2")
#install.packages("writexl")
library(ADMM)
library(cauchypca)
library(cluster)
library(kernlab)
library('nproc')
library(GLDEX)
library(SPRT)
library(stats)
library(mclust)
library("copula")
library("fCopulae")
library("QRM")
library("moments")
library("VineCopula")
library("psych")
library(bayestestR)
library(Metrics)
library(MLmetrics)
library(ggplot2)
library(writexl)
options(scipen = 999)
Y <- read.table('C:/Users/DELL/Downloads/dataset.csv', header = TRUE, sep=',')
Y <- Y$field1
X1 <- Y[1:12210]
#Denoising the matrix using Total variation regularization
xsignal = c(X1)
output = admm.tv(c(X1))
## visualize
opar <- par(no.readonly=TRUE)
plot(xsignal, type="l", main="Total Variance Regularization")
lines(1:12210, output$x, col="red", lwd=2)
par(opar)


X2 <- Y[11701:12000]
#Denoising the matrix using Total variation minimization
xsignal = c(X2)
output = admm.tv(c(X2))
## visualize
opar <- par(no.readonly=TRUE)
plot(xsignal, type="l", main="Total Variance Regularization")
lines(1:300, output$x, col="red", lwd=2)
par(opar)
#Collect both sensor data at Node
X=c(X1, X2)

##########################################Copula and Bayes at Fusion Center#####
val.ln <- as.matrix(X, nrow = 200, ncol = 3)
udata <- prop.table(val.ln)

ro=0.91
ndata <- matrix(udata, nrow=300, ncol=2)
norm.cop <- normalCopula (param = ro, dim=2)
norm.cop
NormCopEst <- fitCopula(norm.cop, ndata, method = "mpl")
NormCopEst
logLik (NormCopEst)
AIC(NormCopEst)
BIC(NormCopEst)



tdata <- matrix(udata, nrow=300, ncol=2)
tCop <- tCopula(param = 0.5, dim = 2, df = 2)
tCop
TCopEst <- fitCopula(tCop, tdata, method = "mpl", estimate.variance = TRUE)
TCopEst
logLik(TCopEst)
AIC(TCopEst)
BIC(TCopEst)


gdata <- matrix(udata, nrow=300, ncol=2)
gumbel <- gumbelCopula(dim = 2, param = 1.1)
gumbCopEst <- fitCopula(gumb.cop0, gdata, method = "mpl")
gumbCopEst
logLik(gumbCopEst)
AIC(gumbCopEst)
BIC(gumbCopEst)

cdata <- matrix(udata, nrow=300, ncol=2)
Clay.cop0 <- claytonCopula(dim = 2, param = ro)
ClayCopEst <- fitCopula(Clay.cop0, cdata, method = "mpl")
ClayCopEst
logLik(ClayCopEst)
AIC(ClayCopEst)
BIC(ClayCopEst)

u <- rCopula(300, gumbel)

if(require(scatterplot3d))
  scatterplot3d(u)

u <- c(rCopula(300, gumbel))

#################Probabilistic Clustering######################################
P_cluster<-Mclust(u)
summary(P_cluster)
plot(P_cluster)
results = data.frame(u,cluster=P_cluster$classification)
clusplot(as.matrix(u), results$cluster, shade = FALSE,labels=5,col.clus="blue",col.p="red",span=FALSE,main="Probabilistic Cluster Mapping",cex=1.2)

v <- results$cluster
c1=0
c2=0
c3=0
c4=0
c5=0
c6=0
for(i in 1:600)
{
  
  if (v[i] == 1)
  {
    c1[i] <- X[i]
    
  }
  if (v[i] == 2)
  {
    c2[i] <- X[i]
    
  }
  if (v[i] == 3)
  {
    c3[i] <- X[i]
    
  }
  if (v[i] == 4)
  {
    c4[i] <- X[i]
    
  }
  if (v[i] == 5)
  {
    c5[i] <- X[i]
    
  }
  if (v[i] == 6)
  {
    c6[i] <- X[i]
    
  }
  
}


c1 <- as.vector(na.omit(c1))
c2 <- as.vector(na.omit(c2))
c3 <- as.vector(na.omit(c3))
c4 <- as.vector(na.omit(c4))
c5 <- as.vector(na.omit(c5))
c6 <- as.vector(na.omit(c6))

c1=fun.zero.omit(c1)
c2=fun.zero.omit(c2)
c3=fun.zero.omit(c3)
c4=fun.zero.omit(c4)
c5=fun.zero.omit(c5)
c6=fun.zero.omit(c6)

if(length(c1)>1)
{
    c1 = matrix(c1)
    ## try different regularization values
    out1 = admm.rpca(c1, lambda=1)
    s1 = out1$L
    r=pchisq(s1[,1], length(c1)-1, lower.tail=FALSE)#gives p value
}
  #out2 = admm.rpca(c1, lambda=0.1)
  #out3 = admm.rpca(c1, lambda=1)
  #opar <- par(no.readonly=TRUE)
  #par(mfrow=c(1,3))
  #image(out1$S, main="lambda=0.01")
  #image(out2$S, main="lambda=0.1")
  #image(out3$S, main="lambda=1")
  #par(opar)

if(length(c2)>1)
{
  c2 = matrix(c2)
  ## try different regularization values
  out2 = admm.rpca(c2, lambda=1)
  s2 = out2$L
  s=pchisq(s2[,1], length(c2)-1, lower.tail=FALSE)#gives p value
}

if(length(c3)>1)
{
    c3 = matrix(c3)
    ## try different regularization values
    out3 = admm.rpca(c3, lambda=1)
    s3 = out3$L
    t=pchisq(s3[,1], length(c3)-1, lower.tail=FALSE)#gives p value
}

  if(length(c4)>1)
  {
    c4 = matrix(c4)
    ## try different regularization values
    out4 = admm.rpca(c4, lambda=1)
    s4 = out4$L
    u=pchisq(s4[,1], length(c4)-1, lower.tail=FALSE)#gives p value
  }
  if(length(c5)>1)
  {
    c5 = matrix(c5)
    ## try different regularization values
    out5 = admm.rpca(c5, lambda=1)
    s5 = out5$L
    v=pchisq(s5[,1], length(c5)-1, lower.tail=FALSE)#gives p value
  }
  if(length(c6)>1)
  {
    c6 = matrix(c6)
    ## try different regularization values
    out6 = admm.rpca(c6, lambda=1)
    s6 = out6$L
    w=pchisq(s6[,1], length(c6)-1, lower.tail=FALSE)#gives p value
  }
###########Sparsity Visualization of robust PCA of different cluster at different regularization parameter lambda###########
#out2 = admm.rpca(c1, lambda=0.1)
#out3 = admm.rpca(c1, lambda=1)
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,3))
image(out1$S, main="lambda=1")
image(out2$S, main="lambda=1")
image(out3$S, main="lambda=1")
image(out4$S, main="lambda=1")
image(out5$S, main="lambda=1")
image(out6$S, main="lambda=1")
#image(out2$S, main="lambda=0.1")
#image(out3$S, main="lambda=1")
par(opar) 
out1$history
out2$history
out3$history
out4$history
out5$history
out6$history

######################Adjusting or correcting p-value through Holm ###
p1=p.adjust(r, method = "holm", n = length(r))
p2=p.adjust(s, method = "holm", n = length(s))
p3=p.adjust(t, method = "holm", n = length(t))
p4=p.adjust(u, method = "holm", n = length(u))
p5=p.adjust(v, method = "holm", n = length(v))
p6=p.adjust(w, method = "holm", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p1,p2,p3,p4,p5,p6)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}
  
T1=c(p1,p2,p3,p4,p5,p6)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T1))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(Holm)/file20.xlsx")


######################Adjusting or correcting p-value through hochberg ###
p7=p.adjust(r, method = "hochberg", n = length(r))
p8=p.adjust(s, method = "hochberg", n = length(s))
p9=p.adjust(t, method = "hochberg", n = length(t))
p10=p.adjust(u, method = "hochberg", n = length(u))
p11=p.adjust(v, method = "hochberg", n = length(v))
p12=p.adjust(w, method = "hochberg", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p7,p8,p9,p10,p11,p12)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}

T2=c(p7,p8,p9,p10,p11,p12)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T2))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(hochberg)/file20.xlsx")

######################Adjusting or correcting p-value through hommel ###
p13=p.adjust(r, method = "hommel", n = length(r))
p14=p.adjust(s, method = "hommel", n = length(s))
p15=p.adjust(t, method = "hommel", n = length(t))
p16=p.adjust(u, method = "hommel", n = length(u))
p17=p.adjust(v, method = "hommel", n = length(v))
p18=p.adjust(w, method = "hommel", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p13,p14,p15,p16,p17,p18)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}

T3=c(p13,p14,p15,p16,p17,p18)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T3))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(hommel)/file20.xlsx")

######################Adjusting or correcting p-value through bonferroni ###
p19=p.adjust(r, method = "bonferroni", n = length(r))
p20=p.adjust(s, method = "bonferroni", n = length(s))
p21=p.adjust(t, method = "bonferroni", n = length(t))
p22=p.adjust(u, method = "bonferroni", n = length(u))
p23=p.adjust(v, method = "bonferroni", n = length(v))
p24=p.adjust(w, method = "bonferroni", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p19,p20,p21,p22,p23,p24)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}

T4=c(p19,p20,p21,p22,p23,p24)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T4))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(bonferroni)/file20.xlsx")

######################Adjusting or correcting p-value through BH ###
p25=p.adjust(r, method = "BH", n = length(r))
p26=p.adjust(s, method = "BH", n = length(s))
p27=p.adjust(t, method = "BH", n = length(t))
p28=p.adjust(u, method = "BH", n = length(u))
p29=p.adjust(v, method = "BH", n = length(v))
p30=p.adjust(w, method = "BH", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p25,p26,p27,p28,p29,p30)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}

T5=c(p25,p26,p27,p28,p29,p30)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T5))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(BH)/file20.xlsx")

######################Adjusting or correcting p-value through BY ###
p31=p.adjust(r, method = "BY", n = length(r))
p32=p.adjust(s, method = "BY", n = length(s))
p33=p.adjust(t, method = "BY", n = length(t))
p34=p.adjust(u, method = "BY", n = length(u))
p35=p.adjust(v, method = "BY", n = length(v))
p36=p.adjust(w, method = "BY", n = length(w))
##############Hypothesis Testing#######################################
minimum_of_p=min(p31,p32,p33,p34,p35,p36)
if(minimum_of_p<=0.05/(length(s1[,1])+length(s2[,1])+length(s3[,1])+length(s4[,1])+length(s5[,1])+length(s6[,1])))
{
  print("Alternate Hypothesis is true")
}else { print("Null Hypothesis is true")}

T6=c(p31,p32,p33,p34,p35,p36)
#plot(sort(rank(T)), type="l", col="blue" , lwd = 0.1)

write_xlsx(data.frame(rank(sort(T6))), "C:/Users/DELL/Desktop/New paper/New Paper 10/paper 10.3 coding(lambda=1)/pgraph(BY)/file20.xlsx")














