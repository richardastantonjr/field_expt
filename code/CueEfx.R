##Analyze changes in community structure before and after 
##addition of predatory bird cues compared to an untreated and procedural controls
##on the same property

## pull in veg and bird survey data; nesting data
PreSurveys <- read.csv("C:/Users/atavistic/Dropbox/Florida/Notes data summaries/Swaziland 2014/Swaziland bird data_proofed.csv")
ExperimentalSurveys <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/ExperimentalSurveys.csv")
ProcControlSurveys <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/ProcControlSurveys.csv")

raMatrixPre <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/raMatrixPre.csv")
PreSurveys[is.na(PreSurveys)] <- 0
raMatrixPre[is.na(raMatrixPre)] <- 0
colnames(raMatrixPre)[1] <- "GPSPt"
included<-c("24",  "26",  "27",  "33",  "34",  "35",  "36",  "40",  "53",  "55",  
  "58",  "60","61",  "62",  "65",  "67",  "68",  "69", "70", "121", "123",
  "125", "126", "128")
raMatrixPre<- raMatrixPre[rownames(raMatrixPre) %in% included,]

raMatrixPost <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/raMatrixPost.csv")

raMatrixProcCont <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/raMatrixProcCont.csv")
raMatrixProcCont[is.na(raMatrixProcCont)] <- 0

PreVeg<-read.csv("C:/Users/atavistic/Dropbox/Florida/Notes data summaries/Swaziland 2014/Primary excel with metadata/Swaziland veg data20150314.csv")
PreVeg<-PreVeg[,1:42]
PreVeg<- PreVeg[PreVeg$GPSpt %in% included,]
PreVeg[is.na(PreVeg)] <- 0

PostVeg <- read.csv("C:/Users/atavistic/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/Cue addition experiment/CueExpRcode/Analysis csv/veg2016.csv")
PostVeg<- PostVeg[PostVeg$GPSpt %in% included,]
NestSearching <- read.csv("C:/Users/atavistic/Dropbox/Florida/Notes data summaries/Swaziland 2015/NestSearching.csv")

## before and after mean shrub cover for each point
preShrub<-(PreVeg$ShrCover1+PreVeg$ShrCover2+PreVeg$ShrCover3)/3
postShrub<-(PostVeg$SH1+PostVeg$SH2+PostVeg$SH3)/3
ShrRatio<-postShrub/preShrub

preGrass<-(PreVeg$GrassCover1+PreVeg$GrassCover2+PreVeg$GrassCover3)/3
postGrass<-(PostVeg$G1+PostVeg$G2.1+PostVeg$G3.1)/3
GrassRatio<-postGrass/preGrass

preTree<-(PreVeg$TreeCover1+PreVeg$TreeCover2+PreVeg$TreeCover3)/3
postTree<-(PostVeg$T1+PostVeg$T2+PostVeg$T3)/3
TreeRatio<-postTree/preTree

## produce species detection-nondetection by site matrices
paMatrixPre<-raMatrixPre
paMatrixPost<-raMatrixPost
paMatrixProcCont<-raMatrixProcCont
paMatrixPre[raMatrixPre >0]<-1    
paMatrixPost[raMatrixPost >0]<-1
paMatrixProcCont[raMatrixProcCont[,]>0]<-1

## abundance by site, pre and post
PreAbund<-rowSums(raMatrixPre)
PostAbund<-rowSums(raMatrixPost)
ProcContAbund<-rowSums(raMatrixProcCont)
AbundRatio<-PreAbund/PostAbund

##species richness by site, pre and post
Prerich<-rowSums(paMatrixPre)                  
Postrich<-rowSums(paMatrixPost)
ProcContrich<-rowSums(paMatrixProcCont)
richRatio<-Prerich/Postrich

##species evenness and Shannon H, pre and post
library(vegan)
ShannonPre<-diversity(raMatrixPre,index = "shannon")
ShannonPost<-diversity(raMatrixPost,index = "shannon")
ShannonProcCont<-diversity(raMatrixProcCont,index="shannon")
ShannonRatio<-ShannonPre/ShannonPost

## Import functional and life-history traits for each species, compute metrics

## manipulate matrices to permit pre-post compositional dissimilarity
(pre_NMDS<-metaMDS(raMatrixPre, distance="jaccard",k=2,trymax=100))
stressplot(pre_NMDS)

par(mfrow=c(1,2))
treat=c(rep("Unthinned",5),rep("Thinned",4))
ordiplot(pre_NMDS,type="n",xlim=c(-0.85,0.9),ylim=c(-0.5,0.8),main="Before")
ordihull(pre_NMDS,groups=treat,draw="polygon",col="grey90",label=T,cex=0.5)

(post_NMDS<-metaMDS(raMatrixPost,distance = "jaccard", k=2, trymax=100))
stressplot(post_NMDS)

(pc_NMDS<-metaMDS(raMatrixProcCont,distance = "jaccard", k=2, trymax=100))
stressplot(pc_NMDS)

treat=c(rep("Control",12),rep("Treatment",12))
## Use ellipses with 95% conf.
par(mfrow=c(1,3))
ordiplot(pre_NMDS,type="n",xlim=c(-1.5,1.5),ylim=c(-0.5,0.8),main="Before cue addition")
ordiellipse(pre_NMDS,groups=treat,draw="polygon",col=c("grey90","blue"),
            label=F,cex=0.5,conf=0.95)
points(pre_NMDS, display = "sites",pch=c(rep(16,5),rep(1,4)))

ordiplot(post_NMDS,type="n",xlim=c(-1.2,1),ylim=c(-0.8,0.8),main="After cue addition")
ordiellipse(post_NMDS,groups=treat,draw="polygon",col="grey90",
            label=F,cex=0.5,conf = 0.95)
points(post_NMDS, display = "sites",pch=c(rep(16,12),rep(1,12)))

treat=c(rep("Control",6),rep("Treatment",6))
ordiplot(pc_NMDS,type="n",xlim=c(-1.2,1),ylim=c(-0.8,0.8),main="After procedural control cue addition")
ordiellipse(pc_NMDS,groups=treat,draw="polygon",col="grey90",
            label=F,cex=0.5,conf = 0.95)
points(pc_NMDS, display = "sites",pch=c(rep(16,5),rep(1,4)))

## test if compositional changes are "sig" using randomization test
(pro <- procrustes(pre_NMDS, post_NMDS,symmetric=T))
protest(pre_NMDS, post_NMDS)
plot(pro)
plot(pro, kind = 2)

treated<-c(1,1,1,1,1,0,0,0,0)
(ef <- envfit(post_NMDS~treated, permu = 999))





Trt<-unique(postThinBirds[,1:2])
data=data.frame(cbind(Trt,PreAbund,PostAbund,Postrich,ShannonPre,ShannonPost,preShrub,postShrub))

library(visreg)
##Abundance
m1<-lm(PostAbund~Treated+PreAbund,data=data)
summary(m1)
m2<-lm(AbundRatio~Treated,data=data)
summary(m2)
visreg(m1)
visreg(m2)

##Richness
m3<-lm(Postrich~Treated+Prerich,data=data)
summary(m3)
m4<-lm(richRatio~Treated,data=data)
summary(m4)
visreg(m3)
visreg(m4)

##Shannon
m5<-lm(ShannonPost~Treated+ShannonPre,data=data)
summary(m5)
m6<-lm(ShannonRatio~Treated,data=data)
summary(m6)
visreg(m5)
visreg(m6)

tiff("v.tif", width = 8, height = 10, units = 'in',res=300, compression = "lzw")
par(mfrow=c(3,2))
visreg(m1)
visreg(m3)
visreg(m5)
dev.off()


tiff("thinefx.tif", width = 8, height = 10, units = 'in',res=300, compression = "lzw")
par(mfrow=c(1,3))
boxplot(AbundRatio~Treated,data,main="Abundance", las=1,cex.lab=1.4,cex.axis=1.3,
        xlab="Thinning treatment", ylab="Initial abundance / post-thin abundance")
boxplot(richRatio~Treated,data,main="Species richness", las=1,cex.lab=1.4,cex.axis=1.3,
        xlab="Thinning treatment", ylab="Initial richness / post-thin richness")
boxplot(ShannonRatio~Treated,data,main="Shannon diversity", las=1,cex.lab=1.4,cex.axis=1.3,
        xlab="Thinning treatment", ylab="Initial Shannon / post-thin Shannon")
dev.off()

abunL<-confint(m2,"Treated")[1,1]
abunU<-confint(m2,"Treated")[1,2]
richL<-confint(m4,"Treated")[1,1]
richU<-confint(m4,"Treated")[1,2]
ShaL<-confint(m6,"Treated")[1,1]
ShaU<-confint(m6,"Treated")[1,2]

studies.ci <- c(-0.1700916,abunL,abunU,-0.3968054,richL,richU,-0.3534447,ShaL,ShaU)
studies.ci <- matrix(studies.ci, ncol=3, byrow=T)
authors <- c("Abundance","Richness", "Shannon diversity")
tiff("coefs.tif", width = 12, height = 8, units = 'in',res=300, compression = "lzw")
forest.plot.or(authors= authors,studies.ci=studies.ci,plot.lim=c(-1.5, 1.5),cex=1.25,
               ci.txt="Coef. (95% CI)",ref.vline.at=0,log.scale=F,tick.lim=F,tick.pch=1,
               study.txt="Response", plot.xaxis=T,standard.or.plot=F)
##title('Shrub thinning and community structure', font.main=1)
dev.off()




##---------------
library(car)
outlierTest(m6) # Bonferonni p-value for most extreme obs
qqPlot(m6, main="QQ Plot") #qq plot for studentized resid 
leveragePlots(m1)
vif(m1)
crPlots(m1)

library(corrplot)
f<-cor(data)
corrplot(f,method = "circle")
corrplot(f,method = "number")
