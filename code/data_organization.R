## Code to prepare data for analysis of a field experiment conducted on the Mbuluzi game
## Reserve in Swaziland, 2014-2016 by Stanton and Subiya. Also preliminary analysis and
## data visulaization.

## We added auditory cues of predator birds across a shrub cover gradient and compared it to
## controls matched for shrub cover. Procedural controls included dummy speakers visited 
## and manipulated the same as broadcast stations, and a subsequent pink noise treatment 
## applied to 1/2 of the control plots, also with a similar range of shrub cover.

## We used repeated 10 min point counts and nest searching to quantify species and community-
##level responses to the treatment. Specifically, we visited each point 4 times between
## Dec 2014 and March 2015 to collect baseline bird occurrence data across an entire rainy 
## season as part of a larger study. We then surveyed birds during the experiment in a manner 
## amenable to analysis using a robust design allowing for temporary immigration and emigration
## of species, in part because we were interested in temporal dynamics but also because the
## immediate response to the treatment was of interest. Surveys therefore consisted of 6 
## replicate periods between 21 October and 14 December 2015, with two surveys in one day per 
## period. We assumed closure in occupancy/abundance within a period but not between periods. 
## we recorded all species detected within 100m and binned the detection bands of 0-50 and 
## 51-100m. Procedural control surveys were completed between 9 and 12 January 2016.

## Nest searching: RAS searched for nests, focusing on species that appeared to be responding
## to the treatments in terms of abundance/detectability- (based on pivottables) He searched
## for 6 hrs per point within 50m of each point, divided into 1-4 hr stints allocated 
## throughout the day and across dates to minimize bias in effort across points and treatments
## He located nests using behavioral cues, e.g. following birds carrying food and nest material,
## listening for calls on nests, observing relief provided when birds switched incubation 
## duties, and observing nests located haphazardly from a distance to determine if they were 
## active. We only recorded active nests. At each active nest we recorded the stage at which 
## it was discovered, i.e. empty, eggs, nestlings, the date, and the species.

## packages----------------------------------------------------------
library(tidyverse)
library(vegan)
library(unmarked)
library(jagsUI)
##-------------------------------------------------------------------

## data--------------------------------------------------------------
## bird surveys
pre_surveys<-read.csv("./data/Swaziland bird data_proofed.csv")
## import 371 sites *5 visits *209 species array
## built from pre_surveys in a separate
## piece of code- written for analyzing the full data set
detectHists <- read.table("./data/abunhists20160919", quote="\"", 
                    comment.char="", stringsAsFactors=FALSE)
post_surveys<-read.csv("./data/ExperimentalSurveys.csv")
noise_control_surveys<-read.csv("./data/ProcControlSurveys.csv")
pre_sampling_covs<-read.csv("./data/SamplingCovs.csv")
##post_sampling_covs<-read.csv("./")                  ######## need to locate and arrange

## vegetation
veg2015<-read.csv("./data/Swaziland veg data.csv")
veg2016<-read.csv("./data/veg2016.csv")

# nests
nests<-read.csv("./data/Nests.csv")
nest_effort<-read.csv("./data/search_efforts.csv")

# traits
traits<-read.csv("./data/TraitData.csv")
##----------------------------------------------------------------------------------
##----------------------------------------------------------------------------------

## SUBSET AND ORGANIZE DATA
## GPSpt,i.e. waypoint numbers,to keep from pre_surveys and detectHists
included<-c("24",  "26",  "27",  "33",  "34",  "35",  "36",  "40",  "53",  "55",  
            "58",  "60","61",  "62",  "65",  "67",  "68",  "69", "70", "121", "123",
            "125", "126", "128")

## subset pre-treatment bird surveys
## put vector back into matrix form and subset
detectHists<-as.vector(detectHists[,1])
dim(detectHists)<-c(371,5,209)
points<-sort(unique(pre_surveys$GPSpt))
## match GPspts to positions in the 1st dimension of the array and drop points not revisited
## during the field experiment
detectHists<-detectHists[which((points) %in% included),,]  
## get unique species names for this data set
pre_surveys<-pre_surveys[,1:5]  ## drop the NumOut and Total columns
SpeciesNames<-sort(unique(pre_surveys$Species)) 
## drop species not detected within 50m; indices from a separate script
SpeciesNames<-SpeciesNames[c(-1,-7,-8,-22,-60,-62,-88,-99,-117,-125,-129,-136,-139,
                                  -143,-165,-172,-175)]
## subset the surveys to remove additional species not detected at any experimental points
pre_surveys<-pre_surveys[pre_surveys$GPSpt %in% included,]
## identify which species slices of detectHists to keep and drop 5th visit column
## all of one of which are NA because not surveyed
detectHists<-detectHists[,1:4,which((SpeciesNames %in% unique(pre_surveys$Species))==TRUE )]
pre_survey_species<-subset(SpeciesNames,SpeciesNames %in% unique(pre_surveys$Species)==TRUE )
#names(detectHists)<-c(rep(NA,24),rep(NA,4),pre_survey_species)

## subset pre-treatment vegetation data
veg2015<-veg2015[,1:42]        ## drop the notes column
veg2015<- veg2015[veg2015$GPSpt %in% included,]
veg2015[is.na(veg2015)] <- 0

## remove BACI study points from the 2016 veg, these were for a concurrent side project
veg2016<- veg2016[veg2016$GPSpt %in% included,]

## before and after % mean shrub, grass, and tree cover for each point
preShrub<-(veg2015$ShrCover1+veg2015$ShrCover2+veg2015$ShrCover3)/3
postShrub<-(veg2016$S1+veg2016$S2.1+veg2016$S3.1)/3
preGrass<-(veg2015$GrassCover1+veg2015$GrassCover2+veg2015$GrassCover3)/3
postGrass<-(veg2016$G1+veg2016$G2.1+veg2016$G3.1)/3
preTree<-(veg2015$TreeCover1+veg2015$TreeCover2+veg2015$TreeCover3)/3
postTree<-(veg2016$T1+veg2016$T2+veg2016$T3)/3


###-----------------------------------------------------------------------------------
## manipulate post_surveys to get a site*period*visit*species array and a table
## of sampling covariates suitable for occupancy modeling
##------------------------------------------------------------------------------------                                   
##create a list of the dates when each location was sampled
GPSptDateMatrix<- post_surveys %>% 
distinct(GPSpt,Date,Time) %>%                                     
unstack(Date~GPSpt) 

## create a vector of number of visits per point
num_visits<-c(rep(NA,24))
for (i in 1:length(GPSptDateMatrix)){
  num_visits[i]<-length(GPSptDateMatrix[[i]])
}
 
## convert from list to data frame 
## this is not doing what I expected at all because the list elements do not seem to be in order
## and the uneven sampling is not captured.
dates_visited_by_site<-data.frame()
for (i in 1:length(GPSptDateMatrix)){
  dates_visited_by_site<-rbind(dates_visited_by_site,GPSptDateMatrix[[i]])
}

rowMax <- max(sapply(GPSptDateMatrix, length))         ## rowMax=5, maximum number of visits to a site in the dataset

SiteVisitDateMatrix<-do.call(rbind, lapply(GPSptDateMatrix, function(x){ 
  length(x) <- rowMax 
  x })) 
SiteVisitDateMatrix<-data.frame(SiteVisitDateMatrix)
SiteVisitDateMatrix$GPSpt<-rownames(GPSptDateMatrix)  

bulbul_surveys<-subset(post_surveys, species="Dark-capped Bulbul")
siteByDate_bulbul <- xtabs(NumIn ~ GPSpt+Date, data=bulbul_surveys) 



for (i in 1:length(GPSptDateMatrix)){
  siteByDate_bulbul[i,]<-subset(siteByDate_bulbul[i,], 
       as.integer(names(siteByDate_bulbul[i,])) %in% as.integer(unlist(GPSptDateMatrix[i]))==TRUE)                                             
}


## create a matrix of values indicating if a GPSpt and Date combination is valid
Junk<-XSurveys
for (i in 1:length(unique(Surveys$Species))){
  for (j in 1:length(unique(Surveys$GPSpt))){
    Junk[j,,i]<-colnames(XSurveys[,,i]) %in% SiteVisitDateMatrix[j,1:5]}}






## work with sampling covariates for the experimental surveys

##-----------------------------------------------------------
## subset the trait data
##-----------------------------------------------------------
traits<-filter(traits,traits$Species %in% unique(pre_surveys$Species)) 
  traits<- traits[,c(1:6,9,19,20,22)]
  traits[is.na(traits)] <- 0 
  
  predator<-rep(NA, nrow(traits))
  for (i in 1:nrow(traits)){
  if((traits[i,8]>0) || (traits[i,9]>0) || (traits[i,10]>0)){
    predator[i]=1
    } else {
    predator[i] =0} 
  }
  traits$predator<-predator
  traits<-traits[,c(1:7,11)]

##-----------------------------------------------------------------------------------------------------------
## Naive occupancy for the pre-surveys
##-----------------------------------------------------------------------------------------------------------
  naive_psi_pre<-data.frame()
  for (i in 1:length(pre_survey_species)){
    tryCatch({
      y<- as.data.frame(detectHists[,,i]) ## create an encounter history for species i
      ## calculate naive occupancy for each species
      naive_psi_pre[i,1]<-pre_survey_species[i]
      naive_psi_pre[i,2]<-sum(apply(y, 1, max, na.rm=TRUE))/24
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  names(naive_psi_pre)<-c("species","naive_psi")
  
detectHists.subset10<-detectHists[,,which(naive_psi_pre$naive_psi>=0.1)]
detectHists.subset20<-detectHists[,,which(naive_psi_pre$naive_psi>=0.2)]
nspec.subset10<-dim(detectHists.subset10)[3]
nspec.subset20<-dim(detectHists.subset20)[3]
  
pre_survey_species10<-subset(naive_psi_pre,naive_psi>=0.1)
pre_survey_species20<-subset(naive_psi_pre,naive_psi>=0.2)
  
##----------------------------------------------------------------------------------------------
## Naive occupancy for the post surveys
##----------------------------------------------------------------------------------------------
naive_psi_post<-data.frame()
for (i in 1:length(post_survey_species)){
  tryCatch({
    y<- as.data.frame(detectHists_post[,,i]) ## create an encounter history for species i
    ## calculate naive occupancy for each species
    naive_psi_post[i,1]<-post_survey_species[i]
    naive_psi_post[i,2]<-sum(apply(y, 1, max, na.rm=TRUE))/24
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
names(naive_psi_post)<-c("species","naive_psi")

detectHists_post.subset10<-detectHists[,,which(naive_psi_post$naive_psi>=0.1)]
detectHists.subset20_post<-detectHists[,,which(naive_psi_post$naive_psi>=0.2)]
nspec_post.subset10<-dim(detectHists.subset10)[3]
nspec_post.subset20<-dim(detectHists.subset20)[3]

post_survey_species10<-subset(naive_psi_post,naive_psi>=0.1)
post_survey_species20<-subset(naive_psi_post,naive_psi>=0.2)
  
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
  
  
  
  