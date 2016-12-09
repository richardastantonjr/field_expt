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
## immediate response to the treatment was of interest. Surveys therefore consisted of 5-6 
## replicate periods between 21 October and 14 December 2015, with two surveys in one day per 
## period. We assumed closure in occupancy/abundance within a period but not between periods. 
## we recorded all species detected within 100m and binned the detection bands of 0-50 and 
## 51-100m. Procedural control surveys were completed between 9 and 12 January 2016.

## Nest searching: RAS searched for nests, focusing on species that appeared to be responding
## to the treatments in terms of abundance/detectability- (based on pivot tables) He searched
## for 6 hrs per point within 50m of each point, divided into 1-4 hr stints allocated 
## throughout the day and across dates to minimize bias in effort across points and treatments
## He located nests using behavioral cues, e.g. following birds carrying food and nest material,
## listening for calls on nests, observing relief provided when birds switched incubation 
## duties, and observing haphazardly-discovered nests from a distance to determine if they were 
## active. We only recorded active nests. At each active nest we recorded the stage at which 
## it was discovered, i.e. empty, eggs, nestlings, the date, and the species.

## The code is intended to be implemented as follows.
## 1. source "data_organization.R" to import and manipulate data. The relevant products are:
##     A. Arrays of detection-nondetection events for pretreatment, treatment, and procedural 
##        control surveys. The pretreatment array is 24 sites * 4 visits * XX species. The treatment 
##        array is 24 sites * 10 surveys * XX species with a separate indicator the surveys are
##        spread over 5 sampling periods with 2 visits/period. The procedural control array is 12 sites * 
##        XX surveys * YY species.
##     B. Data frames of site covariates for the pre-treatment and treatment periods.
##     C. Data frames [?] of sampling covariates for the pretreatment, treatment, and procedural 
##        control surveys.
## packages----------------------------------------------------------
library(tidyverse)
library(vegan)
library(unmarked)
library(jagsUI)
##-------------------------------------------------------------------

## data--------------------------------------------------------------
## bird surveys
pre_surveys<-read.csv("./data/Swaziland bird data_proofed.csv")
## import 371 sites * 5 visits * 209 species array
## built from pre_surveys in a separate
## piece of code- written for analyzing a landscape-scale data set
## This will be subsetted later.
detectHists <- read.table("./data/abunhists20160919", quote="\"", 
                    comment.char="", stringsAsFactors=FALSE)
post_surveys<-read.csv("./data/ExperimentalSurveys.csv")
noise_control_surveys<-read.csv("./data/ProcControlSurveys.csv")

## sampling covariates
pre_sampling_covs<-read.csv("./data/SamplingCovs.csv")
##post_sampling_covs<-read.csv("./")                  ######## need to locate and arrange

## vegetation
veg2015<-read.csv("./data/Swaziland veg data.csv")
veg2016<-read.csv("./data/veg2016.csv")

# nests
nests<-read.csv("./data/Nests.csv")
nest_effort<-read.csv("./data/search_efforts.csv")

## bird foraging traits and masses collected from Hockey et al. 2005, i.e.
## Robert's Birds of Southern Africa, 5th edition.
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
pre_surveys<-pre_surveys[,1:5]  ## drop the NumOut and Total columns so the spatial scale of the 
                                ## bird data and the vegetation measures match
SpeciesNames<-sort(unique(pre_surveys$Species)) 
## drop species not detected within 50m; indices from a separate script
SpeciesNames<-SpeciesNames[c(-1,-7,-8,-22,-60,-62,-88,-99,-117,-125,-129,-136,-139,
                                  -143,-165,-172,-175)]
## subset the surveys to remove additional species not detected at any experimental points
pre_surveys<-pre_surveys[pre_surveys$GPSpt %in% included,]
## identify which species slices of detectHists to keep and drop 5th visit column because
## all but one are NA because not surveyed
detectHists<-detectHists[,1:4,which((SpeciesNames %in% unique(pre_surveys$Species))==TRUE )]
pre_survey_species<-subset(SpeciesNames,SpeciesNames %in% unique(pre_surveys$Species)==TRUE )

## subset pre-treatment vegetation data
veg2015<-veg2015[,1:42]        ## drop the notes column
veg2015<- veg2015[veg2015$GPSpt %in% included,]
veg2015[is.na(veg2015)] <- 0

## subset post-treatment data to exclude locations sampled for a concurrent side project 
## (not in scope of this project)
veg2016<- veg2016[veg2016$GPSpt %in% included,]

## before and after % mean shrub, grass, and tree cover for each point
preShrub<-(veg2015$ShrCover1+veg2015$ShrCover2+veg2015$ShrCover3)/3
postShrub<-(veg2016$S1+veg2016$S2.1+veg2016$S3.1)/3
preGrass<-(veg2015$GrassCover1+veg2015$GrassCover2+veg2015$GrassCover3)/3
postGrass<-(veg2016$G1+veg2016$G2.1+veg2016$G3.1)/3
preTree<-(veg2015$TreeCover1+veg2015$TreeCover2+veg2015$TreeCover3)/3
postTree<-(veg2016$T1+veg2016$T2+veg2016$T3)/3

## subset the trait data
traits<-filter(traits,traits$Species %in% unique(pre_surveys$Species)) 
traits<- traits[,c(1:6,9,19,20,22)]
traits[is.na(traits)] <- 0 

## determine if a species is a predator
predator<-rep(NA, nrow(traits))
for (i in 1:nrow(traits)){
  if((traits[i,8]>0) || (traits[i,9]>0) || (traits[i,10]>0)){
    predator[i]=1
  } else {
    predator[i] =0} 
}
traits$predator<-predator
traits<-traits[,c(1:7,11)]


###-----------------------------------------------------------------------------------
## manipulate post_surveys to get a site*period*visit*species array and a table
## of sampling covariates suitable for occupancy modeling
##------------------------------------------------------------------------------------                                   
##create a list of the dates when each location was sampled
GPSptDateList<- post_surveys %>% 
distinct(GPSpt,Date,Time) %>%                                     
unstack(Date~GPSpt) 

## create a vector of number of visits per point
num_visits<-c(rep(NA,24))
for (i in 1:length(GPSptDateMatrix)){
  num_visits[i]<-length(GPSptDateMatrix[[i]])
}

## convert the list into a data.frame where each cell is a date visited or NA
siteDateMatrix<-do.call(rbind, lapply(GPSptDateList , function(x){ 
  length(x) <- max(num_visits)
  x })) %>% 
  data.frame()


## create empty XX species by 24 sites by 14 visits array to fill with counts as appropriate
post_detections<-array(NA,dim = c(24,14))



## work with sampling covariates for the experimental surveys

##-----------------------------------------------------------------------------------------------------------
## Naive occupancy for the pre-surveys; useful for preliminary analysis and because
## low occurrence can bias occupancy models so we may wish to exclude rare species
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
  
  
  
  