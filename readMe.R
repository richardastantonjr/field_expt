
## TLDR: source the R script named "data.organization.R" to see what was done. Details of why and how are below
## and as comments in the script.

## Code to prepare data for analysis of a field experiment conducted on the Mbuluzi game
## Reserve in Swaziland, 2014-2016 by Stanton and Subiya. Also preliminary analysis and
## data visualization.

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

## I bit off more than I could chew but the content should demonstrate familiarity with tidyr, for loops,
## version control and functions. 
