

########----------------------------------------------------------------------
##### preliminary analysis of treatmeant by shrub cover for some spp.
####--------------------------------------------------------------------------
## total detections treat vs. control
Detections_by_treatment<-select(post_surveys, Treated, Species, NumIn) %>% 
  group_by(Treated) %>% 
  summarise(total_detections = sum(NumIn, na.rm= TRUE))

## total detections for each shrub cover amount
ct_by_shr<-select(post_surveys, Treated, NumIn, GPSpt, meanShrub) %>%
  group_by(GPSpt,meanShrub, Treated) %>% 
  summarise(counts = sum(NumIn, na.rm=TRUE)) %>% 
  arrange(Treated, GPSpt,counts)

## no indication total counts vary 
plot(counts~meanShrub, data=ct_by_shr[1:12,])
points(counts~meanShrub, data=ct_by_shr[13:24,],col='blue', pch=4)

summary(lm(counts~Treated*meanShrub, data = ct_by_shr))
summary(lm(counts~Treated+meanShrub, data = ct_by_shr))
summary(lm(counts~Treated, data = ct_by_shr))
summary(lm(counts~meanShrub, data = ct_by_shr))

## data frame of species counts grouped by treatment/control
veg2015<-mutate(veg2015, meanShrub=(ShrCover1+ShrCover2+ShrCover3)/3)
shrub_covers<-veg2015[,c(1,43)]
post_surveys<-inner_join(post_surveys,shrub_covers)

Detection_sp_treat<-select(post_surveys, Species, Treated, NumIn) %>% 
  group_by(Treated, Species) %>% 
  summarise( total_ct = sum(NumIn, na.rm= TRUE)) %>% 
  ## filter uncommon species
  filter(total_ct > (sum(Detections_by_treatment[,2])/500)) %>% 
  as.data.frame() %>% 
  droplevels() %>% 
  arrange(Species,Treated)



             
ggplot(data = Detection_sp_treat) + 
  geom_boxplot(mapping = aes(x = Treated, y = total_ct, shape = Treated))+
  facet_wrap(~ Species)
  








### incomplete; 
Detection_spp_treat_shr<-select(post_surveys,GPSpt, Species, Treated, NumIn) %>% 
  group_by(GPSpt, Species) %>% 
  summarise( total_ct = sum(NumIn, na.rm= TRUE)) %>% 
  ## filter uncommon species
  filter(total_ct> (sum(Detections_by_treatment[,2])/500)) %>% 
  as.data.frame() %>% 
  ## droplevels(Detection_species_treat) %>% 
  arrange(GPSpt, Species)









