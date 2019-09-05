#load libraries
#install.packages("ggpubr")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(nlme)
library(dplyr)
library(lme4)
library(lmerTest)

# MOLT SIZES
#set up working directory
setwd("C:/Users/emma.reinhardt/Desktop/Research/Size/Data")
#read in file
allMolts <- read.csv("molts.csv", stringsAsFactors = TRUE)
#get rid of empty rows/columns and compile all width columns
#renamed columns .csv file to make this easier
mw <- subset(allMolts, Group != "")
mw <- subset(mw, select = c(Crab, Group, J1, J2, J3, J4, J5, J6))
#turn dataset into long/skinny
mw <- gather(mw, key = Stage, value = moltWidth, J1, J2, J3, J4, J5, J6)
#turn characters into factors
mw$moltWidth <- as.numeric(mw$moltWidth)
#adg in treatments
mw$Treatment <- ""
mw$Treatment[mw$Group == "HA" | mw$Group == "HB" | mw$Group == "HC"] <- "High"
mw$Treatment[mw$Group == "LA" | mw$Group == "LB" | mw$Group == "LC"] <- "Low"

#repeat all above steps for a separate lengths dataframe (fixing the stage names for merging)
ml <- subset(allMolts, Group != "")
ml <- subset(ml, select = c(Crab, Group, J1_L, J2_L, J3_L, J4_L, J5_L, J6_L))
ml <- gather(ml, key = StageL, value = moltLength, J1_L, J2_L, J3_L, J4_L, J5_L, J6_L)
cleanUp <- function(dd){
  dd$Stage = ""
  dd$Stage[dd$StageL == "J1_L"] <- "J1"
  dd$Stage[dd$StageL == "J2_L"] <- "J2"
  dd$Stage[dd$StageL == "J3_L"] <- "J3"
  dd$Stage[dd$StageL == "J4_L"] <- "J4"
  dd$Stage[dd$StageL == "J5_L"] <- "J5"
  dd$Stage[dd$StageL == "J6_L"] <- "J6"
  dd$Stage[dd$StageL == "J7_L"] <- "J7"
  dd <- subset(dd, select = c(Crab, Group, Stage, moltLength))
  return(dd)
}

ml <- cleanUp(ml)

ml$moltLength <- as.numeric(ml$moltLength)
ml$Treatment <- ""
ml$Treatment[ml$Group == "HA" | ml$Group == "HB" | ml$Group == "HC"] <- "High"
ml$Treatment[ml$Group == "LA" | ml$Group == "LB" | ml$Group == "LC"] <- "Low"

#merge the consolidated width and length dataframes
#THIS IS THE SOLE MOLT DF
molt <- merge(ml, mw)


#BODY SIZES
#do the same thing to prep all of the size data
d <- read.csv("sizes.csv", stringsAsFactors = FALSE)
#fix data frame
d <- subset(d, Group != "")
d <- subset(d, select = c(Crab, Group, J1, J2, J3, J4, J5, J6))
#short-wide to long-skinny
dg <- gather(d, key = Stage, value = carapaceWidth, J1, J2, J3, J4, J5, J6)
#convert variable from character to number
dg$carapaceWidth <- as.numeric(dg$carapaceWidth)
#adg Treatment
dg$Treatment <- ""
dg$Treatment[dg$Group == "HA" | dg$Group == "HB" | dg$Group == "HC"] <- "High" 
dg$Treatment[dg$Group == "LA" | dg$Group == "LB" | dg$Group == "LC"] <- "Low"


#combine width data frames
#THIS IS MOLT + LIVE CRAB WIDTHS
#total <- data.frame(mw$Crab, mw$Group, mw$Stage, mw$Treatment, mw$moltWidth, dg$carapaceWidth)
total <- merge(mw, dg)



### MOLT + BODY PLOTS
#pretty title
legendTitle = expression("CO"[2]*" Treatment")
#linear fit of molts vs bodies
ggplot(total, aes(x = moltWidth, y = carapaceWidth)) + 
  geom_point(aes(color = Stage, shape = Treatment)) + 
  theme_bw() + 
  geom_smooth(aes(fill = Treatment), method = lm, se = TRUE) + 
  labs(x = "Molt Width (mm)", y = "Carapace Width (mm)", title = "Molt vs. Carapace Width Across Stages", 
       legend, fill = legendTitle, color = "Stage", shape = legendTitle)
#linear fit of molts vs bodies across stages
ggplot(total, aes(moltWidth, carapaceWidth)) + 
  geom_point(aes(color = Stage, shape = Treatment)) +
  geom_smooth(aes(fill = Treatment), method = lm, se=TRUE) +
  theme_bw() +
  facet_wrap(~Stage, ncol = 2) + 
  labs(x = "Molt Width (mm)", y = "Carapace Width (mm)", title = "Molt vs. Carapace Width Separated By Stage", 
       fill = legendTitle, shape = legendTitle, stage = "Stage")



## MOLT PLOTS
#overview of molt widths across stages
ggplot(mw, aes(x = Stage, y = moltWidth)) + 
  labs(x = "Stage", y = "Molt Widths (mm)", title = "Molt Widths Across Stage", color = legendTitle) + 
  geom_boxplot(aes(colour = Treatment), position= position_dodge(0.9)) + 
  geom_jitter(aes(color = Treatment), trim = FALSE,
              binaxis = 'y', stackdir = 'center', dotsize = 0.8, alpha = 0.3,
              position = position_jitterdodge(0.4)) +
  theme_bw()
#plots separated based on stage
ggplot(mw, aes(x = Stage, y = moltWidth)) + 
  labs(x = "Stage", y = "Molt Widths (mm)", title = "Molt Widths Across Stage") + 
  geom_dotplot(aes(fill = Treatment, color = Treatment), trim = FALSE,
               binaxis = 'y', stackdir = 'center', dotsize = 0.8, alpha = 0.3,
               position = position_dodge(0.8)) + 
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  facet_wrap(~Stage, ncol = 2) + 
  theme_bw(base_size = 10)


#molt length vs. width separated by stage
ggplot(molt, aes(x = moltWidth, y = moltLength)) + 
  labs(x = "Molt Width (mm)", y = "Molt Length (mm)", title = "Length vs. Width of Molts Across Stages", 
       color = legendTitle, fill = legendTitle) + 
  geom_point(aes(color = Treatment), trim = FALSE,
              binaxis = 'y', stackdir = 'center', dotsize = 0.3, alpha = 0.5,
              position = position_jitter(0.2)) + 
  theme_bw() + 
#  facet_wrap(~Stage, ncol=2) + 
  geom_smooth(method = lm, aes(fill = Treatment))




#plots of body sizes vs. molt sizes
#how to get y axis scale to accommodate body & molt measurements?
ggplot(total, aes(x = moltWidth, y = carapaceWidth)) + 
  geom_point(aes(shape = Treatment, color = Stage)) + 
  geom_smooth(aes(shape = Treatment), method=lm, se = TRUE) + 
  theme_bw() + 
  labs(x = "Molt Width (mm)", y = "Living Crab Width (mm)", title = "Living Crab Width vs. Molt Width Across Stages",
       shape = legendTitle, color = "Stage")

ggplot(total, aes(moltWidth, carapaceWidth)) + 
  geom_point(aes(color = Stage, shape = Treatment)) +
  geom_smooth(aes(fill = Treatment), method = lm, se=TRUE) +
  theme_bw() +
  facet_wrap(~Stage, ncol = 2) + 
  labs(title = "Living Crab Width vs. Molt Width By Stage", x = "Molt Width (mm)", y = "Living Crab Width (mm)",
       color = "Stage", fill = legendTitle, shape = legendTitle)



#Function to calc stand errors and CI for error bars
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}  


summaryForErrorBars <- summarySE(data=molt, measurevar = "moltWidth", 
                             groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$moltWidth + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$moltWidth - summaryForErrorBars$ci
moltW <- ggplot(summaryForErrorBars, aes(x = Treatment, y = moltWidth, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=molt, aes(x = Treatment, y = moltWidth, colour = Stage), alpha = 0.6,
              position = position_jitter(0.4)) +
  theme_bw(base_size = 18) +
#  scale_color_manual(values = myColors) +
  labs(x = legendTitle, y = "Molt Width (mm)", title = "Treatment Comparison of Molt Widths") + 
  scale_color_discrete(breaks=c("J6","J5","J4","J3","J2","J1"))

myColors <- c("#f26d54", "#36bddc","#337d40", "#b6471d", "#00c23f", "#0231d1")


#### SAME AS ABOVE BUT FOR MOLT LENGTHS, NOT WIDTHS
summaryForErrorBars <- summarySE(data=molt, measurevar = "moltLength", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$moltLength + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$moltLength - summaryForErrorBars$ci
moltL <- ggplot(summaryForErrorBars, aes(x = Treatment, y = moltLength, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar)) +
  geom_point(size = 5) +
  geom_jitter(data=molt, aes(x = Treatment, y = moltLength, colour = Stage), alpha = 0.4,
              position = position_jitter(0.4)) +
  theme_bw(base_size = 18) +
  labs(x = legendTitle, y = "Molt Length (mm)", title = "Treatment Comparison of Molt Lengths")


ggarrange(moltW, moltL,
          labels = c("A", "B"),
          common.legend=TRUE,
          legend="bottom")





##### MOLT STATS
###helpful functions
#isolate stages to compare each one
isolate <- function(dd, stageName){
  single <- dd[dd$Stage == stageName,]
  return(single)
}

#function for finding live crab width : molt width ratio
ratio <- function(dataName) {
  slant <- dataName$carapaceWidth/dataName$moltWidth
  return(slant)
}

#function to run for each stage
fitComp <- function(dd){
  stat <- lme(Ratio~Stage + Treatment + Stage*Treatment, random = ~1|(Crab), data = dd, na.action = na.omit)
  return(summary(stat))
}

#function to compare each stage for ONLY molt sizes
#for molt widths
statisticW <- function(dataName){
  fitW <- lm(moltWidth ~ Treatment, data = dataName, na.action = na.omit)
  return(summary(fitW))
}

#for molt lengths
statisticL <- function(dataName){
  fitL <- lm(moltLength ~ Treatment, data = dataName, na.action = na.omit)
  return(summary(fitL))
}



### MOLT WIDTHS AND LENGTHS
#overall 
moltStatW <- lmer(moltWidth ~ Treatment + (1|Crab) + (1|Stage), data = molt)
summary(moltStatW)
#p value = 0.00019

moltStatL <- lmer(moltLength ~ Treatment + (1|Crab) + (1|Stage), data = molt)
summary(moltStatL)
#p value = 0.0146


#by stage
j1M <- isolate(molt, "J1")
j2M <- isolate(molt, "J2")
j3M <- isolate(molt, "J3")
j4M <- isolate(molt, "J4")
j5M <- isolate(molt, "J5")
j6M <- isolate(molt, "J6")

#compare j1
statisticW(j1M)
statisticL(j1M)
#compare j2
statisticW(j2M)
statisticL(j2M)
#compare j3
statisticW(j3M) #significant p value = 0.0044
statisticL(j3M)
#compare j4
statisticW(j4M) #significant p value = 0.0023
statisticL(j4M) #significant p value = 0.0025
#compare j5
statisticW(j5M) #significant p value = 0.0003
statisticL(j5M) #significant p value = 0.0005
#compare j6
statisticW(j6M)
statisticL(j6M)



### MOLTS & BODIES STATISTIC --> mostly body/molt ratio

#insert ratios into dataframe
total$Ratio <- ""
all <- ratio(total)
total$Ratio <- all

### AIC (overall)
fitAllT <- lmer(Ratio~Treatment + (1|Crab), data = total, na.action = na.omit)
isSingular(fitAllT, tol=1e-05) #TRUE

fitAllS <- lmer(Ratio~Stage + (1|Crab), data = total, na.action = na.omit)
isSingular(fitAllS, tol=1e-05) #TRUE

fitAllTS <- lmer(Ratio~Stage + Treatment + (1|Crab), data = total, na.action = na.omit)
isSingular(fitAllTS, tol=1e-05) #TRUE

fitAllTSTS <- lmer(Ratio~Stage + Treatment + Stage*Treatment + (1|Crab), data = total, na.action = na.omit)
isSingular(fitAllTSTS, tol=1e-05) #TRUE

AIC(fitAllT, fitAllS, fitAllTS, fitAllTSTS, k=2)
#fitAllTSTS has lowest AIC value
#run with treatment + stage + TR*ST parameter?


#### FOR BODY WIDTHS + MOLT WIDTHS
summaryForErrorBars <- summarySE(data=total, measurevar = "Ratio", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$Ratio + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$Ratio - summaryForErrorBars$ci
ggplot(summaryForErrorBars, aes(x = Treatment, y = Ratio, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar)) +
  geom_point(size = 5) +
  geom_jitter(data=total, aes(x = Treatment, y = Ratio, colour = Stage), alpha = 0.4,
              position = position_jitter(0.4)) +
  theme_bw(base_size = 15) +
  labs(x = legendTitle, y = "Body Width to Molt Width Ratio", title = "Treatment Comparison of Molt Lengths")


### per stage statistics

#sorry, you have to isolate the stages again (for the dataframe with LIVE crab widths, not MOLT LENGTHS)
#another isolation function until I can fix the weird nomenclature in 'total'
isolate2 <- function(dd, stageName){
  single <- dd[dd$Stage == stageName,]
  return(single)
}

j1T <- isolate2(total, "J1")
j2T <- isolate2(total, "J2")
j3T <- isolate2(total, "J3")
j4T <- isolate2(total, "J4")
j5T <- isolate2(total, "J5")
j6T <- isolate2(total, "J6")

#mixed model for each stage
fitComp(j1T)

stat1 <- lme(Ratio ~ Stage + Treatment + Stage*Treatment, random = ~1|Crab, data = j1T, na.action = na.omit)






