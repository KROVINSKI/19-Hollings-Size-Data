#Reading Juv Length data
#install.packages("tidyr")
library(tidyr)
library(ggplot2)
library(nlme)
library(lme4)
library(lmerTest)

#set working directory
setwd("C:/Users/emma.reinhardt/Desktop/Research/Size/Data")
#read the csv file
d <- read.csv("sizes.csv", stringsAsFactors = FALSE)

#get rid of extra (empty) rows and columns
d <- subset(d, Group != "")
d <- subset(d, select = c(Crab, Group, J1_L, J2_L, J3_L, J4_L, J5_L, J6_L, J7_L))

#turn from short-wide to long-skinny dataset
dg <- gather(d, key = StageL, value = carapaceLength, J1_L, J2_L, J3_L, J4_L, J5_L, J6_L, J7_L)

#fix the weird symbols in the dataframe to make it easier
cleanUp <- function(dd){
  dd$Stage = ""
  dd$Stage[dd$StageL == "J1_L"] <- "J1"
  dd$Stage[dd$StageL == "J2_L"] <- "J2"
  dd$Stage[dd$StageL == "J3_L"] <- "J3"
  dd$Stage[dd$StageL == "J4_L"] <- "J4"
  dd$Stage[dd$StageL == "J5_L"] <- "J5"
  dd$Stage[dd$StageL == "J6_L"] <- "J6"
  dd$Stage[dd$StageL == "J7_L"] <- "J7"
  dd <- subset(dd, select = c(Crab, Group, Stage, carapaceLength))
  return(dd)
}

dg <- cleanUp(dg)

#convert variable from character to number
dg$carapaceLength <- as.numeric(dg$carapaceLength)

#add treatment variable
dg$Treatment <- ""
dg$Treatment[dg$Group == "HA" | dg$Group == "HB" | dg$Group == "HC"] <- "High" 
dg$Treatment[dg$Group == "LA" | dg$Group == "LB" | dg$Group == "LC"] <- "Low"

#make a proper legend title
legendTitle = expression("CO"[2]*" Treatment")



### LENGTH STATISTICS

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

dgsub <- subset(dg, Stage != "J7")
summaryForErrorBars <- summarySE(data=dgsub, measurevar = "carapaceLength", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$carapaceLength + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$carapaceLength - summaryForErrorBars$ci
ggplot(summaryForErrorBars, aes(x = Treatment,y = carapaceLength, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=dgsub, aes(x = Treatment,y = carapaceLength, colour = Stage), alpha = 0.4,
              position = position_jitter(0.4)) +
  labs(x = legendTitle, y = "Carapace Length (mm)", title = "Treatment Comparison of Carapace Length") + 
  theme_bw(base_size = 18) + 
  scale_color_discrete(breaks=c("J6","J5","J4","J3","J2","J1"))


#separate confidence intervals by stage
ggplot(summaryForErrorBars, aes(x = Treatment,y = carapaceLength, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=dgsub, aes(x = Treatment,y = carapaceLength, colour = Stage), alpha = 0.4,
              position = position_jitter(0.4)) +
  labs(x = legendTitle, y = "Carapace Length (mm)", title = "Treatment Comparison by Stage") + 
  facet_wrap(~Stage) + 
  theme_bw(base_size = 18)


##### Repeated Measures Anaalysis
#using a mixed model to find significant differences

#overall statistcal significance
fitTreatCrabStage <- lmer(carapaceLength ~ Treatment + (1|Crab) + (1|Stage), data = dg)
summary(fitTreatCrabStage)
#p value = 0.00525, overall treatment is significant
#could also cut out J7 stage


#find significance between stages
#isolate each stage
isolate <- function(stageName){
  single <- dg[dg$Stage == stageName,]
  return(single)
}

#isolate stages to compare each one
j1 <- isolate("J1")
j2 <- isolate("J2")
j3 <- isolate("J3")
j4 <- isolate("J4")
j5 <- isolate("J5")
j6 <- isolate("J6")
j7 <- isolate("J7")

statistic <- function(dataName){
  fit <- lme(carapaceLength ~ Treatment, data = dataName, random = ~1|Crab, na.action = na.omit)
  return(summary(fit))
}

#compare j1
statistic(j1)
#compare j2
statistic(j2)
#compare j3
statistic(j3)
#compare j4
statistic(j4) #significant
#compare j5
statistic(j5) #significant
#compare j6
statistic(j6) #significant
#compare j7
statistic(j7)



#make pretty plots
#lengths across stages
ggplot(dg, aes(x = Stage, y = carapaceLength)) +
  labs(title = "Carapace Length Across Life Stage", y = "Length of Carapace (mm)", x = "Stage", color =  legendTitle) + 
  geom_jitter(aes(color = Treatment), trim = FALSE,
              binaxis = 'y', stackdir = 'center', dotsize = 0.8, alpha = 0.3,
              position = position_jitterdodge(0.3)) +
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  theme_bw()
#separate graphs per stage
ggplot(dg, aes(x = Stage, y = carapaceLength)) + 
  labs(x = "Stage", y = "Length of Carapace (mm)", title = "Carapace Length for Each Stage") +
  geom_jitter(aes(color = Treatment), trim = FALSE,
              stackdir = 'center', binaxis = 'y', dotsize = 0.8, alpha = 0.3,
              position = position_jitterdodge(0.4)) + 
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  theme_bw() + 
  facet_wrap(~Stage, ncol=2)


