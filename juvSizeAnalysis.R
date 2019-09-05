#Reading Juv size data
#install.packages("tidyr")
#install.packages("ggpubr")
#install.packages("nlme")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(plyr)
library(nlme)
library(lmerTest)
library(lme4)

#set working directory
setwd("C:/Users/emma.reinhardt/Desktop/Research/Size/Data")
#read the csv file
d <- read.csv("sizes.csv", stringsAsFactors = FALSE)

#get rid of extra (empty) rows and columns
d <- subset(d, Group != "")
d <- subset(d, select = c(Crab, Group, J1, J2, J3, J4, J5, J6, J7))

#turn from short-wide to long-skinny dataset
dg <- gather(d, key = Stage, value = carapaceWidth, J1, J2, J3, J4, J5, J6, J7)

#convert variable from character to number
dg$carapaceWidth <- as.numeric(dg$carapaceWidth)
#add treatment variable
dg$Treatment <- ""
dg$Treatment[dg$Group == "HA" | dg$Group == "HB" | dg$Group == "HC"] <- "High" 
dg$Treatment[dg$Group == "LA" | dg$Group == "LB" | dg$Group == "LC"] <- "Low"

#legend title - for future use
legendTitle = expression("CO"[2]*" Treatment")



### WIDTH STATISTICS

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
summaryForErrorBars <- summarySE(data=dgsub, measurevar = "carapaceWidth", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$carapaceWidth + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$carapaceWidth - summaryForErrorBars$ci
ggplot(summaryForErrorBars, aes(x = Treatment,y = carapaceWidth, colour = Stage)) + 
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=dgsub, aes(x = Treatment,y = carapaceWidth, colour = Stage), alpha = 0.4,
              position = position_jitter(w = 0.35, h = 0)) +
  theme_bw(base_size = 18) + 
  labs(x = legendTitle, y = "Carapace Width (mm)", title = "Treatment Comparison of Carapace Widths") + 
    scale_color_discrete(breaks=c("J6","J5","J4","J3","J2","J1"))


#separate confidence intervals based on stage
ggplot(summaryForErrorBars, aes(x = Treatment,y = carapaceWidth, colour = Stage)) + 
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=dgsub, aes(x = Treatment,y = carapaceWidth, colour = Stage), alpha = 0.4,
              position = position_jitter(0.3)) +
  theme_bw(base_size = 18) + 
  labs(x = legendTitle, y = "Carapace Width (mm)", title = "Treatment Comparison by Stage") + 
  facet_wrap(~Stage)



### Repeated Measures Analysis

#fitLMtreatment <- lm(carapaceWidth ~ Treatment, data = dg)
#summary(fitLMtreatment)

#fitLMstage <- lm(carapaceWidth ~ Stage, data = dg)
#summary(fitLMstage)
#dg$Crab <- factor(dg$Crab)

#fitTreatmentRMcrab <- lmer(carapaceWidth ~ Treatment + (1|Crab), data = dg)
#summary(fitTreatmentRMcrab)

#fitTreatmentRMcrab <- lmer(carapaceWidth ~ Treatment + (1|Crab), data = dg)
#summary(fitTreatmentRMcrab)

#overall statistcal significance
fitTreatmentRMcrabstage <- lmer(carapaceWidth ~ Treatment + (1|Crab) + (1|Stage), data = dgsub)
summary(fitTreatmentRMcrabstage)
#p value = 0.000103, overall treatment is significant


#finding p values when comparing stages
lmeStat <- function(dd, stageName){
  dd <- dd[dd$Stage == stageName,]
  fit <- lme(carapaceWidth ~ Treatment, data = dd, random = ~1|Crab, na.action = na.omit)
  return(summary(fit))
}


#find p values for each stage comparison
pVAl1 <- lmeStat(dg,"J1")$tTable[2,5]
pVAl2 <- lmeStat(dg,"J2")$tTable[2,5]
pVAl3 <- lmeStat(dg,"J3")$tTable[2,5] #significant begins here
pVAL4 <- lmeStat(dg,"J4")$tTable[2,5] #
pVAL5 <- lmeStat(dg,"J5")$tTable[2,5] #
pVAL6 <- lmeStat(dg,"J6")$tTable[2,5] #signifance ends here
pVAL7 <- lmeStat(dg,"J7")$tTable[2,5] 


#roundabout way of finding p values
#new data frames for each stage 
#just in case these are needed
j1 <- dg[(dg$Stage == "J1"),]
j2 <- dg[(dg$Stage == "J2"),]
j3 <- dg[(dg$Stage == "J3"),]
j4 <- dg[(dg$Stage == "J4"),]
j5 <- dg[(dg$Stage == "J5"),]
j6 <- dg[(dg$Stage == "J6"),]
j7 <- dg[(dg$Stage == "J7"),]



### BODY WIDTH PLOTS
#dotplot of treatments at each stage
ggplot(dg, aes(x = Stage, y = carapaceWidth)) +
  geom_point(aes(color = Treatment), alpha = 0.5) + 
  labs(title = "Carapace Width Across Life Stage", y = "Width of Carapace (mm)", x = "Life Stage", color =  legendTitle) + 
  theme_minimal()
#boxplot + jitter overview
ggplot(dg, aes(x = Stage, y = carapaceWidth, colour = Treatment)) +
  labs(title = "Carapace Width Across Life Stages", y = "Width of Carapace (mm)", x = "Life Stage", color = legendTitle) + 
  geom_jitter(aes(color = Treatment), trim = FALSE,
              binaxis = 'y', stackdir = 'center', dotsize = 0.3, alpha = 0.3,
              position = position_jitterdodge(0.2)) +
  geom_boxplot() + 
  theme_bw()
#boxplot + point at each stage
ggplot(dg, aes(x = Stage, y = carapaceWidth)) + 
  labs(x = "Stage", y = "Carapace Width (mm)", title = "Carapace Width Across Life Stages") + 
  geom_jitter(aes(color = Treatment), trim = FALSE,
              binaxis = 'y', stackdir = 'center', dotsize = 0.8, alpha = 0.3,
              position = position_jitterdodge(0.3)) +
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  #  geom_dotplot(aes(fill = Treatment, color = Treatment), trim = FALSE,
  #             binaxis = 'y', stackdir = 'center', dotsize = 0.8,
  #             position = position_dodge(0.8)) + 
  facet_wrap(~Stage, ncol = 2) + 
  theme_bw(base_size = 10)





##### WEIGHTS
#weight file
w = read.csv("crab_weights.csv")
#add in treatment column
w$Treatment = ""
w$Treatment[w$Group == "HA" | w$Group == "HB" | w$Group == "HC"] <- "HighCO2" 
w$Treatment[w$Group == "LA" | w$Group == "LB" | w$Group == "LC"] <- "LowCO2"
#make everything is numeric
w$Weight = as.numeric(w$Weight)
#merge the weights file with the size file
ws <- merge(w, dg)

#WEIGHT PLOTS
#weight across stages
ggplot(w, aes(x = Stage, y = Weight)) + 
  labs(x = "Stage", y = "Crab Weight (g)", title = "Crab Weights Across Stages") + 
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  geom_jitter(aes(fill = Treatment, color = Treatment),
              binaxis = 'y', stackdir = 'center', dotsize = 0.5, alpha = 0.4,
              position = position_jitterdodge(0.3)) + 
  theme_bw()

#WEIGHT + SIZE PLOTS
ggplot(ws, aes(x = Weight, y = carapaceWidth)) + 
  geom_jitter(aes(color = Stage)) + 
  geom_smooth(aes(color = Stage), method = lm, se = TRUE) + 
  theme_bw()


##confidence interval graph again with weight data

summaryForErrorBars <- summarySE(data=ws, measurevar = "Weight", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$Weight + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$Weight - summaryForErrorBars$ci
ggplot(summaryForErrorBars, aes(x = Treatment, y = Weight, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5, position = position_dodge(0.3)) +
  geom_jitter(data=ws, aes(x = Treatment, y = Weight, colour = Stage), alpha = 0.4,
              position = position_dodge(0.3)) +
  theme_bw(base_size = 14) + 
  labs(x = "Treatment", y = "Weight (g)", title = "Crab Weights Across Stages")











