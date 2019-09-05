#Comparison of width to height over juvenile life stages
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("lme4")
#install.packages("reshape")
library(tidyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(nlme)
library(reshape)
library(lmerTest)

#set working directory
setwd("C:/Users/emma.reinhardt/Desktop/Research/Size/Data")
#read the csv file
size <- read.csv("sizes.csv", stringsAsFactors = FALSE)

#get rid of extra (empty) rows and columns
size <- subset(size, Group != "")

#turn from short-wide to long-skinny dataset
totalW <- gather(size, key = Stage, value = carapaceWidth, J1, J2, J3, J4, J5, J6, J7)
totalL <- gather(size, key = StageL, value = carapaceLength, J1_L, J2_L, 
                              J3_L, J4_L, J5_L, J6_L, J7_L)
#cut out redundant parts of the df's (the individual stages)
totalW <- subset(totalW, select = c(Crab, Group, Stage, carapaceWidth))
totalL <- subset(totalL, select = c(Crab, Group, StageL, carapaceLength))
#fix the excess symbols in totalL dataframe to make life easier when merging
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

totalL <- cleanUp(totalL)

#convert variable from character to number
totalW$carapaceWidth <- as.numeric(totalW$carapaceWidth)
totalL$carapaceLength <- as.numeric(totalL$carapaceLength)

#merge the cleaned length and width dataframes
total <- merge(totalW, totalL)

#add treatment variable
total$Treatment <- ""
total$Treatment[total$Group == "HA" | total$Group == "HB" | total$Group == "HC"] <- "High"
total$Treatment[total$Group == "LA" | total$Group == "LB" | total$Group == "LC"] <- "Low"


#function that determines the average ratio for W/L for each stage
ratio <- function(dataName) {
  slant <- dataName$carapaceWidth/dataName$carapaceLength
  return(slant)
}

#creating a new column for each W/L ratio
total$Ratio <- ""
all <- ratio(total)
total$Ratio <- all

#in case you want to test statistics on the data without J7 since sample size was so small
totalsub <- subset(total, Stage !="J7")

##### STATISTICS
### AIC (overall statistics)
#make different lmer models that randomize differen variables
#use AIC to find the model that best fits the data

#AIC
fitAllT <- lmer(Ratio~Treatment + (1|Crab) + (1|Stage), data = total, na.action = na.omit)
summary(fitAllT)

fitAllS <- lmer(Ratio~Stage + (1|Crab), data = total, na.action = na.omit)

fitAllTS <- lmer(Ratio~Stage + Treatment + (1|Crab), data = total, na.action = na.omit)

fitAllTSTS <- lmer(Ratio~Stage + Treatment + Stage*Treatment + (1|Crab), data = total, na.action = na.omit)

AIC(fitAllT, fitAllS, fitAllTS, fitAllTSTS, k=2)
#fitAllT has smallest AIC value
#use only Treatment parameter with randomizing by Crab and Stage 


### per stage statistics
#new function to isolate the singular stage I want to use 
isolate <- function(stageName){
  single <- total[total$Stage == stageName,]
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

#function for using linear mixed-effects model
findingFit <- function(dd){
  stat <- lme(Ratio~Treatment, random = ~1|Crab, data = dd, na.action = na.omit)
  return(summary(stat))
}

#find significance per stage
findingFit(j1)
findingFit(j2) #significant
findingFit(j3) #
findingFit(j4) #
findingFit(j5) #
findingFit(j6) #not significant
findingFit(j7) #significant




#######SLOPES (probably not useful)
#function to find slope
slope <- function(x, y){
  return(lm(x~y, na.action = na.omit)$coefficients[2])
}


#find slopes for each stage
fit <- lme(carapaceLength~carapaceWidth + (1|Crab), data = total, na.action = na.omit)
slope <- fit$coefficients[2]

#compare each slope of W:L ratio in individual stages
#stage j1
fitj1 <- lm(carapaceLength~carapaceWidth, data = j1)
slopej1 <- fitj1$coefficients[2]
summary(fitj1)
fitj1MIX <- lme(slopej1 ~ Treatment, data = j1, random = ~1|Crab, na.action = na.omit)
#low slope
lowj1 <- j1[!(j1$Treatment == "High CO2"),]
#high slope
highj2 <- j1[!(j1$Treatment == "Low CO2"),]

#stage j2
fitj2 <- lm(carapaceLength~carapaceWidth, data = j2)
slopej2 <- fitj2$coefficients[2]
summary(fitj2)

#stage j3
fitj3 <- lm(carapaceLength~carapaceWidth, data = j3)
slopej1 <- fitj3$coefficients[2]
summary(fitj3)

#stage j4
fitj4 <- lm(carapaceLength~carapaceWidth, data = j4)
slopej4 <- fitj4$coefficients[2]
summary(fitj4)

#stage j5
fitj5 <- lm(carapaceLength~carapaceWidth, data = j5)
slopej5 <- fitj5$coefficients[2]
summary(fitj5)

#stage j6
fitj6 <- lm(carapaceLength~carapaceWidth, data = j6)
slopej6 <- fitj6$coefficients[2]
summary(fitj6)

#stage j7
fitj7 <- lm(carapaceLength~carapaceWidth, data = j7)
slopej7 <- fitj7$coefficients[2]
summary(fitj7)







### PLOTS
#make a proper legend title
legendTitle = expression("CO"[2]*" Treatment")
#scatter plot overview with a regression line
ggplot(total, aes(x = carapaceWidth, y = carapaceLength)) +
  geom_point(aes(shape = Treatment, color = Stage), alpha = 0.5, size = 2.5) +
  labs(title = "Carapace Size Across Stages", x = "Carapace Width (mm)", y = "Carapace Length (mm)", 
       shape = legendTitle, color = "Stage", fill = legendTitle) +
  theme_bw(base_size=18) + 
  geom_smooth(aes(fill = Treatment), method = lm)

#separate data based on stage
ggplot(total, aes(x = carapaceWidth, y = carapaceLength)) +
  geom_point(aes(colour =  Treatment), alpha = 0.25) +
  labs(title = "Carapace Sizes At Each Stage", x = "Carapace Width (mm)", y = "Carapace Length (mm)",
       color = legendTitle) +
  facet_wrap(~Stage, ncol = 2) + theme_bw(base_size = 15) + 
  #add a 95% confidence interval & regression line
  geom_smooth(aes(color = Treatment), method = lm, se = TRUE)

#grouped boxplot + dot plot
ggplot(total, aes(x = carapaceWidth, y = carapaceLength)) + 
  geom_boxplot(aes(color = Treatment), width = 0.5, size = 0.4,
               position = position_dodge(0.8)) + 
  geom_dotplot(aes(fill = Treatment, color = Treatment), trim = FALSE,
               binaxis = 'y', stackdir = 'center', dotsize = 0.8,
               position = position_dodge(0.8)) + 
  facet_wrap(~Stage, ncol = 2) + theme_bw(base_size = 10)

#testing to see the differences in ratio across stages
ggplot(total, aes(x = Treatment, y= Ratio))+
  geom_jitter(aes(color=Stage), alpha=0.3)+
  facet_wrap(~Stage)+
  labs(x= legendTitle, y="Width to Length Carapace Ratio", title="Carapace Ratio Across Stages")+
  theme_bw(base_size=16)


###CONFIDENCE INTERVALS
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

totalsub <- subset(total, Stage != "J7")
summaryForErrorBars <- summarySE(data=totalsub, measurevar = "Ratio", 
                                 groupvars = c("Stage", "Treatment"), na.rm = TRUE)

summaryForErrorBars$upperBar <- summaryForErrorBars$Ratio + summaryForErrorBars$ci
summaryForErrorBars$lowerBar <- summaryForErrorBars$Ratio - summaryForErrorBars$ci
ggplot(summaryForErrorBars, aes(x = Treatment,y = Ratio, colour = Stage)) +
  geom_errorbar(aes(ymin = lowerBar, ymax = upperBar))+
  geom_point(size = 5) +
  geom_jitter(data=totalsub, aes(x = Treatment,y = Ratio, colour = Stage), alpha = 0.4,
              position = position_jitter(w = 0.35, h = 0)) +
  theme_bw(base_size = 18) + 
  labs(x = legendTitle, y = "Width to Length Ratio", title = "Treatment Comparison of Shape") + 
  facet_wrap(~Stage)












