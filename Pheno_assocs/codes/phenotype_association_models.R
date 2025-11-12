library(data.table)
p<-fread('E:/Shuai/Data/olink_data_21259.txt') 
library(reshape2)
p2<-dcast(p, eid ~ protein_id, value.var="result", fun.aggregate = mean)
write.table(p2, file='E:/Shuai/Data/Wide_format_Olink_proteins_UKB21259_Feb2024.txt', quote=F, row.names=F )
library(data.table)
p2<-fread("E:/Shuai/Data/Wide_format_Olink_proteins_UKB21259_Feb2024.txt")

##########################################################################
p2<-as.data.frame(p2)

percentmissing<-NA
### % missing by protein
for(i in 2:2924){
  percentmissing[i-1]<-(sum(is.na(p2[,i]))/53073)*100
  print(i)
}
mean(percentmissing)
max(percentmissing)
hist(percentmissing)
sum(percentmissing>20) #12 proteins have >20% missing
percentmissingb<-which(percentmissing>20)+1
p3<-p2[,-percentmissingb]

percentmissing2<-NA
for(i in 2:2912){
  percentmissing2[i-1]<-(sum(is.na(p3[,i]))/53073)*100
}
sum(percentmissing2>20)
hist(percentmissing2)

#missing proteins per person
missingpercent<-NULL
for(j in 1:53073){
  missingpercent[j]<-(sum(is.na(p3[j, ]))/2911)*100
  print(j)
}
mean(missingpercent)
hist(missingpercent)
sum(missingpercent>50)

missingpercent <- apply(p3, 1, function(x) 100*sum(is.na(x))/2911)
rows_to_keep <- missingpercent <= 50
p4 <- p3[rows_to_keep, ]
fwrite(p4, "E:/Shuai/Data/proteomics_QCd_rowscolumnsexcl.tab")



#################################################
# Bring in UKB PA, covariate, basic phenos
library(data.table)
p4<-fread("E:/Shuai/Data/proteomics_QCd_rowscolumnsexcl.tab")

library(dplyr)
p4 <- rename(p4, f.eid = eid)

bdE <- fread("E:/Shuai/Data/ukb8427.tab")
bdE[1:7, 1:7, with=FALSE]
#names(p)[1] <- "f.eid"
length(intersect(p4$f.eid, bdE$f.eid))
m <- merge(bdE, p4, by="f.eid")
fwrite(m, "E:/Shuai/Data/merged_protein_phenoUKB21259.tab")
#################################################

m<-fread("E:/Shuai/Data/merged_protein_phenoUKB21259.tab")
######make covariates: season, race/ethncity, as.factor(smoking), drinking, TDI, %fat
summary(m$f.55.0.0)
m$Season<-NA
m$Season[which(m$f.55.0.0==12 | m$f.55.0.0==1 | m$f.55.0.0==2 |m$f.55.0.0==3 )]<-0
m$Season[which(m$f.55.0.0==4 | m$f.55.0.0==5 | m$f.55.0.0==10 |m$f.55.0.0==11 )]<-1
m$Season[which(m$f.55.0.0==6 | m$f.55.0.0==7 | m$f.55.0.0==8 |m$f.55.0.0==9 )]<-2
summary(m$Season)
table(m$Season)
table(m$f.21000.0.0)
summary(m$f.21000.0.0)
m$race <- 2
m$race<-ifelse(m$f.21000.0.0%in%c(1,1001, 1002,1003), 1, m$race)
m$race<-ifelse(m$f.21000.0.0%in%c(3,3001,3002, 3003, 3004), 3, m$race)
m$race<-ifelse(m$f.21000.0.0%in%c(4,4001,4002, 4003), 4, m$race)
m$race<-ifelse(m$f.21000.0.0%in%c(5), 5, m$race)
m$race<-ifelse(m$f.21000.0.0%in%c(6), 6, m$race)
m$race<-ifelse(m$f.21000.0.0%in%c(-3,-1,NA), NA, m$race)
table(m$race)
summary(m$race)
table(m$f.20116.0.0)
summary(m$f.20116.0.0)
m$as.factor(smoking)<-NA
m$as.factor(smoking)[which(m$f.20116.0.0==0 )]<-0
m$as.factor(smoking)[which(m$f.20116.0.0==1 )]<-1
m$as.factor(smoking)[which(m$f.20116.0.0==2 )]<-2
table(m$as.factor(smoking))
summary(m$as.factor(smoking))
table(m$f.1558.0.0)
summary(m$f.1558.0.0)
m$drinking<-NA
m$drinking[which(m$f.1558.0.0==1 )]<-6
m$drinking[which(m$f.1558.0.0==2 )]<-3.5
m$drinking[which(m$f.1558.0.0==3 )]<-1.5
m$drinking[which(m$f.1558.0.0==4 )]<-0.5
m$drinking[which(m$f.1558.0.0==5 )]<-0.1
m$drinking[which(m$f.1558.0.0==6 )]<-0
table(m$drinking)
summary(m$drinking) 
table(m$f.23099.0.0)
summary(m$f.23099.0.0)
names(m)[names(m) == "f.23099.0.0"] <- "bodyfat"
table(m$bodyfat)
summary(m$bodyfat) 
hist(m$bodyfat)  
boxplot(m$bodyfat, main="Boxplot of Bodyfat", ylab="Values")
Q1 <- quantile(m$bodyfat, 0.25, na.rm = TRUE)
Q3 <- quantile(m$bodyfat, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1
outliers <- m$bodyfat[(m$bodyfat < (Q1 - 1.5 * IQR)) | (m$bodyfat > (Q3 + 1.5 * IQR))]
z_scores <- scale(m$bodyfat)
outliers <- m$bodyfat[abs(z_scores) > 3]
table(m$f.6138.0.0)
summary(m$f.6138.0.0)
m$Education <- NA
m$Education <- ifelse(m$f.6138.0.0 == -3 , NA, 0)
m$Education[m$f.6138.0.0 == 1] <- 1
table(m$Education)

fwrite(m, "E:/Shuai/Data/covariates_merged_protein_phenoUKB21259.tab")
#################################################





m<-fread("E:/Shuai/Data/covariates_merged_protein_phenoUKB21259.tab")

####################################################################################################################################
#Make PA phenos
#####################################################################################################################################
## Make SSOE variable ########################################################################################
table(m$f.991.0.0)
#assign NA to people who are coded as -3 or -1 for variables f.991.0.0 (Days of SS) and 1001.0.0 (Duration of SS)
m$DaysSS<-ifelse(m$f.991.0.0==-3 | m$f.991.0.0==-1 , NA, m$f.991.0.0)
table(m$DaysSS)

table(m$f.1001.0.0)
m$DurSS<-ifelse(m$f.1001.0.0==-3 | m$f.1001.0.0==-1 , NA, m$f.1001.0.0)
table(m$DurSS)

#Identify those who do >=2 days SS and duration >=15-30 mins #############################
m$TwoplusdaysSS<-ifelse(m$DaysSS>=4 & m$DurSS>1  , 1, 0) #& bd$DurSS>1
table(m$TwoplusdaysSS)
##########################################################################################

#assign NA to people who are coded as -3 or -1 for variables f.3637.0.0 and 3647.0.0
m$DaysOtherEx<-ifelse(m$f.3637.0.0==-3 | m$f.3637.0.0==-1 , NA, m$f.3637.0.0)
table(m$DaysOtherEx)
m$DurOtherEx<-ifelse(m$f.3647.0.0==-3 | m$f.3647.0.0==-1 , NA, m$f.3647.0.0)
table(m$DurOtherEx)

#Identify those who do >=2 days OE and duration >=15-30 mins ##############################
m$TwoplusdaysOther<-ifelse(m$DaysOtherEx>=4 & m$DurOtherEx>1,1,0)

m$StrenorOther<-ifelse(!is.na(m$f.991.0.0) | !is.na(m$f.3637.0.0), NA, 0)
#Controls are those who do not do either strenuous exercise or other exercises
m$StrenorOther[which(m$TwoplusdaysSS==1 | m$TwoplusdaysOther==1)]<-1
table(m$StrenorOther)
#Exclude individuals who can't walk
m$StrenorOther<-ifelse(m$f.864.0.0==-2, NA, m$StrenorOther)

## Make MVPA variable ######################################################################################
############################################################################################################

#Vigorous activity ############
#assign NA to people who are coded as -3 or -1 for variables f.904.0.0 and f.914.0.0, and assign 0 duration for those saying 0 days per week of VPA
table(m$f.904.0.0)
m$DaysperWeekVigPA<-ifelse(m$f.904.0.0==-3 | m$f.904.0.0==-1 , NA, m$f.904.0.0)
m$DurationVigPA<-ifelse(m$f.914.0.0==-3 | m$f.914.0.0==-1, NA, m$f.914.0.0)
m$DurationVigPA<-ifelse(m$DaysperWeekVigPA==0, 0 ,m$DurationVigPA )

#recode people reporting more than 960 minutes of VPA as NA, and those reporting >3 hours of VPA as 3 hours
m$DurationVigPA<-ifelse(m$DurationVigPA>960, NA, m$DurationVigPA)
m$DurationVigPA<-ifelse(m$DurationVigPA>=180, 180, m$DurationVigPA)

#Multiply days/week by duration of Vig PA
m$VigPA_Minsperweek<-m$DaysperWeekVigPA * m$DurationVigPA
#######################################

# Moderate activity ###################
table(m$f.884.0.0)
#assign NA to people who are coded as -3 or -1 for variables f.884.0.0 and f.894.0.0, and assign 0 duration for those saying 0 days per week of MPA
m$DaysperWeekModPA<-ifelse(m$f.884.0.0==-3 | m$f.884.0.0==-1 , NA, m$f.884.0.0)
m$DurationModPA<-ifelse(m$f.894.0.0==-3 | m$f.894.0.0==-1, NA, m$f.894.0.0)
m$DurationModPA<-ifelse(m$DaysperWeekModPA==0, 0,m$DurationModPA )

#recode people reporting more than 960 minutes of MPA as NA, and those reporting >3 hours of MPA as 3 hours
m$DurationModPA<-ifelse(m$DurationModPA>960, NA, m$DurationModPA)
m$DurationModPA<-ifelse(m$DurationModPA>=180, 180, m$DurationModPA)

#Multiply days/week by duration of Mod PA
m$ModPA_Minsperweek<-m$DaysperWeekModPA * m$DurationModPA

#mUlitply 4 for Moderate and 8 for vigorous activity
m$MVPA<-(m$ModPA_Minsperweek*4) + (m$VigPA_Minsperweek*8)

#Exclude individuals who can't walk
m$MVPA<-ifelse(m$f.864.0.0==-2, NA, m$MVPA)

#INVERSE NORMALIZATION
my.invnorm = function(x)
{
  res = rank(x, na.last = "keep")
  res = qnorm(res/(length(res)+0.5))
  return(res)
}
m$InvNorm_MVPA<-my.invnorm(m$MVPA)
hist(m$InvNorm_MVPA)



## Make VPA variable #######################################################################################
#assign NA to people who are coded as -3 or -1 for variables f.904.0.0 and f.914.0.0, and assign 0 duration for those saying 0 days per week of VPA
table(m$f.904.0.0)
m$DaysperWeekVigPA<-ifelse(m$f.904.0.0==-3 | m$f.904.0.0==-1 , NA, m$f.904.0.0)
table(m$DaysperWeekVigPA)

m$DurationVigPA<-ifelse(m$f.914.0.0==-3 | m$f.914.0.0==-1, NA, m$f.914.0.0)
m$DurationVigPA<-ifelse(m$DaysperWeekVigPA==0, 0 ,m$DurationVigPA )

##Multiply days/week by duration of Vig PA, and recode people reporting more than 960 minutes of VPA as NA
m$VPA<-m$DaysperWeekVigPA * m$DurationVigPA
m$VPA<-ifelse(m$DurationVigPA>960, NA, m$VPA)
hist(m$VPA)

#Dichotomoize VPA phenotype #####################################################################################
m$VigPA_Dich2<-NA 
m$VigPA_Dich2[which(m$DaysperWeekVigPA==0)]<-0
m$VigPA_Dich2[which(m$DurationVigPA>=25 & m$DaysperWeekVigPA>=3)]<-1
table(m$VigPA_Dich2)

sum(!is.na(m$DaysperWeekVigPA))
sum(!is.na(m$DurationVigPA))
sum(!is.na(m$DaysperWeekVigPA) & m$GenExcl3==0)
sum(is.na(m$VigPA_Dich2)) 
####MAKE EXCLUSIONS ################################################################################################
#Exclude individuals who can't walk
m$VigPA_Dich2<-ifelse(m$f.864.0.0==-2, NA, m$VigPA_Dich2)
#remove people reporting 16 or more hours/day (16*8=480)
m$VigPA_Dich2<-ifelse(m$DurationVigPA>960, NA, m$VigPA_Dich2)
sum(!is.na(m$DaysperWeekVigPA) & m$GenExcl3==0)
table(m$VigPA_Dich2)
hist(m$VigPA_Dich2)


## Make AveAcc variable ######################################################################################
#############################################################################################################

#get Age ###################################
#first get days between f.90010.0.0 and f.21003.0.0
m$daysfromBasetoAccel<-as.Date(m$f.90010.0.0) - as.Date(m$f.53.0.0)
m$AGE_atAccel<-m$f.21003.0.0 + ((m$daysfromBasetoAccel)/365)
m$AGE_atAccel<-as.numeric(m$AGE_atAccel)
summary(m$AGE_atAccel)
###########################################
#ADD SEASON1 VARIABLE
library(lubridate)
m$f.90010.0.0<-as.Date(m$f.90010.0.0)
m$MonthAccel<-month(m$f.90010.0.0)

m$Season1<-NA
m$Season1[which(m$MonthAccel==12 | m$MonthAccel==1 | m$MonthAccel==2 |m$MonthAccel==3 )]<-0
m$Season1[which(m$MonthAccel==4 | m$MonthAccel==5 | m$MonthAccel==10 |m$MonthAccel==11 )]<-1
m$Season1[which(m$MonthAccel==6 | m$MonthAccel==7 | m$MonthAccel==8 |m$MonthAccel==9 )]<-2
table(m$Season1)


#set NA to f.90143.0.0 for those that are low-quality data ########################################################
m$Accel_Ave_DQexcl<-ifelse(m$f.90015.0.0==0, NA, m$f.90012.0.0)
sum(is.na(m$Accel_Ave_DQexcl))
###################################################################################################################
summary(m$f.90015.0.0)
summary(m$f.90012.0.0)
summary(m$f.864.0.0)
#Exclude individuals who can't walk###################################
m$Accel_Ave_DQexcl<-ifelse(m$f.864.0.0==-2, NA, m$Accel_Ave_DQexcl)
sum(is.na(m$Accel_Ave_DQexcl))

#Set NA to outliers because some major outliers ###############################################################################
summary(m$Accel_Ave_DQexcl)
mean(m$Accel_Ave_DQexcl, na.rm=T)
sd(m$Accel_Ave_DQexcl, na.rm=T)

sum(m$Accel_Ave_DQexcl> (mean(m$Accel_Ave_DQexcl, na.rm=T)+(sd(m$Accel_Ave_DQexcl, na.rm=T) *4)), na.rm=T) 
a<-mean(m$Accel_Ave_DQexcl, na.rm=T)+(sd(m$Accel_Ave_DQexcl, na.rm=T) *4)
m$Accel_Ave_DQexcl_minus46outliers<-ifelse(m$Accel_Ave_DQexcl>a, NA ,m$Accel_Ave_DQexcl)
sum(!is.na(m$Accel_Ave_DQexcl_minus46outliers))
hist(m$Accel_Ave_DQexcl_minus46outliers)
## Make Acc425 variable ######################################################################################
#############################################################################################################
#get Age ###################################
#first get days between f.90010.0.0 and f.21003.0.0
m$daysfromBasetoAccel<-as.Date(m$f.90010.0.0) - as.Date(m$f.53.0.0)
m$AGE_atAccel<-m$f.21003.0.0 + ((m$daysfromBasetoAccel)/365)
m$AGE_atAccel<-as.numeric(m$AGE_atAccel)
summary(m$AGE_atAccel)
###########################################
#ADD SEASON1 VARIABLE
library(lubridate)
m$f.90010.0.0<-as.Date(m$f.90010.0.0)
m$MonthAccel<-month(m$f.90010.0.0)

m$Season1<-NA
m$Season1[which(m$MonthAccel==12 | m$MonthAccel==1 | m$MonthAccel==2 |m$MonthAccel==3 )]<-0
m$Season1[which(m$MonthAccel==4 | m$MonthAccel==5 | m$MonthAccel==10 |m$MonthAccel==11 )]<-1
m$Season1[which(m$MonthAccel==6 | m$MonthAccel==7 | m$MonthAccel==8 |m$MonthAccel==9 )]<-2
table(m$Season1)
#set NA to f.90143.0.0 for those that are low-quality data
m$Accel_425_DQexcl<-ifelse(m$f.90015.0.0==0, NA, m$f.90140.0.0)
sum(is.na(m$Accel_425_DQexcl))
hist(m$f.90140.0.0)
hist(m$Accel_425_DQexcl)
#Set NA to outliers because some major outliers
summary(m$Accel_425_DQexcl)
mean(m$Accel_425_DQexcl, na.rm=T)
sd(m$Accel_425_DQexcl, na.rm=T)
hist(m$Accel_425_DQexcl)
mean(m$f.90140.0.0, na.rm=T)
sum(m$Accel_425_DQexcl<(mean(m$Accel_425_DQexcl, na.rm=T)+(sd(m$Accel_425_DQexcl, na.rm=T))) *4, na.rm=T) 
m$Accel_425_DQexcl_minus465outliers<-ifelse(m$Accel_425_DQexcl<(mean(m$Accel_425_DQexcl, na.rm=T)-(sd(m$Accel_425_DQexcl, na.rm=T)*4)) , NA ,m$Accel_425_DQexcl)

sum(is.na(m$Accel_425_DQexcl_minus465outliers))
hist(m$Accel_425_DQexcl_minus465outliers)
summary(m$Accel_425_DQexcl_minus465outliers)
m$OneMinus_Accel_425_DQexcl_minus465outliers<-1 - m$Accel_425_DQexcl_minus465outliers 
sum(m$Accel_425_DQexcl_minus465outliers==1, na.rm=T)
sum(m$OneMinus_Accel_425_DQexcl_minus465outliers==0, na.rm=T)
sum(!is.na(m$OneMinus_Accel_425_DQexcl_minus465outliers))
summary(m$OneMinus_Accel_425_DQexcl_minus465outliers)
hist(m$OneMinus_Accel_425_DQexcl_minus465outliers)
mean(!is.na(m$OneMinus_Accel_425_DQexcl_minus465outliers))
mean(m$OneMinus_Accel_425_DQexcl_minus465outliers, na.rm = TRUE)
rounded_median <- round(median(m$OneMinus_Accel_425_DQexcl_minus465outliers, na.rm = TRUE), 6)
print(rounded_median)
my.invnorm <- function(x) {
  res <- rank(x, na.last = "keep")
  res <- qnorm(res / (length(res) + 0.5))
  return(res)
}
m$InvNorm_Accel425<-my.invnorm(m$OneMinus_Accel_425_DQexcl_minus465outliers)
hist(m$InvNorm_Accel425)


#########ACCSED
library(data.table)
MachineLearningData <- fread('E:/Shuai/Data/MachineLearningData.csv')
setnames(MachineLearningData, old = "ID_21259", new = "f.eid")
m <- merge(m, MachineLearningData, by = "f.eid", all.x = TRUE, all.y = FALSE)
#get Age ###################################
#first get days between f.90010.0.0 and f.21003.0.0
m$daysfromBasetoAccel<-as.Date(m$f.90010.0.0) - as.Date(m$f.53.0.0)
m$AGE_atAccel<-m$f.21003.0.0 + ((m$daysfromBasetoAccel)/365)
m$AGE_atAccel<-as.numeric(m$AGE_atAccel)
summary(m$AGE_atAccel)
###########################################
#ADD SEASON1 VARIABLE
library(lubridate)
m$f.90010.0.0<-as.Date(m$f.90010.0.0)
m$MonthAccel<-month(m$f.90010.0.0)

m$Season1<-NA
m$Season1[which(m$MonthAccel==12 | m$MonthAccel==1 | m$MonthAccel==2 |m$MonthAccel==3 )]<-0
m$Season1[which(m$MonthAccel==4 | m$MonthAccel==5 | m$MonthAccel==10 |m$MonthAccel==11 )]<-1
m$Season1[which(m$MonthAccel==6 | m$MonthAccel==7 | m$MonthAccel==8 |m$MonthAccel==9 )]<-2
table(m$Season1)

#set NA for those that are low-quality data ########################################################
m$Acc_sed<-ifelse(m$f.90015.0.0==0, NA, m$aveSed)
summary(m$f.90015.0.0)
summary(m$aveSed)

sum(is.na(m$Acc_sed))
###################################################################################################################
summary(m$Acc_sed)
hist(m$Acc_sed)
hist(m$aveSed)
#Exclude individuals who can't walk###################################
m$Acc_sed<-ifelse(m$f.864.0.0==-2, NA, m$Acc_sed)
sum(is.na(m$Acc_sed))

#Set NA to outliers because some major outliers ###############################################################################
summary(m$Acc_sed)
mean(m$Acc_sed, na.rm=T)
sd(m$Acc_sed, na.rm=T)

sum(m$Acc_sed> (mean(m$Acc_sed, na.rm=T)+(sd(m$Acc_sed, na.rm=T) *4)), na.rm=T) 
a<-mean(m$Acc_sed, na.rm=T)+(sd(m$Acc_sed, na.rm=T) *4)
m$Acc_sed_minus46outliers<-ifelse(m$Acc_sed>a, NA ,m$Acc_sed)
sum(!is.na(m$Acc_sed_minus46outliers))

mean_val <- mean(m$Acc_sed, na.rm = TRUE)
sd_val <- sd(m$Acc_sed, na.rm = TRUE)
upper_threshold <- mean_val + (sd_val * 4)
lower_threshold <- mean_val - (sd_val * 4)
m$Acc_sed_minus46outliers <- ifelse(m$Acc_sed >= lower_threshold & m$Acc_sed <= upper_threshold, m$Acc_sed, NA)
sum(!is.na(m$Acc_sed_minus46outliers))
#set NA for people who should be excluded becuase of withdrawal of consent ########################################################################



hist(m$StrenorOther)
hist(m$InvNorm_MVPA)
hist(m$VigPA_Dich2)
hist(m$Accel_Ave_DQexcl_minus46outliers)
hist(m$InvNorm_Accel425)
hist(m$Acc_sed_minus46outliers)

#  Bring in dementia incidence data
dem<-fread("E:/Shuai/Data/biobankDementia.csv")
m2<-merge(m, dem, by='f.eid')
########m2$StrenorOther[which(m2$f.eid%in%excl$V1)]<-NA
excl<-read.table('E:/Shuai/Data/w21259_2023-04-25.csv', header=F)
m2 <- m2[!m2$f.eid %in% excl$V1, ]
sum(is.na(m2$StrenorOther))
sum(excl$V1%in%m2$f.eid)
summary(excl$V1)
class(m2$f.eid)
class(excl$V1)
length(unique(excl$V1))

fwrite(m2, "E:/Shuai/Data/dem_PA_covariates_merged_protein_phenoUKB21259.tab")
library(data.table)
m2<-fread("E:/Shuai/Data/dem_PA_covariates_merged_protein_phenoUKB21259.tab")
hist(m2$f.46.0.0)
hist(m2$f.47.0.0)
m2$gripstrength <- (m2$f.46.0.0 + m2$f.47.0.0) / 2
hist(m2$gripstrength)
hist(m2$f.2020.0.0)
m2$f.2020.0.0[m2$f.2020.0.0 == -1] <- NA
m2$f.2020.0.0[m2$f.2020.0.0 == -3] <- NA

#################################################################################################################################
library("survival")
library("survminer")
###############################################################################################################################################################
# Test association of each of the five PAs with incident dementia in five separate models
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$StrenorOther + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat+ m2$chronic_conditions)
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$InvNorm_MVPA + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat + m2$chronic_conditions)
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$VigPA_Dich2 + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat + m2$chronic_conditions)
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$Accel_Ave_DQexcl_minus46outliers + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat + m2$chronic_conditions)
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$InvNorm_Accel425 + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat + m2$chronic_conditions)
r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$Acc_sed_minus46outliers + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat + m2$chronic_conditions)

#save results for PA association (results row) somewhere
#####################################################################################################################################################################

results_df <- data.frame(
  PA_Variable = character(),
  Coef = numeric(),
  StdErr = numeric(),
  HazardRatio = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

pa_variables <- c("StrenorOther", "InvNorm_MVPA", "VigPA_Dich2", "Accel_Ave_DQexcl_minus46outliers", "InvNorm_Accel425", "Acc_sed_minus46outliers")

for (pa_var in pa_variables) {
  formula <- as.formula(paste("Surv(FU_time_Dem, dementia) ~", pa_var, "+ Season + as.factor(race) + f.21003.0.0 + f.31.0.0 + Education + drinking + as.factor(smoking) + bodyfat", sep = " "))
  model <- coxph(formula, data = m2)
  
  summary_model <- summary(model)
  
  coef <- summary_model$coef[1, "coef"]
  stderr <- summary_model$coef[1, "se(coef)"]
  hazard_ratio <- exp(coef)
  p_value <- summary_model$coef[1, "Pr(>|z|)"]
  ci <- confint(model, level = 0.95)[pa_var, ] 
  ci_lower <- exp(ci[1])  
  ci_upper <- exp(ci[2])  
  
  results_df <- rbind(results_df, data.frame(
    PA_Variable = pa_var,
    Coef = coef,
    StdErr = stderr,
    HazardRatio = hazard_ratio,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    P_Value = p_value
  ))
}
print(results_df)
write.csv(results_df, "E:/Shuai/Data/PA_association_results.csv", row.names = FALSE)



###################################################################################################################################################
# Test association of SSOE with dementia controlling for each protein, and pull out beta, se, 
# and p for SSOE or other PA pheno each time in loop (for each protein) 
####################################################################################################################################################
#get protein info
proteins<-fread("E:/Shuai/Data/coding143.tsv")
summary(proteins$coding)
table(proteins$coding)

index_of_column_1 <- which(names(m2) == "1")
index_of_column_2923 <- which(names(m2) == "2923")
print(index_of_column_1)
print(index_of_column_2923)

my.invnorm <- function(x) {
  res <- rank(x, na.last = "keep")
  res <- qnorm(res / (length(res) + 0.5))
  return(res)
}

OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$StrenorOther + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/SSOE with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')





OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$InvNorm_MVPA + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0 +  m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/MVPA with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')






OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$VigPA_Dich2 + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0  + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/VPA with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')










OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$Accel_Ave_DQexcl_minus46outliers + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0  + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/ACCAVE with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')











OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$InvNorm_Accel425 + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0  + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/ACC425 with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')











OUT <- as.data.frame(matrix(nrow=2911, ncol=8))
names(OUT) <- c("coding", "beta_pa", "se_pa", "pvalue_pa", "beta_protein", "se_protein", "pvalue_protein", "n")
OUT[,1] <- names(m2)[5814:8724] #change me if needed
for (i in 5814:8724) {          #change me here and below if needed
  normalized_protein <- my.invnorm(m2[[i]])
  r <- coxph(formula= Surv(m2$FU_time_Dem, m2$dementia) ~ m2$Acc_sed_minus46outliers + normalized_protein + m2$Season + as.factor(m2$race) + m2$f.21003.0.0 + m2$f.31.0.0  + m2$Education + m2$drinking + as.factor(m2$smoking) + m2$bodyfat)
  OUT[i-5813, 2] <- summary(r)$coeff[1,1]  
  OUT[i-5813, 3] <- summary(r)$coeff[1,2] 
  OUT[i-5813, 4] <- summary(r)$coeff[1,4]
  OUT[i-5813, 5] <- summary(r)$coeff[2,1]  
  OUT[i-5813, 6] <- summary(r)$coeff[2,2] 
  OUT[i-5813, 7] <- summary(r)$coeff[2,4]
  OUT[i-5813, 8] <- length(resid(r))
  print(i)
}
OUT$coding <- as.integer(OUT$coding)
OUT2 <- merge(proteins, OUT, by="coding")
summary(r)$coeff
summary(OUT2[,4])
summary(OUT2[,5])
min(OUT2[,4])
OUT2[which(OUT2$pvalue==0),]
write.table(OUT2, file="E:/Shuai/Data/ACCSED with dementia controlling for each proteinPA_association_results.txt", row.names=FALSE, quote=FALSE, sep='\t')





library(readr)
library(dplyr)
library(purrr)
library(haven)
if (requireNamespace("data.table", quietly = TRUE)) library(data.table)

options(stringsAsFactors = FALSE)

cols_take <- function(DT, cols) {
  if ("data.table" %in% class(DT)) DT[, ..cols] else DT[, cols, drop = FALSE]
}
sanitize <- function(v) gsub("\\.", "_", v)

files_pa <- list(
  MVPA   = "E:/Shuai/Data/8_14_MVPA_Model3.txt",
  VPA    = "E:/Shuai/Data/8_14_VPA_Model3.txt",
  SSOE   = "E:/Shuai/Data/8_14_SSOE_Model3.txt",
  ACC425 = "E:/Shuai/Data/8_14_ACC425_Model3.txt",
  ACCSED = "E:/Shuai/Data/8_14_ACCSED_Model3.txt",
  ACCAVE = "E:/Shuai/Data/8_14_ACCAVE_Model3.txt"
)
file_dem <- "E:/Shuai/Data/protein_to_dementia_scan.txt"     
out_dir  <- "E:/Shuai/Data/stata_bonf"                     
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dem <- read_tsv(file_dem, show_col_types = FALSE)
pcol_dem <- if ("pvalue" %in% names(dem)) "pvalue" else if ("p" %in% names(dem)) "p" else stop("No p column in dem")
ccol_dem <- if ("coding" %in% names(dem)) "coding" else stop("No coding column in dem")

dem2 <- dem %>%
  mutate(
    coding = as.integer(.data[[ccol_dem]]),
    p_dem  = as.numeric(.data[[pcol_dem]])
  ) %>%
  filter(!is.na(coding), !is.na(p_dem)) %>%
  mutate(
    p_bon_dem = p.adjust(p_dem, method = "bonferroni"),
    pass_dem  = p_bon_dem < 0.05
  ) %>%
  select(coding, p_dem, p_bon_dem, pass_dem, any_of(c("HR","beta","se","CI_L","CI_U"))) %>%
  group_by(coding) %>% slice_min(order_by = p_bon_dem, n = 1, with_ties = FALSE) %>% ungroup()

read_one_pa_bon <- function(pa_name, path) {
  tb <- read_tsv(path, show_col_types = FALSE)
  pcol <- if ("pvalue" %in% names(tb)) "pvalue" else if ("p" %in% names(tb)) "p" else stop(paste("No p column in", path))
  ccol <- if ("coding" %in% names(tb)) "coding" else stop(paste("No coding column in", path))
  tb %>%
    mutate(
      coding = as.integer(.data[[ccol]]),
      p_pa   = as.numeric(.data[[pcol]]),
      pa_var = pa_name
    ) %>%
    filter(!is.na(coding), !is.na(p_pa)) %>%
    group_by(pa_var) %>%
    mutate(
      p_bon_pa = p.adjust(p_pa, method = "bonferroni"),
      pass_pa  = p_bon_pa < 0.05
    ) %>%
    ungroup() %>%
    select(coding, pa_var, p_pa, p_bon_pa, pass_pa, any_of(c("beta","se","meaning",
                                                             "n","outlier_count","outlier_pvalue","heteroscedasticity_pvalue"))) %>%
    group_by(coding, pa_var) %>% slice_min(order_by = p_bon_pa, n = 1, with_ties = FALSE) %>% ungroup()
}

pa_long_bon <- purrr::imap_dfr(files_pa, ~ read_one_pa_bon(pa_name = .y, path = .x))

effect_filter <- function(df) {
  if ("HR" %in% names(df)) {
    df %>% mutate(abs_logHR = abs(log(HR))) %>% filter(abs_logHR >= 0.05)
  } else if ("beta" %in% names(df)) {
    df %>% mutate(abs_beta = abs(beta)) %>% filter(abs_beta >= 0.05)
  } else df
}

candidates_bon <- pa_long_bon %>%
  filter(pass_pa) %>%
  inner_join(dem2 %>% filter(pass_dem), by = "coding") %>%
  effect_filter()

summary_counts_bon <- candidates_bon %>%
  count(pa_var, name = "n_candidates") %>%
  arrange(desc(n_candidates))
print(summary_counts_bon)

write_tsv(candidates_bon, file.path(out_dir, "mediator_candidates_by_PA_BON.tsv"))
write_tsv(summary_counts_bon, file.path(out_dir, "mediator_candidates_counts_BON.tsv"))


prot_idx   <- 5814:8724
stopifnot(all(prot_idx %in% seq_along(names(m2))))
prot_names <- names(m2)[prot_idx]

coding_num <- suppressWarnings(as.integer(prot_names))
if (any(is.na(coding_num))) {
  idx_na <- which(is.na(coding_num)); coding_num[idx_na] <- seq_along(coding_num)[idx_na]
  warning("Some protein names are not integer-like; used running indices for those.")
}
protvar <- paste0("p_", coding_num)
if (anyDuplicated(protvar)) protvar <- make.unique(protvar)

missing_p <- setdiff(protvar, names(m2))
if (length(missing_p) > 0) {
  prot_norm <- setNames(
    as.data.frame(lapply(prot_idx, function(j) my.invnorm(m2[[j]])), check.names = FALSE),
    protvar
  )
  bind_cols_idx <- which(!(protvar %in% names(m2)))
  m2 <- dplyr::bind_cols(m2, prot_norm[bind_cols_idx])
}

protein_map <- tibble(coding = coding_num, protvar = protvar) %>% distinct()
if (exists("proteins")) {
  protein_map <- protein_map %>%
    left_join(mutate(proteins, coding = as.integer(coding)), by = "coding") %>%
    tidyr::separate_wider_delim(meaning, delim = ";", names = c("gene_symbol","description"), too_few = "align_start", cols_remove = FALSE)
}
write_tsv(protein_map, file.path(out_dir, "protein_map_with_meaning.tsv"))


base_vars_raw <- c("FU_time_Dem","dementia","Season","race",
                   "f.21003.0.0","f.31.0.0","Education","drinking",
                   "smoking","bodyfat","chronic_conditions")
stopifnot(all(base_vars_raw %in% names(m2)))

base_df <- cols_take(m2, base_vars_raw)
names(base_df) <- sanitize(names(base_df))
base_df$id <- seq_len(nrow(base_df))
base_df <- base_df %>% mutate(
  Season  = as.factor(Season),
  race    = as.factor(race),
  smoking = as.factor(smoking)
)

pa_vars_raw <- c(
  SSOE   = "StrenorOther",
  MVPA   = "InvNorm_MVPA",
  VPA    = "VigPA_Dich2",
  ACCAVE = "Accel_Ave_DQexcl_minus46outliers",
  ACC425 = "InvNorm_Accel425",
  ACCSED = "Acc_sed_minus46outliers"
)
stopifnot(all(unname(pa_vars_raw) %in% names(m2)))

make_stata_dataset_subset <- function(pa_label, pa_raw_name) {
  cand_codes <- candidates_bon %>%
    filter(pa_var == pa_label) %>%
    distinct(coding) %>% pull(coding)
  
  if (length(cand_codes) == 0) {
    message("No Bonferroni candidates for ", pa_label, " — skipping.")
    return(invisible(NULL))
  }
  
  protvars <- protein_map %>%
    filter(coding %in% cand_codes) %>%
    pull(protvar) %>% unique()
  
  missing_in_m2 <- setdiff(protvars, names(m2))
  if (length(missing_in_m2) > 0) stop("Missing p_XXXX in m2: ", paste(missing_in_m2, collapse=", "))
  
  pa_df   <- cols_take(m2, pa_raw_name); names(pa_df) <- sanitize(names(pa_df))
  prot_df <- cols_take(m2, protvars)
  
  out_df <- dplyr::bind_cols(base_df, pa_df, prot_df)
  out_path <- file.path(out_dir, paste0("mediation_input_", pa_label, "_BON_candidates.dta"))
  haven::write_dta(out_df, out_path)
  message("Wrote: ", out_path)
}

purrr::iwalk(pa_vars_raw, ~ make_stata_dataset_subset(pa_label = .y, pa_raw_name = .x))

message("All Bonferroni candidate datasets written to: ", out_dir)
































library(dplyr)
library(survival)
library(readr)
if (requireNamespace("data.table", quietly = TRUE)) library(data.table)

options(stringsAsFactors = FALSE)


cols_take <- function(DT, cols) {
  if ("data.table" %in% class(DT)) DT[, ..cols] else DT[, cols, drop = FALSE]
}

factorize_if_present <- function(df, cols) {
  for (v in cols) if (v %in% names(df)) df[[v]] <- as.factor(df[[v]])
  df
}

fit_cox <- function(df, time_var, event_var, pa_col, covars) {
  fml <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ",
           pa_col, " + ", paste(covars, collapse = " + "))
  )
  fit <- coxph(fml, data = df)
  s   <- summary(fit)
  hr  <- exp(coef(fit)[1])
  ci  <- exp(confint(fit)[1, ])
  data.frame(
    PA = pa_col,
    N = nrow(df),
    Events = s$nevent,
    HR = hr,
    CI_L = ci[1],
    CI_U = ci[2],
    P = s$coefficients[1, "Pr(>|z|)"],
    check.names = FALSE
  )
}


coding_to_cols <- function(m2, codes_integer) {
  nm <- names(m2)
  want_char  <- as.character(codes_integer)
  want_pchar <- paste0("p_", codes_integer)
  unique(c(intersect(want_char, nm), intersect(want_pchar, nm)))
}

prot_idx   <- 5814:8724
stopifnot(all(prot_idx %in% seq_along(names(m2))))

time_var <- "FU_time_Dem"
event_var <- "dementia"
covars <- c("Season","race","f.21003.0.0","f.31.0.0","Education","drinking","smoking","bodyfat","chronic_conditions")
pa_vars <- c(
  SSOE   = "StrenorOther",
  MVPA   = "InvNorm_MVPA",
  VPA    = "VigPA_Dich2",
  ACCAVE = "Accel_Ave_DQexcl_minus46outliers",
  ACC425 = "InvNorm_Accel425",
  ACCSED = "Acc_sed_minus46outliers"
)


candidates_path <- "E:/Shuai/Data/mediator_candidates_by_PA_q05.tsv"
candidates_path <- "E:/Shuai/Data/stata_bonf/mediator_candidates_by_PA_BON.tsv"

candidates_tbl <- read_tsv(candidates_path, show_col_types = FALSE) %>%
  mutate(coding = as.integer(coding))

m2 <- factorize_if_present(m2, c("Season","race","smoking"))

needed <- c(time_var, event_var, covars, unname(pa_vars))
missing_needed <- setdiff(needed, names(m2))
if (length(missing_needed) > 0) {
  stop("nop：", paste(missing_needed, collapse = ", "))
}


have_any_prot  <- rowSums(!is.na(cols_take(m2, prot_idx))) > 0

subset_A_cols <- c(time_var, event_var, covars, unname(pa_vars))
subset_A_base <- cols_take(m2, subset_A_cols)
subset_A_base <- subset_A_base[have_any_prot, , drop = FALSE]

keep_core <- c(time_var, event_var, covars)
subset_A_base <- subset_A_base %>%
  filter(if_all(all_of(keep_core), ~ !is.na(.)))

results_A <- lapply(unname(pa_vars), function(pa) {
  df <- subset_A_base %>% filter(!is.na(.data[[pa]]))
  if (nrow(df) == 0) return(NULL)
  fit_cox(df, time_var, event_var, pa, covars)
}) %>% bind_rows()

if (nrow(results_A) > 0) results_A$Subset <- "Proteomics panel (any protein)"

results_B_list <- list()

for (pa_lab in names(pa_vars)) {
  pa_col <- pa_vars[[pa_lab]]
  
  codes <- candidates_tbl %>%
    filter(pa_var == pa_lab) %>%
    distinct(coding) %>%
    pull(coding) %>% as.integer()
  
  if (length(codes) == 0) {
    message("No candidates for ", pa_lab, " — skip mediation-matched subset.")
    next
  }
  
  prot_cols_for_pa <- coding_to_cols(m2, codes)
  if (length(prot_cols_for_pa) == 0) {
    message("No matching protein columns in m2 for ", pa_lab, ". Check mapping (raw coding vs p_<coding>).")
    next
  }
  
  need_cols <- c(time_var, event_var, covars, pa_col, prot_cols_for_pa)
  dfB <- cols_take(m2, need_cols) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  if (nrow(dfB) == 0) {
    message("No complete cases for ", pa_lab, " after requiring all its candidate proteins.")
    next
  }
  
  resB <- fit_cox(dfB, time_var, event_var, pa_col, covars)
  resB$Subset <- paste0("Mediation-matched (candidates for ", pa_lab, ")")
  results_B_list[[pa_lab]] <- resB
}

results_B <- bind_rows(results_B_list)

results_all <- bind_rows(results_A, results_B) %>%
  mutate(PA_readable = names(pa_vars)[match(PA, unname(pa_vars))]) %>%
  relocate(PA_readable, Subset, N, Events, HR, CI_L, CI_U, P)

out_path <- "E:/Shuai/Data/PA_to_dementia_cox_in_proteomics_subsets.tsv"
write_tsv(results_all, out_path)
print(results_all)
message("Wrote: ", out_path)











