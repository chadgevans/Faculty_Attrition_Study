# Need to create a data frame with
# ID, entry and end date (this allows the calculation of time), age (at entry), treatment, censored (1 means attrits, 0 means not) 

# Institution Level
lev<-paste(levels(df$EDTP), collapse = "','")
lev<-gsub("\\d: ","", lev)
levels(df$EDTP)<-c('NonPostsecondary','NonPostsecondary','2-year','4-year','Medical school/center','University research institute','Other', NA,'NonPostsecondary','2-year','4-year','Medical school/center','University research institute','Other','NonPostsecondary','2-year','4-year','Medical school/center','Other educational institution')

levels(df$EMED)<-c(NA,"No","Yes")
df$EDTP[df$EMED=="No"]<-'NonPostsecondary'  # Using EMED to fill in those working outside academia in the EDTP vector

levels(df$CARNEG)<-c(rep("ResDoc",4), rep("AsscBaMa", 5), rep("Other", 10),NA,NA)
levels(df$PUBPRI)<-c("Public","Private",NA,NA)

#Job Characteristics
levels(df$TENSTA)<-c("Tenured","Tenure-Track","NTT",NA,rep("Non-tenure System/Position",4))
df$SALARY[df$SALARY==9999998]<-NA #also there's a df$SALARY[df$SALARY==9999996]<-NA This is perhaps a max value, but no record of it in codebook
levels(df$WAPRI)<-c("Other",rep("Research",2),rep("Other",9),"Teaching","Other",NA,rep("Research",2), rep("Other",4), NA,"Other") # What o you spend the most hours doing?

#Demographic Characteristics
lev<-paste(levels(df$MARSTA), collapse = "','")
lev<-gsub("\\d: ","", lev)
levels(df$MARSTA)<-c('Married',rep('Not Married',4),'Married',rep('Not Married',4))
levels(df$GENDER)<-c("Female","Male")
levels(df$MINRTY)=c("No","Yes")
levels(df$CHLVIN)=c("No","Yes") # any children living with R (can change over time)

#Sample subjects to save computational time
set.seed(777)
inds<-unique(df$REFID) # Id all the subjects in all years
sel<-sample(inds, 1000)
df<-df[df$REFID %in% sel,]

#Build the survival dataframe

f<-df %>%
  filter(DGRYR>1989) %>%
  select(REFID,YEAR,EDTP,TENSTA) %>%
  gather(Key, Value, -REFID, -YEAR) %>%
  unite(temp, Key, YEAR, sep = ".") %>% 
  spread(temp, Value)

# Determining Entry, Exit and Censoring
sdata<-f["REFID"] # constructor to hold results
d<-f %>%
  select(starts_with("EDTP")) %>%
  `colnames<-`(c(1993,1995,1997,1999,2001,2003,2006,2008,2010,2013)) %>%
  as_data_frame()

# Entry Date
sdata$EntryDate<-rep(NA, nrow(d))
for(j in rev(1:ncol(d))){
  for(i in 1:nrow(d)){
    if(!is.na(d[i,j]) & !(d[i,j] %in% c("NonPostsecondary")) ){
      sdata$EntryDate[i]<-names(d[j])
    }
  }
}

#End date
EDINST<-c("2-year","4-year","Medical school/center","University research institute","Other","Other educational institution")
sdata$EndDate<-rep(NA, nrow(d))
for(i in 1:nrow(d)){
  for(j in 1:ncol(d)){
    if(is.na(sdata$EntryDate[i])){
      print("No Entry Date")
      break
    }
    else if(names(d[i,j])==sdata$EntryDate[i]){
      print(paste("Subject ID",i,"enters study in",sdata$EntryDate[i]), sep=" ")
      if(j==10){
        print("(and exits 2013 as well)")
        sdata$EndDate[i]<-names(d[i,j]) 
        break
      }
      while(d[i,j+1] %in% EDINST){
        print(paste("Study Status in",names(d[i,j+1]), ": In Study"))
        j=j+1
        if(j==10){
          break
        }
      }
      print(paste("Exit Study: ", names(d[i,j])))
      sdata$EndDate[i]<-names(d[i,j]) # This is the last date in the study (not j+1)
    }
    else{
      print("Not in Study")
    }
  }
}

#Censoring
# 1: quits Academia for job elsewhere.
# 0: a) never attrited 
# or b) was lost to follow-up.
sdata$Censor<-rep(NA, nrow(d))
EDINST<-c("2-year","4-year","Medical school/center","University research institute","Other","Other educational institution")

for(i in 1:nrow(d)){
  if(is.na(sdata$EntryDate[i])){
    print(paste("No Entry Date for subject",i,sep = " "))
    next # 
  }
  for(j in 1:ncol(d)){
    if(names(d[i,j])==sdata$EntryDate[i]){
      print(paste("Subject ID",i,"enters study in",sdata$EntryDate[i]), sep=" ")
      if(j==10){
        print("(and exits 2013 as well)")
        sdata$Censor[i]<-0 
        break
      }
      while(d[i,j+1] %in% EDINST){
        print(paste("Subject", i, "in study in", names(d[i,j+1]), sep = " "))
        j=j+1
        if(j==10){
          sdata$Censor[i]<-0
          break
        }
      }
      print(paste("Subject",i, "exits Study in", names(d[i,j]),sep = " "))
      if(j==10){
        break
      }else if(is.na(d[i,j+1])){
        sdata$Censor[i]<-0
      } else{
        sdata$Censor[i]<-1
      }
    }
    else{
      print(paste(names(d[i,j]),"not an entry year", sep = " "))
    }
  }
}


# Continuous Variables
contvars<-c("EntrySALARY","EntryAGE")
for(i in 1:length(contvars)){
  sdata[contvars]<-as.numeric(rep(NA, nrow(sdata)))
}

# Including Time-Constant variables in survival data frame
factvars<-c("EntryTENSTA","EntryMARSTA","EntryCHLVIN","EntryWAPRI","EntryCARNEG","EntryPUBPRI")
temp0<-rep(NA, length(factvars))
for(i in 1:length(factvars)){
  sdata[factvars[i]]<-factor(rep(NA, nrow(sdata)))
  temp0[i]<-gsub("Entry","",factvars[i])
  levels(sdata[,factvars[i]])<-levels(df[,temp0[i]])
}

# I don't know why a double loop doesn't work here below
for(i in 1:nrow(sdata)){
  temp<-df[which(sdata$REFID[i]==df$REFID),]
  sdata$EntryTENSTA[i]<-temp$TENSTA[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryMARSTA[i]<-temp$MARSTA[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryCHLVIN[i]<-temp$CHLVIN[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryWAPRI[i]<-temp$WAPRI[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryCARNEG[i]<-temp$CARNEG[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryPUBPRI[i]<-temp$PUBPRI[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntryAGE[i]<-temp$AGE[sdata$EntryDate[i]==temp$YEAR]
  sdata$EntrySALARY[i]<-temp$SALARY[sdata$EntryDate[i]==temp$YEAR]
}

# Time Constant characteristics
sdata$GENDER<-factor(rep(NA, nrow(sdata))) 
levels(sdata$GENDER)<-levels(df$GENDER)
sdata$MINRTY<-factor(rep(NA, nrow(sdata))) 
levels(sdata$MINRTY)<-levels(df$MINRTY)
for(i in 1:nrow(sdata)){
  temp<-df[which(sdata$REFID[i]==df$REFID),]
  sdata$GENDER[i]<-Mode(temp$GENDER)
  sdata$MINRTY[i]<-Mode(temp$MINRTY)
}


