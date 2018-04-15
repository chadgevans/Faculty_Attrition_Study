Faculty Attrition
================
Chad Evans

Built with R version 3.4.4.

Contents
--------

-   [Configure](#configure)
    -   directories
    -   libraries
    -   data
    -   helper functions
-   [Munge](#munge)
-   [Life Tables](#life-tables)
    -   Life table of all faculty, aggregated
    -   Life table of faculty, by Tenure Status
-   [Graphing Survival Curves of Faculty](#graphing-survival-curves-of-faculty)
    -   Life Table Estimator
    -   Kaplan-Meier Estimator
-   [Cox Proportional Hazards Models](#cox-proportional-hazards-models)
    -   [Naive Model](#naive-model)
    -   [Controls Model](#controls-model)
    -   [Mediators Model](#mediators-model)
    -   [Full Model](#full-model)
-   [Tabulation of All Models](#tabulation-of-all-models)
-   [Conclusions](#conclusions)

Configure
---------

``` r
library(knitr)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(plyr)
library(dplyr)
library(survival)
library(KMsurv)
library(car)
library(RColorBrewer)
library(survminer)
library(icenReg)
library(readxl)
```

Munge
-----

``` r
sdata$EntryDate<-as.numeric(sdata$EntryDate)
sdata$EndDate<-as.numeric(sdata$EndDate)
sdata$lower<-sdata$EndDate-sdata$EntryDate
sdata$upper<-sdata$ExitDate-sdata$EntryDate
sdata$upper[sdata$Censor==0]<-NA
sdata$TIME2<-sdata$EndDate-sdata$EntryDate+1 # we add onen here for the same reason (but time is different calc.)

data<-sdata %>%
  filter(!(EntryTENSTA=="Tenured")) %>% # Remove those few who start off with tenure
  droplevels()
data$TIME<-data$ExitDate-data$EntryDate
data$INTERVAL = floor(data$TIME/2)
  
data$DEG2ENTRY<-as.numeric(data$EntryDate)-data$DGRYR
data$EntryEMTP<-relevel(data$EntryEMTP,ref = "Four-year")
data$EntryWAPRI<-relevel(data$EntryWAPRI,ref = "Teaching")
data$GENDER<-relevel(data$GENDER,ref = "Male")
data$EntryEMTP[data$EntryEMTP=="Other Ed"]<-'Two-year'
data$EntryEMTP<-droplevels(data$EntryEMTP)
vars<-data.frame(model.matrix( ~ EntryTENSTA - 1, data=data )) # Need dummy coded for coxph tt()
data$NTT<-vars$EntryTENSTANTT
data$NTS<-vars$EntryTENSTANon.tenure.System.Position
```

``` r
length(unique(data$REFID)) # The number of individuals who participated in the study
```

Life Tables
-----------

First, lets write a function to create our life tables.

``` r
intervals<-2 # number units per intervals for life table
LifeTableIt <- function(data){
  ltab.data<-data %>% 
  mutate(interval = floor(TIME/intervals)) %>%
  select(interval,Censor) %>%
  group_by(interval) %>%
  dplyr::summarise(attrit=sum(Censor), count=n()) %>%
  mutate(nlost=count-attrit)
  int<-c(ltab.data$interval,NA) #length is 1 + nlost and nevent 
  lifetable<-round(with(ltab.data, lifetab(tis=int, ninit=nrow(data), nlost=nlost, nevent=attrit)),3)
return(lifetable)
}
```

### Life Table of All Faculty, Aggregated

This first Life Table will tabulate the survival probabilities and hazards for all faculty (aggregated).

``` r
ltab<-LifeTableIt(data)
kable(ltab)
```

|       |  nsubs|  nlost|    nrisk|  nevent|   surv|    pdf|  hazard|  se.surv|  se.pdf|  se.hazard|
|-------|------:|------:|--------:|-------:|------:|------:|-------:|--------:|-------:|----------:|
| 1-2   |  21436|   8051|  17410.5|    2025|  1.000|  0.116|   0.123|    0.000|   0.002|      0.003|
| 2-3   |  11360|   3550|   9585.0|     770|  0.884|  0.071|   0.084|    0.002|   0.002|      0.003|
| 3-4   |   7040|   2067|   6006.5|     329|  0.813|  0.045|   0.056|    0.003|   0.002|      0.003|
| 4-5   |   4644|   1304|   3992.0|     139|  0.768|  0.027|   0.035|    0.004|   0.002|      0.003|
| 5-6   |   3201|    365|   3018.5|      72|  0.741|  0.018|   0.024|    0.004|   0.002|      0.003|
| 6-7   |   2764|    742|   2393.0|      48|  0.724|  0.015|   0.020|    0.005|   0.002|      0.003|
| 7-8   |   1974|    630|   1659.0|      20|  0.709|  0.009|   0.012|    0.005|   0.002|      0.003|
| 8-9   |   1324|    470|   1089.0|      12|  0.701|  0.008|   0.011|    0.005|   0.002|      0.003|
| 9-10  |    842|    314|    685.0|       4|  0.693|  0.004|   0.006|    0.006|   0.002|      0.003|
| 10-11 |    524|    269|    389.5|       3|  0.689|  0.005|   0.008|    0.006|   0.003|      0.004|
| 11-NA |    252|    252|    126.0|       0|  0.684|     NA|      NA|    0.007|      NA|         NA|

By the last interval, over 30 percent of faculty have attrited from their careers in academia. The hazard is highest during the early stages of the careers and then plateaus later on.

### Life Table for Faculty, by Tenure Status

Now, let's tabulate lifetables for each of the tenure status groups and extract survival probabilities.

``` r
LT_TT<-data %>% 
  filter(EntryTENSTA=="Tenure-Track") %>%
  LifeTableIt()

LT_NTT<-data %>% 
  filter(EntryTENSTA=="NTT") %>%
  LifeTableIt()

LT_NTS<-data %>% 
  filter(EntryTENSTA=="Non-tenure System/Position") %>%
  LifeTableIt()

LTs<-data.frame(LT_TT[,c("surv","se.surv")], LT_NTT[,c("surv","se.surv")], LT_NTS[,c("surv","se.surv")])
names(LTs)<-c("TT Surv","TT SE","NTT Surv","NTT SE","NTS Surv","NTS SE")
LTs
```

    ##       TT Surv TT SE NTT Surv NTT SE NTS Surv NTS SE
    ## 1-2     1.000 0.000    1.000  0.000    1.000  0.000
    ## 2-3     0.969 0.002    0.878  0.006    0.830  0.004
    ## 3-4     0.935 0.004    0.807  0.008    0.731  0.005
    ## 4-5     0.908 0.005    0.752  0.010    0.675  0.006
    ## 5-6     0.889 0.006    0.725  0.011    0.641  0.007
    ## 6-7     0.878 0.006    0.706  0.011    0.617  0.008
    ## 7-8     0.868 0.007    0.689  0.012    0.599  0.008
    ## 8-9     0.862 0.007    0.677  0.013    0.590  0.009
    ## 9-10    0.857 0.008    0.665  0.014    0.582  0.009
    ## 10-11   0.857 0.008    0.660  0.015    0.574  0.010
    ## 11-NA   0.857 0.008    0.650  0.017    0.565  0.012

``` r
#names(LTs)<-c("Tenure/Track","TT SE","Non_Tenure/Track","NTT SE","No Tenure System","NTS SE")
```

``` r
write.csv(LTs, file.path(Graphs,"Life_table.csv"))
```

Clearly, the survival probabilities cascade much more rapidly for NTT and NTS faculty.

### Graphing Survival Curves of Faculty

#### Life Table (Actuarial) Estimator

The following graph plots the actuarial survival probabilities for each class of faculty.

``` r
TYPE<-c(rep("Tenure-Track",11), rep("Non-Tenure Track",11), rep("Non-Tenure System",11))
SURV<-c(LTs$`TT Surv`,LTs$`NTT Surv`,LTs$`NTS Surv`)
SE<-c(LTs$`TT SE`,LTs$`NTT SE`,LTs$`NTS SE`)
INTERVAL<-rep(c(0:10),3)
TIME<-rep(2*c(0:10),3)
ymin<-SURV-2*SE
ymax<-SURV+2*SE
ltab2<-data.frame(TYPE,SURV,TIME,ymin,ymax)

ggplot(ltab2, aes(x=TIME,y=SURV, group=rev(TYPE), colour=rev(TYPE))) + 
  ylim(.5,1) +
  geom_step(size=1) +
  labs(title="Probability of Attrition from Academia", subtitle= "Life Table (Actuarial) Estimator", x="Years Since First Academic Job", y = "Probability of Remaining in Academia") +
 scale_colour_manual(name="Tenure Status", labels=c("Tenure-Track", "Non-tenure Track","No Tenure System"), values=c("#E7B800", "#2E9FDF", "blue")) +
  theme_bw() +
  geom_rect(aes(xmin = TIME, xmax = TIME+2, ymin = ymin, ymax = ymax), alpha = 0.1)
```

![](graphs/Life_table_graph-1.png)

#### Kaplan-Meier (Product Limit) Estimator

These survival probabilities can also be estimated and graphed using the KM estimator. In this chapter, I will be reporting only the actuarial estimator, as it is better computationally for larger samples.

``` r
km.tenure <- survfit( Surv(TIME, Censor)~ strata(EntryTENSTA), data=data, type="kaplan-meier")
ggsurvplot(km.tenure, data = data, ylim = c(.5,1), size = 1, palette = c("#E7B800", "#2E9FDF", "blue"), conf.int = TRUE, legend.labs = c("Tenure-Track", "Non-Tenure track","No Tenure System"),  ggtheme = theme_bw(),  ncensor.plot.height = 0.25, title="Probability of Remaining in Academia: by Initial Tenure Status", subtitle="Kaplan-Meier Estimator", xlab="Years Since First Academic Job", ylab="Probability of Remaining in Academia", legend="right", legend.title="Tenure Status",ncensor.plot =F )
```

![](graphs/KM_Graph-1.png)

The default confidence interval (above) is log, which calculates intervals based on log(survival).

### Statistical Test for differences

Often it is useful to have a statistical test for whether survival curves of two (or more) groups differ. Let's look at the most common test, the log-rank test (a.k.a. the Mantel-Haenszel test).

``` r
survdiff(Surv(TIME, Censor) ~ EntryTENSTA, data=data,rho=0) # log-rank or Mantel-Haenszel test
```

    ## Call:
    ## survdiff(formula = Surv(TIME, Censor) ~ EntryTENSTA, data = data, 
    ##     rho = 0)
    ## 
    ##                                            N Observed Expected (O-E)^2/E
    ## EntryTENSTA=Tenure-Track                6667      453     1226    487.70
    ## EntryTENSTA=NTT                         3942      653      609      3.12
    ## EntryTENSTA=Non-tenure System/Position 10827     2316     1586    335.74
    ##                                        (O-E)^2/V
    ## EntryTENSTA=Tenure-Track                     808
    ## EntryTENSTA=NTT                                4
    ## EntryTENSTA=Non-tenure System/Position       665
    ## 
    ##  Chisq= 880  on 2 degrees of freedom, p= 0

The null hypothesis for a log-rank test is that the groups have the same survival. In this case, the Chi-square statistic is significant for tenure status of first appointment. The survival curve, therefore, differs based on the tenure status of the PhD recipient's first academic job.

#### Summary of Nonparametric models

Substantially, the actuarial estimator and KM estimator tell similar stories. Attrition is much higher and more rapid for faculty working off the tenure track. The lifetable estimator perhaps shows a slightly higher attrition rate for NTS faculty, compared to the KM estimate.

Regression Models for Survival Data
-----------------------------------

There are two types of regression models in survival analysis, semi-parametric and fully-parametric. During the exploratory phase, this study considered both approaches. However, I ultimately opted for the semiparametric method of cox proportional hazards because of a similar fit and less stringent assumptions.

The final models for this chapter are right censored Cox proportional hazards model. It is true that the data are also interval censored, however, the intervals are very close to being of equal length. As a result, we can treat this as a discrete time analyisis with right censored data. The interval censored functions were poorly developed and did not allow for interactions with time. So this also encouraged the use of a right censoring function which allowed for interactions with time. This interaction effectively allowed me to use the cox proportional hazards model even though the hazards of tenure-track and non-tenure track faculty are not proportional at one or more time points. This is particularly important during the years of tenure review, as the hazard of attrition likely shifts significantly for those faculty on the tenure line. To allow for non-proportional hazards, the final model interacts time with tenure status. As Allison says, "Testing is both the diagnosis and the solution."

### Specification

For specifying our models, we will subset the data and first train models using training data. After we are satisfied with that specification, we'll test all models using the independent test data. This will help prevent overfitting.

``` r
smp_size <- floor(0.60 * nrow(data))
set.seed(777)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
```

Cox Proportional Hazards Models
-------------------------------

First, we fit a simple Cox model predicting attrition from academia from a categorical predictor of tenure status at workforce entry. We employ the efron method of dealing with ties, although other popular methods (e.g., breslow method) are available. The Efron approximation is more accurate when dealing with tied death times, and is as efficient computationally.

The Surv function here frames censoring as follows. It assumes that individuals entered the study halfway through the previous interval. The extra year assumption was necessary to avoid "spontaneous attrition" that would have removed from the analysis any indiviudal who entered and exited during the same interval. As these early attriters are so very important to this study, we assume they entered halfway into the previous interval. This method treats censoring from the "EndDate"--that is, the last year he or she was still employed. Censoring, of course, happens after that point, so an individual who enters in year 0, remains in the study through year 2 in its entirety and subsequently drops out between year 2 and 4 would be coded as "3+." There are two years for the full interval, plus the extra entry year (assumption).

Naive Model
-----------

``` r
#Just the naive model
cox1 <- coxph( Surv(TIME2, Censor) ~ NTT + NTS, data = test,  method="efron", robust = TRUE)
```

Proportional Hazards models and AFT models must be interpreted in different ways. AFT models give the percentage change in survival time. Cox PH models give the percentage change to the hazard at all time points, following this formula: *h*(*t*)=*h*<sub>0</sub>(*t*)*e*<sup>(*β*′*x*)</sup>

In this case, the effect of initial tenure status on time to attrition has an estimated coefficient of 1.0089259 and 1.3237404. Exponentiated, this means that subjects appointed to lower tenure-status jobs multiply their baseline hazards *h*<sub>0</sub>(*t*) by a factor of 2.7426537 and 3.7574493. Their "risk" of attriting from academia is 174.2653671 and 275.744933 percent higher than academics who begin immediately on the tenure-track. Importantly, Cox models state that this is the impact on the subject's hazard at any given time, t. It does not, however, imply an expansion (or contraction) of the lifespan of the subject (as AFT models do).

Controls Model
--------------

Now we introduce controls to the model.

``` r
cox2 <- coxph( Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY +
    EntryWAPRI + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    EntryCTZUSIN, data = test,  method="efron", robust = TRUE)
```

Mediators Model
---------------

Now we introduce mediators to the model.

``` r
#Model with mediators
cox3 <- coxph( Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY +
    EntryWAPRI + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    EntryCTZUSIN + EntryWKTRNI + LogSalary, data = test,  method="efron", robust = TRUE)
```

An important limitation of these models up to this point is that they do not account for the possibility of non-proportional hazards based on tenure status. Let's test whether the proportional hazards assumption may have been violated in this analysis.

``` r
cox.zph <- cox.zph(cox3, transform="km")
round(cox.zph$table,3)
```

    ##                                    rho  chisq     p
    ## NTT                             -0.076  6.799 0.009
    ## NTS                             -0.127 17.073 0.000
    ## DEG2ENTRY                        0.016  0.347 0.556
    ## EntryWAPRIOther                 -0.053  3.161 0.075
    ## EntryWAPRIResearch              -0.073  5.763 0.016
    ## EntryPUBPRIPrivate              -0.049  3.012 0.083
    ## EntryEMTPTwo-year                0.012  0.198 0.656
    ## EntryEMTPMed                     0.139 23.255 0.000
    ## EntryEMTPUni Research Institute -0.010  0.124 0.725
    ## SDRCARNR2                       -0.028  0.952 0.329
    ## SDRCARNDoctorate                -0.018  0.384 0.535
    ## SDRCARNOther                     0.005  0.038 0.846
    ## SDRCARNMedHealth                 0.062  4.203 0.040
    ## EntryAGE                         0.034  1.407 0.236
    ## GENDERFemale                     0.029  1.045 0.307
    ## MINRTYYes                       -0.056  3.881 0.049
    ## EntryMARINDYes                  -0.047  2.796 0.094
    ## EntryCHLVINYes                  -0.034  1.405 0.236
    ## EntryCTZUSINCitizen              0.012  0.170 0.680
    ## EntryWKTRNITraining             -0.010  0.124 0.725
    ## LogSalary                       -0.072  4.287 0.038
    ## GLOBAL                              NA 73.672 0.000

This function tests proportionality by interacting each predictor with log time (km transformation). Rho is the pearson product-moment correlation between the scaled residuals (Schoenfeld) and log(time) for each covariate. The global test jointly tests all of the interactions. Low p-values suggest a violation of the assumption of proportional hazards.

Full Model
----------

Now let's build a final, comprehensive model. Importantly we'll include an interaction between time and our tenure status variables (NTT or NTS dummies) to test the assumption of proportional hazards and also correct for non-proportional hazards. We'll also include the interaction between institution type and control.

``` r
cox4 <- coxph( Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY + EntryWAPRI + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN + EntryWKTRNI + LogSalary + tt(NTT) + tt(NTS) + EntryPUBPRI*EntryEMTP, data = test,  method="efron", robust = TRUE, tt = function(x, t, ...) x * log(t))
```

Even after controlling for background characteristics, there are significant differences. Here, NTT status or working at a college or university without a tenure system impacts the hazard, multiplying the baseline by a factor of 2.7741545 and 3.352443. This is equivalent to saying that each tenure status increases the hazard of attrition by 177.4154456 and 235.2443049 percent, controlling for background characteristics. R output also provides the exponentiated negative coefficient. To my understanding, that just allows you to compare the groups relative to the baseline hazard of the tenure-track group. Robust standard errors were used in this model.

The model reveals other important predictors, including the subject's main job and possible interactions between institution type and public/private status.

Tabulation of All models
------------------------

``` r
Mod1<-rbind(summary(cox1)$coefficients[,c(2,5,6)], fill_spaces(24,3))
Mod2<-rbind(summary(cox2)$coefficients[,c(2,5,6)], fill_spaces(7,3))
Mod3<-rbind(summary(cox3)$coefficients[,c(2,5,6)], fill_spaces(5,3))
Mod4<-summary(cox4)$coefficients[,c(2,5,6)]

Final_Table<-cbind(Mod1, Mod2, Mod3, Mod4)
rownames(Final_Table)<-rownames(summary(cox4)$coef)
colnames(Final_Table)<-c("Mod1_HR","Mod1_z","Mod1_pval","Mod2_HR","Mod2_z","Mod2_pval","Mod3_HR","Mod3_z","Mod3_pval","Mod4_HR","Mod4_z","Mod4_pval")
kable(Final_Table)
```

|                                                    |  Mod1\_HR|   Mod1\_z|  Mod1\_pval|   Mod2\_HR|     Mod2\_z|  Mod2\_pval|   Mod3\_HR|     Mod3\_z|  Mod3\_pval|   Mod4\_HR|     Mod4\_z|  Mod4\_pval|
|----------------------------------------------------|---------:|---------:|-----------:|----------:|-----------:|-----------:|----------:|-----------:|-----------:|----------:|-----------:|-----------:|
| NTT                                                |  2.742654|  10.62914|           0|  2.3565254|   8.6429732|   0.0000000|  2.3000975|   8.3044999|   0.0000000|  2.7741545|   7.2893383|   0.0000000|
| NTS                                                |  3.757449|  16.77793|           0|  2.6151224|  10.9462626|   0.0000000|  2.5368381|  10.3990629|   0.0000000|  3.3524430|   9.7151586|   0.0000000|
| DEG2ENTRY                                          |        NA|        NA|          NA|  0.9780665|  -1.2648273|   0.2059332|  0.9812951|  -1.0704249|   0.2844281|  0.9820066|  -1.0298818|   0.3030655|
| EntryWAPRIOther                                    |        NA|        NA|          NA|  1.6653979|   4.9945942|   0.0000006|  1.6869666|   5.1242022|   0.0000003|  1.6507824|   4.8980905|   0.0000010|
| EntryWAPRIResearch                                 |        NA|        NA|          NA|  1.5743821|   5.2155889|   0.0000002|  1.5867156|   5.2869869|   0.0000001|  1.5585873|   5.0868146|   0.0000004|
| EntryPUBPRIPrivate                                 |        NA|        NA|          NA|  1.0597033|   0.9490728|   0.3425836|  1.0659701|   1.0459951|   0.2955633|  0.9477651|  -0.6118430|   0.5406416|
| EntryEMTPTwo-year                                  |        NA|        NA|          NA|  1.0385268|   0.1853365|   0.8529651|  1.0214819|   0.1035308|   0.9175417|  0.8319744|  -0.7733996|   0.4392859|
| EntryEMTPMed                                       |        NA|        NA|          NA|  1.2152175|   2.6587894|   0.0078422|  1.2146126|   2.6517060|   0.0080086|  1.1397273|   1.3440400|   0.1789354|
| EntryEMTPUni Research Institute                    |        NA|        NA|          NA|  1.1714475|   1.8315566|   0.0670175|  1.1792170|   1.9053315|   0.0567370|  1.0589085|   0.5374091|   0.5909850|
| SDRCARNR2                                          |        NA|        NA|          NA|  1.1470554|   1.4194103|   0.1557794|  1.1477243|   1.4260680|   0.1538487|  1.1498291|   1.4417399|   0.1493758|
| SDRCARNDoctorate                                   |        NA|        NA|          NA|  1.1036693|   0.9943637|   0.3200458|  1.1023669|   0.9817950|   0.3262009|  1.1051577|   1.0194901|   0.3079703|
| SDRCARNOther                                       |        NA|        NA|          NA|  1.2640600|   1.1339058|   0.2568341|  1.2387926|   1.0217128|   0.3069169|  1.2370769|   1.0158110|   0.3097194|
| SDRCARNMedHealth                                   |        NA|        NA|          NA|  1.0764915|   0.5465142|   0.5847126|  1.0707195|   0.5074136|   0.6118647|  1.0728198|   0.4951720|   0.6204787|
| EntryAGE                                           |        NA|        NA|          NA|  0.9956296|  -0.8583280|   0.3907114|  0.9956620|  -0.8501841|   0.3952228|  0.9958946|  -0.8128147|   0.4163243|
| GENDERFemale                                       |        NA|        NA|          NA|  0.9264368|  -1.2869294|   0.1981189|  0.9184368|  -1.4271422|   0.1535389|  0.9181464|  -1.4340915|   0.1515461|
| MINRTYYes                                          |        NA|        NA|          NA|  1.0252689|   0.3220834|   0.7473895|  1.0268595|   0.3425822|   0.7319128|  1.0296921|   0.3813277|   0.7029601|
| EntryMARINDYes                                     |        NA|        NA|          NA|  0.9172141|  -1.2323093|   0.2178336|  0.9156365|  -1.2556640|   0.2092378|  0.9120943|  -1.3077190|   0.1909686|
| EntryCHLVINYes                                     |        NA|        NA|          NA|  0.9905358|  -0.1339133|   0.8934711|  0.9884834|  -0.1629914|   0.8705252|  0.9919569|  -0.1141935|   0.9090844|
| EntryCTZUSINCitizen                                |        NA|        NA|          NA|  1.0570681|   0.8085050|   0.4187999|  1.0555757|   0.7879393|   0.4307322|  1.0587435|   0.8331092|   0.4047832|
| EntryWKTRNITraining                                |        NA|        NA|          NA|         NA|          NA|          NA|  1.0260713|   0.4355781|   0.6631429|  1.0272647|   0.4550056|   0.6491052|
| LogSalary                                          |        NA|        NA|          NA|         NA|          NA|          NA|  0.9291256|  -2.2112916|   0.0270157|  0.9282446|  -2.1343236|   0.0328163|
| tt(NTT)                                            |        NA|        NA|          NA|         NA|          NA|          NA|         NA|          NA|          NA|  0.8229292|  -1.7121377|   0.0868713|
| tt(NTS)                                            |        NA|        NA|          NA|         NA|          NA|          NA|         NA|          NA|          NA|  0.7083724|  -3.5035737|   0.0004591|
| EntryPUBPRIPrivate:EntryEMTPTwo-year               |        NA|        NA|          NA|         NA|          NA|          NA|         NA|          NA|          NA|  2.9717977|   2.2596388|   0.0238437|
| EntryPUBPRIPrivate:EntryEMTPMed                    |        NA|        NA|          NA|         NA|          NA|          NA|         NA|          NA|          NA|  1.1791694|   1.1931475|   0.2328116|
| EntryPUBPRIPrivate:EntryEMTPUni Research Institute |        NA|        NA|          NA|         NA|          NA|          NA|         NA|          NA|          NA|  1.3507137|   1.7709123|   0.0765753|

``` r
write.csv(Final_Table, file.path(Graphs,"Final_Table.csv"))
```

Conclusions
-----------

1.  Tenure status at hiring related to attrition
2.  Attrition more prevalent among NTT and NTS faculty, however, large numbers make a career of it
3.  Attrition is more of a risk for researchers and administrators, not teachers
4.  It might be useful to look into private, two-year institutions
5.  Relationships likely underfit by the model. Acquire more faculty characteristics, particularly time-varying characteristics.
