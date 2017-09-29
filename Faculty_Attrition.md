Faculty Attrition
================
Chad Evans

Built with R version 3.3.2.

Contents
--------

-   [Configure](#configure)
    -   directories
    -   libraries
    -   data
    -   helper functions
-   [Munge](#munge)
-   [Life Tables](#life-table-method)
    -   Aggregated Life table
    -   Life table by Tenure Status
-   [Graphing Survival Curves of Faculty](#graphing-survival-curves-of-faculty)
    -   Life Table Estimator
    -   Kaplan-Meier Estimator
-   [Cox Proportional Hazards Model](#cox-proportional-hazards-model)
    -   [Right Censored Data Model](#right-censored-data-model)
    -   [Interval Censored Data Model](#interval-censored-data-model)
    -   [Comparing Right Censored Model and Interval Censored Model](#comparing-right-censored-model-and-interval-censored-model)
-   [Accelerated Failure Time Models](#accelerated-failure-time-models)
    -   Exponential
    -   Lognormal
    -   [AFT Specification](#aft-specification)
-   [Comparing the Cox model and the Exponential AFT model](#comparing-the-cox-model-and-the-exponential-aft-model)
-   [Extensions](#extensions)
    -   Ridge Regression
    -   Smoothing Splines
    -   Frailty Models

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
sdata$ExitDate2<-sdata$ExitDate
sdata$ExitDate2[sdata$ExitDate==2015]<-NA
data<-sdata %>%
  filter(!(EntryTENSTA=="Tenured")) %>% # Remove those few who start off with tenure
  droplevels() %>%
  mutate(TIME=as.numeric(ExitDate)-as.numeric(EntryDate))
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
# data<-na.omit(data)
length(unique(data$REFID))  # The number of individuals who participated in the study
```

Life Tables
-----------

This is also known as the actuarial method. For the life table method, if observations are censored on the same month (or time unit) that events occurred, they are assumed to be at risk for just half the month.

First let's write a function that takes the data as its argument and generates a lifetable.

``` r
intervals <- 2  # number units per intervals for life table
LifeTableIt <- function(data) {
    ltab.data <- data %>% mutate(interval = floor(TIME/intervals)) %>% select(interval, 
        Censor) %>% group_by(interval) %>% dplyr::summarise(attrit = sum(Censor), 
        count = n()) %>% mutate(nlost = count - attrit)
    int <- c(ltab.data$interval, NA)  #length is 1 + nlost and nevent 
    lifetable <- round(with(ltab.data, lifetab(tis = int, ninit = nrow(data), 
        nlost = nlost, nevent = attrit)), 3)
    return(lifetable)
}
```

### Life Table of Aggregated Faculty

The first Life Table will tabulate the survival probabilities and hazards for all faculty (aggregated).

``` r
ltab <- LifeTableIt(data)
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

### Survival Probabilities for Faculty, by Tenure Status

Now, let's tabulate lifetables for each of the tenure status groups and extract survival probabilities.

``` r
LT_TT <- data %>% filter(EntryTENSTA == "Tenure-Track") %>% LifeTableIt()

LT_NTT <- data %>% filter(EntryTENSTA == "NTT") %>% LifeTableIt()

LT_NTS <- data %>% filter(EntryTENSTA == "Non-tenure System/Position") %>% LifeTableIt()

LTs <- data.frame(LT_TT[, c("surv", "se.surv")], LT_NTT[, c("surv", "se.surv")], 
    LT_NTS[, c("surv", "se.surv")])
names(LTs) <- c("TT Surv", "TT SE", "NTT Surv", "NTT SE", "NTS Surv", "NTS SE")
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

### Graphing Survival Curves of Faculty

#### Life Table (Actuarial) Estimator

``` r
TYPE <- c(rep("Tenure-Track", 11), rep("Non-Tenure Track", 11), rep("Non-Tenure System", 
    11))
SURV <- c(LTs$`TT Surv`, LTs$`NTT Surv`, LTs$`NTS Surv`)
SE <- c(LTs$`TT SE`, LTs$`NTT SE`, LTs$`NTS SE`)
INTERVAL <- rep(c(0:10), 3)
TIME <- rep(2 * c(0:10), 3)
ymin <- SURV - 2 * SE
ymax <- SURV + 2 * SE
ltab2 <- data.frame(TYPE, SURV, TIME, ymin, ymax)

ggplot(ltab2, aes(x = TIME, y = SURV, group = rev(TYPE), colour = rev(TYPE))) + 
    ylim(0.5, 1) + geom_step(size = 1) + labs(title = "Probability of Attrition from Academia", 
    subtitle = "Life Table (Actuarial) Estimator", x = "Years Since First Academic Job", 
    y = "Probability of Remaining in Academia") + scale_colour_manual(name = "Tenure Status", 
    labels = c("Tenure-Track", "Non-tenure Track", "No Tenure System"), values = c("#E7B800", 
        "#2E9FDF", "blue")) + theme_bw() + geom_rect(aes(xmin = TIME, xmax = TIME + 
    2, ymin = ymin, ymax = ymax), alpha = 0.1)
```

![](graphs/Life_table_graph-1.png)

#### Kaplan-Meier (Product Limit) Estimator

``` r
km.tenure <- survfit(Surv(TIME, Censor) ~ strata(EntryTENSTA), data = data, 
    type = "kaplan-meier")
ggsurvplot(km.tenure, data = data, ylim = c(0.5, 1), size = 1, palette = c("#E7B800", 
    "#2E9FDF", "blue"), conf.int = TRUE, legend.labs = c("Tenure-Track", "Non-Tenure track", 
    "No Tenure System"), ggtheme = theme_bw(), ncensor.plot.height = 0.25, title = "Probability of Remaining in Academia: by Initial Tenure Status", 
    subtitle = "Kaplan-Meier Estimator", xlab = "Years Since First Academic Job", 
    ylab = "Probability of Remaining in Academia", legend = "right", legend.title = "Tenure Status", 
    ncensor.plot = F)
```

![](graphs/KM_Graph-1.png)

The default confidence interval (above) is log, which calculates intervals based on log(survival). Other options include plain or log-log. IDRE has code for another method called loghall.

### Statistical Test for differences

Often it is useful to have a statistical test for whether survival curves of two groups differ. The following function implements the G-rho family of Harrington and Fleming (1982), with weights on each death of S(t)^rho, where S is the Kaplan-Meier estimate of survival. With rho = 0 this is the log-rank or Mantel-Haenszel test, and with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test. Allison also discusses other possibilities, including the Wilcoxon test, the Cox test of equality, the Tarone-Ware test of equality, the Peto-Peto-Prentice test of equality and the Generalized Fleming-Harrington test of equality. Let's look at the most common test, the log-rank test.

``` r
survdiff(Surv(TIME, Censor) ~ EntryTENSTA, data = data, rho = 0)  # log-rank or Mantel-Haenszel test
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

The null hypothesis for a log-rank test is that the groups have the same survival. In this case, the Chi-square statistic is significant for both age and tenure status of first appointment. The survival curve, therefore, differs based on the entry age (or tenure status) of the PhD recipient's first academic job). But the KM approach is not able to estimate survival adjusted for covariates. For this reason, we will return to the semi-parametric Cox Proportional hazards model and also parametric survival models shortly.

### Summarizing Distributions

According to Allison, the life table and related graph are sufficent for this part of my analysis. However, it is important to note that these curves can also be summarized in other ways ways. First, let's look at the typical (average) surival times between groups. One downside to this approach is the fact that you need to make an assumption about what the typical career length is.

``` r
km.means.quantiles(km.tenure, 20)
```

    ## [1] "When assuming career span equals 20, the average career duration in years is..."
    ##               strata(EntryTENSTA)=Tenure-Track 
    ##                                         18.271 
    ##                        strata(EntryTENSTA)=NTT 
    ##                                         15.665 
    ## strata(EntryTENSTA)=Non-tenure System/Position 
    ##                                         14.367 
    ## [1] "25 Percent of the sample (or this category) lasts NA or fewer years"
    ## [2] "25 Percent of the sample (or this category) lasts 10 or fewer years"
    ## [3] "25 Percent of the sample (or this category) lasts 6 or fewer years" 
    ## [1] "10 Percent of the sample (or this category) lasts 9 or fewer years"
    ## [2] "10 Percent of the sample (or this category) lasts 3 or fewer years"
    ## [3] "10 Percent of the sample (or this category) lasts 2 or fewer years"

``` r
km.means.quantiles(km.tenure, 35)
```

    ## [1] "When assuming career span equals 35, the average career duration in years is..."
    ##               strata(EntryTENSTA)=Tenure-Track 
    ##                                         31.279 
    ##                        strata(EntryTENSTA)=NTT 
    ##                                         25.852 
    ## strata(EntryTENSTA)=Non-tenure System/Position 
    ##                                         23.356 
    ## [1] "25 Percent of the sample (or this category) lasts NA or fewer years"
    ## [2] "25 Percent of the sample (or this category) lasts 10 or fewer years"
    ## [3] "25 Percent of the sample (or this category) lasts 6 or fewer years" 
    ## [1] "10 Percent of the sample (or this category) lasts 9 or fewer years"
    ## [2] "10 Percent of the sample (or this category) lasts 3 or fewer years"
    ## [3] "10 Percent of the sample (or this category) lasts 2 or fewer years"

Those on the tenure-track tend to have lengthy careers (31.2 years, assuming 35 year-long career). NTT and those outside of tenure systems typically remain in academia between 24 and 26 years. 10 percent of tenure track faculty remain 8 years or fewer. This same proportion lasts two years or fewer in the case of NTT faculty and 1 year or fewer in the case of non-tenure-system faculty.

An important finding here, is that attrition typically pans out for all classes of faculty. In other words, after they have been in the career a handful of years, their probability of staying does not fluctuate that much. This has important implications for human resource policies. Turnover not as concerning for off track faculty after they have remained for a while.

#### Summary of Nonparametric models

Life tables are a great way to summarize distributions and survival curves. Kaplan-Meier curves and log-rank tests are also useful, however, they are most useful when the predictor variable is categorical (e.g., drug vs. placebo), or takes a small number of values (e.g., drug doses 0, 20, 50, and 100 mg/day) that can be treated as categorical. The log-rank test and KM curves don’t work easily with quantitative predictors. For quantitative predictor variables, we turn to alternative methods like the Cox proportional hazards model or Accelerated Failure Time (AFT) models. Such models work also with categorical predictor variables, encoded as {0,1} indicator or dummy variables.

Regression Models for Survival Data
-----------------------------------

There are two types of regression models:

1.  Semi-parametric models, the most common of which is the Cox Proportional Hazards model. Proportional hazard models in general (like the Cox model) assume that the effect of a covariate is to multiply a baseline hazard by some constant. Hazards are “proportional” when the ratio of the hazards for any two individuals is constant, i.e., it does not depend on time.

2.  Fully Parametric AFT models, where it is assumed that log(To) has a specific probability distribution. This class of models assumes that the effect of a covariate is to accelerate or decelerate the life course of a career by some constant.

One case worth noting is that the Weibull distribution (including the exponential distribution as a special case) can be parameterised as either a proportional hazards model or an AFT model. It is the only family of distributions that possesses this property.

### Cox Proportional Hazards Model

The most common semiparametric method is the Cox Proportional Hazards Model. The biggest advantage of the Cox model relates to its flexibilty of functional form. Parametric models require a choice of functional form and often there is no good basis for which to choose. In many instances the choice can be overly restrictive. The Cox model requires no commitment to the functional form of the baseline hazard function. This is a semi-parametric model in that only the effects of covariates are parametrized, and not the hazard function. In other words, we don't make any distributional assumptions about survival times.

Survival analysis is more robust than traditional OLS, in particular because of how it accommodates censoring. Quite commonly, survival data are right censored. This means that for some observations, we do not observe the event of interest. We only know that the event occurs sometime in the future for that observation. In other cases, data are interval censored. This means that we do not directly observe the event times of individuals, we can only say that the event occurs within some specified interval. In this study, we observe several waves of data and we only know the interval during which the even occured. These data, in other words, are interval censored.

R's capabilities are limited with regard to interval censored data. The package most commonly used for cox models is called "coxph()" and this package cannot account for interval censoring--only right censoring. Nevertheless, in this exercise we'll specify a cox model that treats time-to-event as if it were only right censored. After that, we'll turn to another package that is less commonly used in survival analysis, but does handle interval censored data using cox's proportional hazards.

#### Right Censored Data Model

First, we fit a simple Cox model predicting attrition from academia from a categorical predictor of tenure status at workforce entry. We employ the efron method of dealing with ties, although other popular methods (e.g., breslow method) are available. The Efron approximation is more accurate when dealing with tied death times, and is as efficient computationally.

``` r
RC_Mod1 <- coxph(Surv(TIME, Censor) ~ EntryTENSTA, data = data, method = "efron")  # breslow option available
summary(RC_Mod1)
```

    ## Call:
    ## coxph(formula = Surv(TIME, Censor) ~ EntryTENSTA, data = data, 
    ##     method = "efron")
    ## 
    ##   n= 21436, number of events= 3422 
    ## 
    ##                                          coef exp(coef) se(coef)     z
    ## EntryTENSTANTT                        1.09101   2.97727  0.06120 17.83
    ## EntryTENSTANon-tenure System/Position 1.41504   4.11663  0.05152 27.46
    ##                                       Pr(>|z|)    
    ## EntryTENSTANTT                          <2e-16 ***
    ## EntryTENSTANon-tenure System/Position   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                       exp(coef) exp(-coef) lower .95
    ## EntryTENSTANTT                            2.977     0.3359     2.641
    ## EntryTENSTANon-tenure System/Position     4.117     0.2429     3.721
    ##                                       upper .95
    ## EntryTENSTANTT                            3.357
    ## EntryTENSTANon-tenure System/Position     4.554
    ## 
    ## Concordance= 0.639  (se = 0.005 )
    ## Rsquare= 0.045   (max possible= 0.951 )
    ## Likelihood ratio test= 994.6  on 2 df,   p=0
    ## Wald test            = 758.7  on 2 df,   p=0
    ## Score (logrank) test = 876  on 2 df,   p=0

Proportional Hazards models and AFT models must be interpreted in different ways. AFT models give the percentage change in survival time. Cox PH models give the percentage change to the hazard at all time points, following this formula: *h*(*t*)=*h*<sub>0</sub>(*t*)*e*<sup>(*β*′*x*)</sup>

In this case, the effect of initial tenure status on time to attrition has an estimated coefficient of 1.0910072 and 1.4150352. Exponentiated, this means that subjects appointed to lower tenure-status jobs multiply their baseline hazards *h*<sub>0</sub>(*t*) by a factor of 2.9772714 and 4.1166314. Their "risk"" of attriting from academia is 197.7271383 and 311.6631416 percent higher than academics who begin immediately on the tenure-track. Importantly, Cox models state that this is the impact on the subject's hazard at any given time, t. It does not, however, imply an expansion (or contraction) of the lifespan of the subject.

Now let's build a more comprehensive model. Let's specify a complex model using training data and then test that model using independent test data. This will prevent overfitting. Importantly we'll include an interaction between time and our tenure status variable to test the assumption of proportional hazards.

``` r
smp_size <- floor(0.6 * nrow(data))
set.seed(777)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
```

``` r
RC_Mod2 <- coxph(Surv(TIME, Censor) ~ NTT + NTS + DEG2ENTRY + EntryWAPRI + EntryGOVSUP + 
    EntryWKTRNI + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + 
    EntryMARIND + EntryCHLVIN + EntryCTZUSIN + tt(NTT) + tt(NTS), data = test, 
    method = "efron", robust = TRUE, tt = function(x, t, ...) x * log(t))
summary(RC_Mod2)
```

    ## Call:
    ## coxph(formula = Surv(TIME, Censor) ~ NTT + NTS + DEG2ENTRY + 
    ##     EntryWAPRI + EntryGOVSUP + EntryWKTRNI + EntryPUBPRI + EntryEMTP + 
    ##     SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    ##     EntryCTZUSIN + tt(NTT) + tt(NTS), data = test, robust = TRUE, 
    ##     tt = function(x, t, ...) x * log(t), method = "efron")
    ## 
    ##   n= 7866, number of events= 1185 
    ##    (709 observations deleted due to missingness)
    ## 
    ##                                       coef  exp(coef)   se(coef)
    ## NTT                              1.3125751  3.7157299  0.2437654
    ## NTS                              1.5105245  4.5291056  0.2125174
    ## DEG2ENTRY                       -0.0194452  0.9807427  0.0168075
    ## EntryWAPRIOther                  0.5043190  1.6558575  0.1058783
    ## EntryWAPRIResearch               0.4478723  1.5649789  0.0961390
    ## EntryGOVSUPSupport              -0.0003969  0.9996032  0.0675952
    ## EntryWKTRNITraining              0.0260795  1.0264226  0.0594274
    ## EntryPUBPRIPrivate               0.0552103  1.0567629  0.0610538
    ## EntryEMTPTwo-year                0.0466405  1.0477453  0.2059886
    ## EntryEMTPMed                     0.1889138  1.2079368  0.0742695
    ## EntryEMTPUni Research Institute  0.1485878  1.1601947  0.0855925
    ## SDRCARNR2                        0.1187281  1.1260637  0.0974193
    ## SDRCARNDoctorate                 0.0883125  1.0923294  0.0986340
    ## SDRCARNOther                     0.2363906  1.2666690  0.2065230
    ## SDRCARNMedHealth                 0.0791300  1.0823450  0.1408324
    ## EntryAGE                        -0.0052559  0.9947579  0.0052339
    ## GENDERFemale                    -0.0749055  0.9278311  0.0598757
    ## MINRTYYes                        0.0034068  1.0034127  0.0771502
    ## EntryMARINDYes                  -0.0817102  0.9215390  0.0698328
    ## EntryCHLVINYes                  -0.0128549  0.9872273  0.0710196
    ## EntryCTZUSINCitizen              0.0492965  1.0505317  0.0693631
    ## tt(NTT)                         -0.3103874  0.7331629  0.1675865
    ## tt(NTS)                         -0.4034162  0.6680340  0.1419105
    ##                                  robust se      z Pr(>|z|)    
    ## NTT                              0.2402014  5.464 4.64e-08 ***
    ## NTS                              0.2100866  7.190 6.48e-13 ***
    ## DEG2ENTRY                        0.0172884 -1.125  0.26069    
    ## EntryWAPRIOther                  0.1021842  4.935 8.00e-07 ***
    ## EntryWAPRIResearch               0.0909620  4.924 8.49e-07 ***
    ## EntryGOVSUPSupport               0.0674931 -0.006  0.99531    
    ## EntryWKTRNITraining              0.0589793  0.442  0.65836    
    ## EntryPUBPRIPrivate               0.0615966  0.896  0.37008    
    ## EntryEMTPTwo-year                0.2085494  0.224  0.82304    
    ## EntryEMTPMed                     0.0744199  2.538  0.01113 *  
    ## EntryEMTPUni Research Institute  0.0857288  1.733  0.08305 .  
    ## SDRCARNR2                        0.0976316  1.216  0.22395    
    ## SDRCARNDoctorate                 0.0984406  0.897  0.36966    
    ## SDRCARNOther                     0.2084270  1.134  0.25673    
    ## SDRCARNMedHealth                 0.1408488  0.562  0.57425    
    ## EntryAGE                         0.0051335 -1.024  0.30590    
    ## GENDERFemale                     0.0597823 -1.253  0.21022    
    ## MINRTYYes                        0.0769514  0.044  0.96469    
    ## EntryMARINDYes                   0.0704344 -1.160  0.24601    
    ## EntryCHLVINYes                   0.0709031 -0.181  0.85613    
    ## EntryCTZUSINCitizen              0.0688503  0.716  0.47399    
    ## tt(NTT)                          0.1646211 -1.885  0.05937 .  
    ## tt(NTS)                          0.1412603 -2.856  0.00429 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                 exp(coef) exp(-coef) lower .95 upper .95
    ## NTT                                3.7157     0.2691    2.3205    5.9498
    ## NTS                                4.5291     0.2208    3.0005    6.8366
    ## DEG2ENTRY                          0.9807     1.0196    0.9481    1.0145
    ## EntryWAPRIOther                    1.6559     0.6039    1.3553    2.0230
    ## EntryWAPRIResearch                 1.5650     0.6390    1.3094    1.8704
    ## EntryGOVSUPSupport                 0.9996     1.0004    0.8757    1.1410
    ## EntryWKTRNITraining                1.0264     0.9743    0.9144    1.1522
    ## EntryPUBPRIPrivate                 1.0568     0.9463    0.9366    1.1924
    ## EntryEMTPTwo-year                  1.0477     0.9544    0.6962    1.5768
    ## EntryEMTPMed                       1.2079     0.8279    1.0440    1.3976
    ## EntryEMTPUni Research Institute    1.1602     0.8619    0.9808    1.3725
    ## SDRCARNR2                          1.1261     0.8880    0.9299    1.3635
    ## SDRCARNDoctorate                   1.0923     0.9155    0.9007    1.3248
    ## SDRCARNOther                       1.2667     0.7895    0.8419    1.9058
    ## SDRCARNMedHealth                   1.0823     0.9239    0.8212    1.4264
    ## EntryAGE                           0.9948     1.0053    0.9848    1.0048
    ## GENDERFemale                       0.9278     1.0778    0.8252    1.0432
    ## MINRTYYes                          1.0034     0.9966    0.8629    1.1668
    ## EntryMARINDYes                     0.9215     1.0851    0.8027    1.0580
    ## EntryCHLVINYes                     0.9872     1.0129    0.8591    1.1344
    ## EntryCTZUSINCitizen                1.0505     0.9519    0.9179    1.2023
    ## tt(NTT)                            0.7332     1.3640    0.5310    1.0123
    ## tt(NTS)                            0.6680     1.4969    0.5065    0.8811
    ## 
    ## Concordance= 0.66  (se = 0.019 )
    ## Rsquare= 0.047   (max possible= 0.921 )
    ## Likelihood ratio test= 377.3  on 23 df,   p=0
    ## Wald test            = 286.8  on 23 df,   p=0
    ## Score (logrank) test = 335.8  on 23 df,   p=0,   Robust = 384.1  p=0
    ## 
    ##   (Note: the likelihood ratio and score tests assume independence of
    ##      observations within a cluster, the Wald and robust score tests do not).

##### Publication-ready tabulation

``` r
table <- round(as.data.frame(summary(RC_Mod2)$coefficients), 3)
rownames(table) <- c("Non-tenure Track (NTT)", "No Tenure System", "Time between Degree and Job", 
    "Admin/Other", "Researcher", "Gov't Support", "Workplace Training", "Private Control", 
    "Two-year/Other", "Medical", "Research Institute", "PhD Research II", "PhD Doctorate Institution", 
    "PhD Other", "PhD Medical/Health", "Age", "Female", "Minority", "Married", 
    "Children", "Citizen", "Time x NTT", "Time x No Tenure System")
cox_table <- table
kable(cox_table)
```

|                             |    coef|  exp(coef)|  se(coef)|  robust se|       z|  Pr(&gt;|z|)|
|-----------------------------|-------:|----------:|---------:|----------:|-------:|------------:|
| Non-tenure Track (NTT)      |   1.313|      3.716|     0.244|      0.240|   5.464|        0.000|
| No Tenure System            |   1.511|      4.529|     0.213|      0.210|   7.190|        0.000|
| Time between Degree and Job |  -0.019|      0.981|     0.017|      0.017|  -1.125|        0.261|
| Admin/Other                 |   0.504|      1.656|     0.106|      0.102|   4.935|        0.000|
| Researcher                  |   0.448|      1.565|     0.096|      0.091|   4.924|        0.000|
| Gov't Support               |   0.000|      1.000|     0.068|      0.067|  -0.006|        0.995|
| Workplace Training          |   0.026|      1.026|     0.059|      0.059|   0.442|        0.658|
| Private Control             |   0.055|      1.057|     0.061|      0.062|   0.896|        0.370|
| Two-year/Other              |   0.047|      1.048|     0.206|      0.209|   0.224|        0.823|
| Medical                     |   0.189|      1.208|     0.074|      0.074|   2.538|        0.011|
| Research Institute          |   0.149|      1.160|     0.086|      0.086|   1.733|        0.083|
| PhD Research II             |   0.119|      1.126|     0.097|      0.098|   1.216|        0.224|
| PhD Doctorate Institution   |   0.088|      1.092|     0.099|      0.098|   0.897|        0.370|
| PhD Other                   |   0.236|      1.267|     0.207|      0.208|   1.134|        0.257|
| PhD Medical/Health          |   0.079|      1.082|     0.141|      0.141|   0.562|        0.574|
| Age                         |  -0.005|      0.995|     0.005|      0.005|  -1.024|        0.306|
| Female                      |  -0.075|      0.928|     0.060|      0.060|  -1.253|        0.210|
| Minority                    |   0.003|      1.003|     0.077|      0.077|   0.044|        0.965|
| Married                     |  -0.082|      0.922|     0.070|      0.070|  -1.160|        0.246|
| Children                    |  -0.013|      0.987|     0.071|      0.071|  -0.181|        0.856|
| Citizen                     |   0.049|      1.051|     0.069|      0.069|   0.716|        0.474|
| Time x NTT                  |  -0.310|      0.733|     0.168|      0.165|  -1.885|        0.059|
| Time x No Tenure System     |  -0.403|      0.668|     0.142|      0.141|  -2.856|        0.004|

Even after controlling for background characteristics, there are significant differences. Here, NTT status or working at a college or university without a tenure system impacts the hazard, multiplying the baseline by a factor of 3.7157299 and 4.5291056. This is equivalent to saying that each tenure status increases the hazard of attrition by 271.5729917 and 352.9105583 percent, controlling for background characteristics. R output also provides the exponentiated negative coefficient. To my understanding, that just allows you to compare the groups relative to the baseline hazard of the tenure-track group. Robust standard errors were used in this model.

The model reveals several other important predcitors, including the subject main job and institution type. The interaction between time and no tenure system status is somewhat concerning. These issues are being worked through with Allison.

An important assumption of the Cox model is that hazard functions are proportional. We can test each of the variables in the model, as well as test the model as a whole using the cox.zph() function. It tests proportionality by interacting each predictor with log time (km transformation). Rho is the pearson product-moment correlation between the scaled residuals (Schoenfeld) and log(time) for each covariate. The global test jointly tests all of the interactions. Low p-values suggest a violation of the assumption of proportional hazards.

``` r
cox.zph <- cox.zph(RC_Mod2, transform = "km")
round(cox.zph$table, 3)
```

    ##                                    rho   chisq     p
    ## NTT                             -0.523 451.263 0.000
    ## NTS                             -0.601 684.840 0.000
    ## DEG2ENTRY                        0.016   0.336 0.562
    ## EntryWAPRIOther                 -0.025   0.724 0.395
    ## EntryWAPRIResearch              -0.034   1.256 0.262
    ## EntryGOVSUPSupport              -0.012   0.171 0.679
    ## EntryWKTRNITraining             -0.008   0.079 0.778
    ## EntryPUBPRIPrivate              -0.027   0.987 0.320
    ## EntryEMTPTwo-year                0.000   0.000 0.991
    ## EntryEMTPMed                     0.133  21.998 0.000
    ## EntryEMTPUni Research Institute  0.006   0.049 0.825
    ## SDRCARNR2                        0.018   0.416 0.519
    ## SDRCARNDoctorate                 0.010   0.122 0.727
    ## SDRCARNOther                     0.029   1.061 0.303
    ## SDRCARNMedHealth                 0.095  10.949 0.001
    ## EntryAGE                         0.076   7.005 0.008
    ## GENDERFemale                     0.065   5.301 0.021
    ## MINRTYYes                       -0.037   1.693 0.193
    ## EntryMARINDYes                  -0.065   5.412 0.020
    ## EntryCHLVINYes                  -0.017   0.366 0.545
    ## EntryCTZUSINCitizen              0.015   0.287 0.592
    ## tt(NTT)                          0.536 469.680 0.000
    ## tt(NTS)                          0.626 779.592 0.000
    ## GLOBAL                              NA 857.495 0.000

Generally, there is some evidence that the hazards are propoprtional (violating the key assumption). For most covariates, rho is not significant. However, it is concerning that it is significant for our two most important predictors above (NTT and Non-tenure system). Also, the global test is significant. The significant interaction of time and no tenure system status also suggests that hazards may not be proportional across strata.

#### Interval Censored Data Model

The data used in this study are interval censored. We know the interval during when an event took place, but we do not know the particular time of the event. A different package must be used to integrate interval censoring, however, this package (and all others I could find) are not used very commonly by researchers in R.

``` r
IC_Mod1 <- ic_sp(Surv(as.numeric(EndDate), as.numeric(ExitDate2), type = "interval2") ~ 
    NTT + NTS, model = "ph", bs_samples = 10, data = data)
summary(IC_Mod1)
```

    ## 
    ## Model:  Cox PH
    ## Baseline:  semi-parametric 
    ## Call: ic_sp(formula = Surv(as.numeric(EndDate), as.numeric(ExitDate2), 
    ##     type = "interval2") ~ NTT + NTS, data = data, model = "ph", 
    ##     bs_samples = 10)
    ## 
    ##     Estimate Exp(Est) Std.Error z-value p
    ## NTT   0.3824    1.466   0.02495   15.33 0
    ## NTS   0.2483    1.282   0.02115   11.74 0
    ## 
    ## final llk =  -43922.98 
    ## Iterations =  57 
    ## Bootstrap Samples =  10 
    ## WARNING: only  10  bootstrap samples used for standard errors. 
    ## Suggest using more bootstrap samples for inference

``` r
IC_Mod2 <- ic_sp(Surv(as.numeric(EntryDate), as.numeric(ExitDate2), type = "interval2") ~ 
    NTT + NTS + DEG2ENTRY + EntryWAPRI + EntryGOVSUP + EntryWKTRNI + EntryPUBPRI + 
        EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
        EntryCTZUSIN + NTT * (TIME) + NTS * (TIME), model = "ph", bs_samples = 100, 
    data = test)
summary(IC_Mod2)
```

    ## 
    ## Model:  Cox PH
    ## Baseline:  semi-parametric 
    ## Call: ic_sp(formula = Surv(as.numeric(EntryDate), as.numeric(ExitDate2), 
    ##     type = "interval2") ~ NTT + NTS + DEG2ENTRY + EntryWAPRI + 
    ##     EntryGOVSUP + EntryWKTRNI + EntryPUBPRI + EntryEMTP + SDRCARN + 
    ##     EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    ##     EntryCTZUSIN + NTT * (TIME) + NTS * (TIME), data = test, 
    ##     model = "ph", bs_samples = 100)
    ## 
    ##                                  Estimate Exp(Est) Std.Error  z-value
    ## NTT                              0.307700   1.3600  0.070170   4.3850
    ## NTS                             -0.133800   0.8748  0.064960  -2.0590
    ## DEG2ENTRY                       -0.192000   0.8253  0.010620 -18.0800
    ## EntryWAPRIOther                  0.016460   1.0170  0.053900   0.3054
    ## EntryWAPRIResearch               0.094340   1.0990  0.043120   2.1880
    ## EntryGOVSUPSupport              -0.045950   0.9551  0.029950  -1.5340
    ## EntryWKTRNITraining             -0.071850   0.9307  0.031820  -2.2580
    ## EntryPUBPRIPrivate              -0.076370   0.9265  0.027170  -2.8110
    ## EntryEMTPTwo-year                0.029520   1.0300  0.094880   0.3112
    ## EntryEMTPMed                     0.181500   1.1990  0.041880   4.3340
    ## EntryEMTPUni Research Institute -0.071930   0.9306  0.043580  -1.6510
    ## SDRCARNR2                       -0.080180   0.9229  0.050640  -1.5830
    ## SDRCARNDoctorate                -0.125400   0.8822  0.050730  -2.4710
    ## SDRCARNOther                    -0.346800   0.7069  0.104900  -3.3060
    ## SDRCARNMedHealth                 0.060390   1.0620  0.081720   0.7390
    ## EntryAGE                         0.004273   1.0040  0.002338   1.8280
    ## GENDERFemale                    -0.103700   0.9015  0.030570  -3.3940
    ## MINRTYYes                       -0.289300   0.7488  0.039580  -7.3090
    ## EntryMARINDYes                  -0.071050   0.9314  0.036400  -1.9520
    ## EntryCHLVINYes                   0.094040   1.0990  0.034390   2.7340
    ## EntryCTZUSINCitizen              0.235000   1.2650  0.033790   6.9570
    ## TIME                            -0.086010   0.9176  0.008583 -10.0200
    ## NTT:TIME                        -0.016210   0.9839  0.011010  -1.4720
    ## NTS:TIME                        -0.008626   0.9914  0.010670  -0.8083
    ##                                         p
    ## NTT                             1.160e-05
    ## NTS                             3.950e-02
    ## DEG2ENTRY                       0.000e+00
    ## EntryWAPRIOther                 7.601e-01
    ## EntryWAPRIResearch              2.867e-02
    ## EntryGOVSUPSupport              1.250e-01
    ## EntryWKTRNITraining             2.393e-02
    ## EntryPUBPRIPrivate              4.946e-03
    ## EntryEMTPTwo-year               7.557e-01
    ## EntryEMTPMed                    1.461e-05
    ## EntryEMTPUni Research Institute 9.880e-02
    ## SDRCARNR2                       1.134e-01
    ## SDRCARNDoctorate                1.347e-02
    ## SDRCARNOther                    9.472e-04
    ## SDRCARNMedHealth                4.599e-01
    ## EntryAGE                        6.761e-02
    ## GENDERFemale                    6.892e-04
    ## MINRTYYes                       2.687e-13
    ## EntryMARINDYes                  5.091e-02
    ## EntryCHLVINYes                  6.251e-03
    ## EntryCTZUSINCitizen             3.477e-12
    ## TIME                            0.000e+00
    ## NTT:TIME                        1.410e-01
    ## NTS:TIME                        4.189e-01
    ## 
    ## final llk =  -11872.16 
    ## Iterations =  37 
    ## Bootstrap Samples =  100

``` r
# Need to hack to extract the table
tab <- read_excel("/Users/chadgevans/Research/Projects/Faculty_Attrition/cache/fit_ph.xlsx")
tab <- select(tab, -X__1)
names(tab) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
rownames(tab) <- c("Non-tenure Track (NTT)", "No Tenure System", "Time between Degree and Job", 
    "Admin/Other", "Researcher", "Gov't Support", "Workplace Training", "Private Control", 
    "Two-year/Other", "Medical", "Research Institute", "PhD Research II", "PhD Doctorate Institution", 
    "PhD Other", "PhD Medical/Health", "Age", "Female", "Minority", "Married", 
    "Children", "Citizen", "Time", "Time x NTT", "Time x No Tenure System")
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
cox_table2 <- tab
kable(cox_table2)
```

|                             |       coef|  exp(coef)|  se(coef)|         z|          p|
|-----------------------------|----------:|----------:|---------:|---------:|----------:|
| Non-tenure Track (NTT)      |   0.307700|     1.3600|  0.069350|    4.4360|  0.0000091|
| No Tenure System            |  -0.133800|     0.8748|  0.069430|   -1.9260|  0.0540500|
| Time between Degree and Job |  -0.192000|     0.8253|  0.010530|  -18.2400|  0.0000000|
| Admin/Other                 |   0.016460|     1.0170|  0.049900|    0.3299|  0.7415000|
| Researcher                  |   0.094340|     1.0990|  0.041040|    2.2990|  0.0215300|
| Gov't Support               |  -0.045950|     0.9551|  0.033210|   -1.3840|  0.1665000|
| Workplace Training          |  -0.071850|     0.9307|  0.031320|   -2.2950|  0.0217600|
| Private Control             |  -0.076370|     0.9265|  0.032070|   -2.3820|  0.0172400|
| Two-year/Other              |   0.029520|     1.0300|  0.103800|    0.2845|  0.7760000|
| Medical                     |   0.181500|     1.1990|  0.038010|    4.7760|  0.0000018|
| Research Institute          |  -0.071930|     0.9306|  0.050220|   -1.4320|  0.1521000|
| PhD Research II             |  -0.080180|     0.9229|  0.045670|   -1.7560|  0.0791600|
| PhD Doctorate Institution   |  -0.125400|     0.8822|  0.045580|   -2.7500|  0.0059590|
| PhD Other                   |  -0.346800|     0.7069|  0.107200|   -3.2340|  0.0012210|
| PhD Medical/Health          |   0.060390|     1.0620|  0.073720|    0.8191|  0.4127000|
| Age                         |   0.004273|     1.0040|  0.002380|    1.7950|  0.0726200|
| Female                      |  -0.103700|     0.9015|  0.028620|   -3.6250|  0.0002891|
| Minority                    |  -0.289300|     0.7488|  0.032440|   -8.9180|  0.0000000|
| Married                     |  -0.071050|     0.9314|  0.034950|   -2.0330|  0.0420300|
| Children                    |   0.094040|     1.0990|  0.034380|    2.7360|  0.0062240|
| Citizen                     |   0.235000|     1.2650|  0.033890|    6.9360|  0.0000000|
| Time                        |  -0.086010|     0.9176|  0.007699|  -11.1700|  0.0000000|
| Time x NTT                  |  -0.016210|     0.9839|  0.011370|   -1.4260|  0.1538000|
| Time x No Tenure System     |  -0.008626|     0.9914|  0.011110|   -0.7765|  0.4375000|

#### Comparing Right Censored Model and Interval Censored Model

``` r
cox_compare <- data.frame(cbind(rbind(cox_table[, c(2, 6)], rep(NA, 2)), cox_table2[, 
    c(2, 5)]))
rownames(cox_compare) <- rownames(cox_table2)
colnames(cox_compare) <- c("RC_coef", "RC_pvalue", "IC_coef", "IC_pvalue")
kable(cox_compare)
```

|                             |  RC\_coef|  RC\_pvalue|  IC\_coef|  IC\_pvalue|
|-----------------------------|---------:|-----------:|---------:|-----------:|
| Non-tenure Track (NTT)      |     3.716|       0.000|    1.3600|   0.0000091|
| No Tenure System            |     4.529|       0.000|    0.8748|   0.0540500|
| Time between Degree and Job |     0.981|       0.261|    0.8253|   0.0000000|
| Admin/Other                 |     1.656|       0.000|    1.0170|   0.7415000|
| Researcher                  |     1.565|       0.000|    1.0990|   0.0215300|
| Gov't Support               |     1.000|       0.995|    0.9551|   0.1665000|
| Workplace Training          |     1.026|       0.658|    0.9307|   0.0217600|
| Private Control             |     1.057|       0.370|    0.9265|   0.0172400|
| Two-year/Other              |     1.048|       0.823|    1.0300|   0.7760000|
| Medical                     |     1.208|       0.011|    1.1990|   0.0000018|
| Research Institute          |     1.160|       0.083|    0.9306|   0.1521000|
| PhD Research II             |     1.126|       0.224|    0.9229|   0.0791600|
| PhD Doctorate Institution   |     1.092|       0.370|    0.8822|   0.0059590|
| PhD Other                   |     1.267|       0.257|    0.7069|   0.0012210|
| PhD Medical/Health          |     1.082|       0.574|    1.0620|   0.4127000|
| Age                         |     0.995|       0.306|    1.0040|   0.0726200|
| Female                      |     0.928|       0.210|    0.9015|   0.0002891|
| Minority                    |     1.003|       0.965|    0.7488|   0.0000000|
| Married                     |     0.922|       0.246|    0.9314|   0.0420300|
| Children                    |     0.987|       0.856|    1.0990|   0.0062240|
| Citizen                     |     1.051|       0.474|    1.2650|   0.0000000|
| Time                        |     0.733|       0.059|    0.9176|   0.0000000|
| Time x NTT                  |     0.668|       0.004|    0.9839|   0.1538000|
| Time x No Tenure System     |        NA|          NA|    0.9914|   0.4375000|

#### Summary of Cox Proportional Hazard Models

If the proportional hazards assumption holds, then it is possible to estimate the effect of parameter(s) without any consideration of the baseline hazard function. As mentioned, this is in contrast to parametric models--the focus of our next section.

### Accelerated Failure Time Models

Next, we fit a parametric survival regression model. These are location-scale models for an arbitrary transform of the time variable; the most common cases use a log transformation, leading to accelerated failure time (AFT) models. First, we assume the outcome has an exponential distribution--a good baseline distribution to start with (simplifies calculations). The exponential distribution implies a constant hazard rate. I also model with the log-logistic transformation. This is one of the more popular transformations because, unlike the Weibull distribution, it can exhibit a non-monotonic hazard function which increases at early times and decreases at later times. It also has a closed form solution that speeds up computation (important because of the consequences of censoring). The advantage of the Weibull (and by extention the exponential), of course, is that it can be parameterised as a PH model or an AFT model. In other words, the Weibull family can be interpreted as affecting the risk of event occurance or the duration of the lifespan. Other potential functions include log normal, gamma and inverse gaussian functions.

``` r
exp.mod <- survreg(Surv(TIME, Censor) ~ EntryTENSTA + EntryAGE, data = data, 
    dist = "exponential")
summary(exp.mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(TIME, Censor) ~ EntryTENSTA + EntryAGE, 
    ##     data = data, dist = "exponential")
    ##                                          Value Std. Error      z         p
    ## (Intercept)                            4.27689     0.1065  40.17  0.00e+00
    ## EntryTENSTANTT                        -1.11364     0.0612 -18.20  4.93e-74
    ## EntryTENSTANon-tenure System/Position -1.46464     0.0514 -28.50 1.28e-178
    ## EntryAGE                               0.00872     0.0027   3.23  1.25e-03
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -14986.1   Loglik(intercept only)= -15536.4
    ##  Chisq= 1100.68 on 3 degrees of freedom, p= 0 
    ## Number of Newton-Raphson Iterations: 6 
    ## n= 21436

For this first model, we parameterized log(t) using the exponential distribution. In the case of the exponential distribution, there is one extra parameter that allows *ϵ* to take on one of the extreme value distributions:

*f*(*ϵ*)=*e*<sup>(*ϵ* − *e*<sup>*ϵ*</sup>)</sup>

In the case of AFT models, covariates are interpreted as having a multiplicative effect on the survival time (the expected life span) of an individual. So, NTT status has a slope of -1.1136357 in the simple model above. This means that beginning in a non-tenure track position causes time to attrition to change by -67.1637026 percent. The expected time spent in academia before leaving, in other words, is lower for those starting in non-tenure track. Faculty beginning in no-tenure systems tend to spend even less time in academia -76.8839867 percent, compared to their tenure-track peers. In both cases, the life course is "accelerated."

We can also examine the coefficient for job-entry age. For a one unit increase in age, we expect a 0.8761806 percent change in survival time. Because the coefficient is so small, you can actually just multiply the coefficient by 100 to find an approximation of the percentage change. It "decelerates" the life course by about .9 percent, thereby increasing the total time span spent in academia.

Because the exponential distribution is a special case of the Weibull family, we can interpret the coefficients as changes in the hazard (as we did for cox models). While this is a good baseline, the log logistic transformation is more commonly used, thanks to its computational efficiently and (relative) flexibility of functional form. For the log-logistic parameterization, the errors are scaled as follows:

Now let's estimate a log logistic model. We'll use the same specification, but allow epsilon to take on a slightly different density:
*f*(*ϵ*)=*e*<sup>*ϵ*</sup>/(1 + *e*<sup>*ϵ*</sup>)<sup>2</sup>

``` r
ll.mod <- survreg(Surv(TIME, Censor) ~ EntryTENSTA + EntryAGE, data = data, 
    dist = "loglogistic")
summary(ll.mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(TIME, Censor) ~ EntryTENSTA + EntryAGE, 
    ##     data = data, dist = "loglogistic")
    ##                                          Value Std. Error     z         p
    ## (Intercept)                            3.58158    0.09028  39.7  0.00e+00
    ## EntryTENSTANTT                        -0.97540    0.05007 -19.5  1.56e-84
    ## EntryTENSTANon-tenure System/Position -1.29172    0.04173 -31.0 2.46e-210
    ## EntryAGE                               0.00841    0.00228   3.7  2.19e-04
    ## Log(scale)                            -0.31827    0.01321 -24.1 3.30e-128
    ## 
    ## Scale= 0.727 
    ## 
    ## Log logistic distribution
    ## Loglik(model)= -14718.6   Loglik(intercept only)= -15345.9
    ##  Chisq= 1254.65 on 3 degrees of freedom, p= 0 
    ## Number of Newton-Raphson Iterations: 5 
    ## n= 21436

In this case, the log logistic performs similarly to the exponential distribution, both in terms of fit and coefficient estimates. The simple curvature of the survival curves in this study makes both of these distributions similar. I'll opt for the exponential because it fits well, it is simpler and it can be interpreted in the same paradigm as before. We'll also fit the full model with robust standard errors.

#### Mathmatics

1.  Construct the likelihood function
    -   uncensored cases contribute their density probability at t
    -   censored cases contribute their survival probability at t

2.  Simplify the expression
3.  Take the log to simplify the math
4.  Calculate the derivative and set it equal to zero
5.  Solve using, typically, the Newton-Raphson algorithm

``` r
mod <- survreg(Surv(TIME, Censor) ~ EntryTENSTA + DEG2ENTRY + EntryWAPRI + EntryGOVSUP + 
    EntryWKTRNI + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + 
    EntryMARIND + EntryCHLVIN + EntryCTZUSIN + EntryTENSTA * (TIME), data = data, 
    dist = "exponential", robust = TRUE)
summary(mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(TIME, Censor) ~ EntryTENSTA + DEG2ENTRY + 
    ##     EntryWAPRI + EntryGOVSUP + EntryWKTRNI + EntryPUBPRI + EntryEMTP + 
    ##     SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    ##     EntryCTZUSIN + EntryTENSTA * (TIME), data = data, dist = "exponential", 
    ##     robust = TRUE)
    ##                                               Value Std. Err (Naive SE)
    ## (Intercept)                                 2.76153  0.14005    0.15368
    ## EntryTENSTANTT                             -0.72667  0.10331    0.11250
    ## EntryTENSTANon-tenure System/Position      -0.92211  0.09216    0.10018
    ## DEG2ENTRY                                   0.02820  0.00947    0.01024
    ## EntryWAPRIOther                            -0.45294  0.05856    0.06536
    ## EntryWAPRIResearch                         -0.33731  0.05498    0.06050
    ## EntryGOVSUPSupport                         -0.01231  0.03791    0.04220
    ## EntryWKTRNITraining                        -0.05121  0.03329    0.03745
    ## EntryPUBPRIPrivate                         -0.03143  0.03397    0.03823
    ## EntryEMTPTwo-year                          -0.14879  0.10773    0.11922
    ## EntryEMTPMed                               -0.09138  0.04119    0.04685
    ## EntryEMTPUni Research Institute            -0.11804  0.04669    0.05286
    ## SDRCARNR2                                  -0.02080  0.05593    0.06239
    ## SDRCARNDoctorate                           -0.02488  0.05529    0.06178
    ## SDRCARNOther                               -0.17264  0.11020    0.12479
    ## SDRCARNMedHealth                            0.03573  0.07973    0.08913
    ## EntryAGE                                    0.00424  0.00292    0.00327
    ## GENDERFemale                                0.12758  0.03341    0.03760
    ## MINRTYYes                                   0.07100  0.04502    0.05002
    ## EntryMARINDYes                              0.02196  0.03909    0.04362
    ## EntryCHLVINYes                              0.00335  0.03985    0.04437
    ## EntryCTZUSINCitizen                        -0.22482  0.03899    0.04390
    ## TIME                                        0.26510  0.01357    0.01476
    ## EntryTENSTANTT:TIME                         0.00383  0.01885    0.02028
    ## EntryTENSTANon-tenure System/Position:TIME  0.03916  0.01640    0.01748
    ##                                                   z        p
    ## (Intercept)                                 19.7183 1.50e-86
    ## EntryTENSTANTT                              -7.0338 2.01e-12
    ## EntryTENSTANon-tenure System/Position      -10.0050 1.45e-23
    ## DEG2ENTRY                                    2.9789 2.89e-03
    ## EntryWAPRIOther                             -7.7350 1.03e-14
    ## EntryWAPRIResearch                          -6.1356 8.49e-10
    ## EntryGOVSUPSupport                          -0.3248 7.45e-01
    ## EntryWKTRNITraining                         -1.5381 1.24e-01
    ## EntryPUBPRIPrivate                          -0.9252 3.55e-01
    ## EntryEMTPTwo-year                           -1.3812 1.67e-01
    ## EntryEMTPMed                                -2.2188 2.65e-02
    ## EntryEMTPUni Research Institute             -2.5280 1.15e-02
    ## SDRCARNR2                                   -0.3719 7.10e-01
    ## SDRCARNDoctorate                            -0.4501 6.53e-01
    ## SDRCARNOther                                -1.5666 1.17e-01
    ## SDRCARNMedHealth                             0.4482 6.54e-01
    ## EntryAGE                                     1.4542 1.46e-01
    ## GENDERFemale                                 3.8179 1.35e-04
    ## MINRTYYes                                    1.5772 1.15e-01
    ## EntryMARINDYes                               0.5617 5.74e-01
    ## EntryCHLVINYes                               0.0841 9.33e-01
    ## EntryCTZUSINCitizen                         -5.7662 8.11e-09
    ## TIME                                        19.5359 5.44e-85
    ## EntryTENSTANTT:TIME                          0.2034 8.39e-01
    ## EntryTENSTANon-tenure System/Position:TIME   2.3874 1.70e-02
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -11652.7   Loglik(intercept only)= -13835
    ##  Chisq= 4364.53 on 24 degrees of freedom, p= 0 
    ## (Loglikelihood assumes independent observations)
    ## Number of Newton-Raphson Iterations: 7 
    ## n=19693 (1743 observations deleted due to missingness)

### AFT Specification

``` r
table <- as.data.frame(summary(mod)$table)
rownames(table) <- c("Intercept", "Non-tenure Track", "No Tenure System", "Time between Degree and Job", 
    "Administration/Other", "Research", "Gov't Support", "Workplace Training", 
    "Private Control", "Two-year", "Medical", "Research Institute", "PhD Research II", 
    "PhD Doctorate Institution", "PhD Other", "PhD Medical/Health", "Age", "Female", 
    "Minority", "Married", "Parent", "Citizen", "Time", "Time x NTT", "Time x No Tenure System")
table$expCoef <- exp(table$Value)
names(table)[names(table) == "Std. Err"] <- "Robust SE"
aft_table <- round(table[, c("Value", "expCoef", "Robust SE", "z", "p")], 3)
kable(aft_table)
```

|                             |   Value|  expCoef|  Robust SE|        z|      p|
|-----------------------------|-------:|--------:|----------:|--------:|------:|
| Intercept                   |   2.762|   15.824|      0.140|   19.718|  0.000|
| Non-tenure Track            |  -0.727|    0.484|      0.103|   -7.034|  0.000|
| No Tenure System            |  -0.922|    0.398|      0.092|  -10.005|  0.000|
| Time between Degree and Job |   0.028|    1.029|      0.009|    2.979|  0.003|
| Administration/Other        |  -0.453|    0.636|      0.059|   -7.735|  0.000|
| Research                    |  -0.337|    0.714|      0.055|   -6.136|  0.000|
| Gov't Support               |  -0.012|    0.988|      0.038|   -0.325|  0.745|
| Workplace Training          |  -0.051|    0.950|      0.033|   -1.538|  0.124|
| Private Control             |  -0.031|    0.969|      0.034|   -0.925|  0.355|
| Two-year                    |  -0.149|    0.862|      0.108|   -1.381|  0.167|
| Medical                     |  -0.091|    0.913|      0.041|   -2.219|  0.027|
| Research Institute          |  -0.118|    0.889|      0.047|   -2.528|  0.011|
| PhD Research II             |  -0.021|    0.979|      0.056|   -0.372|  0.710|
| PhD Doctorate Institution   |  -0.025|    0.975|      0.055|   -0.450|  0.653|
| PhD Other                   |  -0.173|    0.841|      0.110|   -1.567|  0.117|
| PhD Medical/Health          |   0.036|    1.036|      0.080|    0.448|  0.654|
| Age                         |   0.004|    1.004|      0.003|    1.454|  0.146|
| Female                      |   0.128|    1.136|      0.033|    3.818|  0.000|
| Minority                    |   0.071|    1.074|      0.045|    1.577|  0.115|
| Married                     |   0.022|    1.022|      0.039|    0.562|  0.574|
| Parent                      |   0.003|    1.003|      0.040|    0.084|  0.933|
| Citizen                     |  -0.225|    0.799|      0.039|   -5.766|  0.000|
| Time                        |   0.265|    1.304|      0.014|   19.536|  0.000|
| Time x NTT                  |   0.004|    1.004|      0.019|    0.203|  0.839|
| Time x No Tenure System     |   0.039|    1.040|      0.016|    2.387|  0.017|

The exponential distribution has a constant hazard *λ*(*t*)=*λ* and thus a survival function of *S*(*t*)=*e*<sup>−*λ*(*t*)</sup> and density of *f*(*t*)=*λ* × *e*<sup>−*λ*(*t*)</sup>. An interesting occurance is that the expected survival time for this distribution is *E*(*t*)=1/*λ* and its variance is *E*(*t*)=1/*λ*<sup>2</sup>. This makes the mean survival time equal to e^intercept (15.8241007). It's inverse (0.0631947) is the MLE of the (constant) hazard rate.

AFT models are typcally interpreted in a way that covariates have a multiplicative effect on the expected survival time. So, with regard to tenure status, taking your first job as NTT accelerates the time to attrition by a factor of 0.4835179 (0.4835179 times shorter survival time compared to the baseline survival). Beginning an academic career in a non-tenure system accelerates the time to attrition by a factor of 0.3976797. The life course for these states is -51.6482143 and -60.2320289 percent shorter, respectively.

The Weibull family of distributions (of which the exponential is a sub-class) has the advantage that covariates can also be interpreted as an impact on the hazard ratios. For this famiily of distributions, the coefficient is multiplied by -1 and then multiplied by a shape parameter (1/scale parameter). In the case of the exponential distributuion, the shape parameter is simply 1/1. So in our case, the hazard ratio comparing NTT to tenure-track positions is 2.0681759. The risk of attrition increases by a factor of 2.0681759 when one begins an academic career in a NTT position. Faculty with a first job at a non-teure system institution increases their risk by a factor of 2.5145864.

Comparing the Cox model and the Exponential AFT model
-----------------------------------------------------

``` r
cox_table2 <- rbind(rep(NA, 6), cox_table)
aft_table <- aft_table[-23, ]
ctable <- cbind(cox_table2[, c(2, 6)], aft_table[, c(2, 5)])
rownames(ctable) <- rownames(aft_table)
ctable$expCoef <- round(exp(aft_table$Value * -1 * 1/1), 3)  # convert to hazard ratios
colnames(ctable) <- c("Cox:HR", "Cox p-val", "AFT:HR", "AFT p-val")
kable(ctable)
```

|                             |  Cox:HR|  Cox p-val|  AFT:HR|  AFT p-val|
|-----------------------------|-------:|----------:|-------:|----------:|
| Intercept                   |      NA|         NA|   0.063|      0.000|
| Non-tenure Track            |   3.716|      0.000|   2.069|      0.000|
| No Tenure System            |   4.529|      0.000|   2.514|      0.000|
| Time between Degree and Job |   0.981|      0.261|   0.972|      0.003|
| Administration/Other        |   1.656|      0.000|   1.573|      0.000|
| Research                    |   1.565|      0.000|   1.401|      0.000|
| Gov't Support               |   1.000|      0.995|   1.012|      0.745|
| Workplace Training          |   1.026|      0.658|   1.052|      0.124|
| Private Control             |   1.057|      0.370|   1.031|      0.355|
| Two-year                    |   1.048|      0.823|   1.161|      0.167|
| Medical                     |   1.208|      0.011|   1.095|      0.027|
| Research Institute          |   1.160|      0.083|   1.125|      0.011|
| PhD Research II             |   1.126|      0.224|   1.021|      0.710|
| PhD Doctorate Institution   |   1.092|      0.370|   1.025|      0.653|
| PhD Other                   |   1.267|      0.257|   1.189|      0.117|
| PhD Medical/Health          |   1.082|      0.574|   0.965|      0.654|
| Age                         |   0.995|      0.306|   0.996|      0.146|
| Female                      |   0.928|      0.210|   0.880|      0.000|
| Minority                    |   1.003|      0.965|   0.931|      0.115|
| Married                     |   0.922|      0.246|   0.978|      0.574|
| Parent                      |   0.987|      0.856|   0.997|      0.933|
| Citizen                     |   1.051|      0.474|   1.252|      0.000|
| Time x NTT                  |   0.733|      0.059|   0.996|      0.839|
| Time x No Tenure System     |   0.668|      0.004|   0.962|      0.017|

Whether choosing to model as a Cox model or an exponential AFT model, the results are comparable. There are some minor differences in the hazard ratios associated with the demographic variables and institution type. But they are pretty minor.

Extensions
----------

Overfitting is always a risk when developing statistical models. There are several extensions for dealing with overfitting (Penalized Regression). Methods include ridge regression (ridge), smoothing splines (pspline), and frailty models (frailty).

### Ridge regression

We can also fit a model with different baseline survival shapes for each of the three tenure groups (i.e., two different scale parameters). Here, I also robustified the standard errors using the robust argument. I also included a ridge penalty on both predictors to prevent overfitting.

When employing ridge regression, the likelihood is penalised by theta/2 times the linear combination of the squared coefficients. When scaling, penalties are calculated based on rescaling covariates to have unit variance. for coefficients based on rescaling the predictors to have unit variance. If df is specified then theta is chosen based on an approximate degrees of freedom

``` r
rmod <- survreg(Surv(TIME, Censor) ~ EntryTENSTA + EntryWAPRI + EntryGOVSUP + 
    EntryWKTRNI + EntryCARNEG + EntryPUBPRI + EntryEMTP + SDRCARN + EntryAGE + 
    GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN + ridge(EntryTENSTA, 
    EntryWAPRI, EntryGOVSUP, EntryWKTRNI, EntryCARNEG, EntryPUBPRI, EntryEMTP, 
    SDRCARN, EntryAGE, GENDER, MINRTY, EntryMARIND, EntryCHLVIN, EntryCTZUSIN, 
    theta = 5, scale = T), data = data, dist = "exponential")
summary(rmod)
```

``` r
rtable <- as.data.frame(summary(rmod)$table)
rownames(rtable) <- c("Intercept", "Non-tenure Track", "No Tenure System", "Researcher", 
    "Teacher", "Gov't Support", "Workplace Training", "Assoc/BA/MA Institution", 
    "Other Institution", "Private Control", "Four-year", "Medical", "Research Institute", 
    "Other Educ", "Research II", "Doctorate Institution", "Other", "Medical/Health", 
    "Age", "Male", "Minority", "Married", "Parent", "Citizen", seq(1:14))
rtable$expCoef <- exp(rtable$Value)
names(rtable)[names(rtable) == "Std. Error"] <- "Ridge SE"
rtable <- round(rtable[1:24, ], 3)
rtable <- rtable[, c("Value", "expCoef", "Ridge SE", "z", "p")]
kable(rtable)
```

``` r
ctable <- cbind(table, rtable[1:24, ])
kable(ctable)
```

### Frailty models

These models allow you to add a simple random effects term to a survival model (Cox or survreg).
