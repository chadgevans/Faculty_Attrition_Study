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
-   [Life Tables](#life-tables)
    -   Life table of all faculty, aggregated
    -   Life table of faculty, by Tenure Status
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
-   [Comparing the Exponential AFT Model and the Cox Models](#comparing-the-exponential-aft-model-and-the-cox-models)
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
length(unique(data$REFID))  # The number of individuals who participated in the study
```

Life Tables
-----------

The life table method or estimator is also known as the actuarial method/estimator. For the life table method, if observations are censored on the same month (or time unit) that events occurred, they are assumed to be at risk for just half the month. This is in contrast to the Kaplan-Meier estimator which assumes they are at risk for the whole month. The KM estimator is generally better for smaller datasets with exact times of event occurance.

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

### Life Table of All Faculty, Aggregated

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

By the last interval, over 30 percent of faculty have attrited from their careers in academia. The hazard is highest during the early stages of the careers and then plateaus later on.

### Life Table for Faculty, by Tenure Status

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

Clearly, the survival probabilities cascade much more rapidly for NTT and NTS faculty.

### Graphing Survival Curves of Faculty

#### Life Table (Actuarial) Estimator

Survival probabilities are sometimes communicated most efficiently in graphical form. The following graph plots the actuarial survival probabilities for each class of faculty.

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

These survival probabilities can also be estimated and graphed using the KM estimator. Allison recommends reporting only the actuarial estimator, probably because it is more robust for larger samples.

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

The null hypothesis for a log-rank test is that the groups have the same survival. In this case, the Chi-square statistic is significant for tenure status of first appointment. The survival curve, therefore, differs based on the tenure status of the PhD recipient's first academic job. But the KM approach is not able to estimate survival adjusted for covariates. For this reason, we will return to the semi-parametric Cox Proportional hazards model and also parametric survival models shortly.

#### Summary of Nonparametric models

Life tables are a great way to summarize distributions and survival curves. Kaplan-Meier curves and log-rank tests are also useful, however, they are most useful when the predictor variable is categorical (e.g., drug vs. placebo), or takes a small number of values (e.g., drug doses 0, 20, 50, and 100 mg/day) that can be treated as categorical. The log-rank test and KM curves don’t work easily with quantitative predictors. For quantitative predictor variables, we turn to alternative methods like the Cox proportional hazards model or Accelerated Failure Time (AFT) models. Such models work also with categorical predictor variables, encoded as {0,1} indicator or dummy variables.

Substantially, the actuarial estimator and KM estimator tell similar stories. Attrition is much higher and more rapid for faculty working off the tenure track. The lifetable estimator perhaps shows a slightly higher attrition rate for NTS faculty, comparted to the KM estimate.

Regression Models for Survival Data
-----------------------------------

There are two types of regression models in survival analysis:

1.  Semi-parametric models, the most common of which is the Cox Proportional Hazards model. Proportional hazard models in general (like the Cox model) assume that the effect of a covariate is to multiply a baseline hazard by some constant. Hazards are “proportional” when the ratio of the hazards for any two individuals is constant, i.e., it does not depend on time.

2.  Fully Parametric AFT models, where it is assumed that log(To) has a specific probability distribution. This class of models assumes that the effect of a covariate is to accelerate or decelerate the life course of a career by some constant.

One case worth noting is that the Weibull distribution (including the exponential distribution as a special case) can be parameterised as either a proportional hazards model or an AFT model. It is the only family of distributions that possesses this property.

### Cox Proportional Hazards Model

The most common semiparametric method is the Cox Proportional Hazards Model. The biggest advantage of the Cox model relates to its flexibilty of functional form. Parametric models require a choice of functional form and often there is no good basis for which to choose. In many instances the choice can be overly restrictive. The Cox model requires no commitment to the functional form of the baseline hazard function. This is a semi-parametric model in that only the effects of covariates are parametrized, and not the hazard function. In other words, we don't make any distributional assumptions about survival times.

Cox models work mathematically because, when hazards are proportional, it is possible to factor the Likelihood function into a part with betas and a part with betas and the baseline hazard function. We sacrifice the information of part two and use the standard likelihood approach to estimate the betas in part one (which is consistent and asymptotically normal). This approach is not fully efficient, but we gain flexibilty with functional form when we do this (making estimates more robust). It turns out that estimating part one only requires knowing the order in which events took place.

There are some implicit assumptions, however. For one, we must assume that all characteristics impacting the hazard are included in the model. Also, censoring must be noninformative, observations must be independent and variable must also be measured without error.

Survival analysis is more robust than traditional OLS, in particular because of how it accommodates censoring. Quite commonly, survival data are right censored. This means that for some observations, we do not observe the event of interest. We only know that the event occurs sometime in the future for that observation. In other cases, data are interval censored. This means that we do not directly observe the event times of individuals, we can only say that the event occurs within some specified interval. Of course, sometimes data are both right censored (for some observations) and interval censored. This in the case for this study. We observe several waves of data and we know the interval during which the event occured for most observations. Sometimes we only know that the event occurs some point after the last interval.

Unfortunately, R's capabilities are somewhat limited with regard to interval censored data. coxph() in the "survival" package (the function most commonly used for cox models) does not account for interval censoring--it only handles right censoring. Nevertheless, in this exercise we'll first specify a cox model treating time-to-event as if it were entirely right censored. After that, we'll turn to another package that is less commonly used in survival analysis, but does handle both right censoring and interval censored data.

For specifying our model, we will subset the data and first train a model using training data. After we are satisfied with that specification, we'll test using the independent test data. This will help prevent overfitting.

``` r
smp_size <- floor(0.6 * nrow(data))
set.seed(777)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
```

#### Right Censored Data Model

First, we fit a simple Cox model predicting attrition from academia from a categorical predictor of tenure status at workforce entry. We assume here that we know exact event times (not the interval). We employ the efron method of dealing with ties, although other popular methods (e.g., breslow method) are available. The Efron approximation is more accurate when dealing with tied death times, and is as efficient computationally.

These data are also censored (partially) on the left. Specifically, we know that individuals entered academia between two study waves, but we do not know the particular entry year. So we make an assumption for all of our models that time of entry was one year prior to the first wave when we observe their participation in academia. In other words, we assume that they entered the position midway through the previous interval. This will be about right on average.

The Surv function here frames censoring as follows. It assumes that individuals entered the study halfway through the previous interval. It treats censoring from the "EndDate"--that is, the last year he or she was still employed. Censoring, of course, happens after that point, so an individual who enters in year 0, remains in the study through year 2 in its entirely and subsequently drops out between year 2 and 4 would be coded as "3+." There are two years for the full interval, plus the extra entry year.

``` r
RC_Mod1 <- coxph(Surv(TIME2, Censor) ~ NTT + NTS, data = test, method = "efron")  # breslow option available
summary(RC_Mod1)
```

    ## Call:
    ## coxph(formula = Surv(TIME2, Censor) ~ NTT + NTS, data = test, 
    ##     method = "efron")
    ## 
    ##   n= 8575, number of events= 1362 
    ## 
    ##        coef exp(coef) se(coef)     z Pr(>|z|)    
    ## NTT 1.00893   2.74265  0.09565 10.55   <2e-16 ***
    ## NTS 1.32374   3.75745  0.07982 16.59   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##     exp(coef) exp(-coef) lower .95 upper .95
    ## NTT     2.743     0.3646     2.274     3.308
    ## NTS     3.757     0.2661     3.213     4.394
    ## 
    ## Concordance= 0.637  (se = 0.008 )
    ## Rsquare= 0.04   (max possible= 0.935 )
    ## Likelihood ratio test= 353.4  on 2 df,   p=0
    ## Wald test            = 276.6  on 2 df,   p=0
    ## Score (logrank) test = 313.6  on 2 df,   p=0

Proportional Hazards models and AFT models must be interpreted in different ways. AFT models give the percentage change in survival time. Cox PH models give the percentage change to the hazard at all time points, following this formula: *h*(*t*)=*h*<sub>0</sub>(*t*)*e*<sup>(*β*′*x*)</sup>

In this case, the effect of initial tenure status on time to attrition has an estimated coefficient of 1.0089259 and 1.3237404. Exponentiated, this means that subjects appointed to lower tenure-status jobs multiply their baseline hazards *h*<sub>0</sub>(*t*) by a factor of 2.7426537 and 3.7574493. Their "risk"" of attriting from academia is 174.2653671 and 275.744933 percent higher than academics who begin immediately on the tenure-track. Importantly, Cox models state that this is the impact on the subject's hazard at any given time, t. It does not, however, imply an expansion (or contraction) of the lifespan of the subject.

Now let's build a more comprehensive model. Let's specify a complex model using training data and then test that model using independent test data. This will prevent overfitting. Importantly we'll include an interaction between time and our tenure status variables (NTT or NTS dummies) to test the assumption of proportional hazards. Here we use the TIME2 variable which simple adds an extra year to everyone's event time. You are basically assuming then that their interval began halfway through the previous interval.

``` r
RC_Mod2 <- coxph(Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY + EntryWAPRI + 
    EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + SDRCARN + 
    EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN + 
    tt(NTT) + tt(NTS), data = test, method = "efron", robust = TRUE, tt = function(x, 
    t, ...) x * log(t))
summary(RC_Mod2)
```

    ## Call:
    ## coxph(formula = Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY + 
    ##     EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * 
    ##     EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + 
    ##     EntryCHLVIN + EntryCTZUSIN + tt(NTT) + tt(NTS), data = test, 
    ##     robust = TRUE, tt = function(x, t, ...) x * log(t), method = "efron")
    ## 
    ##   n= 7905, number of events= 1190 
    ##    (670 observations deleted due to missingness)
    ## 
    ##                                                         coef exp(coef)
    ## NTT                                                 1.050756  2.859812
    ## NTS                                                 1.242949  3.465820
    ## DEG2ENTRY                                          -0.021551  0.978680
    ## EntryWAPRIOther                                     0.486410  1.626466
    ## EntryWAPRIResearch                                  0.436570  1.547390
    ## EntryWKTRNITraining                                 0.021265  1.021492
    ## EntryPUBPRIPrivate                                 -0.060413  0.941376
    ## EntryEMTPTwo-year                                  -0.168894  0.844598
    ## EntryEMTPMed                                        0.127429  1.135904
    ## EntryEMTPUni Research Institute                     0.051707  1.053067
    ## SDRCARNR2                                           0.138130  1.148125
    ## SDRCARNDoctorate                                    0.100120  1.105303
    ## SDRCARNOther                                        0.229712  1.258237
    ## SDRCARNMedHealth                                    0.076845  1.079875
    ## EntryAGE                                           -0.004249  0.995760
    ## GENDERFemale                                       -0.078362  0.924630
    ## MINRTYYes                                           0.026342  1.026692
    ## EntryMARINDYes                                     -0.090126  0.913816
    ## EntryCHLVINYes                                     -0.006429  0.993591
    ## EntryCTZUSINCitizen                                 0.058315  1.060049
    ## tt(NTT)                                            -0.199849  0.818854
    ## tt(NTS)                                            -0.345253  0.708041
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                1.084622  2.958321
    ## EntryPUBPRIPrivate:EntryEMTPMed                     0.170195  1.185537
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute  0.297362  1.346303
    ##                                                     se(coef) robust se
    ## NTT                                                 0.142225  0.139420
    ## NTS                                                 0.127150  0.123364
    ## DEG2ENTRY                                           0.016889  0.017496
    ## EntryWAPRIOther                                     0.105574  0.102241
    ## EntryWAPRIResearch                                  0.091283  0.087047
    ## EntryWKTRNITraining                                 0.059346  0.059029
    ## EntryPUBPRIPrivate                                  0.087221  0.087645
    ## EntryEMTPTwo-year                                   0.234258  0.236871
    ## EntryEMTPMed                                        0.097365  0.097210
    ## EntryEMTPUni Research Institute                     0.107171  0.106534
    ## SDRCARNR2                                           0.096694  0.096822
    ## SDRCARNDoctorate                                    0.098305  0.098070
    ## SDRCARNOther                                        0.206401  0.208361
    ## SDRCARNMedHealth                                    0.140849  0.141829
    ## EntryAGE                                            0.005181  0.005068
    ## GENDERFemale                                        0.059659  0.059419
    ## MINRTYYes                                           0.076835  0.076643
    ## EntryMARINDYes                                      0.069745  0.070267
    ## EntryCHLVINYes                                      0.070845  0.070647
    ## EntryCTZUSINCitizen                                 0.069132  0.068481
    ## tt(NTT)                                             0.116969  0.113868
    ## tt(NTS)                                             0.099629  0.098430
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                0.476484  0.481549
    ## EntryPUBPRIPrivate:EntryEMTPMed                     0.137364  0.138026
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute  0.169968  0.169697
    ##                                                         z Pr(>|z|)    
    ## NTT                                                 7.537 4.82e-14 ***
    ## NTS                                                10.075  < 2e-16 ***
    ## DEG2ENTRY                                          -1.232 0.218043    
    ## EntryWAPRIOther                                     4.757 1.96e-06 ***
    ## EntryWAPRIResearch                                  5.015 5.29e-07 ***
    ## EntryWKTRNITraining                                 0.360 0.718669    
    ## EntryPUBPRIPrivate                                 -0.689 0.490641    
    ## EntryEMTPTwo-year                                  -0.713 0.475832    
    ## EntryEMTPMed                                        1.311 0.189903    
    ## EntryEMTPUni Research Institute                     0.485 0.627425    
    ## SDRCARNR2                                           1.427 0.153684    
    ## SDRCARNDoctorate                                    1.021 0.307302    
    ## SDRCARNOther                                        1.102 0.270258    
    ## SDRCARNMedHealth                                    0.542 0.587945    
    ## EntryAGE                                           -0.838 0.401872    
    ## GENDERFemale                                       -1.319 0.187233    
    ## MINRTYYes                                           0.344 0.731076    
    ## EntryMARINDYes                                     -1.283 0.199625    
    ## EntryCHLVINYes                                     -0.091 0.927486    
    ## EntryCTZUSINCitizen                                 0.852 0.394461    
    ## tt(NTT)                                            -1.755 0.079242 .  
    ## tt(NTS)                                            -3.508 0.000452 ***
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                2.252 0.024300 *  
    ## EntryPUBPRIPrivate:EntryEMTPMed                     1.233 0.217550    
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute  1.752 0.079720 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                                    exp(coef) exp(-coef)
    ## NTT                                                   2.8598     0.3497
    ## NTS                                                   3.4658     0.2885
    ## DEG2ENTRY                                             0.9787     1.0218
    ## EntryWAPRIOther                                       1.6265     0.6148
    ## EntryWAPRIResearch                                    1.5474     0.6462
    ## EntryWKTRNITraining                                   1.0215     0.9790
    ## EntryPUBPRIPrivate                                    0.9414     1.0623
    ## EntryEMTPTwo-year                                     0.8446     1.1840
    ## EntryEMTPMed                                          1.1359     0.8804
    ## EntryEMTPUni Research Institute                       1.0531     0.9496
    ## SDRCARNR2                                             1.1481     0.8710
    ## SDRCARNDoctorate                                      1.1053     0.9047
    ## SDRCARNOther                                          1.2582     0.7948
    ## SDRCARNMedHealth                                      1.0799     0.9260
    ## EntryAGE                                              0.9958     1.0043
    ## GENDERFemale                                          0.9246     1.0815
    ## MINRTYYes                                             1.0267     0.9740
    ## EntryMARINDYes                                        0.9138     1.0943
    ## EntryCHLVINYes                                        0.9936     1.0065
    ## EntryCTZUSINCitizen                                   1.0600     0.9434
    ## tt(NTT)                                               0.8189     1.2212
    ## tt(NTS)                                               0.7080     1.4123
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                  2.9583     0.3380
    ## EntryPUBPRIPrivate:EntryEMTPMed                       1.1855     0.8435
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute    1.3463     0.7428
    ##                                                    lower .95 upper .95
    ## NTT                                                   2.1760    3.7585
    ## NTS                                                   2.7214    4.4138
    ## DEG2ENTRY                                             0.9457    1.0128
    ## EntryWAPRIOther                                       1.3311    1.9873
    ## EntryWAPRIResearch                                    1.3047    1.8352
    ## EntryWKTRNITraining                                   0.9099    1.1468
    ## EntryPUBPRIPrivate                                    0.7928    1.1178
    ## EntryEMTPTwo-year                                     0.5309    1.3436
    ## EntryEMTPMed                                          0.9389    1.3743
    ## EntryEMTPUni Research Institute                       0.8546    1.2976
    ## SDRCARNR2                                             0.9497    1.3880
    ## SDRCARNDoctorate                                      0.9120    1.3395
    ## SDRCARNOther                                          0.8364    1.8929
    ## SDRCARNMedHealth                                      0.8178    1.4259
    ## EntryAGE                                              0.9859    1.0057
    ## GENDERFemale                                          0.8230    1.0388
    ## MINRTYYes                                             0.8835    1.1931
    ## EntryMARINDYes                                        0.7962    1.0487
    ## EntryCHLVINYes                                        0.8651    1.1411
    ## EntryCTZUSINCitizen                                   0.9269    1.2123
    ## tt(NTT)                                               0.6551    1.0236
    ## tt(NTS)                                               0.5838    0.8587
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                  1.1512    7.6022
    ## EntryPUBPRIPrivate:EntryEMTPMed                       0.9045    1.5538
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute    0.9654    1.8775
    ## 
    ## Concordance= 0.665  (se = 0.016 )
    ## Rsquare= 0.047   (max possible= 0.923 )
    ## Likelihood ratio test= 382.2  on 25 df,   p=0
    ## Wald test            = 299.5  on 25 df,   p=0
    ## Score (logrank) test = 343.6  on 25 df,   p=0,   Robust = 386.6  p=0
    ## 
    ##   (Note: the likelihood ratio and score tests assume independence of
    ##      observations within a cluster, the Wald and robust score tests do not).

Even after controlling for background characteristics, there are significant differences. Here, NTT status or working at a college or university without a tenure system impacts the hazard, multiplying the baseline by a factor of 2.8598117 and 3.4658201. This is equivalent to saying that each tenure status increases the hazard of attrition by 185.9811676 and 246.5820096 percent, controlling for background characteristics. R output also provides the exponentiated negative coefficient. To my understanding, that just allows you to compare the groups relative to the baseline hazard of the tenure-track group. Robust standard errors were used in this model.

The model reveals other important predcitors, including the subject's main job and possible interactions between institution type and public/private status. The interaction between time and no tenure system status is somewhat concerning. These issues are being worked through with Allison.

An important assumption of the Cox model is that hazard functions are proportional. We can test each of the variables in the model, as well as test the model as a whole using the cox.zph() function. It tests proportionality by interacting each predictor with log time (km transformation). Rho is the pearson product-moment correlation between the scaled residuals (Schoenfeld) and log(time) for each covariate. The global test jointly tests all of the interactions. Low p-values suggest a violation of the assumption of proportional hazards.

``` r
cox.zph <- cox.zph(RC_Mod2, transform = "km")
round(cox.zph$table, 3)
```

    ##                                                       rho    chisq     p
    ## NTT                                                -0.519  438.477 0.000
    ## NTS                                                -0.595  634.254 0.000
    ## DEG2ENTRY                                           0.012    0.196 0.658
    ## EntryWAPRIOther                                    -0.016    0.303 0.582
    ## EntryWAPRIResearch                                 -0.032    1.188 0.276
    ## EntryWKTRNITraining                                 0.016    0.335 0.563
    ## EntryPUBPRIPrivate                                 -0.051    3.397 0.065
    ## EntryEMTPTwo-year                                   0.014    0.237 0.627
    ## EntryEMTPMed                                        0.104   13.156 0.000
    ## EntryEMTPUni Research Institute                     0.002    0.003 0.956
    ## SDRCARNR2                                           0.006    0.049 0.824
    ## SDRCARNDoctorate                                   -0.001    0.001 0.970
    ## SDRCARNOther                                        0.030    1.198 0.274
    ## SDRCARNMedHealth                                    0.111   15.238 0.000
    ## EntryAGE                                            0.063    4.779 0.029
    ## GENDERFemale                                        0.056    3.918 0.048
    ## MINRTYYes                                          -0.052    3.420 0.064
    ## EntryMARINDYes                                     -0.052    3.480 0.062
    ## EntryCHLVINYes                                      0.000    0.000 0.988
    ## EntryCTZUSINCitizen                                 0.028    0.961 0.327
    ## tt(NTT)                                             0.606  658.787 0.000
    ## tt(NTS)                                             0.698 1117.724 0.000
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                0.025    0.814 0.367
    ## EntryPUBPRIPrivate:EntryEMTPMed                     0.030    1.116 0.291
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute -0.002    0.003 0.954
    ## GLOBAL                                                 NA 1212.760 0.000

Generally, there is some evidence that the hazards are not propoprtional (violating the key assumption). For most covariates, rho is not significant. However, it is concerning that it is significant for our two most important predictors above (NTT and Non-tenure system). Also, the global test is significant. The significant interaction of time and no tenure system status also suggests that hazards may not be proportional across strata.

#### Interval Censored Data Model

The data used in this study not only have right censoring, they are interval censored. We know the interval during when an event took place, but we do not know the particular time of the event. Unfortunately, the well-developed and supported function coxph() from the survival package does not handle interval censoring. I searched for a different R package capable of handling interval censoring and found the icenreg package. There is some documentation for this package, however, it is considerably less developed than function in the suvival package.

We'll follow the same model development approach as in the last specification. We'll try to develop a complex model using training data and then apply this model to test data.

Using thsi method does not require the assumption that the subject enters the study midway through the previous wave. It only depends on the interval itself. I suspect this is related to the fact that Cox models do not require the time of occurance, only the order in which events took place.

``` r
IC_Mod1 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS, model = "ph", 
    bs_samples = 100, data = data)
summary(IC_Mod1)
```

    ## 
    ## Model:  Cox PH
    ## Baseline:  semi-parametric 
    ## Call: ic_sp(formula = Surv(lower, upper, type = "interval2") ~ NTT + 
    ##     NTS, data = data, model = "ph", bs_samples = 100)
    ## 
    ##     Estimate Exp(Est) Std.Error z-value p
    ## NTT    1.156    3.177   0.05737   20.15 0
    ## NTS    1.516    4.555   0.05209   29.11 0
    ## 
    ## final llk =  -10106.08 
    ## Iterations =  27 
    ## Bootstrap Samples =  100

``` r
IC_Mod2 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS + DEG2ENTRY + 
    EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + 
    SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN, 
    model = "ph", bs_samples = 100, data = test)
summary(IC_Mod2)
```

    ## 
    ## Model:  Cox PH
    ## Baseline:  semi-parametric 
    ## Call: ic_sp(formula = Surv(lower, upper, type = "interval2") ~ NTT + 
    ##     NTS + DEG2ENTRY + EntryWAPRI + EntryWKTRNI + EntryPUBPRI + 
    ##     EntryEMTP + EntryPUBPRI * EntryEMTP + SDRCARN + EntryAGE + 
    ##     GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN, 
    ##     data = test, model = "ph", bs_samples = 100)
    ## 
    ##                                                      Estimate Exp(Est)
    ## NTT                                                 0.9620000   2.6170
    ## NTS                                                 1.1070000   3.0250
    ## DEG2ENTRY                                          -0.0007627   0.9992
    ## EntryWAPRIOther                                     0.4588000   1.5820
    ## EntryWAPRIResearch                                  0.3826000   1.4660
    ## EntryWKTRNITraining                                 0.0086660   1.0090
    ## EntryPUBPRIPrivate                                 -0.0678400   0.9344
    ## EntryEMTPTwo-year                                  -0.2166000   0.8052
    ## EntryEMTPMed                                        0.1372000   1.1470
    ## EntryEMTPUni Research Institute                     0.0459500   1.0470
    ## SDRCARNR2                                           0.1606000   1.1740
    ## SDRCARNDoctorate                                    0.1026000   1.1080
    ## SDRCARNOther                                        0.2810000   1.3240
    ## SDRCARNMedHealth                                    0.0885100   1.0930
    ## EntryAGE                                           -0.0071180   0.9929
    ## GENDERFemale                                       -0.0341300   0.9664
    ## MINRTYYes                                           0.0074100   1.0070
    ## EntryMARINDYes                                     -0.0946700   0.9097
    ## EntryCHLVINYes                                     -0.0250200   0.9753
    ## EntryCTZUSINCitizen                                -0.0776900   0.9253
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                1.1910000   3.2920
    ## EntryPUBPRIPrivate:EntryEMTPMed                     0.1644000   1.1790
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute  0.2798000   1.3230
    ##                                                    Std.Error  z-value
    ## NTT                                                 0.090890 10.58000
    ## NTS                                                 0.087350 12.67000
    ## DEG2ENTRY                                           0.018310 -0.04166
    ## EntryWAPRIOther                                     0.097340  4.71400
    ## EntryWAPRIResearch                                  0.086160  4.44100
    ## EntryWKTRNITraining                                 0.058350  0.14850
    ## EntryPUBPRIPrivate                                  0.091520 -0.74130
    ## EntryEMTPTwo-year                                   0.233400 -0.92800
    ## EntryEMTPMed                                        0.086250  1.59100
    ## EntryEMTPUni Research Institute                     0.103000  0.44620
    ## SDRCARNR2                                           0.106800  1.50300
    ## SDRCARNDoctorate                                    0.088850  1.15500
    ## SDRCARNOther                                        0.232800  1.20700
    ## SDRCARNMedHealth                                    0.132500  0.66790
    ## EntryAGE                                            0.005012 -1.42000
    ## GENDERFemale                                        0.065940 -0.51750
    ## MINRTYYes                                           0.080430  0.09214
    ## EntryMARINDYes                                      0.068410 -1.38400
    ## EntryCHLVINYes                                      0.065880 -0.37980
    ## EntryCTZUSINCitizen                                 0.070860 -1.09600
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                1.528000  0.77960
    ## EntryPUBPRIPrivate:EntryEMTPMed                     0.142900  1.15000
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute  0.169800  1.64800
    ##                                                            p
    ## NTT                                                0.000e+00
    ## NTS                                                0.000e+00
    ## DEG2ENTRY                                          9.668e-01
    ## EntryWAPRIOther                                    2.430e-06
    ## EntryWAPRIResearch                                 8.943e-06
    ## EntryWKTRNITraining                                8.819e-01
    ## EntryPUBPRIPrivate                                 4.585e-01
    ## EntryEMTPTwo-year                                  3.534e-01
    ## EntryEMTPMed                                       1.116e-01
    ## EntryEMTPUni Research Institute                    6.554e-01
    ## SDRCARNR2                                          1.327e-01
    ## SDRCARNDoctorate                                   2.482e-01
    ## SDRCARNOther                                       2.274e-01
    ## SDRCARNMedHealth                                   5.042e-01
    ## EntryAGE                                           1.555e-01
    ## GENDERFemale                                       6.048e-01
    ## MINRTYYes                                          9.266e-01
    ## EntryMARINDYes                                     1.664e-01
    ## EntryCHLVINYes                                     7.041e-01
    ## EntryCTZUSINCitizen                                2.729e-01
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year               4.356e-01
    ## EntryPUBPRIPrivate:EntryEMTPMed                    2.501e-01
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute 9.938e-02
    ## 
    ## final llk =  -3639.088 
    ## Iterations =  25 
    ## Bootstrap Samples =  100

``` r
# This is time-intensive at bs_samples=10, maybe a minute per sample.
IC_Mod3 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS + DEG2ENTRY + 
    EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + 
    EntryPUBPRI * EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + 
    EntryCHLVIN + EntryCTZUSIN + TIME * NTT + TIME * NTS, model = "ph", bs_samples = 10, 
    data = data)
summary(IC_Mod3)
# Also, this model is wrong because there should not be a coefficient
# produced for time.  Furthermore, that coefficient sops up all the
# variation in event times
```

Now let's compare the results of the first model treating the data as right censored and the results when treating the data as interval censored (and also right censored).

#### Comparing Right Censored Model and Interval Censored Model

``` r
space <- c("", "")
temp <- rbind(IC_cox_table[1:20, c(2, 5)], space, space, IC_cox_table[21:23, 
    c(2, 5)])
cox_compare <- data.frame(cbind(RC_cox_table[, c(2, 6)], temp))
rownames(cox_compare) <- rownames(RC_cox_table)
colnames(cox_compare) <- c("RC Cox: HR", "RC Cox pval", "IC Cox: HR", "IC Cox pval")
kable(cox_compare)
```

|                              |  RC Cox: HR|  RC Cox pval| IC Cox: HR | IC Cox pval |
|------------------------------|-----------:|------------:|:-----------|:------------|
| Non-tenure Track (NTT)       |       2.860|        0.000| 2.617      | 0           |
| No Tenure System             |       3.466|        0.000| 3.025      | 0           |
| Time between Degree and Job  |       0.979|        0.218| 0.999      | 0.963       |
| Admin/Other                  |       1.626|        0.000| 1.582      | 0           |
| Researcher                   |       1.547|        0.000| 1.466      | 0           |
| Workplace Training           |       1.021|        0.719| 1.009      | 0.876       |
| Private Control              |       0.941|        0.491| 0.934      | 0.48        |
| Two-year/Other               |       0.845|        0.476| 0.805      | 0.346       |
| Medical                      |       1.136|        0.190| 1.147      | 0.175       |
| Research Institute           |       1.053|        0.627| 1.047      | 0.682       |
| PhD Research II              |       1.148|        0.154| 1.174      | 0.12        |
| PhD Doctorate Institution    |       1.105|        0.307| 1.108      | 0.336       |
| PhD Other                    |       1.258|        0.270| 1.324      | 0.22        |
| PhD Medical/Health           |       1.080|        0.588| 1.093      | 0.492       |
| Age                          |       0.996|        0.402| 0.993      | 0.177       |
| Female                       |       0.925|        0.187| 0.966      | 0.571       |
| Minority                     |       1.027|        0.731| 1.007      | 0.926       |
| Married                      |       0.914|        0.200| 0.91       | 0.248       |
| Children                     |       0.994|        0.927| 0.975      | 0.704       |
| Citizen                      |       1.060|        0.394| 0.925      | 0.313       |
| Time x NTT                   |       0.819|        0.079|            |             |
| Time x No Tenure System      |       0.708|        0.000|            |             |
| Private x Two-Year           |       2.958|        0.024| 3.292      | 0.007       |
| Private x Medical            |       1.186|        0.218| 1.179      | 0.255       |
| Private x Research Institute |       1.346|        0.080| 1.323      | 0.12        |

Again, these models differ slightly. The Right censored (RC) model includes the time interaction with tenure status, but it fails to deal with the interval censoring of the data and assume every individual entered his or her job midway through the previous interval. The Interval censored (IC) model effectivey handles the interval censoring, but it does not allow for an interaction between time and tenure status. Nevertheless, the results are comparable.

The fit statistics of these two models is difficult to reconcile. The RC model has an Rsquare of 0.047. The log likelihood is listed as -1.013048210^{4}, -9939.400298. I'm not sure what two numbers mean. Perhaps one is only for the intercept? The IC model reports a log likelihood of -3639.0884452, which is considerably different from the RC model.

#### Summary of Cox Proportional Hazard Models

If the proportional hazards assumption holds, then it is possible to estimate the effect of parameter(s) without any consideration of the baseline hazard function. As mentioned, this is in contrast to parametric models--the focus of our next section.

### Accelerated Failure Time Models

Next, we fit a parametric survival regression model. These are location-scale models for an arbitrary transform of the time variable; the most common cases use a log transformation, leading to accelerated failure time (AFT) models. First, we assume the outcome has an exponential distribution--a good baseline distribution to start with (simplifies calculations). The exponential distribution implies a constant hazard rate. I also model with the log-logistic transformation. This is one of the more popular transformations because, unlike the Weibull distribution, it can exhibit a non-monotonic hazard function which increases at early times and decreases at later times. It also has a closed form solution that speeds up computation (important because of the consequences of censoring). The advantage of the Weibull (and by extention the exponential), of course, is that it can be parameterised as a PH model or an AFT model. In other words, the Weibull family can be interpreted as affecting the risk of event occurance or the duration of the lifespan. Other potential functions include log normal, gamma and inverse gaussian functions.

For fully parameterized models, the timing of an event matters (unlike the cox models that only require the order). So here we must also make an assumption that subject begin their academic work midway through the previous interval. So we add unity to the upper and lower interval bounds to adjust for this.

``` r
exp.mod <- survreg(Surv(lower + 1, upper + 1, type = "interval2") ~ EntryTENSTA + 
    EntryAGE, data = data, dist = "exponential")
summary(exp.mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(lower + 1, upper + 1, type = "interval2") ~ 
    ##     EntryTENSTA + EntryAGE, data = data, dist = "exponential")
    ##                                          Value Std. Error      z         p
    ## (Intercept)                            4.14264     0.1064  38.92  0.00e+00
    ## EntryTENSTANTT                        -1.13216     0.0612 -18.50  1.95e-76
    ## EntryTENSTANon-tenure System/Position -1.50406     0.0514 -29.26 3.21e-188
    ## EntryAGE                               0.00767     0.0027   2.84  4.53e-03
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -11584.1   Loglik(intercept only)= -12165
    ##  Chisq= 1161.78 on 3 degrees of freedom, p= 0 
    ## Number of Newton-Raphson Iterations: 6 
    ## n= 21436

For this first model, we parameterized log(t) using the exponential distribution. In the case of the exponential distribution, there is one extra parameter that allows *ϵ* to take on one of the extreme value distributions:

*f*(*ϵ*)=*e*<sup>(*ϵ* − *e*<sup>*ϵ*</sup>)</sup>

In the case of AFT models, covariates are interpreted as having a multiplicative effect on the survival time (the expected life span) of an individual. So, NTT status has a slope of -1.1321606 in the simple model above. This means that beginning in a non-tenure track position causes time to attrition to change by -67.7663932 percent. The expected time spent in academia before leaving, in other words, is lower for those starting in non-tenure track. Faculty beginning in no-tenure systems tend to spend even less time in academia -77.7773762 percent, compared to their tenure-track peers. In both cases, the life course is "accelerated."

We can also examine the coefficient for job-entry age. For a one unit increase in age, we expect a 0.7704377 percent change in survival time. Because the coefficient is so small, you can actually just multiply the coefficient by 100 to find an approximation of the percentage change. An increase in age "decelerates" the life course, thereby increasing the total expected time span spent in academia.

Because the exponential distribution is a special case of the Weibull family, we can interpret the coefficients as changes in the hazard (as we did for cox models). While this is a good baseline, the log logistic transformation is sometimes more commonly used, thanks to its computational efficiently and (relative) flexibility of functional form. For the log-logistic parameterization, the errors are scaled as follows:

Now let's estimate a log logistic model. We'll use the same specification, but allow epsilon to take on a slightly different density:
*f*(*ϵ*)=*e*<sup>*ϵ*</sup>/(1 + *e*<sup>*ϵ*</sup>)<sup>2</sup>

``` r
ll.mod <- survreg(Surv(lower + 1, upper + 1, type = "interval2") ~ EntryTENSTA + 
    EntryAGE, data = data, dist = "loglogistic")
summary(ll.mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(lower + 1, upper + 1, type = "interval2") ~ 
    ##     EntryTENSTA + EntryAGE, data = data, dist = "loglogistic")
    ##                                          Value Std. Error      z         p
    ## (Intercept)                            3.45151    0.09169  37.65  0.00e+00
    ## EntryTENSTANTT                        -1.00622    0.05076 -19.82  1.93e-87
    ## EntryTENSTANon-tenure System/Position -1.34887    0.04215 -32.00 1.10e-224
    ## EntryAGE                               0.00775    0.00233   3.33  8.71e-04
    ## Log(scale)                            -0.32328    0.01323 -24.44 6.12e-132
    ## 
    ## Scale= 0.724 
    ## 
    ## Log logistic distribution
    ## Loglik(model)= -11296.4   Loglik(intercept only)= -11959.6
    ##  Chisq= 1326.38 on 3 degrees of freedom, p= 0 
    ## Number of Newton-Raphson Iterations: 5 
    ## n= 21436

In this case, the log logistic performs similarly to the exponential distribution, both in terms of fit and coefficient estimates. The simple curvature of the survival curves in this study makes both of these distributions similar. I'll opt for the exponential because it fits well, it is simpler and it can be interpreted in the same paradigm as before. We'll also fit the full model with robust standard errors.

#### Mathematics

1.  Construct the likelihood function
    -   uncensored cases contribute their density probability at t
    -   censored cases contribute their survival probability at t

2.  Simplify the expression
3.  Take the log to simplify the math
4.  Calculate the derivative and set it equal to zero
5.  Solve using, typically, the Newton-Raphson algorithm

``` r
AFT_mod <- survreg(Surv(lower + 1, upper + 1, type = "interval2") ~ NTT + NTS + 
    DEG2ENTRY + EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * 
    EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    EntryCTZUSIN, data = test, dist = "exponential", robust = TRUE)
summary(AFT_mod)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(lower + 1, upper + 1, type = "interval2") ~ 
    ##     NTT + NTS + DEG2ENTRY + EntryWAPRI + EntryWKTRNI + EntryPUBPRI + 
    ##         EntryEMTP + EntryPUBPRI * EntryEMTP + SDRCARN + EntryAGE + 
    ##         GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN, 
    ##     data = test, dist = "exponential", robust = TRUE)
    ##                                                        Value Std. Err
    ## (Intercept)                                         4.319204  0.21092
    ## NTT                                                -0.936627  0.10232
    ## NTS                                                -1.085575  0.09090
    ## DEG2ENTRY                                          -0.000536  0.01686
    ## EntryWAPRIOther                                    -0.470917  0.10534
    ## EntryWAPRIResearch                                 -0.389984  0.08999
    ## EntryWKTRNITraining                                -0.026677  0.06067
    ## EntryPUBPRIPrivate                                  0.059996  0.09029
    ## EntryEMTPTwo-year                                   0.214057  0.23861
    ## EntryEMTPMed                                       -0.169671  0.09833
    ## EntryEMTPUni Research Institute                    -0.066439  0.11016
    ## SDRCARNR2                                          -0.159451  0.09874
    ## SDRCARNDoctorate                                   -0.114469  0.10154
    ## SDRCARNOther                                       -0.256709  0.21150
    ## SDRCARNMedHealth                                   -0.104528  0.13836
    ## EntryAGE                                            0.005156  0.00528
    ## GENDERFemale                                        0.034937  0.06110
    ## MINRTYYes                                          -0.019952  0.07917
    ## EntryMARINDYes                                      0.100431  0.07207
    ## EntryCHLVINYes                                      0.016524  0.07298
    ## EntryCTZUSINCitizen                                 0.065453  0.06982
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year               -1.349057  0.42357
    ## EntryPUBPRIPrivate:EntryEMTPMed                    -0.164325  0.13993
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute -0.290719  0.17782
    ##                                                    (Naive SE)        z
    ## (Intercept)                                           0.21062  20.4776
    ## NTT                                                   0.10267  -9.1543
    ## NTS                                                   0.09300 -11.9422
    ## DEG2ENTRY                                             0.01656  -0.0318
    ## EntryWAPRIOther                                       0.10588  -4.4706
    ## EntryWAPRIResearch                                    0.09184  -4.3334
    ## EntryWKTRNITraining                                   0.05943  -0.4397
    ## EntryPUBPRIPrivate                                    0.08730   0.6645
    ## EntryEMTPTwo-year                                     0.23464   0.8971
    ## EntryEMTPMed                                          0.09738  -1.7255
    ## EntryEMTPUni Research Institute                       0.10734  -0.6031
    ## SDRCARNR2                                             0.09698  -1.6148
    ## SDRCARNDoctorate                                      0.09840  -1.1274
    ## SDRCARNOther                                          0.20720  -1.2138
    ## SDRCARNMedHealth                                      0.14194  -0.7555
    ## EntryAGE                                              0.00523   0.9761
    ## GENDERFemale                                          0.05983   0.5718
    ## MINRTYYes                                             0.07685  -0.2520
    ## EntryMARINDYes                                        0.07001   1.3936
    ## EntryCHLVINYes                                        0.07131   0.2264
    ## EntryCTZUSINCitizen                                   0.06909   0.9374
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year                  0.47780  -3.1850
    ## EntryPUBPRIPrivate:EntryEMTPMed                       0.13752  -1.1744
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute    0.17041  -1.6349
    ##                                                           p
    ## (Intercept)                                        3.41e-93
    ## NTT                                                5.47e-20
    ## NTS                                                7.13e-33
    ## DEG2ENTRY                                          9.75e-01
    ## EntryWAPRIOther                                    7.80e-06
    ## EntryWAPRIResearch                                 1.47e-05
    ## EntryWKTRNITraining                                6.60e-01
    ## EntryPUBPRIPrivate                                 5.06e-01
    ## EntryEMTPTwo-year                                  3.70e-01
    ## EntryEMTPMed                                       8.44e-02
    ## EntryEMTPUni Research Institute                    5.46e-01
    ## SDRCARNR2                                          1.06e-01
    ## SDRCARNDoctorate                                   2.60e-01
    ## SDRCARNOther                                       2.25e-01
    ## SDRCARNMedHealth                                   4.50e-01
    ## EntryAGE                                           3.29e-01
    ## GENDERFemale                                       5.67e-01
    ## MINRTYYes                                          8.01e-01
    ## EntryMARINDYes                                     1.63e-01
    ## EntryCHLVINYes                                     8.21e-01
    ## EntryCTZUSINCitizen                                3.49e-01
    ## EntryPUBPRIPrivate:EntryEMTPTwo-year               1.45e-03
    ## EntryPUBPRIPrivate:EntryEMTPMed                    2.40e-01
    ## EntryPUBPRIPrivate:EntryEMTPUni Research Institute 1.02e-01
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -4127.6   Loglik(intercept only)= -4344.7
    ##  Chisq= 434.2 on 23 degrees of freedom, p= 0 
    ## (Loglikelihood assumes independent observations)
    ## Number of Newton-Raphson Iterations: 6 
    ## n=7905 (670 observations deleted due to missingness)

### AFT Specification

``` r
table <- as.data.frame(summary(AFT_mod)$table)
rownames(table) <- c("Intercept", "Non-tenure Track", "No Tenure System", "Time between Degree and Job", 
    "Administration/Other", "Research", "Workplace Training", "Private Control", 
    "Two-year", "Medical", "Research Institute", "PhD Research II", "PhD Doctorate Institution", 
    "PhD Other", "PhD Medical/Health", "Age", "Female", "Minority", "Married", 
    "Parent", "Citizen", "Private x Two-Year", "Private x Medical", "Private x Research Institute")
table$expCoef <- exp(table$Value)
names(table)[names(table) == "Std. Err"] <- "Robust SE"
aft_table <- round(table[, c("Value", "expCoef", "Robust SE", "z", "p")], 3)
kable(aft_table)
```

|                              |   Value|  expCoef|  Robust SE|        z|      p|
|------------------------------|-------:|--------:|----------:|--------:|------:|
| Intercept                    |   4.319|   75.129|      0.211|   20.478|  0.000|
| Non-tenure Track             |  -0.937|    0.392|      0.102|   -9.154|  0.000|
| No Tenure System             |  -1.086|    0.338|      0.091|  -11.942|  0.000|
| Time between Degree and Job  |  -0.001|    0.999|      0.017|   -0.032|  0.975|
| Administration/Other         |  -0.471|    0.624|      0.105|   -4.471|  0.000|
| Research                     |  -0.390|    0.677|      0.090|   -4.333|  0.000|
| Workplace Training           |  -0.027|    0.974|      0.061|   -0.440|  0.660|
| Private Control              |   0.060|    1.062|      0.090|    0.664|  0.506|
| Two-year                     |   0.214|    1.239|      0.239|    0.897|  0.370|
| Medical                      |  -0.170|    0.844|      0.098|   -1.726|  0.084|
| Research Institute           |  -0.066|    0.936|      0.110|   -0.603|  0.546|
| PhD Research II              |  -0.159|    0.853|      0.099|   -1.615|  0.106|
| PhD Doctorate Institution    |  -0.114|    0.892|      0.102|   -1.127|  0.260|
| PhD Other                    |  -0.257|    0.774|      0.211|   -1.214|  0.225|
| PhD Medical/Health           |  -0.105|    0.901|      0.138|   -0.755|  0.450|
| Age                          |   0.005|    1.005|      0.005|    0.976|  0.329|
| Female                       |   0.035|    1.036|      0.061|    0.572|  0.567|
| Minority                     |  -0.020|    0.980|      0.079|   -0.252|  0.801|
| Married                      |   0.100|    1.106|      0.072|    1.394|  0.163|
| Parent                       |   0.017|    1.017|      0.073|    0.226|  0.821|
| Citizen                      |   0.065|    1.068|      0.070|    0.937|  0.349|
| Private x Two-Year           |  -1.349|    0.259|      0.424|   -3.185|  0.001|
| Private x Medical            |  -0.164|    0.848|      0.140|   -1.174|  0.240|
| Private x Research Institute |  -0.291|    0.748|      0.178|   -1.635|  0.102|

The exponential distribution has a constant hazard *λ*(*t*)=*λ* and thus a survival function of *S*(*t*)=*e*<sup>−*λ*(*t*)</sup> and density of *f*(*t*)=*λ* × *e*<sup>−*λ*(*t*)</sup>. An interesting occurance is that the expected survival time for this distribution is *E*(*t*)=1/*λ* and its variance is *E*(*t*)=1/*λ*<sup>2</sup>. This makes the mean survival time equal to e^intercept (75.1287755). It's inverse (0.0133105) is the MLE of the (constant) hazard rate. Of course, the model performs poorly extrapolating to such an extreme timepoint.

AFT models are typcally interpreted in a way that covariates have a multiplicative effect on the expected survival time. So, with regard to tenure status, taking your first job as NTT accelerates the time to attrition by a factor of 0.3919478 (0.3919478 times shorter survival time compared to the baseline survival). Beginning an academic career in a non-tenure system accelerates the time to attrition by a factor of 0.3377075. The life course for these states is -60.8052173 and -66.229245 percent shorter, respectively.

The Weibull family of distributions (of which the exponential is a sub-class) has the advantage that covariates can also be interpreted as an impact on the hazard ratios. For this famiily of distributions, the coefficient is multiplied by -1 and then multiplied by a shape parameter (1/scale parameter). In the case of the exponential distributuion, the shape parameter is simply 1/1. So in our case, the hazard ratio comparing NTT to tenure-track positions is 2.55136. The risk of attrition increases by a factor of 2.55136 when one begins an academic career in a NTT position. Faculty with a first job at a non-teure system institution increases their risk by a factor of 2.961142.

Comparing the Exponential AFT Model and the Cox Models
------------------------------------------------------

``` r
aft_table2 <- rbind(aft_table[1:21, c(1, 5)], space, space, aft_table[22:24, 
    c(1, 5)])
RC_cox_table2 <- rbind(space, RC_cox_table[, c(2, 6)])
IC_cox_table2 <- rbind(space, IC_cox_table[1:20, c(2, 5)], space, space, IC_cox_table[21:23, 
    c(2, 5)])
ctable <- cbind(aft_table2, RC_cox_table2, IC_cox_table2)
rownames(ctable) <- rownames(RC_cox_table2)
rownames(ctable)[1] <- "Intercept"
ctable$Value <- round(exp(as.numeric(aft_table2$Value) * -1 * 1/1), 3)  # convert to hazard ratios
colnames(ctable) <- c("AFT:HR", "AFT p-val", "RC Cox:HR", "RC Cox p-val", "IC Cox:HR", 
    "IC Cox p-val")
kable(ctable)
```

|                              |  AFT:HR| AFT p-val | RC Cox:HR | RC Cox p-val | IC Cox:HR | IC Cox p-val |
|------------------------------|-------:|:----------|:----------|:-------------|:----------|:-------------|
| Intercept                    |   0.013| 0         |           |              |           |              |
| Non-tenure Track (NTT)       |   2.552| 0         | 2.86      | 0            | 2.617     | 0            |
| No Tenure System             |   2.962| 0         | 3.466     | 0            | 3.025     | 0            |
| Time between Degree and Job  |   1.001| 0.975     | 0.979     | 0.218        | 0.999     | 0.963        |
| Admin/Other                  |   1.602| 0         | 1.626     | 0            | 1.582     | 0            |
| Researcher                   |   1.477| 0         | 1.547     | 0            | 1.466     | 0            |
| Workplace Training           |   1.027| 0.66      | 1.021     | 0.719        | 1.009     | 0.876        |
| Private Control              |   0.942| 0.506     | 0.941     | 0.491        | 0.934     | 0.48         |
| Two-year/Other               |   0.807| 0.37      | 0.845     | 0.476        | 0.805     | 0.346        |
| Medical                      |   1.185| 0.084     | 1.136     | 0.19         | 1.147     | 0.175        |
| Research Institute           |   1.068| 0.546     | 1.053     | 0.627        | 1.047     | 0.682        |
| PhD Research II              |   1.172| 0.106     | 1.148     | 0.154        | 1.174     | 0.12         |
| PhD Doctorate Institution    |   1.121| 0.26      | 1.105     | 0.307        | 1.108     | 0.336        |
| PhD Other                    |   1.293| 0.225     | 1.258     | 0.27         | 1.324     | 0.22         |
| PhD Medical/Health           |   1.111| 0.45      | 1.08      | 0.588        | 1.093     | 0.492        |
| Age                          |   0.995| 0.329     | 0.996     | 0.402        | 0.993     | 0.177        |
| Female                       |   0.966| 0.567     | 0.925     | 0.187        | 0.966     | 0.571        |
| Minority                     |   1.020| 0.801     | 1.027     | 0.731        | 1.007     | 0.926        |
| Married                      |   0.905| 0.163     | 0.914     | 0.2          | 0.91      | 0.248        |
| Children                     |   0.983| 0.821     | 0.994     | 0.927        | 0.975     | 0.704        |
| Citizen                      |   0.937| 0.349     | 1.06      | 0.394        | 0.925     | 0.313        |
| Time x NTT                   |      NA|           | 0.819     | 0.079        |           |              |
| Time x No Tenure System      |      NA|           | 0.708     | 0            |           |              |
| Private x Two-Year           |   3.854| 0.001     | 2.958     | 0.024        | 3.292     | 0.007        |
| Private x Medical            |   1.178| 0.24      | 1.186     | 0.218        | 1.179     | 0.255        |
| Private x Research Institute |   1.338| 0.102     | 1.346     | 0.08         | 1.323     | 0.12         |

The tenure status coefficients associated with the exponential AFT model are a little bit less, but generally the results across these models are comparable.

Conclusions
-----------

1.  Tenure status at hiring related to attrition
2.  Attrition more prevalent among NTT and NTS faculty, however, large numbers make a career of it
3.  Attrition is more of a risk for researchers and administrators, not teachers
4.  Relationships likely underfit by the model. Acquire more faculty characteristics, particularly time-varying characteristics.
