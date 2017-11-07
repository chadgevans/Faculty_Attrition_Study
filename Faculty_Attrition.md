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
-   [Accelerated Failure Time Models](#accelerated-failure-time-models)
    -   Exponential
    -   Lognormal
    -   [AFT Specification](#aft-specification)
-   [Comparing the Exponential AFT Model and the Cox Models](#comparing-the-exponential-aft-model-and-the-cox-models)
-   [Final Model](#final-model)
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

First, lets write a function to create our life tables.

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

This first Life Table will tabulate the survival probabilities and hazards for all faculty (aggregated).

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

The following graph plots the actuarial survival probabilities for each class of faculty.

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

These survival probabilities can also be estimated and graphed using the KM estimator. In this chapter, I will be reporting only the actuarial estimator, as it is more robust for larger samples.

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

The default confidence interval (above) is log, which calculates intervals based on log(survival).

### Statistical Test for differences

Often it is useful to have a statistical test for whether survival curves of two (or more) groups differ. Let's look at the most common test, the log-rank test (a.k.a. the Mantel-Haenszel test).

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

The null hypothesis for a log-rank test is that the groups have the same survival. In this case, the Chi-square statistic is significant for tenure status of first appointment. The survival curve, therefore, differs based on the tenure status of the PhD recipient's first academic job.

#### Summary of Nonparametric models

Substantially, the actuarial estimator and KM estimator tell similar stories. Attrition is much higher and more rapid for faculty working off the tenure track. The lifetable estimator perhaps shows a slightly higher attrition rate for NTS faculty, comparted to the KM estimate.

Regression Models for Survival Data
-----------------------------------

There are two types of regression models in survival analysis, semi-parametric and fully-parametric. In this study, we will consider both approaches. For semi-parametric analysis, I'll utilize the cox proportional hazards model. I will first treat the censoring as right censoring. Then I will treat it as interval censored. Technically, these data are interval censored, however, there is not a good package yet for R that allows for interval censoring and also allows for interactions between tenure status and time. This interaction is crucial, because it allows for the cox model to handle non-proportional hazards. As Allison says, "Testing is both the diagnosis and the solution."

For the parametric analysis (AFT models), we'll try different functional forms. At the end, I compare the results of the cox models with the results of the exponential form.

For specifying our models, we will subset the data and first train models using training data. After we are satisfied with that specification, we'll test all models using the independent test data. This will help prevent overfitting.

``` r
smp_size <- floor(0.6 * nrow(data))
set.seed(777)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
```

### Cox Proportional Hazards Model

#### Right Censored Data Model

First, we fit a simple Cox model predicting attrition from academia from a categorical predictor of tenure status at workforce entry. We employ the efron method of dealing with ties, although other popular methods (e.g., breslow method) are available. The Efron approximation is more accurate when dealing with tied death times, and is as efficient computationally.

The Surv function here frames censoring as follows. It assumes that individuals entered the study halfway through the previous interval. The extra year assumption was necessary to avoid "spontaneous attrition" that would have removed from the analysis any indiviudal who entered and exited during the same interval. As these early attriters are so very important to this study, we assume they entered halfway into the previous interval. This method treats censoring from the "EndDate"--that is, the last year he or she was still employed. Censoring, of course, happens after that point, so an individual who enters in year 0, remains in the study through year 2 in its entirety and subsequently drops out between year 2 and 4 would be coded as "3+." There are two years for the full interval, plus the extra entry year (assumption).

``` r
RC_Mod1 <- coxph(Surv(TIME2, Censor) ~ NTT + NTS, data = test, method = "efron")  # breslow option available
summary(RC_Mod1)
```

Proportional Hazards models and AFT models must be interpreted in different ways. AFT models give the percentage change in survival time. Cox PH models give the percentage change to the hazard at all time points, following this formula: *h*(*t*)=*h*<sub>0</sub>(*t*)*e*<sup>(*β*′*x*)</sup>

In this case, the effect of initial tenure status on time to attrition has an estimated coefficient of 1.0089259 and 1.3237404. Exponentiated, this means that subjects appointed to lower tenure-status jobs multiply their baseline hazards *h*<sub>0</sub>(*t*) by a factor of 2.7426537 and 3.7574493. Their "risk"" of attriting from academia is 174.2653671 and 275.744933 percent higher than academics who begin immediately on the tenure-track. Importantly, Cox models state that this is the impact on the subject's hazard at any given time, t. It does not, however, imply an expansion (or contraction) of the lifespan of the subject.

An important limitation of this model is that it does not account for the possibility of non-proportional hazards based on tenure status. Let's test whether the proportional hazards assumption may have been violated in this analysis.

``` r
cox.zph <- cox.zph(RC_Mod1, transform = "km")
round(cox.zph$table, 3)
```

    ##           rho  chisq     p
    ## NTT    -0.052  3.716 0.054
    ## NTS    -0.103 14.291 0.000
    ## GLOBAL     NA 15.170 0.001

This function tests proportionality by interacting each predictor with log time (km transformation). Rho is the pearson product-moment correlation between the scaled residuals (Schoenfeld) and log(time) for each covariate. The global test jointly tests all of the interactions. Low p-values suggest a violation of the assumption of proportional hazards.

Now let's build a more comprehensive model. Importantly we'll include an interaction between time and our tenure status variables (NTT or NTS dummies) to test the assumption of proportional hazards and also correct for non-proportional hazards. Here we use the TIME2 variable which simply adds an extra year to everyone's event time.

``` r
RC_Mod2 <- coxph(Surv(TIME2, Censor) ~ NTT + NTS + DEG2ENTRY + EntryWAPRI + 
    EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + SDRCARN + 
    EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN + 
    tt(NTT) + tt(NTS), data = test, method = "efron", robust = TRUE, tt = function(x, 
    t, ...) x * log(t))
summary(RC_Mod2)
```

Even after controlling for background characteristics, there are significant differences. Here, NTT status or working at a college or university without a tenure system impacts the hazard, multiplying the baseline by a factor of 2.8598117 and 3.4658201. This is equivalent to saying that each tenure status increases the hazard of attrition by 185.9811676 and 246.5820096 percent, controlling for background characteristics. R output also provides the exponentiated negative coefficient. To my understanding, that just allows you to compare the groups relative to the baseline hazard of the tenure-track group. Robust standard errors were used in this model.

The model reveals other important predcitors, including the subject's main job and possible interactions between institution type and public/private status.

#### Interval Censored Data Model

The data used in this study not only have right censoring, they are interval censored. We know the interval during when an event took place, but we do not know the particular time of the event. Unfortunately, the well-developed and supported function coxph() from the survival package does not handle interval censoring. I searched for a different R package capable of handling interval censoring and found the icenreg package. There is some documentation for this package, however, it is considerably less developed than functions in the suvival package.

We'll follow the same model development approach as in the last specification. We'll try to develop a complex model using training data and then apply this model to test data.

Using this method does not require the assumption that the subject enters the study midway through the previous wave. It only depends on the interval itself.

``` r
IC_Mod1 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS, model = "ph", 
    bs_samples = 100, data = data)
summary(IC_Mod1)
```

``` r
IC_Mod1 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS + DEG2ENTRY + 
    EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + 
    SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + EntryCTZUSIN, 
    model = "ph", bs_samples = 100, data = test)
summary(IC_Mod1)
```

``` r
# This is time-intensive at bs_samples=10, maybe a minute per sample.
IC_Mod2 <- ic_sp(Surv(lower, upper, type = "interval2") ~ NTT + NTS + DEG2ENTRY + 
    EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * EntryEMTP + 
    EntryPUBPRI * EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + 
    EntryCHLVIN + EntryCTZUSIN + TIME * NTT + TIME * NTS, model = "ph", bs_samples = 10, 
    data = data)
summary(IC_Mod2)
# Also, this model is wrong because there should not be a coefficient
# produced for time.  Furthermore, that coefficient sops up all the
# variation in event times
```

#### Summary of Cox Proportional Hazard Models

If the proportional hazards assumption holds, then it is possible to estimate the effect of parameter(s) without any consideration of the baseline hazard function. As mentioned, this is in contrast to parametric models--the focus of our next section. Even when proportional hazards are not truly proportional, we can include an interaction between the variable and time to correct for it.

### Accelerated Failure Time Models

First, we assume the outcome has an exponential distribution--a good baseline distribution to start with (simplifies calculations). The exponential distribution implies a constant hazard rate. The advantage of the Weibull (and by extention the exponential), is that it can be parameterised as a PH model or an AFT model. In other words, the Weibull family can be interpreted as affecting the risk of event occurance or the duration of the lifespan.

For fully parameterized models, the timing of an event matters (unlike the cox models that only require the order). So here we must also make an assumption that subject begin their academic work midway through the previous interval. So we add unity to the upper and lower interval bounds to adjust for this.

#### Mathematics

1.  Construct the likelihood function
    -   uncensored cases contribute their density probability at t
    -   censored cases contribute their survival probability at t

2.  Simplify the expression
3.  Take the log to simplify the math
4.  Calculate the derivative and set it equal to zero
5.  Solve using, typically, the Newton-Raphson algorithm

``` r
Exp_mod <- survreg(Surv(lower + 1, upper + 1, type = "interval2") ~ NTT + NTS + 
    DEG2ENTRY + EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * 
    EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    EntryCTZUSIN, data = test, dist = "exponential", robust = TRUE)
```

``` r
loglog_mod <- survreg(Surv(lower + 1, upper + 1, type = "interval2") ~ NTT + 
    NTS + DEG2ENTRY + EntryWAPRI + EntryWKTRNI + EntryPUBPRI + EntryEMTP + EntryPUBPRI * 
    EntryEMTP + SDRCARN + EntryAGE + GENDER + MINRTY + EntryMARIND + EntryCHLVIN + 
    EntryCTZUSIN, data = test, dist = "loglogistic", robust = TRUE)
```

``` r
summary(Exp_mod)$loglik[2]
```

    ## [1] -4127.559

``` r
summary(loglog_mod)$loglik[2]
```

    ## [1] -4039.913

Because the log-likelihood is higher (less negative), the log logistic model actually fits the data better. The higher logliklihood simply means that that the probability of the data is marginally closer to 1 (certitude). However, the exponential distribution fits the data similarly. In addition, parameterizing log(time) as an exponential function allows me to convert coefficients to hazards ratios. I can therefore compare the results to the findings of the cox models. As the final analysis will be a cox model, I use the exponential AFT for validation purposes.

### AFT Specification

``` r
table <- as.data.frame(summary(Exp_mod)$table)
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

The exponential distribution has a constant hazard *λ*(*t*)=*λ* and thus a survival function of *S*(*t*)=*e*<sup>−*λ*(*t*)</sup> and density of *f*(*t*)=*λ* × *e*<sup>−*λ*(*t*)</sup>. An interesting occurance is that the expected survival time for this distribution is *E*(*t*)=1/*λ* and its variance is *E*(*t*)=1/*λ*<sup>2</sup>. This makes the mean survival time equal to e^intercept (75.1287755). It's inverse (0.0133105) is the MLE of the (constant) hazard rate. Of course, the model performs poorly extrapolating to such an extreme timepoint.

AFT models are typcally interpreted in a way that covariates have a multiplicative effect on the expected survival time. So, with regard to tenure status, taking your first job as NTT accelerates the time to attrition by a factor of 0.3919478 (0.3919478 times shorter survival time compared to the baseline survival). Beginning an academic career in a non-tenure system accelerates the time to attrition by a factor of 0.3377075. The life course for these states is -60.8052173 and -66.229245 percent shorter, respectively.

The Weibull family of distributions (of which the exponential is a sub-class) has the advantage that covariates can also be interpreted as an impact on the hazard ratios. For this famiily of distributions, the coefficient is multiplied by -1 and then multiplied by a shape parameter (1/scale parameter). In the case of the exponential distributuion, the shape parameter is simply 1/1. So in our case, the hazard ratio comparing NTT to tenure-track positions is 2.55136. The risk of attrition increases by a factor of 2.55136 when one begins an academic career in a NTT position. Faculty with a first job at a non-teure system institution increases their risk by a factor of 2.961142.

Comparing the Exponential AFT Model and the Cox Models
------------------------------------------------------

``` r
space <- c("", "")
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

The results of the AFT model (exponential) look similar to the results from the cox models. This helps to validate the findings.

Comparing the RC and IC cox models, we see that the models differ slightly. The Right censored (RC) model includes the time interaction with tenure status, but it fails to deal with the interval censoring of the data and assume every individual entered his or her job midway through the previous interval. The Interval censored (IC) model effectivey handles the interval censoring, but it does not allow for an interaction between time and tenure status. Nevertheless, the results are comparable.

The fit statistics of these two models is difficult to reconcile. The RC model has an Rsquare of 0.047. The log likelihood is listed as -1.013048210^{4}, -9939.400298. I'm not sure what these two numbers mean. Perhaps one is only for the intercept? The IC model reports a log likelihood of -3639.0884452, which is considerably different from the RC model.

The tenure status coefficients associated with the exponential AFT model are a little bit less, but generally the results across these models are comparable.

Final Model
-----------

The final model for this chapter is a right censored Cox proportional hazards model. It is true that the data are also interval censored, however, the intervals are very close to being of equal length. As a result, we can treat this as a discrete time analyisis with right censored data. Importantly, the right censored model allows us to include a timme interaction. This effectively allows us to use the cox proportional hazards model even though the hazards of tenure-track and non-tenure track faculty are not proportional at one or more time points. This is particularly important during the years of tenure review, as the hazard of attrition likely shifts significantly for those faculty on the tenure line.

``` r
kable(RC_cox_table)
```

|                              |    coef|  exp(coef)|  se(coef)|  robust se|       z|  Pr(&gt;|z|)|
|------------------------------|-------:|----------:|---------:|----------:|-------:|------------:|
| Non-tenure Track (NTT)       |   1.051|      2.860|     0.142|      0.139|   7.537|        0.000|
| No Tenure System             |   1.243|      3.466|     0.127|      0.123|  10.075|        0.000|
| Time between Degree and Job  |  -0.022|      0.979|     0.017|      0.017|  -1.232|        0.218|
| Admin/Other                  |   0.486|      1.626|     0.106|      0.102|   4.757|        0.000|
| Researcher                   |   0.437|      1.547|     0.091|      0.087|   5.015|        0.000|
| Workplace Training           |   0.021|      1.021|     0.059|      0.059|   0.360|        0.719|
| Private Control              |  -0.060|      0.941|     0.087|      0.088|  -0.689|        0.491|
| Two-year/Other               |  -0.169|      0.845|     0.234|      0.237|  -0.713|        0.476|
| Medical                      |   0.127|      1.136|     0.097|      0.097|   1.311|        0.190|
| Research Institute           |   0.052|      1.053|     0.107|      0.107|   0.485|        0.627|
| PhD Research II              |   0.138|      1.148|     0.097|      0.097|   1.427|        0.154|
| PhD Doctorate Institution    |   0.100|      1.105|     0.098|      0.098|   1.021|        0.307|
| PhD Other                    |   0.230|      1.258|     0.206|      0.208|   1.102|        0.270|
| PhD Medical/Health           |   0.077|      1.080|     0.141|      0.142|   0.542|        0.588|
| Age                          |  -0.004|      0.996|     0.005|      0.005|  -0.838|        0.402|
| Female                       |  -0.078|      0.925|     0.060|      0.059|  -1.319|        0.187|
| Minority                     |   0.026|      1.027|     0.077|      0.077|   0.344|        0.731|
| Married                      |  -0.090|      0.914|     0.070|      0.070|  -1.283|        0.200|
| Children                     |  -0.006|      0.994|     0.071|      0.071|  -0.091|        0.927|
| Citizen                      |   0.058|      1.060|     0.069|      0.068|   0.852|        0.394|
| Time x NTT                   |  -0.200|      0.819|     0.117|      0.114|  -1.755|        0.079|
| Time x No Tenure System      |  -0.345|      0.708|     0.100|      0.098|  -3.508|        0.000|
| Private x Two-Year           |   1.085|      2.958|     0.476|      0.482|   2.252|        0.024|
| Private x Medical            |   0.170|      1.186|     0.137|      0.138|   1.233|        0.218|
| Private x Research Institute |   0.297|      1.346|     0.170|      0.170|   1.752|        0.080|

``` r
summary(RC_Mod2)$rsq  # Fit statistic
```

    ##        rsq     maxrsq 
    ## 0.04719445 0.92293120

``` r
write.csv(RC_cox_table, file.path(Graphs, "RC_cox_table.csv"))
```

R output also reports the concordance, likelihood ratio test, wald test and score (logrank) test.

Conclusions
-----------

1.  Tenure status at hiring related to attrition
2.  Attrition more prevalent among NTT and NTS faculty, however, large numbers make a career of it
3.  Attrition is more of a risk for researchers and administrators, not teachers
4.  It might be useful to look into private, two-year institutions
5.  Relationships likely underfit by the model. Acquire more faculty characteristics, particularly time-varying characteristics.
