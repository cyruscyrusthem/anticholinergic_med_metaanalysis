---
title: "Anticholinergic Meta analysis"
output: 
  html_document:
    keep_md: true
---

### Load packages

```r
library(meta)
library(metafor) #As a technical note, with exception of the DerSimonian–Laird and the Paule–Mandel methods the rma.uni function of R package metafor is called internally in the metagen function. Thus, it is a good idea to install R package metafor to make all estimation methods available.
library(tidyverse)
```

### Data import
Import data (.csv) 

```r
library(readr)
Anticholinergic_Medication_R_data <- read_csv("Z:/PROJECTS/2018_Anticholinergic_Med_SR/Analysis/Anticholinergic_Medication_R_data.csv")
```

Assign data to `dat`

```r
dat <- Anticholinergic_Medication_R_data %>%
  filter(include == 1) #%>% #this currently filters out the 2 Stevenson outcomes where px were no longer on medication. This decision does need to be checked, though. 
  #filter(quality > 4) # this removes all studies with quality below the median (5)
```

### Meta-analysis
#### Step 1: Average across studies
This step is necessary to ensure independence of each effect size entered into the meta-analysis. In the current form, most studies report multiple outcomes (time-points, medication, cognitive tests), and these results are dependent on each other, as they come from they come from overlapping samples of participants. 

Data are averaged **within studies** so that each study represents a single effect size.
In all cases of studies which included outcomes for different medications, all medications were either low or high potency, so outcomes were able averaged together (within the study).

Averaged results are assigned to `dat_study`


```r
dat_study <- dat %>%
  group_by(study) %>% #group by study so further calculations are done within each study
  mutate(g = mean(g, na.rm = T),
          st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study, removing NAs from calculation
  filter(row_number()==1) %>% #make each study only appear on one row
  select(study, author, year, Potency, g, st_err) #select relevant data for analysis
```

#### Step 2: Run random effects model of all data
**Model 1: meta package**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
meta_all_mod1 <- metagen(g,
                    st_err,
                    data = dat_study,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") #calculate SMD
```


Code for forest plot:

```r
meta::forest(meta_all_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightcols = c("TE", "ci"),
             rightlabs = c("SMD", "95% CI"),
             text.random = "Overall effect") #found this here: https://rdrr.io/cran/meta/man/forest.html
```

![](meta_analysis_files/figure-html/forest mod1-1.png)<!-- -->

**Model 2: meta package**

* Use the DerSimonian-Laird method (default) to estimate tau
* Not using Knapp-Hartung method


```r
meta_all_mod2 <- metagen(g,
                          st_err,
                          data = dat_study,
                          studlab = paste(author, year),
                          comb.fixed = F,
                          comb.random = T,
                          hakn = F, # not using the Knapp-Hartung method
                          prediction = T, #True = print prediction interval for future studies based on present evidence
                          sm = "SMD")
```

Code for forest plot (not displayed):

```r
meta::forest(meta_all_mod2,
             sortvar = TE,
             xlim = c(-1.5, 1.5),
             leftcols = "studlab") #found this here: https://rdrr.io/cran/meta/man/forest.html
```

**Model 3: meta package**

* Use the Restricted Maximum-Likelihood method to estimate tau
* Not using Knapp-Hartung method


```r
meta_all_mod3 <- metagen(g,
                          st_err,
                          data = dat_study,
                          studlab = paste(author, year),
                          comb.fixed = F,
                          comb.random = T,
                          method.tau = "REML",
                          hakn = F, # not using the Knapp-Hartung method
                          prediction = T, #True = print prediction interval for future studies based on present evidence
                          sm = "SMD")
```

Code for forest plot (not displayed):

```r
meta::forest(meta_all_mod2,
             sortvar = TE,
             xlim = c(-1.5, 1.5),
             leftcols = "studlab") #found this here: https://rdrr.io/cran/meta/man/forest.html
```

Compare the models:

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.0345117 |-0.035378 | 0.1044014 | 0.9939728 | 0.3254361
2 | 0.051671 | -0.0147278 | 0.1180698 | 1.5252269 | 0.1272025
3 | 0.0431682 | -0.0329675 | 0.119304 | 1.1112806 | 0.2664476


#### Step 3: Run subgroup analyses
##### *Potency (low/high)*
Run sub-analyses by medication potency (low/high) for each model


```r
#Model 1
potency_subgroup_mod1 <- update.meta(meta_all_mod1, 
                             byvar=Potency, 
                             comb.random = TRUE, 
                             comb.fixed = FALSE)
#Model 2
potency_subgroup_mod2 <- update.meta(meta_all_mod2, 
                                   byvar=Potency, 
                                   comb.random = TRUE, 
                                   comb.fixed = FALSE)
```
Does not seem to be able to run sub-analysis using REML method:
`Error = Error in rma.uni(yi = TE[sel], sei = seTE[sel], method = method.tau, control = control) :Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.`

Code for forest plot (not shown):

```r
meta::forest(potency_subgroup_mod1,
             sortvar = TE,
             leftcols = "studlab")
```

Potency = **low**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 37 | 0.0271271 | -0.0483943 | 0.1026485 | 0.4710284 | 0.0284255 | 28.0191867
2 | 37 | 0.028565 | -0.0443065 | 0.1014365 | 0.4423155 | 0 | 28.0191867

Potency = **high**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 10 | 0.0645153 | -0.1481803 | 0.2772108 | 0.5098995 | 0.0319191 | 12.5042251
2 | 10 | 0.0646297 | -0.1550929 | 0.2843523 | 0.5642716 | 0.0318397 | 12.5042251

##### *Cognitive domain*

Create data set (`dat_study_domain`) which has the average effect size by cognitive domain within each study.

Cognitive domains were defined using domains outlined by Lezak (2012): 

* Memory
* Executive function
* Attention
* Psychomotor functioning
* Concept formation & reasoning
* Language
* Intelligence
* Perception

Effect sizes based on cognitive composite scores (reported by Stevenson, Robles, Operto) were not included in the cognitive domain sub analysis. 


```r
dat_study_domain <- dat %>%
  group_by(study, cog_domain_lezak) %>% #group by study and cognitive domains (based on Lezak) within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  filter(cog_domain_lezak != "Not Subdomain") %>% #remove outcomes based on cognitive composite scores
  select(study, author, year, Potency, cog_domain_lezak, g, st_err) #select relevant data
```

Then run random effects model by subgroup (cognitive domain).

**Model 1**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
meta_domain_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_domain,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
domain_subgroup_mod1 <- update.meta(meta_domain_mod1, 
                                byvar=cog_domain_lezak, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | 38 | 0.0322832 | -0.0470373 | 0.1116038 | 0.4148524
Psychomotor Functioning | 17 | -0.1170568 | -0.3365868 | 0.1024733 | 0.2749827
Concept Formation & Reasoning | 13 | 0.1021407 | -0.060555 | 0.2648364 | 0.1964251
Perception | 3 | 0.2500101 | -0.900842 | 1.4008621 | 0.4486141
Memory | 16 | 0.0429447 | -0.0739579 | 0.1598474 | 0.4458151
Executive Function | 15 | -0.0198762 | -0.2863981 | 0.2466458 | 0.8752055
General | 14 | 0.0791189 | -0.174362 | 0.3325999 | 0.5119244
Language | 6 | 0.082419 | -0.1117334 | 0.2765713 | 0.3249427


Code for forest plot (not shown):

```r
meta::forest(domain_subgroup_mod1,
             sortvar = TE,
             leftcols = "studlab")
```

**Model 2**


```r
meta_domain_mod2 <- metagen(g,
                       st_err,
                       data = dat_study_domain,
                       studlab = paste(author, year),
                       comb.fixed = F,
                       comb.random = T,
                       hakn = F, # not using the Knapp-Hartung method
                       prediction = T, #True = print prediction interval for future studies based on present evidence
                       sm = "SMD") # says we want to calculate SMD
domain_subgroup_mod2 <- update.meta(meta_domain_mod2, 
                               byvar=cog_domain_lezak, 
                               comb.random = TRUE, 
                               comb.fixed = FALSE)
```

Code for forest plot (not shown):

```r
meta::forest(domain_subgroup_mod2,
             sortvar = TE,
             leftcols = "studlab")
```

Domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | 38 | 0.0457063 | -0.0310838 | 0.1224964 | 0.2433751
Psychomotor Functioning | 17 | -0.0996516 | -0.2919519 | 0.0926487 | 0.3097868
Concept Formation & Reasoning | 13 | 0.1364551 | -0.0098462 | 0.2827565 | 0.0675411
Perception | 3 | 0.2436719 | -0.2697523 | 0.7570961 | 0.3522666
Memory | 16 | 0.0393785 | -0.0750444 | 0.1538013 | 0.4999816
Executive Function | 15 | -0.0036663 | -0.2055026 | 0.1981701 | 0.9715999
General | 14 | 0.0834596 | -0.1449501 | 0.3118694 | 0.4738926
Language | 6 | 0.1141191 | -0.0311372 | 0.2593755 | 0.123603


### Assess publication bias
Steps (based on Amit's paper):

1. Visually inspect funnel plots (effect size vs standard error) for bias
2. **If 10 or more studies**, use Egger's Test to test for small-study effect
    + If statistically significant asymmetry found (1-sided P<0.1), use Duval and Tweedie's Trima nd Fill method to quantify magnitude of bias
3. **If <10 studies**, locate outliers in funnel plot and recalculate effect sizes after their removal

To run Egger's, load Egger's function (found [here](https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/eggers.test.R))


Note: it does not matter which meta-analysis model the funnel plot or Egger's test are conducted based on, as only the individual study effect sizes matter, not the pooled effect size.
HOWEVER, models must be run separately for Trim and Fill.

#### Whole meta-analysis

```r
funnel.meta(meta_all_mod1)
```

![](meta_analysis_files/figure-html/funnel meta all-1.png)<!-- -->


```r
eggers.test(x = meta_all_mod1)
```

```
##              Intercept ConfidenceInterval      t     p
## Egger's test    -0.478        -1.066-0.11 -1.476 0.147
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_all_mod1_trim <- trimfill.meta(meta_all_mod1)
meta_all_mod2_trim <- trimfill.meta(meta_all_mod2)
```


```r
funnel.meta(meta_all_mod1_trim)
```

![](meta_analysis_files/figure-html/funnel trim full meta-1.png)<!-- -->

```r
funnel.meta(meta_all_mod2_trim)
```

![](meta_analysis_files/figure-html/funnel trim full meta-2.png)<!-- -->

Results of trim and fill analysis (with 12 studies imputed)

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.1101953 |0.0290584 | 0.1913322 | 2.7206831 | 0.0086631
2 | 0.1100384 | 0.0339308 | 0.1861461 | 2.8337681 | 0.0046003

#### Potency
Does not seem possible to plot (using funnel plots) the individual parts of the sub-analyses separately. So, first I will run the meta-analyses separately for each subgroup (code not shown). This will provide the same results as previously found, it just requires extra data wrangling (filtering by subgroup), and more code. 



**Low potency**

```r
funnel.meta(meta_low_mod1)
```

![](meta_analysis_files/figure-html/funnel low pot-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_low_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.084         -0.868-0.7 -0.236 0.81469
```

P>0.1 so nothing further required, no indication of substantial publication bias.



**High potency**

```r
funnel.meta(meta_high_mod1)
```

![](meta_analysis_files/figure-html/funnel high pot-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_high_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.734       -2.91--0.558 -3.136 0.01388
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_high_mod1_trim <- trimfill.meta(meta_high_mod1)
meta_high_mod2_trim <- trimfill.meta(meta_high_mod2)
```


```r
funnel.meta(meta_high_mod1_trim)
```

![](meta_analysis_files/figure-html/funnel trim high pot-1.png)<!-- -->

```r
funnel.meta(meta_high_mod2_trim)
```

![](meta_analysis_files/figure-html/funnel trim high pot-2.png)<!-- -->

Results of trim and fill analysis (with **XX** studies imputed)

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.2650632 |0.0056982 | 0.5244283 | 2.191907 | 0.0457893
2 | 0.2684451 | 0.0424758 | 0.4944144 | 2.3283814 | 0.0198919



#### Cognitive domains
Again, run subgroup analyses individually first (code not displayed).



**Attention**

```r
funnel.meta(meta_att_mod1)
```

![](meta_analysis_files/figure-html/funnel att-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_att_mod1)
```

```
##              Intercept ConfidenceInterval      t      p
## Egger's test    -0.317       -1.101-0.467 -0.866 0.3923
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Memory**

```r
funnel.meta(meta_mem_mod1)
```

![](meta_analysis_files/figure-html/funnel mem-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_mem_mod1)
```

```
##              Intercept ConfidenceInterval     t       p
## Egger's test     0.146        -0.638-0.93 0.393 0.70029
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Executive function**

```r
funnel.meta(meta_exec_mod1)
```

![](meta_analysis_files/figure-html/funnel exec-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_exec_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.152        -2.72-0.416 -1.358 0.19745
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Psychomotor functioning**

```r
funnel.meta(meta_psych_mod1)
```

![](meta_analysis_files/figure-html/funnel psych-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_psych_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -3.083      -4.651--1.515 -3.654 0.00235
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_psych_mod1_trim <- trimfill.meta(meta_psych_mod1)
meta_psych_mod2_trim <- trimfill.meta(meta_psych_mod2)
```


```r
funnel.meta(meta_psych_mod1_trim)
```

![](meta_analysis_files/figure-html/funnel trim psych-1.png)<!-- -->

```r
funnel.meta(meta_psych_mod2_trim)
```

![](meta_analysis_files/figure-html/funnel trim psych-2.png)<!-- -->

Results of trim and fill analysis (with **XX** studies imputed)

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.1146913 |-0.1651765 | 0.394559 | 0.8498842 | 0.404544
2 | 0.1147141 | -0.0907476 | 0.3201757 | 1.0942938 | 0.2738261

**DECIDE WHETHER TO PRINT THE RESULTS PRE TRIM AND FILL NEXT TO THIS AS WELL**



**Concept formation**

```r
funnel.meta(meta_conc_mod1)
```

![](meta_analysis_files/figure-html/funnel conc-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_conc_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.386      -2.366--0.406 -2.566 0.02626
```

**FURTHER ACTION REQUIRED**



**Language**

```r
funnel.meta(meta_lang_mod1)
```

![](meta_analysis_files/figure-html/funnel lang-1.png)<!-- -->

<10 studies so Egger's not done

**FURTHER ACTION - assess visually for outliers**



**General**

```r
funnel.meta(meta_gen_mod1)
```

![](meta_analysis_files/figure-html/funnel gen-1.png)<!-- -->

Egger's:

```r
eggers.test(x = meta_gen_mod1)
```

```
## Warning in metabias.meta(x, k.min = 3, method = "linreg"): 1 observation(s)
## dropped due to missing values
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -3.024      -5.572--0.476 -2.346 0.03701
```

**FURTHER ACTION REQUIRED**



**Perception**

```r
funnel.meta(meta_perc_mod1)
```

![](meta_analysis_files/figure-html/funnel perc-1.png)<!-- -->

<10 studies so Egger's not done

**FURTHER ACTION - assess visually for outliers**

### Models in metafor
To contrast with analyses in the `meta` package, we will now also run analyses in the `metafor` package. Upon initial inspection, I noticed the I^2 value was different early in conducting these analyses. I want to ensure we choose the correct package.

Notes from the book: First, the amount of (residual) heterogeneity (i.e., τ2) is estimated with one of the various estimators that have been suggested in the literaturen (e.g. Hunter-Schmidt, Hedges estimator, DerSimonian-Laird estimator, Sidik-Jonkman estimator, maximum-likelihood or restricted maximum-likelihood estimator, empirical Bayes estimator)

#### Step 1: Load the package

```r
library(metafor)
```

#### Step 2: Run the random-effects model
`metafor` asks for the sampling variance via the `vi` argument. However, we have the standard error (square root of the variances), which can be supplied via the `sei` argument. When specifying the data in this way, we **must** set `measure = "GEN"` (Which is the default).

REML estimator is the default t^2 estimator. The various (residual) heterogeneity measures that can be specified via the `method` argument are the:

* "HS" = Hunter-Schmidt estimator
* "HE" = Hedges estimator
* "DL" = DerSimonian-Laird estimator
* "SJ" = Sidik-Jonkman estimator
* "ML" = Maximum-likelihood estimator
* "REML" = Restricted maximum-likelihood estimator
* "EB" = Empirical Bayes estimator

*Knapp and Hartung adjustment*

By default, the test statistics of the individual coefficients in the model (and the corresponding confidence intervals) are based on the normal distribution, while the omnibus test is based on a χ2 distribution with m degrees of freedom (m being the number of coefficients tested). The Knapp and Hartung (2003) method (`knha = TRUE`) is an adjustment to the standard errors of the estimated coefficients, which helps to account for the uncertainty in the estimate of τ2 and leads to different reference distributions. Individual coefficients and confidence intervals are then based on the t-distribution with k − p degrees of freedom, while the omnibus test statistic then uses an F-distribution with m and k − p degrees of freedom (p being the total number of model coefficients including the intercept if it is present). The Knapp and Hartung adjustment is only meant to be used in the context of random- or mixed-effects model.

**Model 1**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
metafor_whole_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "SJ",
                          knha = TRUE)
```

**Model 2**

* Use the DerSimonian-Laird estimator method to estimate tau
* Not using the Knapp-Hartung method


```r
metafor_whole_mod2 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "DL",
                          knha = FALSE)
```

**Model 3**

* Use the DerSimonian-Laird estimator method to estimate tau
* Use the Knapp-Hartung method


```r
metafor_whole_mod3 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "DL",
                          knha = TRUE)
```

**Model 4**

* Use the Restricted maximum-likelihood estimator method to estimate tau
* Use the Knapp-Hartung method


```r
metafor_whole_mod4 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "REML",
                          knha = TRUE)
```
**meta**

Model | k | SMD | LL | UL | z or t val | p | I^2
------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ 
1 | 47 |0.0345117 |-0.035378 | 0.1044014 | 0.9939728 | 0.3254361 | 0
2 | 47 | 0.051671 | -0.0147278 | 0.1180698 | 1.5252269 | 0.1272025 | 0
3 | 47 | 0.0431682 | -0.0329675 | 0.119304 | 1.1112806 | 0.2664476 | 0



**metafor**

Model | k | estimate | LL | UL | z or tval | p | I^2
------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ 
1 | 47 | 0.0345117 |  -0.035378 | 0.1044014 | 0.9939728 | 0.3254361 | 33.9942459
2 | 47 | 0.051671 |  -0.0147278 | 0.1180698 | 1.5252269 | 0.1272025 | 0
3 | 47 | 0.051671 |  -0.0141052 | 0.1174472 | 1.5812446 | 0.1206745 | 0
4 | 47 | 0.0431682 |  -0.0247359 | 0.1110724 | 1.2796443 | 0.2070869 | 14.1995619

##### Subgroup: Potency
**Low**

```r
metafor_lowpot_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "SJ",
                          knha = TRUE,
                          subset = (Potency=="Low"))
```

##### Publication bias (metafor)

```r
metafor::funnel(metafor_whole_mod1)
```

![](meta_analysis_files/figure-html/funnel metafor whole-1.png)<!-- -->

regression test (Eggers? but gave different results to one done through meta: Conducting meta-analyses in R with the metafor package)

```r
regtest(metafor_whole_mod1)
```

```
## 
## Regression Test for Funnel Plot Asymmetry
## 
## model:     mixed-effects meta-regression model
## predictor: standard error
## 
## test for funnel plot asymmetry: t = -1.3705, df = 45, p = 0.1773
```

Test trim and fill method

```r
trim <- metafor::trimfill(metafor_whole_mod1, estimator = "R0", side = NULL)
trim
```

```
## 
## Estimated number of missing studies on the right side: 0 (SE = 1.4142)
## Test of H0: no missing studies on the right side:      p-val = 0.5000
## 
## Random-Effects Model (k = 47; tau^2 estimator: SJ)
## 
## tau^2 (estimated amount of total heterogeneity): 0.0286 (SE = 0.0108)
## tau (square root of estimated tau^2 value):      0.1690
## I^2 (total heterogeneity / total variability):   33.99%
## H^2 (total variability / sampling variability):  1.52
## 
## Test for Heterogeneity:
## Q(df = 46) = 42.7985, p-val = 0.6071
## 
## Model Results:
## 
## estimate      se    tval    pval    ci.lb   ci.ub 
##   0.0345  0.0347  0.9940  0.3254  -0.0354  0.1044    
## 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
