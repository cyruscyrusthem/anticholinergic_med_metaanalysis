---
title: "Anticholinergic Meta analysis"
output: 
  html_document:
    keep_md: true
---

### Load packages

```r
library(meta)
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
**Model 1**

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

**Model 2**

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

**Model 3**

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
1 | 0.0379123 |-0.0322658 | 0.1080904 | 1.087425 | 0.2825141
2 | 0.0607255 | -0.0067365 | 0.1281875 | 1.7642485 | 0.0776902
3 | 0.0481821 | -0.0288177 | 0.125182 | 1.2264337 | 0.2200355


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
1 | 37 | 0.031282 | -0.0447059 | 0.10727 | 0.4092767 | 0.0285213 | 28.3623326
2 | 37 | 0.0387776 | -0.0354397 | 0.1129949 | 0.3058087 | 0 | 28.3623326

Potency = **high**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 10 | 0.0650287 | -0.147939 | 0.2779964 | 0.5071516 | 0.0314259 | 12.4438199
2 | 10 | 0.0646321 | -0.1559363 | 0.2852004 | 0.5657522 | 0.0316974 | 12.4438199

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
Attention | 38 | 0.0331916 | -0.0461727 | 0.1125559 | 0.4022254
Psychomotor Functioning | 17 | -0.1135276 | -0.3307731 | 0.1037179 | 0.2843183
Concept Formation & Reasoning | 13 | 0.1062523 | -0.0569431 | 0.2694477 | 0.1814708
Perception | 3 | 0.2500101 | -0.900842 | 1.4008621 | 0.4486141
Memory | 16 | 0.0688595 | -0.0530282 | 0.1907472 | 0.247191
Executive Function | 15 | -0.0198762 | -0.2863981 | 0.2466458 | 0.8752055
General | 14 | 0.0673139 | -0.1844095 | 0.3190373 | 0.5733312
Language | 6 | 0.0778642 | -0.1204409 | 0.2761693 | 0.359135


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
Attention | 38 | 0.0487378 | -0.0286072 | 0.1260827 | 0.2168151
Psychomotor Functioning | 17 | -0.0965445 | -0.2877927 | 0.0947038 | 0.3224607
Concept Formation & Reasoning | 13 | 0.1438694 | -0.0021101 | 0.2898488 | 0.0534047
Perception | 3 | 0.2436719 | -0.2697523 | 0.7570961 | 0.3522666
Memory | 16 | 0.0884569 | -0.0319845 | 0.2088983 | 0.150016
Executive Function | 15 | -0.0036663 | -0.2055026 | 0.1981701 | 0.9715999
General | 14 | 0.0679172 | -0.1679626 | 0.3037969 | 0.5725261
Language | 6 | 0.1119348 | -0.0387519 | 0.2626215 | 0.1454141


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
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.621      -1.209--0.033 -1.891 0.06503
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
1 | 0.1264678 |0.0447317 | 0.208204 | 3.0971962 | 0.003009
2 | 0.1269536 | 0.050128 | 0.2037793 | 3.2388227 | 0.0012002

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
##              Intercept ConfidenceInterval      t     p
## Egger's test    -0.237       -1.021-0.547 -0.639 0.527
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
## Egger's test    -1.704       -2.88--0.528 -3.076 0.01521
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
1 | 0.2663537 |0.0069972 | 0.5257102 | 2.2026512 | 0.0448767
2 | 0.2697652 | 0.0428978 | 0.4966327 | 2.3305683 | 0.0197761



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
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.364        -1.148-0.42 -0.988 0.32997
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
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.209       -0.993-0.575 -0.487 0.63406
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
## Egger's test    -3.037      -4.605--1.469 -3.612 0.00256
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
1 | 0.1146356 |-0.1620447 | 0.391316 | 0.8592577 | 0.3994646
2 | 0.1146823 | -0.0894536 | 0.3188182 | 1.1010958 | 0.270855

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
## Egger's test    -1.465      -2.445--0.485 -2.791 0.01755
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
## Egger's test    -3.355      -6.099--0.611 -2.396 0.03376
```

**FURTHER ACTION REQUIRED**



**Perception**

```r
funnel.meta(meta_perc_mod1)
```

![](meta_analysis_files/figure-html/funnel perc-1.png)<!-- -->

<10 studies so Egger's not done

**FURTHER ACTION - assess visually for outliers**
