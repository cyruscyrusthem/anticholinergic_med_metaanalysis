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
  filter(!is.na(study)) #a blank row keeps being added in the .csv file. This removes it
```

### Meta-analysis
#### Step 1: Average across studies
This step is necessary to ensure independence of each effect size entered into the meta-analysis. In the current form, most studies report multiple outcomes (time-points, medication, cognitive tests), and these results are dependent on each other, as they come from they come from the same, or overlapping samples of participants. 

Data are averaged **within studies** so that each study represents a single effect size.
In all cases of studies which included outcomes for different medications, all medications were either low or high potency, so outcomes were able to be averaged together (within the study).

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

* Using the Paule-Mandel method


```r
meta_all_mod1 <- metagen(g,
                    st_err,
                    data = dat_study,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    #prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") #calculate SMD
```


Code for forest plot:

```r
png(file='fig1_meta_forest_whole_1.png', width = 8, height = 11, units = "in", res = 300)
meta::forest(meta_all_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI", "Weight"),
             text.random = "Overall effect",
             overall = TRUE,
             smlab = "",
             weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE)
```

```
## Warning in forest.meta(meta_all_mod1, sortvar = TE, xlim = c(-1.5, 1.5), :
## Deprecated argument 'weight' ignored as argument 'weight.study' is also
## provided.
```

```r
dev.off()
```

```
## png 
##   2
```

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


Compare the models:

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.0473952 |-0.0197967 | 0.1145871 | 1.3825013 | 0.1668178
2 | 0.0473952 | -0.0197967 | 0.1145871 | 1.3825013 | 0.1668178
3 | 0.0376968 | -0.0399904 | 0.115384 | 0.9510491 | 0.3415794


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

Code for forest plot (for Model 1):

```r
png(file='fig2_meta_forest_potency_2.png', width = 8, height = 11.5, units = "in", res = 300)
meta::forest(potency_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             #label.test.effect.subgroup.random = "x",
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

```
## png 
##   2
```

Compare the models:

Potency = **low**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 36 | 0.0227174 | -0.0512063 | 0.0966412 | 0.5469645 | 0 | 27.3954939
2 | 36 | 0.0227174 | -0.0512063 | 0.0966412 | 0.5469645 | 0 | 27.3954939

Potency = **high**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 10 | 0.1083688 | -0.0803732 | 0.2971107 | 0.2604449 | 0.0116812 | 12.5042251
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

Effect sizes based on cognitive composite scores (reported in 3 studies: Stevenson, Robles, Operto) were not included in the cognitive domain sub analysis. 


```r
dat_study_domain <- dat %>%
  group_by(study, cog_domain_lezak) %>% #group by study and cognitive domains (based on Lezak) within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  filter(cog_domain_lezak != "Not Subdomain") %>% #remove outcomes based on cognitive composite scores
  select(study, author, year, Potency, cog_domain_lezak, g, st_err) #select relevant data
```

Then, run random effects model by subgroup (cognitive domain).

**Model 1**

* Use the Paule-Mandel method to estimate tau


```r
meta_domain_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_domain,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    #hakn = T, #using the Knapp-Hartung method
                    #prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
domain_subgroup_mod1 <- update.meta(meta_domain_mod1, 
                                byvar=cog_domain_lezak, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Display the results of analysis:

Domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | 37 | 0.0398772 | -0.0381468 | 0.1179012 | 0.3164811
Psychomotor Functioning | 17 | -0.1034579 | -0.3016644 | 0.0947485 | 0.3062871
Concept Formation & Reasoning | 13 | 0.1382341 | -0.0030341 | 0.2795023 | 0.0551276
Perception | 3 | 0.2456687 | -0.2769132 | 0.7682505 | 0.3568469
Memory | 16 | 0.0393785 | -0.0750444 | 0.1538013 | 0.4999816
Executive Function | 15 | -0.0137945 | -0.2495899 | 0.2220009 | 0.9087134
General | 14 | 0.0761191 | -0.1544082 | 0.3066464 | 0.517521
Language | 6 | 0.1131588 | -0.0317128 | 0.2580304 | 0.1257886


Code for forest plot (not shown):

```r
png(file='fig3_meta_forest_cog_1.png', width = 8, height = 31, units = "in", res = 300)
meta::forest(domain_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             #label.test.effect.subgroup.random = "x",
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

```
## png 
##   2
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

Display results of analysis:

Domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | 37 | 0.0398772 | -0.0381468 | 0.1179012 | 0.3164811
Psychomotor Functioning | 17 | -0.0996516 | -0.2919519 | 0.0926487 | 0.3097868
Concept Formation & Reasoning | 13 | 0.1286619 | -0.0203716 | 0.2776953 | 0.0906356
Perception | 3 | 0.2436719 | -0.2697523 | 0.7570961 | 0.3522666
Memory | 16 | 0.0393785 | -0.0750444 | 0.1538013 | 0.4999816
Executive Function | 15 | -0.0036663 | -0.2055026 | 0.1981701 | 0.9715999
General | 14 | 0.0770246 | -0.1515051 | 0.3055543 | 0.5088726
Language | 6 | 0.1131588 | -0.0317128 | 0.2580304 | 0.1257886

##### *Drug class*

Create data set (`dat_study_class`) which has the average effect size by drug class within each study.

Drug class was defined by ATC group into:

* Urological
* Opioid analgesic
* Antiepileptic
* Antiparkinson
* Antipsychotic
* Antidepressant
* Respiratory


```r
dat_study_class <- dat %>%
  group_by(study, drug_class) %>% #group by study and drug class within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  select(study, author, year, drug_class, g, st_err) #select relevant data
```

Then, run random effects model by subgroup (drug class).

**Model 1**

* Use the Paule-Mandel method to estimate tau


```r
meta_class_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_class,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    #hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
class_subgroup_mod1 <- update.meta(meta_class_mod1, 
                                byvar=drug_class, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Class | k | g | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Antiepileptic | 14 | -0.0322294 | -0.1612011 | 0.0967422 | 0.6242852
Antipsychotic | 14 | 0.0609384 | -0.0531475 | 0.1750244 | 0.2951434
Antidepressant | 7 | 0.2391687 | 0.0582267 | 0.4201107 | 0.0095788
Antiparkinson | 1 | 0.0390292 | -1.0830858 | 1.1611442 | 0.9456494
Respiratory | 5 | 0.0211953 | -0.1756152 | 0.2180058 | 0.8328279
Opioid analgesic | 3 | -0.1777412 | -0.4704684 | 0.1149861 | 0.2340184
Urological | 2 | -0.1258015 | -0.7415743 | 0.4899713 | 0.6888489

Forest plot:

```r
png(file='fig4_meta_forest_class_1.png', width = 8, height = 15, units = "in", res = 300)
meta::forest(class_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             #label.test.effect.subgroup.random = "x",
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

```
## png 
##   2
```

##### *Length of administration*

Create data set (`dat_study_duration`) which has the average effect size by length of administration within each study.

Length of administration was categorised as 3 groups:

* Current administration <=1 month (Acute)
* Current administration >1 month (Long-term)
* Historical administration (Historical)


```r
dat_study_duration <- dat %>%
  group_by(study, duration) %>% #group by study and drug class within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  select(study, author, year, duration, g, st_err) #select relevant data
```

Then, run random effects model by subgroup (duration).

**Model 1**

* Use the Paule-Mandel method to estimate tau


```r
meta_duration_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_duration,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    #hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
duration_subgroup_mod1 <- update.meta(meta_duration_mod1, 
                                byvar=duration, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Duration | k | g | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Long-term | 29 | 0.0671079 | -0.0296005 | 0.1638163 | 0.1738115
Acute | 19 | 0.0372 | -0.092257 | 0.1666571 | 0.573296
Historical | 3 | -0.1777412 | -0.4704684 | 0.1149861 | 0.2340184

Forest plot:

```r
png(file='fig4_meta_forest_duration_1.png', width = 8, height = 13, units = "in", res = 300)
meta::forest(duration_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             #label.test.effect.subgroup.random = "x",
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             prediction = F,
             overall = F)
dev.off()
```

```
## png 
##   2
```

### Assess publication bias
Steps (based on Amit's paper):

1. Visually inspect funnel plots (effect size vs standard error) for bias
2. **If 10 or more studies**, use Egger's Test to test for small-study effect
    + If statistically significant asymmetry found (1-sided P<0.1), use Duval and Tweedie's Trim and Fill method to quantify magnitude of bias
3. **If <10 studies**, locate outliers in funnel plot and recalculate effect sizes after their removal

To run Egger's, load Egger's function (found [here](https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/eggers.test.R))


Note: it does not matter which meta-analysis model the funnel plot or Egger's test are conducted based on, as only the individual study effect sizes matter, not the pooled effect size.
HOWEVER, models must be run separately for Trim and Fill.

#### Whole meta-analysis

```r
png(file='figS1_meta_funnel_whole_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_all_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```


```r
eggers_whole <- eggers.test(x = meta_all_mod1)
eggers_whole
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.495       -1.083-0.093 -1.516 0.13657
```

p>0.1, so no further action required.

#### Potency
Does not seem possible to plot (using funnel plots) the individual parts of the sub-analyses separately. So, first I will run the meta-analyses separately for each subgroup (code not shown). This will provide the same results as previously found, it just requires extra data wrangling (filtering by subgroup), and more code. 



**Low potency**

```r
png(file='figS2_meta_funnel_low_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_low_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_low_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.107       -0.891-0.677 -0.297 0.76844
```

```r
eggers_meta_low_mod1 <- eggers.test(x = meta_low_mod1)
```

P>0.1 so nothing further required, no indication of substantial publication bias.



**High potency**

```r
png(file='figS3_meta_funnel_high_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_high_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_high_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.734       -2.91--0.558 -3.136 0.01388
```

```r
eggers_meta_high_mod1 <- eggers.test(x = meta_high_mod1)
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_high_mod1_trim <- trimfill.meta(meta_high_mod1)
```


```r
png(file='figS4_meta_trimfunnel_high_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_high_mod1_trim, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Results of trim and fill analysis (with **XX** studies imputed)

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.2675877 |0.0348041 | 0.5003713 | 2.2530032 | 0.0242589



#### Cognitive domains
Again, run subgroup analyses individually first (code not displayed).



**Attention**

```r
png(file='figS5_meta_funnel_att_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_att_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_att_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.322       -1.106-0.462 -0.869 0.39056
```

```r
eggers_att <- eggers.test(x = meta_att_mod1)
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Memory**

```r
png(file='figS12_meta_funnel_mem_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_mem_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_mem_mod1)
```

```
##              Intercept ConfidenceInterval     t       p
## Egger's test     0.146        -0.638-0.93 0.393 0.70029
```

```r
eggers_mem <- eggers.test(x = meta_mem_mod1)
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Executive function**

```r
png(file='figS13_meta_funnel_exec_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_exec_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_exec_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.152        -2.72-0.416 -1.358 0.19745
```

```r
eggers_exec <- eggers.test(x = meta_exec_mod1)
```

P>0.1 so nothing further required, no indication of substantial publication bias



**Psychomotor functioning**

```r
png(file='figS6_meta_funnel_psych_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_psych_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_psych_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -3.083      -4.651--1.515 -3.654 0.00235
```

```r
eggers_psych <- eggers.test(x = meta_psych_mod1)
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_psych_mod1_trim <- trimfill.meta(meta_psych_mod1)
meta_psych_mod2_trim <- trimfill.meta(meta_psych_mod2)
```


```r
png(file='figS7_meta_trimfunnel_psych_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_psych_mod1_trim, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Results of trim and fill analysis (with **XX** studies imputed)

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.1146221 |-0.1472288 | 0.376473 | 0.8579509 | 0.3909196
2 | 0.1147141 | -0.0907476 | 0.3201757 | 1.0942938 | 0.2738261



**Concept formation**

```r
png(file='figS8_meta_funnel_conc_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_conc_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_conc_mod1)
```

```
##              Intercept ConfidenceInterval      t     p
## Egger's test     -1.49        -2.47--0.51 -2.777 0.018
```

```r
eggers_conc <- eggers.test(x = meta_conc_mod1)
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.



```r
meta_conc_mod1_trim <- trimfill.meta(meta_conc_mod1)

png(file='figS9_meta_trimfunnel_conc_2.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_conc_mod1_trim, xlab = "Hedges' g")
dev.off()
```



**Language**

```r
png(file='figS14_meta_funnel_lang_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_lang_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

<10 studies so Egger's not done. Inspected visually for outliers, didn't seem to be any indication of substantial outliers.



**General**

```r
png(file='figS15_meta_funnel_gen_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_gen_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

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
## Egger's test    -3.233      -5.781--0.685 -2.513 0.02725
```

```r
eggers_gen <- eggers.test(x = meta_gen_mod1)
```

```
## Warning in metabias.meta(x, k.min = 3, method = "linreg"): 1 observation(s)
## dropped due to missing values
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_gen_mod1_trim <- trimfill.meta(meta_gen_mod1)
```

```
## Warning in trimfill.meta(meta_gen_mod1): 1 observation(s) dropped due to missing
## values
```

```r
png(file='figS16_meta_trimfunnel_gen_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_gen_mod1_trim, xlab = "Hedges' g")
dev.off()
```

**Perception**

```r
png(file='figS10_meta_funnel_perc_2.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_perc_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

k<10 so visual inspection for outliers undertaken. Outlier (Aman) identified. Sensitivity analysis conducted with outlier removed.


```r
dat_dom_perc_sens <- dat_study_domain %>%
  filter(cog_domain_lezak == "Perception") %>%
  filter(author != "Aman")
```


```r
meta_perc_sens <- metagen(g,
                    st_err,
                    data = dat_dom_perc_sens,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "PM",
                    #hakn = T,
                    prediction = T, 
                    sm = "SMD")
```


```r
png(file='figS11_meta_funnel_perc_sens_2.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_perc_sens, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```


#### Drug class

Again, run subgroup analyses individually first (code not displayed).



**Antiepileptic**

```r
png(file='figS17_meta_funnel_antiepi_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_antiepi_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_antiepi_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -1.924        -3.1--0.748 -2.967 0.01176
```

```r
eggers_antiepi <- eggers.test(x = meta_antiepi_mod1)
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_antiepi_mod1_trim <- trimfill.meta(meta_antiepi_mod1)

png(file='figS18_meta_trimfunnel_antiepi_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_antiepi_mod1_trim, xlab = "Hedges' g")
dev.off()
```

**Antipsychotic**

```r
png(file='figS19_meta_funnel_antipsych_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_antipsych_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_antipsych_mod1)
```

```
##              Intercept ConfidenceInterval     t       p
## Egger's test     0.533       -0.251-1.317 1.418 0.18167
```

```r
eggers_antipsych <- eggers.test(x = meta_antipsych_mod1)
```

p>0.1, nothing further required.

**Antidepressant**

```r
png(file='figS20_meta_funnel_antidep_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_antidep_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

**Respiratory**

```r
png(file='figS21_meta_funnel_resp_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_resp_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

**Opioid analgesic**

```r
png(file='figS22_meta_funnel_opi_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_opi_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

**Urological**

```r
png(file='figS22_meta_funnel_uro_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_uro_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```


#### Length of administration

Again, run subgroup analyses individually first (code not displayed).



**Long-term**

```r
png(file='figS23_meta_funnel_long_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_long_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_long_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -0.497       -1.477-0.483 -1.064 0.29685
```

```r
eggers_long <- eggers.test(x = meta_long_mod1)
```

p>0.1, nothing further required.

**Acute**

```r
png(file='figS24_meta_funnel_acute_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_acute_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

Egger's:

```r
eggers.test(x = meta_acute_mod1)
```

```
##              Intercept ConfidenceInterval     t       p
## Egger's test     0.268       -0.712-1.248 0.587 0.56462
```

```r
eggers_acute <- eggers.test(x = meta_acute_mod1)
```

p>0.1, nothing further required.

**Historical**

```r
png(file='figS25_meta_funnel_hist_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_hist_mod1, xlab = "Hedges' g")
dev.off()
```

```
## png 
##   2
```

### Forest plots for subgroups

#### Cognitive domains

**Attention**

```r
#png(file='fig3_meta_forest_cog_1.png', width = 8, height = 31, units = "in", res = 300)
meta::forest(meta_att_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             #rightcols = c("TE", "ci"),
             rightlabs = c("Hedges' g", "95% CI", "Weight"),
             text.random = "Overall effect",
             #weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             #label.test.effect.subgroup.random = "x",
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             #test.overall = F,
             #overall = F,
             prediction = F)
```

![](meta_analysis_files/figure-html/meta forest att domain mod1-1.png)<!-- -->

```r
#dev.off()
```

### Models in metafor
To contrast with analyses in the `meta` package, we will now also run analyses in the `metafor` package. Upon initial inspection, I noticed the I^2^ value was different early in conducting these analyses. I want to ensure we choose the correct package.

Notes from the book: First, the amount of (residual) heterogeneity (i.e., τ^2^) is estimated with one of the various estimators that have been suggested in the literaturen (e.g. Hunter-Schmidt, Hedges estimator, DerSimonian-Laird estimator, Sidik-Jonkman estimator, maximum-likelihood or restricted maximum-likelihood estimator, empirical Bayes estimator)

#### Step 1: Load the package

```r
library(metafor)
```

#### Step 2: Run the random-effects model
`metafor` asks for the sampling variance via the `vi` argument. However, we have the standard error (square root of the variances), which can be supplied via the `sei` argument. When specifying the data in this way, we **must** set `measure = "GEN"` (Which is the default).

REML estimator is the default t^2^ estimator. The various (residual) heterogeneity measures that can be specified via the `method` argument are the:

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
1 | 46 |0.0473952 |-0.0197967 | 0.1145871 | 1.3825013 | 0.1668178 | 0
2 | 46 | 0.0473952 | -0.0197967 | 0.1145871 | 1.3825013 | 0.1668178 | 0
3 | 46 | 0.0376968 | -0.0399904 | 0.115384 | 0.9510491 | 0.3415794 | 0



**metafor**

Model | k | estimate | LL | UL | z or tval | p | I^2
------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ 
1 | 46 | 0.0288211 |  -0.0424049 | 0.1000471 | 0.814994 | 0.4193683 | 34.2740508
2 | 46 | 0.0473952 |  -0.0197967 | 0.1145871 | 1.3825013 | 0.1668178 | 0
3 | 46 | 0.0473952 |  -0.0195985 | 0.1143889 | 1.4248917 | 0.1610876 | 0
4 | 46 | 0.0376968 |  -0.0315756 | 0.1069692 | 1.096038 | 0.2788933 | 14.991096

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

**High**

```r
metafor_highpot_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study,
                          method = "SJ",
                          knha = TRUE,
                          subset = (Potency=="High"))
```

##### Subgroup: Cognitive domains
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

Run Model 1 for each cognitive domain

```r
metafor_att_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Attention"))
metafor_mem_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Memory"))
metafor_exec_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Executive Function"))
metafor_psych_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Psychomotor Functioning"))
metafor_conc_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Concept Formation & Reasoning"))
metafor_lang_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Language"))
metafor_gen_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="General"))
```

```
## Warning in rma(yi = g, sei = st_err, data = dat_study_domain, method = "SJ", :
## Studies with NAs omitted from model fitting.
```

```r
metafor_perc_mod1 <- rma(yi = g, 
                          sei = st_err, 
                          data = dat_study_domain,
                          method = "SJ",
                          knha = TRUE,
                          subset = (cog_domain_lezak=="Perception"))
```

#### Step 3: Publication bias

```r
metafor::funnel(metafor_whole_mod1)
```

![](meta_analysis_files/figure-html/metafor funnel whole-1.png)<!-- -->

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
## test for funnel plot asymmetry: t = -1.3797, df = 44, p = 0.1747
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
## Random-Effects Model (k = 46; tau^2 estimator: SJ)
## 
## tau^2 (estimated amount of total heterogeneity): 0.0290 (SE = 0.0110)
## tau (square root of estimated tau^2 value):      0.1704
## I^2 (total heterogeneity / total variability):   34.27%
## H^2 (total variability / sampling variability):  1.52
## 
## Test for Heterogeneity:
## Q(df = 45) = 42.3623, p-val = 0.5843
## 
## Model Results:
## 
## estimate      se    tval    pval    ci.lb   ci.ub 
##   0.0288  0.0354  0.8150  0.4194  -0.0424  0.1000    
## 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


```r
options(digits = 2)
```


### **Results section**
Based on meta PM model

*Overall cognition*

Overall, 46 studies which reported a cumulative 536 effect sizes, were available for analysis. The pooled effect size of the difference between cognition on- and off-medication across the 46 studies was negligible and non-significant (*g* = 0.05, 95% confidence interval (CI): -0.02 to 0.11, *p* = 0.17; see Figure 1), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot  did not reveal significant asymmetry (Egger’s intercept = -0.5, *p* = 0.14; see Table S1).

*Potency* 

**Low potency.** The pooled effect size across low potency medication outcomes was negligible and statistically non-significant (*k* = 36, *g* = 0.02, 95% CI: -0.05 to 0.1, *p* = 0.55; see Figure 2), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.11, *p* = 0.77; see Table S1).

**High potency.** The pooled effect size across high potency medication outcomes was negligible and statistically non-significant (*k* = 10, *g* = 0.11, 95% CI: -0.08 to 0.3, *p* = 0.26; see Figure 2), with low heterogeneity between studies (*tau^2^* = 0.01, *I^2^* = 28.02%). The funnel plot revealed significant asymmetry, indicating a potential small-study effect (Egger’s intercept = -1.73, *p* = 0.01; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

*Cognitive domains* 

**Attention.** The pooled effect size across attention outcomes was negligible and non-significant (*k* = 37, *g* = 0.04, 95% CI: -0.04 to 0.12, *p* = 0.32; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.32, *p* = 0.39; see Table S1).

**Psychomotor Functioning.** The pooled effect size across psychomotor functioning outcomes was negligible and non-significant (*k* = 17, *g* = -0.1, 95% CI: -0.3 to 0.09, *p* = 0.31; see Figure 3), with moderate heterogeneity between studies (*tau^2^* = 0.1, *I^2^* = 63.24%). The funnel plot revealed significant asymmetry (Egger’s intercept = -3.08, *p* = 0; see Table S1). Trim and fill estimation led to an increase in effect size and a change in effect direction (Table S1).

**Concept Formation & Reasoning.** The pooled effect size across concept formation and reasoning outcomes was negligible and non-significant (*k* = 13, *g* = 0.14, 95% CI: 0 to 0.28, *p* = 0.06; see Figure 3), with negligible heterogeneity between studies (*tau^2^* = 0.01, *I^2^* = 15.96%). The funnel plot revealed significant asymmetry (Egger’s intercept = -1.49, *p* = 0.02; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

**Perception.** The pooled effect size across perception outcomes was small and non-significant (*k* = 3, *g* = 0.25, 95% CI: -0.28 to 0.77, *p* = 0.36; see Figure 3), with moderate heterogeneity between studies (*tau^2^* = 0.11, *I^2^* = 50.18%). Visual inspection (*k* < 10) of the funnel plot indicated potential small-study effect (see Table S1). A sensitivity analysis was performed removing the outlier (Aman 2008), which resulted in a decrease in effect size (Table S1).

**Memory.** The pooled effect size across memory outcomes was negligible and non-significant (*k* = 16, *g* = 0.04, 95% CI: -0.08 to 0.15, *p* = 0.5; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.15, *p* = 0.7; see Table S1).

**Executive Function.** The pooled effect size across executive function outcomes was negligible and non-significant (*k* = 15, *g* = -0.01, 95% CI: -0.25 to 0.22, *p* = 0.91; see Figure 3), with low heterogeneity between studies (*tau^2^* = 0.12, *I^2^* = 48.5%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -1.15, *p* = 0.2; see Table S1).

**Language.** The pooled effect size across language outcomes was negligible and non-significant (*k* = 6, *g* = 0.11, 95% CI: -0.03 to 0.26, *p* = 0.13; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**General.** The pooled effect size across intelligence outcomes was negligible and non-significant (*k* = 14, *g* = 0.08, 95% CI: -0.15 to 0.31, *p* = 0.52; see Figure 3), with high heterogeneity between studies (*tau^2^* = 0.13, *I^2^* = 76.23%). The funnel plot revealed significant asymmetry (Egger’s intercept = -3.23, *p* = 0.03; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

*Drug class* 

**Antiepileptic.** The pooled effect size across antiepileptic medications were negligible and non-significant (*k* = 14, *g* = -0.03, 95% CI: -0.16 to 0.1, *p* = 0.62; see Figure 4), with negligible heterogeneity between studies (*tau^2^* = 0, *I^2^* = 9.67%). The funnel plot revealed significant asymmetry (Egger’s intercept = -1.92, *p* = 0.01; see Table S1). Trim and fill estimation led to an decrease in effect size (Table S1).

**Antipsychotic.** The pooled effect size for cognitive outcomes on antipsychotic medications was negligible and non-significant (*k* = 14, *g* = 0.06, 95% CI: -0.05 to 0.18, *p* = 0.3; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.53, *p* = 0.18; see Table S1).

**Antidepressant.** The pooled effect size for cognitive outcomes on antidepressant medications was small and statistically significant (*k* = 7, *g* = 0.24, 95% CI: 0.06 to 0.42, *p* = 0.01; see Figure 4), with negligible heterogeneity between studies (*tau^2^* = 0, *I^2^* = 13.22%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Respiratory.** The pooled effect size for cognitive outcomes on respiratory medications was negligible and non-significant (*k* = 5, *g* = 0.02, 95% CI: -0.18 to 0.22, *p* = 0.83; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Opioid analgesic.** The pooled effect size for cognitive outcomes on opioid analgesic medications was negligible and non-significant (*k* = 3, *g* = -0.18, 95% CI: -0.47 to 0.11, *p* = 0.23; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Urological.** The pooled effect size for cognitive outcomes on urological medications was negligible and non-significant (*k* = 2, *g* = -0.13, 95% CI: -0.74 to 0.49, *p* = 0.69; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

*Length of administration* 

**Long-term.** The pooled effect size across cognitive outcomes following long-term medication use (> 1-month) was negligible and non-significant (*k* = 29, *g* = 0.07, 95% CI: -0.03 to 0.16, *p* = 0.17; see Figure 5), with low heterogeneity between studies (*tau^2^* = 0.01, *I^2^* = 24.17%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.5, *p* = 0.3; see Table S1).

**Acute.** The pooled effect size across cognitive outcomes following acute medication use (<= 1-month) was negligible and non-significant (*k* = 19, *g* = 0.04, 95% CI: -0.09 to 0.17, *p* = 0.57; see Figure 5), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.27, *p* = 0.56; see Table S1).

**Historical.** The pooled effect size across cognitive outcomes following historical medication use was negligible and non-significant (*k* = 3, *g* = -0.18, 95% CI: -0.47 to 0.11, *p* = 0.23; see Figure 5), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Whole analysis**

k | g | LL | UL | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ----
46 |0.05 |-0.02 | 0.11 | 0.17 | 0 | 0%


**Potency sub-analysis**

Potency | k | g | LL | UL | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
Low | 36 | 0.02 | -0.05 | 0.1 | 0.55 | 0 | 0%
High | 10 | 0.11 | -0.08 | 0.3 | 0.26 | 0.01 | 28.02%


**Cognitive domain sub-analysis**

Domain | k | g | LL | UL | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
Attention | 37 | 0.04 | -0.04 | 0.12 | 0.32 | 0 | 0%
Psychomotor Functioning | 17 | -0.1 | -0.3 | 0.09 | 0.31 | 0.1 | 63.24%
Concept Formation & Reasoning | 13 | 0.14 | 0 | 0.28 | 0.06 | 0.01 | 15.96%
Perception | 3 | 0.25 | -0.28 | 0.77 | 0.36 | 0.11 | 50.18%
Memory | 16 | 0.04 | -0.08 | 0.15 | 0.5 | 0 | 0%
Executive Function | 15 | -0.01 | -0.25 | 0.22 | 0.91 | 0.12 | 48.5%
General | 14 | 0.08 | -0.15 | 0.31 | 0.52 | 0.13 | 76.23%
Language | 6 | 0.11 | -0.03 | 0.26 | 0.13 | 0 | 0%


**Sub-analysis table

Sub-analysis | k | g | LL | UL | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
Antiepileptic | 14 | -0.03 | -0.16 | 0.1 | 0.62 | 0 | 9.67
Antipsychotic | 14 | 0.06 | -0.05 | 0.18 | 0.3 | 0 | 0
Antidepressant | 7 | 0.24 | 0.06 | 0.42 | 0.01 | 0 | 13.22
Antiparkinson | 1 | 0.04 | -1.08 | 1.16 | 0.95 | NA | NA
Respiratory | 5 | 0.02 | -0.18 | 0.22 | 0.83 | 0 | 0
Opioid analgesic | 3 | -0.18 | -0.47 | 0.11 | 0.23 | 0 | 0
Urological | 2 | -0.13 | -0.74 | 0.49 | 0.69 | 0 | 0
Low | 36 | 0.02 | -0.05 | 0.1 | 0.55 | 0 | 0
High | 10 | 0.11 | -0.08 | 0.3 | 0.26 | 0.01 | 28.02
Long-term | 29 | 0.07 | -0.03 | 0.16 | 0.17 | 0.01 | 24.17
Acute | 19 | 0.04 | -0.09 | 0.17 | 0.57 | 0 | 0
Historical | 3 | -0.18 | -0.47 | 0.11 | 0.23 | 0 | 0
Attention | 37 | 0.04 | -0.04 | 0.12 | 0.32 | 0 | 0
Psychomotor Functioning | 17 | -0.1 | -0.3 | 0.09 | 0.31 | 0.1 | 63.24
Concept Formation & Reasoning | 13 | 0.14 | 0 | 0.28 | 0.06 | 0.01 | 15.96
Perception | 3 | 0.25 | -0.28 | 0.77 | 0.36 | 0.11 | 50.18
Memory | 16 | 0.04 | -0.08 | 0.15 | 0.5 | 0 | 0
Executive Function | 15 | -0.01 | -0.25 | 0.22 | 0.91 | 0.12 | 48.5
General | 14 | 0.08 | -0.15 | 0.31 | 0.52 | 0.13 | 76.23
Language | 6 | 0.11 | -0.03 | 0.26 | 0.13 | 0 | 0


**Sub-analysis table: 95% CI in one column

Sub-analysis | k | g | 95% CI | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ---- 
Antiepileptic | 14 | -0.03 | -0.16 - 0.1 | 0.62 | 0 | 9.67
Antipsychotic | 14 | 0.06 | -0.05 - 0.18 | 0.3 | 0 | 0
Antidepressant | 7 | 0.24 | 0.06 - 0.42 | 0.01 | 0 | 13.22
Antiparkinson | 1 | 0.04 | -1.08 - 1.16 | 0.95 | NA | NA
Respiratory | 5 | 0.02 | -0.18 - 0.22 | 0.83 | 0 | 0
Opioid analgesic | 3 | -0.18 | -0.47 - 0.11 | 0.23 | 0 | 0
Urological | 2 | -0.13 | -0.74 - 0.49 | 0.69 | 0 | 0
Low | 36 | 0.02 | -0.05 - 0.1 | 0.55 | 0 | 0
High | 10 | 0.11 | -0.08 - 0.3 | 0.26 | 0.01 | 28.02
Long-term | 29 | 0.07 | -0.03 - 0.16 | 0.17 | 0.01 | 24.17
Acute | 19 | 0.04 | -0.09 - 0.17 | 0.57 | 0 | 0
Historical | 3 | -0.18 | -0.47 - 0.11 | 0.23 | 0 | 0
Attention | 37 | 0.04 | -0.04 - 0.12 | 0.32 | 0 | 0
Psychomotor Functioning | 17 | -0.1 | -0.3 - 0.09 | 0.31 | 0.1 | 63.24
Concept Formation & Reasoning | 13 | 0.14 | 0 - 0.28 | 0.06 | 0.01 | 15.96
Perception | 3 | 0.25 | -0.28 - 0.77 | 0.36 | 0.11 | 50.18
Memory | 16 | 0.04 | -0.08 - 0.15 | 0.5 | 0 | 0
Executive Function | 15 | -0.01 | -0.25 - 0.22 | 0.91 | 0.12 | 48.5
General | 14 | 0.08 | -0.15 - 0.31 | 0.52 | 0.13 | 76.23
Language | 6 | 0.11 | -0.03 - 0.26 | 0.13 | 0 | 0
