---
title: "Meta analysis"
output: 
  html_document:
    keep_md: true
---

### Load packages

```r
library(meta)
library(metafor)
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
dat <- Anticholinergic_Medication_R_data
```

### Whole meta-analysis
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
  select(study, Potency, g, st_err) #select relevant data for analysis
```

#### Step 2: Run random effects model of all data
**Model 1**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
meta_all_mod1 <- metagen(g,
                    st_err,
                    data = dat_study,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
```

**Model 2**

* Use the DerSimonian-Laird method (default) to estimate tau
* Do not use Knapp-Hartung method


```r
meta_all_mod2 <- metagen(g,
                          st_err,
                          data = dat_study,
                          studlab = paste(study),
                          comb.fixed = F,
                          comb.random = T,
                          hakn = F, # not using the Knapp-Hartung method
                          prediction = T, #True = print prediction interval for future studies based on present evidence
                          sm = "SMD")
```

Compare the models:

Model | SMD | LL | UL | t (mod1)/z (mod 2) | p-value
------ | ------ | ------ | ------ | ------ | ------
1 | 0.0399374 |-0.0299174 | 0.1097921 | 1.1508125 | 0.2557565
2 | 0.0604699 | -0.0074851 | 0.1284249 | 1.7440775 | 0.0811456


#### Step 3: Run subgroup analyses
##### Potency (low/high)
Run sub-analyses by medication potency (low/high) for each model


```r
potency_subgroup_mod1 <- update.meta(meta_all_mod1, 
                             byvar=Potency, 
                             comb.random = TRUE, 
                             comb.fixed = FALSE)
potency_subgroup_mod2 <- update.meta(meta_all_mod2, 
                                   byvar=Potency, 
                                   comb.random = TRUE, 
                                   comb.fixed = FALSE)
```

Potency = **low**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 37 | 0.0268455 | -0.0491769 | 0.1028679 | 0.4785041 | 0.0270949 | 28.0324573
2 | 37 | 0.0326629 | -0.0424046 | 0.1077303 | 0.3937664 | 0 | 28.0324573

Potency = **high**

Model | k | SMD | LL | UL | p-value | tau^2 | Q
------ | ------ | ------ | ------ | ------ | ------| ------ | ------
1 | 10 | 0.0942258 | -0.1087058 | 0.2971574 | 0.3209322 | 0.0299876 | 11.1168868
2 | 10 | 0.1148865 | -0.0840702 | 0.3138432 | 0.2577307 | 0.0187628 | 11.1168868

##### Cognitive domain analyses

Create data set (`dat_study_domain`) which has the average effect size by cognitive domain within each study.


```r
dat_study_domain <- dat %>%
  group_by(study, cog_domain_lezak) %>% #group by study and cognitive domains (based on Lezak) within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  select(study, Potency, cog_domain_lezak, g, st_err) #select relevant data
```

Then run random effects model by subgroup (cognitive domain).

**Model 1**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
meta_domain_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_domain,
                    studlab = paste(study),
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
Attention | 38 | 0.0463823 | -0.0348533 | 0.1276178 | 0.2547381
Psychomotor Functioning | 17 | -0.1126038 | -0.3399906 | 0.114783 | 0.3094116
Concept Formation & Reasoning | 13 | 0.1135531 | -0.0506499 | 0.2777562 | 0.1577423
Perception | 3 | 0.2500101 | -0.900842 | 1.4008621 | 0.4486141
Memory | 16 | 0.0610615 | -0.0604692 | 0.1825923 | 0.3011356
Executive Function | 15 | -0.0255449 | -0.2899779 | 0.2388882 | 0.8388438
General | 15 | 0.0744954 | -0.1600963 | 0.3090871 | 0.5069277
Language | 6 | 0.1144948 | -0.0524442 | 0.2814338 | 0.138183

**Model 2**


```r
meta_domain_mod2 <- metagen(g,
                       st_err,
                       data = dat_study_domain,
                       studlab = paste(study),
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

Domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | 38 | 0.0602792 | -0.0185533 | 0.1391118 | 0.1339554
Psychomotor Functioning | 17 | -0.0975502 | -0.2987649 | 0.1036645 | 0.3420089
Concept Formation & Reasoning | 13 | 0.1533006 | 0.0078964 | 0.2987048 | 0.0387905
Perception | 3 | 0.2436719 | -0.2697523 | 0.7570961 | 0.3522666
Memory | 16 | 0.0824155 | -0.0411404 | 0.2059715 | 0.1910921
Executive Function | 15 | -0.0094805 | -0.2121213 | 0.1931603 | 0.9269392
General | 15 | 0.0737153 | -0.152454 | 0.2998847 | 0.5229461
Language | 6 | 0.1339417 | -0.0198242 | 0.2877075 | 0.0877706
