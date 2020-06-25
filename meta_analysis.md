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
Copy_of_Copy_of_Anticholinergic_Medication_R_17062020 <- read_csv("Z:/PROJECTS/2018_Anticholinergic_Med_SR/Analysis/Copy of Copy of Anticholinergic Medication_R_17062020.csv")
```

Assign data to `dat`

```r
dat <- Copy_of_Copy_of_Anticholinergic_Medication_R_17062020
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

#### Step 3: Run subgroups analysis
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

### Sub-group meta-analysis
#### Cognitive domain analyses
##### Step 1: Average across domains within studies
Create data set (`dat_study_domain`) which has the average effect size by cognitive domain within each study.


```r
dat_study_domain <- dat %>%
  group_by(study, cog_domain_lezak) %>% #group by study and cognitive domains (based on Lezak) within studies
  mutate(g = mean(g, na.rm = T),
         st_err = mean(st_err, na.rm = T)) %>% #obtain average g and st err within each study/domain
  filter(row_number()==1) %>% #make each study/domain only appear on one row
  select(study, Potency, cog_domain_lezak, g, st_err) #select relevant data
```

Then, create a data set for each cognitive domain (+ general/global cognitive tests)


```r
dat_att <- dat_study_domain %>%
  filter(cog_domain_lezak == "Attention")
dat_psych <-  dat_study_domain %>%
  filter(cog_domain_lezak == "Psychomotor Functioning")
dat_concept <- dat_study_domain %>%
  filter(cog_domain_lezak == "Concept Formation & Reasoning")
dat_mem <- dat_study_domain %>%
  filter(cog_domain_lezak == "Memory")
dat_exec <- dat_study_domain %>%
  filter(cog_domain_lezak == "Executive Function")
dat_lang <- dat_study_domain %>%
  filter(cog_domain_lezak == "Language")
dat_perc <- dat_study_domain %>%
  filter(cog_domain_lezak == "Perception")
dat_general <- dat_study_domain %>%
  filter(cog_domain_lezak == "General")
```

##### Step 2: Run random effects model for each cognitive domain
**Model 1**

* Use the Sidik-Jonkman method to estimate tau
* Use the Knapp-Hartung method


```r
meta_att_mod1 <- metagen(g,
                    st_err,
                    data = dat_att,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_psych_mod1 <- metagen(g,
                    st_err,
                    data = dat_psych,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_concept_mod1 <- metagen(g,
                    st_err,
                    data = dat_concept,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_mem_mod1 <- metagen(g,
                    st_err,
                    data = dat_mem,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_exec_mod1 <- metagen(g,
                    st_err,
                    data = dat_exec,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_lang_mod1 <- metagen(g,
                    st_err,
                    data = dat_lang,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_perc_mod1 <- metagen(g,
                    st_err,
                    data = dat_perc,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD

meta_general_mod1 <- metagen(g,
                    st_err,
                    data = dat_general,
                    studlab = paste(study),
                    comb.fixed = F,
                    comb.random = T,
                    method.tau = "SJ", #use Sidik-Jonkman method
                    hakn = T, #using the Knapp-Hartung method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
```

Compare the models:

**Model 1**

Cognitive domain | k | SMD | LL | UL | p
---- | ---- | ---- | ---- | ---- | ---- 
Attention | meta_att_mod1[["k"]] | meta_att_mod1[["TE.random"]] | meta_att_mod1[["lower.random"]] | meta_att_mod1[["upper.random"]] | meta_att_mod1[["pval.random"]]
Psychomotor Functioning | meta_psych_mod1[["k"]] | meta_psych_mod1[["TE.random"]] | meta_psych_mod1[["lower.random"]] | meta_psych_mod1[["upper.random"]] | meta_psych_mod1[["pval.random"]]
Concept Formation & Reasoning | meta_concept_mod1[["k"]] | meta_concept_mod1[["TE.random"]] | meta_concept_mod1[["lower.random"]] | meta_concept_mod1[["upper.random"]] | meta_concept_mod1[["pval.random"]]
Memory | meta_mem_mod1[["k"]] | meta_mem_mod1[["TE.random"]] | meta_mem_mod1[["lower.random"]] | meta_mem_mod1[["upper.random"]] | meta_mem_mod1[["pval.random"]]
Executive Function | meta_exec_mod1[["k"]] | meta_exec_mod1[["TE.random"]] | meta_exec_mod1[["lower.random"]] | meta_exec_mod1[["upper.random"]] | meta_exec_mod1[["pval.random"]]
Language | meta_lang_mod1[["k"]] | meta_lang_mod1[["TE.random"]] | meta_lang_mod1[["lower.random"]] | meta_lang_mod1[["upper.random"]] | meta_lang_mod1[["pval.random"]]
Perception | meta_perc_mod1[["k"]] | meta_perc_mod1[["TE.random"]] | meta_perc_mod1[["lower.random"]] | meta_perc_mod1[["upper.random"]] | meta_perc_mod1[["pval.random"]]
General | meta_general_mod1[["k"]] | meta_general_mod1[["TE.random"]] | meta_general_mod1[["lower.random"]] | meta_general_mod1[["upper.random"]] | meta_general_mod1[["pval.random"]]
