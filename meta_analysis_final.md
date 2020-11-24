---
title: "Anticholinergic Meta analysis"
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

___

### Data import
Import data (.csv) and assign to `dat`

```r
library(readr)
Anticholinergic_Medication_R_data <- read_csv("Z:/PROJECTS/2018_Anticholinergic_Med_SR/Analysis/Final analysis/anticholinergic_med_metaanalysis/dat_antichol.csv")
```

Assign data to `dat`

```r
dat <- Anticholinergic_Medication_R_data %>%
  filter(!is.na(study)) #a blank row keeps being added in the .csv file. This removes it
```

____

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
  select(study, author, year, potency, g, st_err) #select relevant data for analysis
```

___

#### Step 2: Run random effects model of all data

* Using the Paule-Mandel method


```r
meta_all_mod1 <- metagen(g,
                    st_err,
                    data = dat_study,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    comb.random = T,
                    hakn = T,
                    method.tau = "PM", #use Paule-Mandel method
                    sm = "SMD") #calculate SMD
```

Code for forest plot:

```r
tiff(file='fig2_meta_forest_whole.tiff', width = 7, height = 10.75, units = "in", res = 600, compression = "lzw")
meta::forest(meta_all_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI", "Weight"),
             text.random = "Overall effect",
             overall = TRUE,
             smlab = "",
             weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE)
dev.off()
```

___

#### Step 3: Run subgroup analyses
##### *Potency (low/high)*
Run sub-analyses by medication potency (low/high) for each model.


```r
potency_subgroup_mod1 <- update.meta(meta_all_mod1, 
                             byvar=potency, 
                             comb.random = TRUE, 
                             comb.fixed = FALSE)
```

Code for forest plot:

```r
tiff(file='fig4_meta_forest_potency.tiff', width = 6, height = 11.2, units = "in", res = 600, compression = "lzw")
meta::forest(potency_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

___

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
  select(study, author, year, potency, cog_domain_lezak, g, st_err) #select relevant data
```

Then, run random effects model by subgroup (cognitive domain):

* Use the Paule-Mandel method to estimate tau


```r
meta_domain_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_domain,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    hakn = T,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    sm = "SMD") # says we want to calculate SMD
domain_subgroup_mod1 <- update.meta(meta_domain_mod1, 
                                byvar=cog_domain_lezak, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Code for forest plot:

```r
meta::forest(domain_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
```


```r
#individual cog plots
dat_study_domain_group1 <- dat_study_domain %>%
  filter(cog_domain_lezak == "Attention" | cog_domain_lezak =="Psychomotor Functioning" | cog_domain_lezak =="Concept Formation & Reasoning" | cog_domain_lezak =="Perception")
dat_study_domain_group2 <- dat_study_domain %>%
  filter(cog_domain_lezak == "Memory" | cog_domain_lezak =="Executive Function" | cog_domain_lezak =="Intelligence" | cog_domain_lezak =="Language")

meta_domain_group1 <- metagen(g,
                            st_err,
                            data = dat_study_domain_group1,
                            studlab = paste(author, year),
                            comb.fixed = F,
                            hakn = T,
                            comb.random = T,
                            method.tau = "PM", #use Paule-Mandel method
                            sm = "SMD") # says we want to calculate SMD
domain_subgroup_group1 <- update.meta(meta_domain_group1, 
                                      byvar=cog_domain_lezak, 
                                      comb.random = TRUE, 
                                      comb.fixed = FALSE)
meta_domain_group2 <- metagen(g,
                              st_err,
                              data = dat_study_domain_group2,
                              studlab = paste(author, year),
                              comb.fixed = F,
                              hakn = T,
                              comb.random = T,
                              method.tau = "PM", #use Paule-Mandel method
                              sm = "SMD") # says we want to calculate SMD
domain_subgroup_group2 <- update.meta(meta_domain_group2, 
                                      byvar=cog_domain_lezak, 
                                      comb.random = TRUE, 
                                      comb.fixed = FALSE)
```


```r
tiff(file='fig6a_meta_forest_cog.tiff', width = 6, height = 17.25, units = "in", res = 600, compression = "lzw")
meta::forest(domain_subgroup_group1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```


```r
tiff(file='fig6b_meta_forest_cog.tiff', width = 6, height = 17.25, units = "in", res = 600, compression = "lzw")
meta::forest(domain_subgroup_group2,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

Now, combine the two cognition plots into a single plot:

```r
library(multipanelfigure)
figure1 <- multi_panel_figure(width = 12, height = 17.25, unit = "in", columns = 2, rows = 1, row_spacing = 0) %>%
  fill_panel("fig6a_meta_forest_cog.tiff", column = 1) %>%
  fill_panel("fig6b_meta_forest_cog.tiff", column = 2)
tiff(file='fig6_meta_forest_cog.tiff', width = 12, height = 17.25, units = "in", res = 600, compression = "lzw")
figure1
dev.off()
```
____

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

Then, run random effects model by subgroup (drug class):

* Use the Paule-Mandel method to estimate tau


```r
meta_class_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_class,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    hakn = T,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    sm = "SMD") # says we want to calculate SMD
class_subgroup_mod1 <- update.meta(meta_class_mod1, 
                                byvar=drug_class, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Forest plot:

```r
tiff(file='fig3_meta_forest_class.tiff', width = 6, height = 14.5, units = "in", res = 600, compression = "lzw")
meta::forest(class_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             overall = F,
             prediction = F)
dev.off()
```

_____

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

Then, run random effects model by subgroup (duration):

* Use the Paule-Mandel method to estimate tau


```r
meta_duration_mod1 <- metagen(g,
                    st_err,
                    data = dat_study_duration,
                    studlab = paste(author, year),
                    comb.fixed = F,
                    hakn = T,
                    comb.random = T,
                    method.tau = "PM", #use Paule-Mandel method
                    prediction = T, #True = print prediction interval for future studies based on present evidence
                    sm = "SMD") # says we want to calculate SMD
duration_subgroup_mod1 <- update.meta(meta_duration_mod1, 
                                byvar=duration, 
                                comb.random = TRUE, 
                                comb.fixed = FALSE)
```

Forest plot:

```r
tiff(file='fig5_meta_forest_duration.tiff', width = 6, height = 13, units = "in", res = 600, compression = "lzw")
meta::forest(duration_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             text.overall.random = TRUE,
             overall.hetstat = FALSE,
             print.byvar = F, #this removed "potency =" and just left subgroup label
             smlab = "",
             hetstat = FALSE,
             test.overall = F,
             prediction = F,
             overall = F)
dev.off()
```

____

### Assess publication bias
Steps:

1. Visually inspect funnel plots (effect size vs standard error) for bias
2. **If 10 or more studies**, use Egger's Test to test for small-study effect
    + If statistically significant asymmetry found (1-sided P<0.1), use Duval and Tweedie's Trim and Fill method to quantify magnitude of bias
3. **If <10 studies**, locate outliers in funnel plot and recalculate effect sizes after their removal

To run Egger's, load Egger's function (found [here](https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/eggers.test.R))


____

#### Whole meta-analysis

```r
png(file='figS1_meta_funnel_whole_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_all_mod1, xlab = "Hedges' g")
dev.off()
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

____

### Results

**Whole analysis**

k | g | LL | UL | p | tau^2^ | I^2^ | Q
---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
46 |0.0473952 |-0.0195985 | 0.1143889 | 0.1610876 | 0 | 0% | 42.3623376

**Sub-analysis table (combining with random-effects)**


Sub-analysis | k | g | 95% CI | p | tau^2^ | I^2^ | Q
---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- 
Antiepileptic | 14 | -0.0322294 | -0.1743892 - 0.1099303 | 0.6324484 | 0.0048268 | 9.6687088 | 14.3914693 
Antipsychotic | 14 | 0.0609384 | -0.0342013 - 0.1560782 | 0.1897344 | 0 | 0 | 7.441203 
Antidepressant | 7 | 0.2391687 | 0.0132144 - 0.465123 | 0.0412112 | 0.0045757 | 13.2206831 | 6.9140899 
Antiparkinson | 1 | 0.0390292 | -1.0830858 - 1.1611442 | 0.9456494 | NA | NA | 0 
Respiratory | 5 | 0.0211953 | -0.1483864 - 0.1907769 | 0.7460674 | 0 | 0  | 1.4799235 
Opioid analgesic | 3 | -0.1777412 | -0.7939699 - 0.4384876 | 0.3404144 | 0 | 0 | 1.8391246 
Urological | 2 | -0.1258015 | -1.8281093 - 1.5765064 | 0.5200227 | 0 | 0 | 0.1818443 
Low | 36 | 0.0227174 | -0.0450249 - 0.0904598 | 0.5004789 | 0 | 0 | 27.3954939
High | 10 | 0.1083688 | -0.1095259 - 0.3262634 | 0.2896666 | 0.0116812 | 28.0243284 | 12.5042251
Current + long-term | 29 | 0.0662493 | -0.0338582 - 0.1663568 | 0.1860613 | 0.0125514 | 23.391868 | 36.5496446
Current + acute | 20 | 0.0502412 | -0.0385721 - 0.1390544 | 0.25101 | 0 | 0 | 8.0617008
Historical | 3 | -0.1777412 | -0.7939699 - 0.4384876 | 0.3404144 | 0 | 0 | 1.8391246
Attention | 37 | 0.0398772 | -0.0402848 - 0.1200391 | 0.3197601 | 0 | 0 | 35.4897181
Psychomotor Functioning | 17 | -0.1034579 | -0.3178389 - 0.110923 | 0.3215072 | 0.1022348 | 63.2378901 | 43.5230732
Concept Formation & Reasoning | 13 | 0.1382341 | -0.0188303 - 0.2952985 | 0.0792704 | 0.0068826 | 15.9645404 | 14.2796863
Perception | 3 | 0.2456687 | -0.9015413 - 1.3928786 | 0.4541163 | 0.110928 | 50.1845504 | 4.0148187
Memory | 16 | 0.0393785 | -0.0573014 - 0.1360583 | 0.3989859 | 0 | 0 | 9.0548962
Executive Function | 15 | -0.0137945 | -0.2718257 - 0.2442367 | 0.9103418 | 0.1199335 | 48.5043879 | 27.1867824
Intelligence | 14 | 0.0761191 | -0.1779935 - 0.3302317 | 0.5288061 | 0.1341982 | 76.2331004 | 54.6979211
Language | 6 | 0.1131588 | -0.0679091 - 0.2942267 | 0.1690754 | 0 | 0 | 4.5406659

**Sub-analysis table (combining with fixed-effects)**


Sub-analysis | k | g | 95% CI | p | tau^2^ | I^2^ | Q
---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- 
Antiepileptic | 14 | -0.023635 | -0.1450149 - 0.0977449 | 0.702726 | 0.0048268 | 9.6687088 | 14.3914693 
Antipsychotic | 14 | 0.0609384 | -0.0531475 - 0.1750244 | 0.2951434 | 0 | 0 | 7.441203 
Antidepressant | 7 | 0.2578004 | 0.0908379 - 0.4247628 | 0.0024756 | 0.0045757 | 13.2206831 | 6.9140899 
Antiparkinson | 1 | 0.0390292 | -1.0830858 - 1.1611442 | 0.9456494 | NA | NA | 0 
Respiratory | 5 | 0.0211953 | -0.1756152 - 0.2180058 | 0.8328279 | 0 | 0  | 1.4799235 
Opioid analgesic | 3 | -0.1777412 | -0.4704684 - 0.1149861 | 0.2340184 | 0 | 0 | 1.8391246 
Urological | 2 | -0.1258015 | -0.7415743 - 0.4899713 | 0.6888489 | 0 | 0 | 0.1818443 
Low | 36 | 0.0227174 | -0.0512063 - 0.0966412 | 0.5469645 | 0 | 0 | 27.3954939
High | 10 | 0.1646766 | 0.0035209 - 0.3258323 | 0.0452004 | 0.0116812 | 28.0243284 | 12.5042251
Current + long-term | 29 | 0.0771692 | -0.0010683 - 0.1554067 | 0.0532115 | 0.0125514 | 23.391868 | 36.5496446
Current + acute | 20 | 0.0502412 | -0.0774366 - 0.1779189 | 0.4405616 | 0 | 0 | 8.0617008
Historical | 3 | -0.1777412 | -0.4704684 - 0.1149861 | 0.2340184 | 0 | 0 | 1.8391246
Attention | 37 | 0.0398772 | -0.0381468 - 0.1179012 | 0.3164811 | 0 | 0 | 35.4897181
Psychomotor Functioning | 17 | 0.0103703 | -0.0979011 - 0.1186418 | 0.8510908 | 0.1022348 | 63.2378901 | 43.5230732
Concept Formation & Reasoning | 13 | 0.1590038 | 0.0321123 - 0.2858954 | 0.0140505 | 0.0068826 | 15.9645404 | 14.2796863
Perception | 3 | 0.1793376 | -0.1659284 - 0.5246036 | 0.308658 | 0.110928 | 50.1845504 | 4.0148187
Memory | 16 | 0.0393785 | -0.0750444 - 0.1538013 | 0.4999816 | 0 | 0 | 9.0548962
Executive Function | 15 | 0.057334 | -0.0725232 - 0.1871911 | 0.3868441 | 0.1199335 | 48.5043879 | 27.1867824
Intelligence | 14 | 0.2058602 | 0.1021236 - 0.3095969 | 1.0047122\times 10^{-4} | 0.1341982 | 76.2331004 | 54.6979211
Language | 6 | 0.1131588 | -0.0317128 - 0.2580304 | 0.1257886 | 0 | 0 | 4.5406659


**Sub-analysis: Test of subgroup differences**

Sub-analysis | Q (random) | p (random) | Q (fixed) | p (fixed)
--- | --- | --- |  --- | ---
Class | 9.9761328 | 0.1256607 | 10.1146831 | 0.1199048 
Potency | 0.7059872 | 0.4007792 | 2.4626187 | 0.1165843 
Duration | 2.6221744 | 0.2695269 | 2.7360294 | 0.2546119 
Domain | 5.5882051 | 0.5885671 | 10.8420528 | 0.1456585 
