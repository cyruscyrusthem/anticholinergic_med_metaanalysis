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
             rightlabs = c("Hedges' g", "95% CI", "Weight"),
             text.random = "Overall effect",
             overall = TRUE,
             smlab = "",
             weight = FALSE,
             text.overall.random = TRUE,
             overall.hetstat = FALSE)
dev.off()
```


#### Step 3: Run subgroup analyses
##### *Potency (low/high)*
Run sub-analyses by medication potency (low/high) for each model.


```r
#Model 1
potency_subgroup_mod1 <- update.meta(meta_all_mod1, 
                             byvar=Potency, 
                             comb.random = TRUE, 
                             comb.fixed = FALSE)
```

Code for forest plot:

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

Then, run random effects model by subgroup (cognitive domain):

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

Code for forest plot:

```r
png(file='fig3_meta_forest_cog_1.png', width = 8, height = 31, units = "in", res = 300)
meta::forest(domain_subgroup_mod1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
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
                            comb.random = T,
                            method.tau = "PM", #use Paule-Mandel method
                            #hakn = T, #using the Knapp-Hartung method
                            #prediction = T, #True = print prediction interval for future studies based on present evidence
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
                              comb.random = T,
                              method.tau = "PM", #use Paule-Mandel method
                              #hakn = T, #using the Knapp-Hartung method
                              #prediction = T, #True = print prediction interval for future studies based on present evidence
                              sm = "SMD") # says we want to calculate SMD
domain_subgroup_group2 <- update.meta(meta_domain_group2, 
                                      byvar=cog_domain_lezak, 
                                      comb.random = TRUE, 
                                      comb.fixed = FALSE)
```


```r
png(file='fig3a_meta_forest_cog_1.png', width = 8, height = 18, units = "in", res = 300)
meta::forest(domain_subgroup_group1,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
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
png(file='fig3b_meta_forest_cog_1.png', width = 8, height = 18, units = "in", res = 300)
meta::forest(domain_subgroup_group2,
             sortvar = TE, #sort by effect size
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightlabs = c("Hedges' g", "95% CI"),
             text.random = "Overall effect",
             #weight = FALSE,
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

Forest plot:

```r
png(file='fig4_meta_forest_duration_1.png', width = 8, height = 13, units = "in", res = 300)
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

### Assess publication bias
Steps:

1. Visually inspect funnel plots (effect size vs standard error) for bias
2. **If 10 or more studies**, use Egger's Test to test for small-study effect
    + If statistically significant asymmetry found (1-sided P<0.1), use Duval and Tweedie's Trim and Fill method to quantify magnitude of bias
3. **If <10 studies**, locate outliers in funnel plot and recalculate effect sizes after their removal

To run Egger's, load Egger's function (found [here](https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/eggers.test.R))


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

#### Cognitive domains
Again, run subgroup analyses individually first (code not displayed).



**Attention**

```r
png(file='figS5_meta_funnel_att_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_att_mod1, xlab = "Hedges' g")
dev.off()
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
```


```r
png(file='figS7_meta_trimfunnel_psych_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_psych_mod1_trim, xlab = "Hedges' g")
dev.off()
```

**Concept formation**

```r
png(file='figS8_meta_funnel_conc_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_conc_mod1, xlab = "Hedges' g")
dev.off()
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

<10 studies so Egger's not done. Inspected visually for outliers, didn't seem to be any indication of substantial outliers.

**Intelligence**

```r
png(file='figS15_meta_funnel_gen_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_gen_mod1, xlab = "Hedges' g")
dev.off()
```

Egger's:

```r
eggers.test(x = meta_gen_mod1)
```

```
##              Intercept ConfidenceInterval      t       p
## Egger's test    -3.233      -5.781--0.685 -2.513 0.02725
```

```r
eggers_gen <- eggers.test(x = meta_gen_mod1)
```

p<0.1, therefore Duval and Tweedie's trim and fill method was used to quantify the magnitude of bias.


```r
meta_gen_mod1_trim <- trimfill.meta(meta_gen_mod1)

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

k<10 so visual inspection for outliers undertaken. Outlier (Aman) identified. Sensitivity analysis conducted with outlier removed:


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


#### Drug class

Again, run subgroup analyses individually first (code not displayed).



**Antiepileptic**

```r
png(file='figS17_meta_funnel_antiepi_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_antiepi_mod1, xlab = "Hedges' g")
dev.off()
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

**Respiratory**

```r
png(file='figS21_meta_funnel_resp_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_resp_mod1, xlab = "Hedges' g")
dev.off()
```

**Opioid analgesic**

```r
png(file='figS22_meta_funnel_opi_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_opi_mod1, xlab = "Hedges' g")
dev.off()
```

**Urological**

```r
png(file='figS22_meta_funnel_uro_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_uro_mod1, xlab = "Hedges' g")
dev.off()
```


#### Length of administration

Again, run subgroup analyses individually first (code not displayed).



**Long-term**

```r
png(file='figS23_meta_funnel_long_1.png', width = 8, height = 6, units = "in", res = 300)
funnel.meta(meta_long_mod1, xlab = "Hedges' g")
dev.off()
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

### **Results section**
Based on meta PM model.

*Overall cognition*

Overall, 46 studies which reported a cumulative 536 effect sizes, were available for analysis. The pooled effect size of the difference between cognition on- and off-medication across the 46 studies was negligible and non-significant (*g* = 0.0473952, 95% confidence interval (CI): -0.0197967 to 0.1145871, *p* = 0.1668178; see Figure 1), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot  did not reveal significant asymmetry (Egger’s intercept = -0.495, *p* = 0.13657; see Table S1).

*Potency* 

**Low potency.** The pooled effect size across low potency medication outcomes was negligible and statistically non-significant (*k* = 36, *g* = 0.0227174, 95% CI: -0.0512063 to 0.0966412, *p* = 0.5469645; see Figure 2), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.107, *p* = 0.76844; see Table S1).

**High potency.** The pooled effect size across high potency medication outcomes was negligible and statistically non-significant (*k* = 10, *g* = 0.1083688, 95% CI: -0.0803732 to 0.2971107, *p* = 0.2604449; see Figure 2), with low heterogeneity between studies (*tau^2^* = 0.0116812, *I^2^* = 28.0243284%). The funnel plot revealed significant asymmetry, indicating a potential small-study effect (Egger’s intercept = -1.734, *p* = 0.01388; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

*Cognitive domains* 

**Attention.** The pooled effect size across attention outcomes was negligible and non-significant (*k* = 37, *g* = 0.0398772, 95% CI: -0.0381468 to 0.1179012, *p* = 0.3164811; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.322, *p* = 0.39056; see Table S1).

**Psychomotor Functioning.** The pooled effect size across psychomotor functioning outcomes was negligible and non-significant (*k* = 17, *g* = -0.1034579, 95% CI: -0.3016644 to 0.0947485, *p* = 0.3062871; see Figure 3), with moderate heterogeneity between studies (*tau^2^* = 0.1022348, *I^2^* = 63.2378901%). The funnel plot revealed significant asymmetry (Egger’s intercept = -3.083, *p* = 0.00235; see Table S1). Trim and fill estimation led to an increase in effect size and a change in effect direction (Table S1).

**Concept Formation & Reasoning.** The pooled effect size across concept formation and reasoning outcomes was negligible and non-significant (*k* = 13, *g* = 0.1382341, 95% CI: -0.0030341 to 0.2795023, *p* = 0.0551276; see Figure 3), with negligible heterogeneity between studies (*tau^2^* = 0.0068826, *I^2^* = 15.9645404%). The funnel plot revealed significant asymmetry (Egger’s intercept = -1.49, *p* = 0.018; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

**Perception.** The pooled effect size across perception outcomes was small and non-significant (*k* = 3, *g* = 0.2456687, 95% CI: -0.2769132 to 0.7682505, *p* = 0.3568469; see Figure 3), with moderate heterogeneity between studies (*tau^2^* = 0.110928, *I^2^* = 50.1845504%). Visual inspection (*k* < 10) of the funnel plot indicated potential small-study effect (see Table S1). A sensitivity analysis was performed removing the outlier (Aman 2008), which resulted in a decrease in effect size (Table S1).

**Memory.** The pooled effect size across memory outcomes was negligible and non-significant (*k* = 16, *g* = 0.0393785, 95% CI: -0.0750444 to 0.1538013, *p* = 0.4999816; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.146, *p* = 0.70029; see Table S1).

**Executive Function.** The pooled effect size across executive function outcomes was negligible and non-significant (*k* = 15, *g* = -0.0137945, 95% CI: -0.2495899 to 0.2220009, *p* = 0.9087134; see Figure 3), with low heterogeneity between studies (*tau^2^* = 0.1199335, *I^2^* = 48.5043879%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -1.152, *p* = 0.19745; see Table S1).

**Language.** The pooled effect size across language outcomes was negligible and non-significant (*k* = 6, *g* = 0.1131588, 95% CI: -0.0317128 to 0.2580304, *p* = 0.1257886; see Figure 3), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Intelligence.** The pooled effect size across intelligence outcomes was negligible and non-significant (*k* = 14, *g* = 0.0761191, 95% CI: -0.1544082 to 0.3066464, *p* = 0.517521; see Figure 3), with high heterogeneity between studies (*tau^2^* = 0.1341982, *I^2^* = 76.2331004%). The funnel plot revealed significant asymmetry (Egger’s intercept = -3.233, *p* = 0.02725; see Table S1). Trim and fill estimation led to an increase in effect size (Table S1).

*Drug class* 

**Antiepileptic.** The pooled effect size across antiepileptic medications were negligible and non-significant (*k* = 14, *g* = -0.0322294, 95% CI: -0.1612011 to 0.0967422, *p* = 0.6242852; see Figure 4), with negligible heterogeneity between studies (*tau^2^* = 0.0048268, *I^2^* = 9.6687088%). The funnel plot revealed significant asymmetry (Egger’s intercept = -1.924, *p* = 0.01176; see Table S1). Trim and fill estimation led to an decrease in effect size (Table S1).

**Antipsychotic.** The pooled effect size for cognitive outcomes on antipsychotic medications was negligible and non-significant (*k* = 14, *g* = 0.0609384, 95% CI: -0.0531475 to 0.1750244, *p* = 0.2951434; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.533, *p* = 0.18167; see Table S1).

**Antidepressant.** The pooled effect size for cognitive outcomes on antidepressant medications was small and statistically significant (*k* = 7, *g* = 0.2391687, 95% CI: 0.0582267 to 0.4201107, *p* = 0.0095788; see Figure 4), with negligible heterogeneity between studies (*tau^2^* = 0.0045757, *I^2^* = 13.2206831%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Respiratory.** The pooled effect size for cognitive outcomes on respiratory medications was negligible and non-significant (*k* = 5, *g* = 0.0211953, 95% CI: -0.1756152 to 0.2180058, *p* = 0.8328279; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Opioid analgesic.** The pooled effect size for cognitive outcomes on opioid analgesic medications was negligible and non-significant (*k* = 3, *g* = -0.1777412, 95% CI: -0.4704684 to 0.1149861, *p* = 0.2340184; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Urological.** The pooled effect size for cognitive outcomes on urological medications was negligible and non-significant (*k* = 2, *g* = -0.1258015, 95% CI: -0.7415743 to 0.4899713, *p* = 0.6888489; see Figure 4), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

*Length of administration* 

**Current + long-term.** The pooled effect size across cognitive outcomes following long-term medication use (> 1-month) was negligible and non-significant (*k* = 29, *g* = 0.0671079, 95% CI: -0.0296005 to 0.1638163, *p* = 0.1738115; see Figure 5), with low heterogeneity between studies (*tau^2^* = 0.0134227, *I^2^* = 24.168451%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = -0.497, *p* = 0.29685; see Table S1).

**Current + acute.** The pooled effect size across cognitive outcomes following acute medication use (<= 1-month) was negligible and non-significant (*k* = 19, *g* = 0.0372, 95% CI: -0.092257 to 0.1666571, *p* = 0.573296; see Figure 5), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). The funnel plot did not reveal significant asymmetry (Egger’s intercept = 0.268, *p* = 0.56462; see Table S1).

**Historical.** The pooled effect size across cognitive outcomes following historical medication use was negligible and non-significant (*k* = 3, *g* = -0.1777412, 95% CI: -0.4704684 to 0.1149861, *p* = 0.2340184; see Figure 5), with null heterogeneity between studies (*tau^2^* = 0, *I^2^* = 0%). Visual inspection (*k* < 10) of the funnel plot did not indicate a substantial small-study effect (see Table S1).

**Whole analysis**

k | g | LL | UL | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ----
46 |0.0473952 |-0.0197967 | 0.1145871 | 0.1668178 | 0 | 0%


**Sub-analysis table: 95% CI in one column


Sub-analysis | k | g | 95% CI | p | tau^2^ | I^2^
---- | ---- | ---- | ---- | ---- | ---- | ---- 
Antiepileptic | 14 | -0.0322294 | -0.1612011 - 0.0967422 | 0.6242852 | 0.0048268 | 9.6687088
Antipsychotic | 14 | 0.0609384 | -0.0531475 - 0.1750244 | 0.2951434 | 0 | 0
Antidepressant | 7 | 0.2391687 | 0.0582267 - 0.4201107 | 0.0095788 | 0.0045757 | 13.2206831
Antiparkinson | 1 | 0.0390292 | -1.0830858 - 1.1611442 | 0.9456494 | NA | NA
Respiratory | 5 | 0.0211953 | -0.1756152 - 0.2180058 | 0.8328279 | 0 | 0
Opioid analgesic | 3 | -0.1777412 | -0.4704684 - 0.1149861 | 0.2340184 | 0 | 0
Urological | 2 | -0.1258015 | -0.7415743 - 0.4899713 | 0.6888489 | 0 | 0
Low | 36 | 0.0227174 | -0.0512063 - 0.0966412 | 0.5469645 | 0 | 0
High | 10 | 0.1083688 | -0.0803732 - 0.2971107 | 0.2604449 | 0.0116812 | 28.0243284
Current + long-term | 29 | 0.0671079 | -0.0296005 - 0.1638163 | 0.1738115 | 0.0134227 | 24.168451
Current + acute | 19 | 0.0372 | -0.092257 - 0.1666571 | 0.573296 | 0 | 0
Historical | 3 | -0.1777412 | -0.4704684 - 0.1149861 | 0.2340184 | 0 | 0
Attention | 37 | 0.0398772 | -0.0381468 - 0.1179012 | 0.3164811 | 0 | 0
Psychomotor Functioning | 17 | -0.1034579 | -0.3016644 - 0.0947485 | 0.3062871 | 0.1022348 | 63.2378901
Concept Formation & Reasoning | 13 | 0.1382341 | -0.0030341 - 0.2795023 | 0.0551276 | 0.0068826 | 15.9645404
Perception | 3 | 0.2456687 | -0.2769132 - 0.7682505 | 0.3568469 | 0.110928 | 50.1845504
Memory | 16 | 0.0393785 | -0.0750444 - 0.1538013 | 0.4999816 | 0 | 0
Executive Function | 15 | -0.0137945 | -0.2495899 - 0.2220009 | 0.9087134 | 0.1199335 | 48.5043879
Intelligence | 14 | 0.0761191 | -0.1544082 - 0.3066464 | 0.517521 | 0.1341982 | 76.2331004
Language | 6 | 0.1131588 | -0.0317128 - 0.2580304 | 0.1257886 | 0 | 0
