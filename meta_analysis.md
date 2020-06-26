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
dat <- Anticholinergic_Medication_R_data
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
                    sm = "SMD") # says we want to calculate SMD
```

Code for forest plot:

```r
meta::forest(meta_all_mod1,
             sortvar = TE,
             xlim = c(-1.5, 1.5),
             leftcols = "studlab",
             rightcols = c("TE", "ci"),
             rightlabs = c("SMD", "95% CI"),
             text.random = "Overall effect") #found this here: https://rdrr.io/cran/meta/man/forest.html
```

![](meta_analysis_files/figure-html/forest mod1-1.png)<!-- -->

**Model 2**

* Use the DerSimonian-Laird method (default) to estimate tau
* Do not use Knapp-Hartung method


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
1 | 0.0399374 |-0.0299174 | 0.1097921 | 1.1508125 | 0.2557565
2 | 0.0604699 | -0.0074851 | 0.1284249 | 1.7440775 | 0.0811456
3 | 0.0484235 | -0.0297162 | 0.1265631 | 1.2145971 | 0.2245198


#### Step 3: Run subgroup analyses
*##### Potency (low/high)*
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
Does not seem to be able to run sub-analysis using REML method:
Error = Error in rma.uni(yi = TE[sel], sei = seTE[sel], method = method.tau, control = control) : Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.

Code for forest plot (not shown):

```r
meta::forest(potency_subgroup_mod1,
             sortvar = TE,
             leftcols = "studlab")
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

*##### Cognitive domain*

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
Attention | 38 | 0.0463823 | -0.0348533 | 0.1276178 | 0.2547381
Psychomotor Functioning | 17 | -0.1126038 | -0.3399906 | 0.114783 | 0.3094116
Concept Formation & Reasoning | 13 | 0.1135531 | -0.0506499 | 0.2777562 | 0.1577423
Perception | 3 | 0.2500101 | -0.900842 | 1.4008621 | 0.4486141
Memory | 16 | 0.0610615 | -0.0604692 | 0.1825923 | 0.3011356
Executive Function | 15 | -0.0255449 | -0.2899779 | 0.2388882 | 0.8388438
General | 15 | 0.0744954 | -0.1600963 | 0.3090871 | 0.5069277
Language | 6 | 0.1144948 | -0.0524442 | 0.2814338 | 0.138183


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
Attention | 38 | 0.0602792 | -0.0185533 | 0.1391118 | 0.1339554
Psychomotor Functioning | 17 | -0.0975502 | -0.2987649 | 0.1036645 | 0.3420089
Concept Formation & Reasoning | 13 | 0.1533006 | 0.0078964 | 0.2987048 | 0.0387905
Perception | 3 | 0.2436719 | -0.2697523 | 0.7570961 | 0.3522666
Memory | 16 | 0.0824155 | -0.0411404 | 0.2059715 | 0.1910921
Executive Function | 15 | -0.0094805 | -0.2121213 | 0.1931603 | 0.9269392
General | 15 | 0.0737153 | -0.152454 | 0.2998847 | 0.5229461
Language | 6 | 0.1339417 | -0.0198242 | 0.2877075 | 0.0877706


### Assess publication bias
Steps (based on Amit's paper):

1. Visually inspect funnel plots (effect size vs standard error) for bias
2. **If 10 or more studies**, use Egger's Test to test for small-study effect
    + If statistically significant asymmetry found (1-sided P<0.1), use Duval and Tweedie's Trima nd Fill method to quantify magnitude of bias
3. **If <10 studies**, locate outliers in funnel plot and recalculate effect sizes after their removal

To run Egger's, load Egger's function (found [here](https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/eggers.test.R))


Note: it does not matter which meta-analysis model the funnel plot or Egger's test are conducted based on, as only the individual study effect sizes matter, not the pooled effect size.
HOWEVER, models must be called separately for Trim and Fill.

#### Whole meta-analysis

```r
funnel.meta(meta_all_mod1)
```

![](meta_analysis_files/figure-html/funnel meta all-1.png)<!-- -->


```r
eggers.test(x = meta_all_mod1)
```

```
##              Intercept        ConfidenceInterval      t       p
## Egger's test    -0.546 -1.134-0.0419999999999999 -1.698 0.09642
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
1 | 0.1230909 |0.0428239 | 0.203358 | 3.0696693 | 0.0032574
2 | 0.1227195 | 0.0464688 | 0.1989702 | 3.1544073 | 0.0016082

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
## Egger's test    -0.152       -0.936-0.632 -0.426 0.67243
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
## Egger's test    -1.839      -2.819--0.859 -3.943 0.00428
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
1 | 0.3446468 |0.0791996 | 0.6100941 | 2.7673945 | 0.0143741
2 | 0.3489399 | 0.127657 | 0.5702228 | 3.0906574 | 0.0019971



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
## Egger's test    -0.314        -1.098-0.47 -0.876 0.38702
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
## Egger's test     -0.17       -0.954-0.614 -0.417 0.68335
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
## Egger's test    -1.146       -2.714-0.422 -1.402 0.18445
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
## Egger's test    -2.393      -3.961--0.825 -2.848 0.01222
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
1 | 0.1343854 |-0.1521207 | 0.4208915 | 0.9727479 | 0.3412534
2 | 0.1340769 | -0.0809688 | 0.3491225 | 1.2220003 | 0.2217075

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
## Egger's test    -1.457      -2.437--0.477 -2.821 0.01664
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
##              Intercept ConfidenceInterval      t       p
## Egger's test    -2.924      -5.472--0.376 -2.333 0.03636
```

**FURTHER ACTION REQUIRED**



**Perception**

```r
funnel.meta(meta_perc_mod1)
```

![](meta_analysis_files/figure-html/funnel perc-1.png)<!-- -->

<10 studies so Egger's not done

**FURTHER ACTION - assess visually for outliers**
