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
This step is necessary to ensure independence of each effect size entered into the meta-analysis. 

Data are averaged:

* Within studies
  + In all cases of studies which included outcomes for different medications, all medications were either low or high potency, so outcomes were averaged together.
  
Saved into new data frame: `whole`


Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
