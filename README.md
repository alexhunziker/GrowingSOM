# gsom_r
Implementation of Growing Self Organizing Maps in R

Mainly based on:
__Damminda Alahakoon, Saman K. Halgamuge (2000)__: Dynamic Self-Organizing Maps with Controlled Growth for Knowledge Discovery. IEEE TRANSACTIONS ON NEURAL NETWORKS, VOL. 11.

## Functionality

### train.gsom
Generates and trains a new gsom model according to a numeric data.frame or matrix.

### train_xy.gsom
Supervised learning version of train.gsom

### map.gsom
Maps testdata onto a trained gsom_model.

### predict.gsom
Provides a summary of a gsom_model.

### plot.gsom
Prints a gsom_model.

### summary.gsom 
Summarizes gsom 

### print.gsom
Prints gsom

### 

## Installation
In order to use the functionality, the file ``main.r`` has to be loaded.

Build R package

```bash
	R CMD build gsom-pkg
```

Compile / Install up in R:
```r
    remove.packages("GrowingSOM", lib="~/R/x86_64-pc-linux-gnu-	library/3.2")
    install.packages(repos=NULL,"GrowingSOM_0.1.tar.gz")
    detach("package:GrowingSOM", unload=TRUE)
    library(GrowingSOM)
```

## Remarks
Software is provided as is. No guarantee for functionality and or correctness of the code.

As the programm is under development, not all the functions have been implemented yet. 