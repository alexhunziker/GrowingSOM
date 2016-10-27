# gsom_r
Implementation of Growing Self Organizing Maps in R

Mainly based on:
__Damminda Alahakoon, Saman K. Halgamuge (2000)__: Dynamic Self-Organizing Maps with Controlled Growth for Knowledge Discovery. IEEE TRANSACTIONS ON NEURAL NETWORKS, VOL. 11.

## Functionality

### train()
Generates and trains a new gsom model according to a numeric data.frame or matrix.

### map()
Maps testdata onto a trained gsom_model.

### summary()
Provides a summary of a gsom_model.

### plot()
Prints a gsom_model.

### cont_train()
Train an existing model with continuously new data.

### predict()
Missing feature

### train_dept()
Missing feature

### 

## Technical Stuff
In order to use the functionality, the file ``main.r`` has to be loaded.

Current code includes a c file that needs to be compiled in order to run the code. In order to compile the file execute

```bash
	R CMD SHLIB trainloop.r
```

Required R Libraries:
```r
	install.packages("fields")
	install.packages("plotrix")
```

## Remarks
Software is provided as is. No guarantee for functionality and or correctness of the code.

As the programm is under development, not all the functions have been implemented yet. 