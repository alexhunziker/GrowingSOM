#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
#setwd("~/gsom")
source("main.r")

# Load Data
#load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
load("/media/SWW/Alex/Validation/SBR_raw.RData")
testdf <- testdata$n_07_06[1:2000,3:ncol(testdata$n_07_06)]
traindata <- traindata[1:10000,3:ncol(traindata)]

# Generate gsom model
gsom_model <- gsom.train(traindata, 
                         keepdata = FALSE, iterations = 25, 
                         spreadFactor = 0.9, alpha = 0.5)

# Map data
#mapped_testdata <- gsom.map(testdf, gsom_model)

# Plot
gsom.plot(gsom_model)
gsom.plot(gsom_model, type="property")
gsom.plot(gsom_model, type="property", dim=8, main="Testplot")
gsom.plot(gsom_model, type="training")
gsom.plot(gsom_model, type="dist")

#gsom.plot(mapped_realdata)
#gsom.plot(mapped_testdata)



# Summary
gsom.summary(gsom_model)
