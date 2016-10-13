#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
source("main.r")

# Load Data
load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
testdf <- testdata$n_07_06[1:2000,3:ncol(testdata$n_07_06)]
traindata <- traindata[1:1000,3:ncol(traindata)]

# Generate gsom model
gsom_model <- gsom.train(traindata, 
                         keepdata = FALSE, iterations = 40, 
                         spreadFactor = 0.5, alpha = 0.5)

# Map data
mapped_testdata <- gsom.map(testdf, gsom_model)
mapped_realdata <- gsom.map(traindata, gsom_model)

# Plot
gsom.plot(gsom_model)
gsom.plot(gsom_model, type="property")
gsom.plot(gsom_model, type="property", dim=8, main="Testplot")
gsom.plot(gsom_model, type="training")

gsom.plot(mapped_testdata)
print_crude(mapped_testdata)



# Summary
gsom.summary(gsom_model)
