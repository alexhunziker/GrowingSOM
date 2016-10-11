#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

source("main.r")

# Load Data
load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
testdf <- testdata$k_1617_06

# Generate gsom model
gsom_model <- gsom.train(traindata, 
                         retaintestdata = FALSE, iterations = 100, 
                         spreadFactor = 0.5, alpha = 0.5)

# Map data
mapped_testdata <- gsom.map(testdf, gsom_model)

# Plot

# Summary
gsom.summary(gsom_model)