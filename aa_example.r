#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

#setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
setwd("~/gsom")
source("main.r")

# Load Data
#load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
load("/media/SWW/Alex/Validation/SBR_raw.RData")
testdf <- testdata$n_07_06[,3:ncol(testdata$n_07_06)]
testdf2 <- testdata$k_171819_05[1:6000,3:ncol(testdata$n_07_06)]
testdf2[7] <- 0
traindata <- traindata[1:50000,3:ncol(traindata)]

# Generate gsom model
gsom_model <- train.gsom(traindata, 
                         keepdata = FALSE, iterations = 50, 
                         spreadFactor = 0.9, alpha = 0.5)

# Map data
mapped_realdata <- gsom.map(traindata, gsom_model)
mapped_testdata <- gsom.map(testdf, gsom_model)
mapped_testdata2 <- gsom.map(testdf2, gsom_model)

# Plot
plot(gsom_model)
plot(gsom_model, type="property")
plot(gsom_model, type="property", dim=8, main="Testplot")
plot(gsom_model, type="training")
plot(gsom_model, type="dist")

predicted <- gsom.predict(traindata, gsom_model, retaindata=TRUE)

plot(mapped_realdata)
plot(mapped_testdata)
plot(mapped_testdata2)


# Summary
summary(gsom_model)
