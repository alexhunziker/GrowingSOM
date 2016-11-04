#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

#setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
setwd("~/gsom")

install.packages(repos=NULL,"GrowingSOM_0.1.tar.gz")
detach("package:GrowingSOM", unload=TRUE)
library(GrowingSOM)

# Load Data
#load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
load("/media/SWW/Alex/Validation/SBR_raw.RData")
#data("complex_sampledata")

# for loading wines...
require(kohonen)
data("wines")
training.wines <- sample(nrow(wines), 120)
test.wines <- wines[-training,]

testdf <- testdata$n_07_06[,3:ncol(testdata$n_07_06)]
testdf2 <- testdata$k_171819_05[1:6000,3:ncol(testdata$n_07_06)]
testdf2[7] <- 0
traindata <- traindata[1:5000,3:ncol(traindata)]

# Generate gsom model
gsom_model <- train.gsom(traindata, 
                         keepdata = FALSE, iterations = 50, 
                         spreadFactor = 0.9, alpha = 0.5)

supervised_gsom <- train_xy.gsom(training.wines, wine.classes[training.wines])

# Map data
mapped_realdata <- map(gsom_model, traindata)
mapped_testdata <- map(gsom_model, testdf, retaindata = TRUE)
mapped_testdata2 <- map(gsom_model, testdf2)

# Plot
plot(gsom_model)
plot(gsom_model, type="property")
plot(gsom_model, type="property", dim=8, main="Testplot")
plot(gsom_model, type="training")
plot(gsom_model, type="dist")

predicted <- predict(gsom_model, traindata, retaindata=TRUE)

plot(mapped_realdata)
plot(mapped_testdata)
plot(mapped_testdata2)


# Summary
summary(gsom_model)
