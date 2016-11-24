#######################################
#GSOM - Growing Self Organizing Maps
#example.r
#11/10/16 - Alex Hunziker
#######################################

# This is a example file to test and demonstrate the 
# functionality of the gsom functions.

#setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
setwd("~/gsom")

remove.packages("GrowingSOM", lib="~/R/x86_64-pc-linux-gnu-library/3.2")
install.packages(repos=NULL,"GrowingSOM_0.1.tar.gz")
try(detach("package:GrowingSOM", unload=TRUE), silent = TRUE)
library(GrowingSOM)

# Load Data
#load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
load("/media/SWW/Alex/Validation/SBR_raw.RData")
#data("complex_sampledata")

# for loading wines...
require(kohonen)
data("wines")
training.wines <- sample(nrow(wines), 120)
test.wines <- wines[-training.wines,] 

data(iris)
iris[,5] <- as.numeric(iris[,5])
iris[,6] <- iris[,5]^2   #to test 2 dependent variables.

testdf <- testdata$n_07_06[,3:ncol(testdata$n_07_06)]
testdf2 <- testdata$k_171819_05[1:6000,3:ncol(testdata$n_07_06)]
testdf2[7] <- 0
traindata <- traindata[1:5000,3:ncol(traindata)]

# Generate gsom model
gsom_model <- train.gsom(traindata, 
                         keepdata = T, iterations = 50, 
                         spreadFactor = 0.9, alpha = 0.5)

supervised_gsom <- train_xy.gsom(wines[training.wines,], wine.classes[training.wines], spreadFactor = 0.95)
supervised_gsom2 <- train_xy.gsom(iris[,1:4], iris[,5:6], spreadFactor = 0.95)

gsom_model_hex <- train.gsom(traindata, 
                             keepdata = FALSE, iterations = 50, 
                             spreadFactor = 0.9, alpha = 0.5, nhood = "hex")

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

plot(gsom_model_hex)
plot(gsom_model_hex, type="property")
plot(gsom_model_hex, type="dist")

plot(supervised_gsom, type="property")
plot(supervised_gsom, type="predict")
predicted <- predict(supervised_gsom, test.wines, retaindata=TRUE)
compare <- cbind(wine.classes[-training.wines], (predicted$prediction$predict)*(supervised_gsom$norm_param$maxy[1:1]-supervised_gsom$norm_param$miny[1:1])+supervised_gsom$norm_param$miny[1:1])
print(compare)

plot(supervised_gsom2, type="property")
plot(supervised_gsom2, type="predict")
predicted2 <- predict(supervised_gsom2, iris[,5:6], retaindata=TRUE)
compare <- cbind(iris[,5:6], predicted2$prediction[,3:4])
print(compare)

plot(mapped_realdata)
plot(mapped_testdata)
plot(mapped_testdata2)


# Summary
summary(gsom_model)
summary(supervised_gsom)
summary(mapped_testdata)
summary(predicted)

print(gsom_model_hex)
print(mapped_testdata2)
print(supervised_gsom)
