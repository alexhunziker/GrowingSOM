#################################
# Check Training
# 02.11.16 - Alex Hunziker
#################################

# This file provides some sanity checks for the training functionality of the gsom model.

#setwd("Q:/Abteilungsprojekte/eng/SWWData/Alex/gsom")
setwd("~/gsom")
source("main.r")

load("~/gsom/checkdfs.RData")

# TEST 1:
# Must return 1 node.
gsom_model <- gsom.train(checkdfs[[1]], keepdata = FALSE, iterations = 50, spreadFactor = 0.99, alpha = 0.5)
if(nrow(gsom_model$nodes$position) != 1) error("Test 1 failed. Wrong number of nodes found.")
if(gsom_model$training$meandist[50] > 0.0001) warning("Distance is not negligible for very easy Test 1")

# TEST 2:
# Must return 2 nodes.
gsom_model <- gsom.train(checkdfs[[2]], keepdata = FALSE, iterations = 50, spreadFactor = 0.99, alpha = 0.5)
if(nrow(gsom_model$nodes$position) != 2) error("Test 1 failed. Wrong number of nodes found.")
if(gsom_model$training$meandist[50] > 0.0001) warning("Distance is not negligible for very easy Test 2")

# TEST 3:
# Must return 2 clusters, each with a few nodes
gsom_model <- gsom.train(checkdfs[[3]], keepdata = FALSE, iterations = 50, spreadFactor = 0.95, alpha = 0.5)
print("Test 3: There should be 2 visually seperable clusters")
gsom.plot(gsom_model, type="property", main="Test3")
if(gsom_model$training$meandist[50] > 0.005) warning(paste0("Distance for Test 3 is ", gsom_model$training$meandist[50], ", which is rather high for this kind of data"))
Sys.sleep(2)

# TEST 4:
# Must return 3-4 clusters, each with a few nodes, property 5 is random
gsom_model <- gsom.train(checkdfs[[4]], keepdata = FALSE, iterations = 50, spreadFactor = 0.90, alpha = 0.5)
print("Test 4: There should be 3-4 clusters for Properties 1-4. Property 5 should be somewhat arranged.")
gsom.plot(gsom_model, type="property", main="Test4")
if(gsom_model$training$meandist[50] > 0.01) warning(paste0("Distance for Test 4 is ", gsom_model$training$meandist[50], ", which is rather high for this kind of data"))
Sys.sleep(5)

# TEST 5:
# Must return 5 nodes.
gsom_model <- gsom.train(checkdfs[[5]], keepdata = FALSE, iterations = 100, spreadFactor = 0.99, alpha = 0.5)
gsom.plot(gsom_model, type="property", main="Test5")
if(nrow(gsom_model$nodes$position) != 5) error("Test 5 failed. Wrong number of nodes found.")
if(gsom_model$training$meandist[100] > 0.0001) warning("Distance is not negligible for very easy Test 5")
Sys.sleep(1)

# TEST 6a:
# Should return 5 nodes. Noise in data included.
gsom_model <- gsom.train(checkdfs[[6]], keepdata = FALSE, iterations = 50, spreadFactor = 0.05, alpha = 0.5)
gsom.plot(gsom_model, type="property", main="Test6a")
if(nrow(gsom_model$nodes$position) != 5) warning(paste("Test 6a returned unexpected amount of nodes (5 were expected): ", nrow(gsom_model$nodes$position)))
Sys.sleep(1)

# TEST 6b:
# Should return 5 clusters whith a few nodes each. Noise in data included.
gsom_model <- gsom.train(checkdfs[[6]], keepdata = FALSE, iterations = 100, spreadFactor = 0.92, alpha = 0.5)
print("Test 6b: 5 Seperate Clusters should be visible.")
gsom.plot(gsom_model, type="property", main="Test6b")
if(gsom_model$training$meandist[100] > 0.01) warning(paste0("Distance for Test 6b is ", gsom_model$training$meandist[50], ", which is rather high for this kind of data"))
Sys.sleep(5)

# TEST 7:
# Random data. Should return some kind of order for some of the dimensions
gsom_model <- gsom.train(checkdfs[[7]], keepdata = FALSE, iterations = 50, spreadFactor = 0.90, alpha = 0.5)
print("Test 7: Random data. Should return some kind of order for some of the dimensions")
gsom.plot(gsom_model, type="property", main="Test7")
if(gsom_model$training$meandist[50] > 0.01) warning(paste0("Distance for Test 7 is ", gsom_model$training$meandist[50], ", which is rather high for this kind of data"))

rm(gsom_model, checkdfs)

print("test run complete.")
