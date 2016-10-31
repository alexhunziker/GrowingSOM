#######################################
#GSOM - Growing Self Organizing Maps
#map.r
#11/10/16 - Alex Hunziker
#######################################

# Loads all functionality for the GSOM

library(plotrix)
library(fields)

source("map.r")
source("summary.r")
source("plot.r")
source("train.r")
source("predict.r")

dyn.load("trainloop.so")
dyn.load("map.so")

print("GSOM loaded.")