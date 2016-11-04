setwd("~/gsom")
df = iris[,1:4]

dyn.load("trainloop.so")

rep <- 20;
lentr <- 100

df = as.matrix(df)

min <- apply(df, 2, function(x){min(x)})
max <- apply(df, 2, function(x){max(x)})
df <- t(apply(df, 1, function(x){(x-min)/(max-min)}))

weights <- matrix(0, nrow=100, ncol=4)
weights[1:4,] <- c(0.1, 0.2, 0.3, 0.5, 0.6, 0.3, 0.1, 0.9, 0.7, 0.8, 0.9, 0.8, 0.7, 0.6, 0.1, 0.5)

distnd <- rep(0, 100) #Error per node
freq <- rep(0, 100) #Frequ of nodes

gt = -ncol(df) * log(0.3)

npos <- matrix(0, nrow=100, ncol=2)
npos[1:4,] <- c(0, 1, 1, 0, 1, 0, 1, 0)

if(rep > lentr) error("Max nr of iterations exceeded.")

currtrain <- matrix(0, nrow=100, ncol=5)

outc = .C("som_train_loop",
   df = as.double(df),
   weights = as.double(weights),
   distnd = as.double(distnd),
   prep = as.integer(100), #repetitions
   plendf = as.integer(nrow(iris)),
   plennd = as.integer(4),
   plrinit = as.double(0.8), #initial learning rate
   freq = as.double(freq),
   alpha = as.double(0.99), #for lr depreciation
   pdim = as.integer(ncol(df)),
   gt = as.double(gt), # Growth Rate
   npos = as.double(npos),
   pradius = as.integer(2), #Missing feature
   plentn = as.integer(100), #Max num of nodes
   plentd = as.integer(nrow(iris)),
   currtrain = as.double(currtrain),
   plentr = as.integer(lentr) #Max num of iterations
)

training <- matrix(outc$currtrain, ncol=5)
training <- training[1:21,] #Workaround +1
weights <- matrix(outc$weights, ncol=4)
weights <- weights[1:outc$plennd,]
npos <- matrix(outc$npos, ncol=2)
npos <- npos[1:outc$plennd,]
freq <- outc$freq[1:outc$plennd]
distnd <- outc$distnd

plot(npos, type="n", col="white")
text(npos[,1], npos[,2], freq)
