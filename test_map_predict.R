require(kohonen)
require(GrowingSOM)

setwd("~/gsom")
remove.packages("GrowingSOM", lib="~/R/x86_64-pc-linux-gnu-library/3.2")
install.packages(repos=NULL,"GrowingSOM_0.1.tar.gz")
try(detach("package:GrowingSOM", unload=TRUE), silent = TRUE)
library(GrowingSOM)

data(iris)

iris_p = as.numeric(iris[,5])
iris_d = iris[,1:4]

minx <- apply(iris_d, 2, function(x){min(x)})
maxx <- apply(iris_d, 2, function(x){max(x)})
iris_d = t(apply(iris_d, 1, function(x){(x-minx)/ifelse(maxx==minx,1,(maxx-minx))}))

iris_p = (iris_p - min(iris_p))/(max(iris_p)-min(iris_p))

testsample = sample(1:150, 50)

iris_p_tr = iris_p[-testsample]
iris_p_te = iris_p[testsample]

iris_d_tr = as.matrix(iris_d[-testsample,])
iris_d_te = as.matrix(iris_d[testsample,])

a=vector()
b=vector()
c=vector()

for (i in 1:100){
  grid = somgrid(xdim=10, ydim=10)
  Ksom = som(data = iris_d_tr, grid=grid, rlen=50)

  Fsom = train.gsom(iris_d_tr, keepdata = T, iterations = 50, gridsize = 10)

  Gsom = train.gsom(iris_d_tr, keepdata = T, iterations = 50, spreadFactor = 0.99)

  Kmap = kohonen::map(Ksom, iris_d_te)
  Fmap = GrowingSOM::map.gsom(Fsom, iris_d_te)
  Gmap = GrowingSOM::map.gsom(Gsom, iris_d_te)


  kmat = Ksom$codes[Kmap$unit.classif,]
  kdist = iris_d_te - kmat
  a[i] = sum((kdist)^2)

  fmat = Fsom$nodes$codes[Fmap$mapped$bmn,]
  fdist = iris_d_te - fmat
  b[i] = sum((fdist)^2)

  gmat = Gsom$nodes$codes[Gmap$mapped$bmn,]
  gdist = iris_d_te - gmat
  c[i] = sum((gdist)^2)
}
cat("Mapping of 4 Dimensional Data. Sum of all absolute distances: Traditional Kohonen =", mean(a), ", 
    GrowingSOM =", mean(b), mean(c), "(Fixed and Growing respectively)\n")

for(i in 1:100){
  Ksom = xyf(data = iris_d_tr, Y = iris_p_tr, grid=grid, rlen=50)
  Fsom = train_xy.gsom(iris_d_tr, iris_p_tr, iterations=50, gridsize = 10, initrad = 6)
  Gsom = train_xy.gsom(iris_d_tr, iris_p_tr, iterations=50, spreadFactor = 0.95)

  Kpredict = predict(Ksom, iris_d_te)
  Fpredict = predict.gsom(Fsom, iris_d_te)
  Gpredict = predict.gsom(Fsom, iris_d_te)

  kdist = Kpredict$prediction - iris_p_te
  a[i] = sum(kdist^2)

  fdist = Fpredict$prediction$predict - iris_p_te
  b[i] = sum(fdist^2)

  gdist = Gpredict$prediction$predict - iris_p_te
  b[i] = sum(gdist^2)
}
cat("Prediction of 1 Variable, 4 Dimensional Input data. Sum of all absolute distances: Traditional Kohonen =", mean(a), ", 
    GrowingSOM =", mean(b), mean(c), "(Fixed and Growing respectively)\n")
