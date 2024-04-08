library(ineq)


# Self-coding of Gini Coefficient
gini.coef = function(x, freq=NA, norm=F){
  
  # the increasing index
  x.sort.ix = order(x)
  x.sort = sort(x)
  n = length(x)
  
  if (is.numeric(freq)) {
    freq.sort = freq[x.sort.ix] # the increasing index
    u = c(0, freq.sort)         # x-axis
    u = cumsum(u) / sum(u)
    v = c(0, x.sort*freq.sort)  # y-axis
    v = cumsum(v) / sum(v)
  }
  else {
    u = 0:n
    u = u/n                # x-axis
    v = c(0, x.sort)
    v = cumsum(v) / sum(v) # y-axis
  }
  
  gini = 0.
  for (i in 2:(n+1)){ 
    gini = gini + v[i-1] + v[i] 
  }
  gini = 1 - (1/n)* gini
  
  # standardized Gini coefficient
  if (norm == T){ 
    gini = n / (n-1) * gini 
  }
  
  res = list(gini=gini, v=v, u=u)
  return(res)
}


# Example of computing Gini coefficient
x <- c(20, 14, 59,  9, 36, 23,  3)

gini.coef(x, norm=F)$gini
gini.coef(x, norm=T)$gini
Gini(x, corr=F)
Gini(x, corr=T)



# Example of computing Lorenz curve
x <- c(20, 14, 59,  9, 36, 23,  3)
f <- c(10, 20, 35, 45, 10,  5, 25)

lc.res = Lc(x, n=f, plot=TRUE)
res = gini.coef(x, freq=f)

sum(abs(res$u - lc.res$p))
sum(abs(res$v - lc.res$L))


# Example 2
x <- c(1, 2, 3, 4, 5)  # income
f <- c(1, 2, 3, 2, 1)  # population

f <- c(1, 1, 1, 1, 1)  # equal income
f <- c(0, 0, 0, 0, 1)  # equal income

gini.coef(x, freq=f, norm=F)$gini
gini.coef(x, freq=f, norm=TRUE)$gini
Lc(x, n=f, plot=TRUE)



# Example of computing Lorenz curve
Lasym(x, n=f, interval=FALSE)

x = c(5, 10, 15, 15, 15, 20, 25)
Lasym(x, interval=FALSE)
Lasym(x, interval=TRUE)
Lc(x, plot=TRUE)






