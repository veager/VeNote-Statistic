


# how to make a QQ plot in R
x = rnorm(100, mean=0, sd=1)
y = rt(100, df=10)


qqplot.t = function(x1, x2){
  x = unique(c(x1, x2))
  
  p.x1 = ecdf(x1)(x)
  y = quantile(x, )
}




x1 = rnorm(10, mean=0, sd=1)
x2 = rt(1000, df=10)

x1.n = length(x1)
x2.n = length(x2)

if (x1.n >= x2.n){
  ref = x2
} else if (x1.n < x2.n){
  ref = x1
}

ref.p = ecdf(ref)(ref)


x1.q = quantile(x1, ref.p)
x2.q = quantile(x2, ref.p)

qqplot(x1, x2)
points(x1.q, x2.q, col='red')



#
#
# QQ-Norm
# 
# =============================


qqnorm.points <- function(x){
  n = length(x)
  x = sort(x)
  x.p = ecdf(x)(x)
  x.p = x.p - 1/(2*n)
  norm.q = qnorm(x.p, mean=0, sd=1)
  points(norm.q, x, col='red')
}



qqnorm.lines <- function(x){
  n = length(x)
  x = sort(x)
  x.p = ecdf(x)(x)
  x.p = x.p - 1/(2*n)
  
  mu = mean(x)
  sd = sd(x)
  
  y1 = min(x)
  y2 = max(x)
  x1 = (y1 - mu) / sd
  x2 = (y2 - mu) / sd
  
  lines(c(x1, x2), c(y1, y2), col='blue')
}


x = rnorm(200, mean=10, sd=2)  # Ñש±¾
qqnorm(x)
qqline(x, col='green')

qqnorm.points(x)
qqnorm.lines(x)



