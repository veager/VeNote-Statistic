#
# chi-square goodness-of-fit
# 
# =======================================================

x = c(315 ,  108,  101,   32)
p = c(9/16, 3/16, 3/16, 1/16)
chisq.test(x, p=p, rescale.p=F)

qchisq(df=3, p=0.95)


no = 0:10
x = c(15, 47, 70, 81, 52, 25, 16, 9, 4, 1, 0)
p = dpois(no, lambda=3, log=FALSE) 
chisq.test(x, p=p, rescale.p=T)

qchisq(df=length(x)-1-1, p=0.95)


no = 0:10
x = c(57, 203, 383, 525, 532, 408, 273, 139, 45, 27, 10)
p = dpois(no, lambda=3.870, log=FALSE) 
chisq.test(x, p=p, rescale.p=T)

sum(x) * p



#
# k-s goodness-of-fit
# 
# =======================================================

ks.test.norm.D <- function(y) {
  ecdf.y = ecdf(y)    # 样本经验分布函数，function 类型
  
  y1 = knots(ecdf.y)  # 样本点
  p = ecdf.y(y1)      # 样本点对应的概率
  n = length(p)       
  
  norm.p = pnorm(y, mean=0, sd=1)  # 标准正态分布
  
  d1 = max(abs(norm.p - p))
  d2 = max(abs(norm.p[2:n] - p[1:(n-1)]))
  D = max(d1, d2)
  return (D)
}


ks.test.D <- function(y1, y2) {
  ecdf.y1 = ecdf(y1)
  ecdf.y2 = ecdf(y2)
  
  x = c(y1, y2)
  x = unique(sort(x))
  
  n = length(x)
  
  p1 = ecdf.y1(x)
  p2 = ecdf.y2(x)
  
  d1 = max(abs(p2 - p1))
  d2 = max(abs(p2[2:n] - p1[1:n-1]))
  D = max(d1, d2)
  
  return(D)
}




y1 = rnorm(20, mean=0.1, sd=1.1)
y2 = rnorm(40, mean=-0.1, sd=1.2)

# 单样本检验
ks.test(y1, "pnorm", mean=0, sd=1, alternative='two.sided')
ks.test.norm.D(y1)


# 双样本检验
ks.test(y1, y2)
ks.test.D(y1, y2)

library(NSM3)

pKolSmirn(y1, y2, method='Asymptotic')
pKolSmirn(y1, y2, method='Exact')

qKolSmirnLSA(0.006489)
qKolSmirn(0.006489)

