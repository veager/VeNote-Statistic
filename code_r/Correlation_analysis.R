



## Example

# Example from "Introduction to Statistics and Data Analysis" Table 4.7 pp 83 
# for Pearson 
# non-tied samples
x1 = c(10.85, 10.44, 10.50, 10.89, 10.62)
x2 = c( 7.84,  7.96,  7.81,  7.47,  7.74)


# From "Nonparametric statistical methods (Third edition)" (2014) Example 8.5, pp. 430
# for Spearman rank 
# tied samples
x1 = c(7.1, 7.1, 7.2, 8.3, 9.4, 10.5, 11.4)
x2 = c(2.8, 2.9, 2.8, 2.6, 3.5,  4.6,  5.0)


# From "Nonparametric statistical methods (Third edition)" (2014) Example 8.1, pp. 397
# for Kendall rank
# non-tied samples
x1 = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
x2 = c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)




states <- state.x77
states = as.data.frame(states)
states = na.omit(states)
# names(states)
# [1] "Population" "Income"     "Illiteracy" "Life Exp"  
# [5] "Murder"     "HS Grad"    "Frost"      "Area"  
x1 = states$Income
x2 = states$'Life Exp'




##
## Pearson Correlation Coefficient
## 

cor.pearson = function(x, y) {
  x.mean = mean(x)
  y.mean = mean(y)
  r = sum((x - x.mean) * (y - y.mean)) / sqrt(sum((x - x.mean)^2) * sum((y - y.mean)^2))
  return (r)
}

cor.pearson.test.t = function(x, y) {
  n = length(x)
  r = cor.pearson(x, y)
  Sr = sqrt((1 - r^2) / (n-2))
  t = abs(r) / Sr
  return (t)
}

cor.pearson.test.z = function(x, y, rho=0) {
  n = length(x)
  r = cor.pearson(x, y)
  z1 = atanh(r)
  mu1 = atanh(rho)
  sigma1 = 1 / sqrt(n-3)
  z = (z1 - mu1) / sigma1
  return (z)
}


# Example from "Introduction to Statistics and Data Analysis" Table 4.7 pp 83 
# non-tied samples
x1 = c(10.85, 10.44, 10.50, 10.89, 10.62)
x2 = c( 7.84,  7.96,  7.81,  7.47,  7.74)


# Pearson Correlation Coefficient
cor(x1, x2, method='pearson')
cor.pearson(x1, x2)


# t-test : using cor.test()
cor.test(x1, x2, method="pearson", conf.level=0.95, alternative="two.sided")$statistic
cor.pearson.test.t(x1, x2)


# t-test : using corr.test()
states <- state.x77
states = as.data.frame(states)
states = na.omit(states)

library(psych)
corr.test(states, method="pearson", use="pairwise", alpha=0.05)


# Z-test : self codes 
# arctanh function, z-transform
# x = 0.2  # [0, 1]
# y1 = 1/2 * log((1+x)/(1-x)) # ln
# y2 = atanh(x)
x1 = states$Population
x2 = states$Income 

cor.pearson.test.z(x1, x2, rho=0.2)
cor.pearson.test.z(x1, x2, rho=0.7)
cor.pearson.test.z(x1, x2, rho=-0.7)




# ========================================================
##
## Spearman Correlation Coefficient
## 

cor.spearman.1 = function(x, y) {
  # compute by definition
  x.rank = rank(x)
  y.rank = rank(y)
  r = cor.pearson(x.rank, y.rank)
  return (r)
}

cor.spearman.2 = function(x, y) {
  # computationally efficient formula
  n = length(x)
  x.rank = rank(x)
  y.rank = rank(y)
  r = 1 - 6 * sum((x.rank - y.rank)^2) / (n * (n^2 - 1))
  return (r)
}

cor.spearman.test.z = function(x, y) {
  n = length(x)
  S = sum((rank(x) - rank(y))^2)
  
  r = cor.spearman.1(x, y)
  z = r / sqrt(1 / (n-1))
  
  res = list(z=z, S=S)
  return (res)
}



# Example from "Introduction to Statistics and Data Analysis" Table 4.7 pp 83 
# non-tied samples
x1 = c(10.85, 10.44, 10.50, 10.89, 10.62)
x2 = c( 7.84,  7.96,  7.81,  7.47,  7.74)

# Spearman Correlation Coefficient
cor(x1, x2, method='spearman')
cor.spearman.1(x1, x2)
cor.spearman.2(x1, x2)

# Hypothesis tesing: using self-codes
cor.spearman.test.z(x1, x2)
z = cor.spearman.test.z(x1, x2)$z
pnorm(-abs(z), lower.tail=T)*2

# Hypothesis tesing: using cor.test()
cor.test(x1, x2, method="spearman", conf.level=0.95)

# Hypothesis tesing: using pSpearman()
r = cor(x1, x2, method="spearman")
library(SuppDists)
pSpearman(q=-abs(r), r=length(x1), lower.tail=T)*2


# - - - - - - - - - - 
# Testing for tied samples
# From "Nonparametric statistical methods (Third edition)" (2014) Example 8.5, pp. 430
x1 = c(7.1, 7.1, 7.2, 8.3, 9.4, 10.5, 11.4)
x2 = c(2.8, 2.9, 2.8, 2.6, 3.5,  4.6,  5.0)

# Spearman Correlation Coefficient
cor(x1, x2, method='spearman')
cor.spearman.1(x1, x2)
cor.spearman.2(x1, x2)

# Hypothesis tesing: using self-codes
cor.spearman.test.z(x1, x2)
z = cor.spearman.test.z(x1, x2)$z
pnorm(-abs(z), lower.tail=T)*2


# Hypothesis tesing: using cor.test for tied samples with raise a warning
cor.test(x1, x2, method="spearman", conf.level=0.95)

# Hypothesis tesing: using pSpearman()
r = cor(x1, x2, method="spearman")
library(SuppDists)
pSpearman(q=-abs(r), r=length(x1), lower.tail=T)*2



# ========================================================
##
## Kendall Correlation Coefficient
## 

# tau_A, only used for non-tied samples
cor.kendall.1 = function(x, y) {
  n = length(x)
  S = 0
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      S = S + sign(x[i]-x[j]) * sign(y[i]-y[j])
    }
  }
  r = 2 / (n * (n-1)) * S
  return (r)
}


# tau_B, used for tied and non-tied samples. When sample is non-tied, tau_B = tau_A
cor.kendall.2 = function(x, y) {
  n = length(x)
  S = c(0, 0, 0, 0)
  # names(S) = c('nc', 'nd', 'n1', 'n2')
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      sign1 = sign(x[i]-x[j])
      sign2 = sign(y[i]-y[j])
      if (sign1 * sign2 > 0) {
        S[1] = S[1] + 1
      }
      else if (sign1 * sign2 < 0) {
        S[2] = S[2] + 1
      }
      else if (sign1 * sign2 == 0) {
        if (sign1 != 0) {
          S[4] = S[4] + 1
        } else if (sign2 != 0) {
          S[3] = S[3] + 1
        }
      }
    }
  }
  r = (S[1] - S[2]) / sqrt((n*(n-1) / 2 - S[3]) * (n*(n-1) / 2 - S[4]))
  return (r)
}


cor.kendall.test.z = function(x, y){
  n = length(x)
  r = cor.kendall.1(x, y)
  z = r / sqrt( 2 * (2*n+5) / (9 * n * (n-1)))
  
  S = 0   # the number of concordant pairs
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      s1 = sign(x[i]-x[j]) * sign(y[i]-y[j])
      if (s1 > 0) {
        S = S + s1
      }
    }
  }
  res = list("S"=S, "z"=z)
  return (res)
}


# Example: non-tied samples
# From "Nonparametric statistical methods (Third edition)" (2014) Example 8.1, pp. 397
x1 = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
x2 = c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)

# Kendall Correlation Coefficient
cor(x1, x2, method='kendall')
cor.kendall.2(x1, x2)
cor.kendall.1(x1, x2)

# Hypothesis tesing: using self-code
cor.kendall.test.z(x1, x2)
z = cor.kendall.test.z(x1, x2)$z
pnorm(-abs(z), lower.tail=T)*2  # two-sided test

# Hypothesis tesing: using cor.test()
cor.test(x1, x2, method='kendall', conf.level=0.95, alternative='two.sided')

# Hypothesis tesing: using pKendall()
r = cor(x1, x2, method='kendall')
library(SuppDists)
pKendall(-abs(r), N=length(x1), lower.tail=T)*2


# Example: tied samples
x1 = c(7.1, 7.1, 7.2, 8.3, 9.4, 10.5, 11.4)
x2 = c(2.8, 2.9, 2.8, 2.6, 3.5,  4.6,  5.0)

# Kendall Correlation Coefficient
cor(x1, x2, method='kendall')
cor.kendall.2(x1, x2)
cor.kendall.1(x1, x2)

# Hypothesis tesing: using cor.test() for tied samples with raise a warning
cor.kendall.test.z(x1, x2)
cor.test(x1, x2, method='kendall', conf.level=0.95, alternative='two.sided')

# Hypothesis tesing: using pKendall()
r = cor(x1, x2, method='kendall')
library(SuppDists)
pKendall(-abs(r), N=length(x1), lower.tail=T)*2











cor.1 = function(x, y, method="pearson"){
  if (method == "pearson"){
    r = cor.pearson(x, y)
  } 
  else if (method == "spearman") {
    r = cor.spearman.1(x, y)
  } 
  else if (method == "kendall") {
    r = cor.kendall.2(x, y)
  } 
  return (r)
}









