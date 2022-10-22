# 1. Engineering Experiments

- **Terminology**
  - The levels of the factor are sometimes called **treatments**
  
  - Each treatment has six **observations** or **replicates**
  
  - The runs are run in **random** order
  
- Every experiment involves a sequence of activities:

  - **Conjecture (推测)** – Original hypothesis that motivates the experiment.

  - **Experiment** – Test performance to investigate the conjecture.

  - **Analysis** – Statistical analysis of data from the experiment.

  - **Conclusion** – What has been learned about the original conjecture from the experiment. Often the experiment will lead to a revised conjecture, and a new experiment, and so forth.

## 1.1 Typical data

Typical data arrangement for a single factor experiment

- Suppose there are $a$ different levels of a single factor that we wish to compare.

$$
\begin{array}{c|cccc|cc}
\hline \text { Trt. } 
  & \text {Obs.} & \text {Obs.} & \text {Obs.} & \text {Obs.}
  & \text { Tot. } 
  & \text { Ave. } 
\\
\hline 1 & y_{11} & y_{12} & \cdots & y_{1 n} & y_1 \cdot & \bar{y}_1 \cdot \\
2 & y_{21} & y_{22} & \cdots & y_{2 n} & y_2 . & \bar{y}_2 . \\
\vdots & \vdots & \vdots & \vdots \ \vdots \ \vdots & \vdots & \vdots & \vdots \\
a & y_{a 1} & y_{a 2} & \cdots & y_{a n} & y_a . & \bar{y}_a . \\
\hline & & & & & y_{. .} & \bar{y} \cdot . \\
\hline
\end{array}
$$

where:
$$
\begin{array}{rcl}
y_{i \cdot} & = & \sum \limits_{j=1}^n y_{i j}, \quad i=1,2, \cdots, a \\
\bar{y}_{i \cdot} & = & \dfrac{y_{i \cdot}}{n}, \quad i=1,2, \cdots, a \\
y_{\cdot \cdot} & = & \sum \limits_{i=1}^a \sum \limits_{j=1}^n y_{i j} \\
\bar{y}_{\cdot \cdot} & = & \dfrac{y_{\cdot \cdot}}{N} \\ 
N & = & a \times n
\end{array}
$$

# 2. Fixed Effect Model

## 2.1 Formulation

### 2.1.1 Probabilistic model (i.e., linear statistical model)

We can formulate the observations in the table by the linear statistical model:
$$
Y_{i j} = \mu + \tau_i + \epsilon_{i j},  \qquad
\left\{\begin{array}{l}
i=1,2, \ldots, a \\
i=1,2, \ldots, n
\end{array}\right.
$$

where $\mu$ is the **overall mean**; $\tau_i$ is $i$th treatment effect; $\epsilon_{i j}$ is random error term and assume $\epsilon_{i j} \overset{i.i.d.}{\sim} \mathcal{N}(\mu, \sigma), \, \forall i,j$ (i.e., independently and identically distributed)

The model could be written as:
$$
Y_{i j} = \mu_i + \epsilon_{i j},  \qquad
\left\{\begin{array}{l}
i=1,2, \ldots, a \\
j=1,2, \ldots, n
\end{array}\right.
$$
where $\mu_i = \mu + \tau_i$

### 2.1.2 Assumptions

For the fixed effects models, the treatment effects are usually defined as deviations from the overall mean so that:
$$
\frac{\sum \limits_{i=1} ^a \mu_i}{a} 
  = \frac{\sum \limits_{i=1}^a (\mu+\tau_i )}{a} = \mu 
\quad \Rightarrow \quad
\sum_{i=1}^a \tau_i = 0
$$

The random error term $\epsilon_{i j}$ is assumed $\epsilon_{i j} \overset{i.i.d.}{\sim} \mathcal{N}(\mu, \sigma), \, \forall i,j$ (i.e., independently and identically distributed)

### 2.1.3 Fixed-effects model vs. random effects model

- **Fixed-effects**:

  - First, the $a$ treatments could have been **specifically** chosen by the experimenter. 

  - In this situation, we wish to test hypotheses about treatment means, and our conclusions will **apply only to the factor levels considered in the analysis**. 

  - The conclusion **cannot** be extended to similar treatments that were not explicitly considered. 

  - **We may also wish to estimate the model parameters** $(\mu, \tau_i, \sigma^2)$

- **Random-effects**

  - Alternatively, the a treatments could be a **random sample** from a larger population of treatments.
  
  - In this situation, we should like to be able to extend the conclusions (which are based on the sample treatments) to all treatments in the population, whether or not they were explicitly considered in the analysis.

  - Here, the $\tau_i$ is a **random variable**, and knowledge about the particular ones investigated is relatively useless. 

  - Instead, we test hypotheses about the variability of the and try to estimate this variability

### 2.1.4 Hypotheses

We wish to test the hypotheses:
$$
\begin{aligned}
& H_0 : \mu_1=\mu_2=\cdots=\mu_a \\
& H_1 : \mu_i \neq \mu_j \text{ for at least one pair } (i, j)
\end{aligned}
$$

An equivalent statement for the fixed-effect models:

$$
\begin{aligned}
& H_0 : \tau_1 = \tau_2 = \cdots = \tau_a = 0 \\
& H_1 : \tau_i \neq 0 \text { for at least one } i
\end{aligned}
$$

## 2.2 Variance Decomposition

### 2.2.1 Formulation

Total variability:
$$
SS_T = \sum_{i=1}^a \sum_{j=1}^n\left(y_{i j} - \bar{y}_{\cdot \cdot}\right)^2
$$

Decomposition of the total sum of squares $SS_T$:
$$
\sum_{i=1}^a \sum_{j=1}^n \left(y_{i j} -\bar{y}_{\cdot}\right)^2 = 
  n \sum_{i=1}^a \left(\bar{y}_{i \cdot}-\bar{y}_{\cdot}\right)^2 
  + \sum_{i=1}^a \sum_{j=1}^n\left(y_{i j}-\bar{y}_{i \cdot}\right)^2
$$

- Where $SS_{\text{Treatments}}$ is the sum of squares due to **treatment** (i.e, between treatments); 

  $SS_E$ is the sum of squares due to **error** (i.e., within treatments).

Then will have:

$$
SS_T = SS_{\text{Treatments}} + SS_E
$$

### 2.2.2 Number of degrees of freedom

Number of degrees of freedom for $SS_T$, $SS_{\text{Treatments}}$ and $SS_E$

The number of degrees of freedom of a sum of squares is equal to **the number of independent elements in that sum of squares**.

- $SS_T$ has $N-1$ degrees of freedom, where $N = a \times n$

- $SS_{\text{Treatments}}$ has $a-1$ degrees of freedom

- $SS_E$ has $a(n-1)$ degrees of freedom $[N-1-(a-1)] =N-a$

## 2.3 Sum Square

### 2.3.1 Analysis of $SS_E$

Note that:
$$
S S_E = \sum_{i=1}^a \sum_{j=1}^n (y_{i j}-\bar{y}_{i \cdot} )^2
= \sum_{i=1}^a \left[\sum_{j=1}^n \left(y_{i j}-\bar{y}_{i \cdot} \right)^2 \right]
$$

- **Sample variance** in the $i$th treatment
$$
S_i^2 = \frac{ \sum \limits_{j=1}^n \left(y_{i j}-\bar{y}_{i \cdot} \right)^2}{n-1}, 
\qquad i=1,2, \cdots, a
$$

- **Pooled estimate** of the common variance $\sigma^2$
$$
\frac{(n-1) S_1^2+(n-1) S_2^2 + \cdots+(n-1) S_a^2}{(n-1)+(n-1)+\cdots+(n-1)}
= \frac{ \sum \limits_{i=1}^a \left[ \sum \limits_{j=1}^n\left(y_{i j} - \bar{y}_{i \cdot}\right)^2 \right]}{\sum \limits_{i=1}^a(n-1)}
= \frac{S S_E}{N-a}
$$

- **Mean squares**:
$$
M S_E \triangleq \frac{S S_E}{N-a}
$$

- **Inference 1:** If there are **no differences** between the $a$ treatment means, the pooled estimate $MS_E$ can be used to estimate $\sigma^2$

- **Expected values** of $MS_E$
$$
\mathbb{E} \left[ MS_E \right] = \mathbb{E} \left[ \frac{S S_E}{N-a} \right]  = \sigma^2
$$

#### Proof of $MS_E = \sigma^2$

$$
\begin{aligned}
\mathbb{E} \left[ MS_E \right] &= \mathbb{E} \left[ \frac{S S_E}{N-a} \right] 
  = \frac{1}{N-a} \cdot \mathbb{E} \left[\sum_{i=1}^a \sum_{j=1}^n \left(y_{i j} -\bar{y}_{i \cdot}\right)^2 \right] \\
  & = \frac{1}{N-a} \cdot \mathbb{E} \left[
    \sum_{i=1}^a \sum_{j=1}^n \left(y_{i j}^2-2 y_{i j} \bar{y}_{i \cdot} + \bar{y}_{i \cdot}^2 \right) \right] \\
  & = \frac{1}{N-a} \cdot \mathbb{E} \left[ 
    \sum_{i=1}^a \sum_{j=1}^n y_{i j}^2 
    - 2 n \sum_{i=1}^a \bar{y}_{i \cdot}^2
    + n \sum_{i=1}^a \bar{y}_{i \cdot}^2 \right] \\
  & = \frac{1}{N-a} \cdot \mathbb{E} \left[
    \sum_{i=1}^a \sum_{j=1}^n y_{i j}^2 
    - \frac{1}{n} \sum_{i=1}^a y_{i \cdot}^2 \right] \\
  & =\frac{1}{N-a} \cdot \mathbb{E} \left[ 
    \sum_{i=1}^a \sum_{j=1}^n \left(\mu+\tau_i+\varepsilon_{i j}\right)^2 
    - \frac{1}{n} \sum_{i=1}^a \left(\sum_{j=1}^n \left(\mu + \tau_i + \varepsilon_{i j} \right) \right)^2 \right] \\
  & = \frac{1}{N-a} \cdot \mathbb{E} \left[n \sum_{i=1}^a \left(\mu + \tau_i \right)^2 + N \sigma^2 - n \sum_{i=1}^a \left(\mu + \tau_i \right)^2 - a \sigma^2 \right] \\
  & = \sigma^2
\end{aligned}
$$

### 2.3.2 Analysis of $SS_{\text{Treatments}}$

Note that
$$
SS_{\text {Treatments }} 
  = n \sum_{i=1}^a\left(\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot}\right)^2
$$

- Mean squares of treatments, $MS_{\text{Treatments}}$

$$
M S_{\text{Treatment}} \triangleq \frac{S S_{\text{Treatment}}}{N-a}
$$

- **Inference 2:** If there are **no differences** between the $a$ treatment means, we could use the variation of the treatment average from the grand average to estimate $\sigma^2$, i.e.,
$$
\frac{S S_{\text {Treatments }}}{a-1}
  = \frac{n \sum \limits_{i=1}^a \left(\bar{y}_{i \cdot} -\bar{y}_{\cdot \cdot}\right)^2}{a-1}
$$

- Expected values of $M S_{\text{Treatment}}$
$$
\mathbb{E} \left[ M S_{\text {Treatment }} \right]
 = \sigma^2 + \frac{n}{a-1}  \sum \limits_{i=1}^a \tau_i^2
$$

#### Proof of $\mathbb{E} \left[ M S_{\text {Treatment }} \right] = \sigma^2 + \frac{n}{a-1}  \sum_{i=1}^a \tau_i^2$

$$
\begin{split}
MS_{\text {Treatments }} &= \frac{n}{a-1} \sum_{i=1}^a\left(\bar{y}_{i \cdot} - \bar{y}_{\cdot \cdot}\right)^2 \\
  &= \frac{n}{a-1} \left[ \sum_{i=1}^a \bar{y}_{i \cdot}^2
                          - 2 \sum_{i=1}^a \bar{y}_{i \cdot} \, \bar{y}_{\cdot \cdot}
                          + \sum_{i=1}^a \bar{y}_{\cdot \cdot}^2 \right] \\
  &=\frac{n}{a-1} \left[ \sum_{i=1}^a \bar{y}_{i \cdot}^2-  a \bar{y}_{\cdot \cdot}^{2} \right]
\end{split}
$$

$$
\begin{split}
\mathbb{E} \left[ M S_{\text {Treatment }} \right] 
  &= \frac{n}{a-1} \ \mathbb{E} \left[ \sum_{i=1}^a \bar{y}_{i \cdot}^2-a \bar{y}_{\cdot \cdot}^2 \right]  \\
  &= \frac{n}{a-1} \left( \sum_{i=1}^a \mathbb{E} \left[ \bar{y}_{i \cdot}^2 \right] 
      - a \, \mathbb{E} \left[ \bar{y}_{\cdot \cdot}^2 \right] \right) \\
  &= \frac{n}{a-1} \left[ \sum_{i=1}^a \left(
          \mathbb{E} \left[\bar{y}_{i \cdot}\right]^2 + \mathbb{D} \left[ \bar{y}_{i \cdot} \right] 
      \right) 
      - a \left( \mathbb{E} \left[ \bar{y}_{\cdot \cdot}\right]^2 + \mathbb{D} \left[ \bar{y}_{\cdot \cdot} \right] \right)
      \right] \\
  & = \frac{n}{a-1} \left\{ 
    \sum_{i=1}^a \left[ \left(\mu+\tau_i\right)^2 + \frac{\sigma^2}{n} \right] 
    - a \left[ \left( \mu + \frac{ \sum_{i=1}^a \tau_i}{a}\right)^2 + \frac{\sigma^2}{a \cdot n} \right] 
      \right\} \\
  & = \sigma^2 + \frac{n}{a-1} \left[ \sum_{i=1}^a \left( \mu+\tau_i\right)^2-a \mu^2 \right] \\
  & = \sigma^2 + \frac{n}{a-1} \left( \sum_{i=1}^a \tau_i^2+2 \mu \sum_{i=1}^a \tau_i \right) \\
  & = \sigma^2 + \frac{n}{a-1} \sum \limits_{i=1}^a \tau_i^2
\end{split}
$$

- Deduction:
$$
\mathbb{E} \left[ M S_{\text {Treatment }} \right] \ge \mathbb{E} \left[ M S_{E} \right]
$$

### 2.4 Statistic Distribution

Two chi-square distributions and one F-distribution

- $SS_E/ \sigma^2$ follows chi-square with $N-a$ degrees of freedom

$$
\dfrac{SS_E}{\sigma^2} \sim \chi^2(N-a)
$$

- If the **null hypothesis** $H_0:\tau_1=\cdots=\tau_a=0$ is **true**, we can show that $SS_{\text{Treatment}}/ \sigma^2$ follows chi-square with $a-1$ degrees of freedom

$$
\dfrac{SS_{\text{Treatment}}}{\sigma^2} \sim \chi^2(a-1)
$$

- If the **null hypothesis** $H_0:\tau_1=\cdots=\tau_a=0$ is **true**, we can show that $F_0 = M S_{E} / M S_{\text {Treatment }}$ is distributed as F-distribution with $a-1$ and $N-a$ degrees of freedom:

$$
F_0 = \frac{ MS_{\text {Treatment }} }{ MS_{E} } 
  = \frac{SS_{\text {Treatments }} / (a-1)}{SS_E /(N-a)} 
  \sim F(a-1, N-a)
$$

### 2.5 Summary

- Hypotheses test

$$
\begin{aligned}
H_0 &: \mu_1 = \mu_2 = \cdots = \mu_a \\
H_1 &: \mu_i \neq \mu_j \text{ for at least one pair } (i, j)
\end{aligned}
$$

- The appropriate test statistic is
$$
F_0 = \frac{ MS_{\text {Treatment }} }{ MS_{E} } 
  = \frac{SS_{\text {Treatments }} / (a-1)}{SS_E /(N-a)} 
$$

- Rejection region of $H_0$ if
$$
f_0 > f_{\alpha,a-1,N-a}
$$

- Analysis of variance table
$$
\begin{array}{l|l|l|l}
\hline 
\text{ Source of Variation} 
  & \text{ Sum of Squares }
  & \text{ DoF }
  & \text{ Mean Square }
  & F_0
\\ \hline 
\text { Between treatments } 
  & S S_{\text {Treatments }} = n \sum \limits_{i=1}^a \left(\bar{y}_{i \cdot} - \bar{y}_{\cdot \cdot} \right)^2   
  & a-1 
  & MS_{\text {Treatments }} 
  & F_0 = \frac{M S_{\text {Treatments }}}{M S_E} \\
\text { Error (within treatments) } & S S_E=S S_T-S S_{\text {Treatments }}
  & N-a & MS_E 
\\
\text{ Total}
  & SS_T = \sum \limits_{i=1}^a \sum \limits_{j=1}^n \left(y_{i j}-\bar{y}_{\cdot \cdot} \right)^2 & N-1 & \\
\hline
\end{array}
$$

## 2.6 Estimation of model parameters

The linear statistical model
$$
Y_{ij} = \mu + \tau_i + \epsilon_{ij} = \mu_i + \epsilon_{ij}
$$

### 2.6.1 Parameters Estimation

- Estimate of overall mean $\mu$ and treatment effects $\tau_i$
$$
\hat{\mu} = \bar{y}_{\cdot \cdot} \ ,
\qquad 
\hat{\tau_i} = \bar{y}_{i \cdot} - \bar{y}_{\cdot \cdot}
$$

- Estimates of treatment mean
$$
\hat{\mu_i} = \hat{\mu} + \hat{\tau_i} = \bar{y}_{i \cdot}
$$

- Estimate of variance $\sigma^2$
$$
\hat{\sigma}^2 = MS_E
$$

### 2.6.2 Confidence interval

Confidence interval **on the mean of the $i$th treatment** $\mu_i$

- $\bar{Y}_{i \cdot} = \dfrac{1}{n} \sum \limits_{j=1}^{n}Y_{ij}$ is normal distributed random variable 
$$
\bar{Y}_{i \cdot} = \frac{\sum \limits_{j=1}^{n}Y_{ij}}{n} \ \sim \ \mathcal{N} \left(\mu_i, \frac{\sigma^2}{n} \right)
$$

- $T = \dfrac{\bar{Y}_{i \cdot}-\mu_i}{\sqrt{MS_E/n}}$ is a $t$-distribution with $N-a$ degrees of freedom
$$
T = \frac{\bar{Y}_{i \cdot}-\mu_i}{\sqrt{MS_E/n}} \ \sim \ t(N-a)
$$

- A $100(1-a)$ percent confidence interval on the $i$th treatment mean $\mu_i$
$$
\bar{y}_{i \cdot} - t_{ \frac{\alpha}{2}, N-a} \ \sqrt{\frac{M S_E}{n}} 
\leq \mu_i 
\leq \bar{y}_{i \cdot} + t_{\frac{\alpha}{2}, N-a} \ \sqrt{\frac{M S_E}{n}}
$$

Confidence interval on **the difference in two treatment means** $\mu_i-\mu_j$

- $\bar{Y}_{i \cdot} - \bar{Y}_{j \cdot}$ is normal distributed random variable:
$$
\bar{Y}_{i \cdot} - \bar{Y}_{j \cdot} \sim \ \mathcal{N} \left(\mu_i - \mu_j, \frac{2 \sigma^2}{n} \right)
$$

- $T = \dfrac{\bar{Y}_{i \cdot} - \bar{Y}_{j \cdot} - (\mu_i - \mu_j)}{\sqrt{2MS_E/n}}$ has a $t$-distribution with $N-a$ degrees of freedom
$$
T = \dfrac{ \left( \bar{Y}_{i \cdot} - \bar{Y}_{j \cdot} \right) - \left( \mu_i - \mu_j \right)}{\sqrt{2MS_E/n}} 
\ \sim \ t(N-a)
$$

- A $100(1-a)$ percent confidence interval on the $i$th treatment mean $\mu_i-\mu_j$
$$
\left( \bar{y}_{i \cdot} - \bar{y}_{j \cdot} \right) - t_{ \frac{\alpha}{2}, N-a} \ \sqrt{\frac{2 \, MS_E}{n}} 
\leq \mu_i 
\leq \left( \bar{y}_{i \cdot} - \bar{y}_{j \cdot} \right) + t_{\frac{\alpha}{2}, N-a} \ \sqrt{\frac{2 \, MS_E}{n}}
$$

## 2.7 Unbalanced experiments

$$
\begin{array}{rcl}
SS_T &=& \sum \limits_{i=1}^a \sum \limits_{j=1}^{n_i} y_{i j}^2 - \dfrac{y_{\cdot \cdot}^2}{N} \\
SS_{\text {Treatments }} &=& \sum \limits_{i=1}^a \dfrac{y_{\cdot i}^2}{n_i} - \dfrac{y_{\cdot \cdot}^2}{N} \\
SS_E &=& SS_T-SS_{\text {Treatments }}
\end{array}
$$
where
$$
N = \sum_{i=1}^{a}n_i
$$


## 2.8 Multiple comparisons following the ANOVA

When the null hypothesis $H_0: \mu_1 = \mu_2 = \cdots = \mu_a$ is **rejected** in the ANOVA, we need methods that can identify which means are different, called **multiple comparisons methods**.

**Fisher's least significant difference (LSD)** method

- Hypothesis: $H_0: \mu_i = \mu_j; \ H_1: \mu_i \neq \mu_j$

- Test statistic: $t_0 = \dfrac{\bar{y}_{i \cdot} - \bar{y}_{j \cdot}}{\sqrt{2 \, MS_E / n}}$

- $LSD = t_{\alpha/2, \ a(n-1)} \sqrt{2 \, MS_E / n}$. If $|\bar{y}_{i \cdot} - \bar{y}_{j \cdot}|> LSD$, the pair of means $\mu_i$ and $\mu_j$ and would be declared significant different.

- How does it work? : $T = \frac{ \left( \bar{Y}_{i \cdot} - \bar{Y}_{j \cdot} \right) - \left( \mu_i - \mu_j \right)}{\sqrt{2MS_E/n}} 
\ \sim \ t(N-a)$

- If the sample sizes are different in each treatment (i.e., **unbalanced experiments**):
$$
LSD = t_{\frac{\alpha}{2}, \ a(n-1)} \sqrt{MS_E \left( \frac{1}{n_i} + \frac{1}{n_j} \right)}
$$

## 2.9 Model Adequacy Checking

The linear statistical model
$$
Y_{ij} = \mu + \tau_i + \epsilon_{ij} = \mu_i + \epsilon_{ij}, 
\ \ i=1,\cdots,a; \ j=1,\cdots,n;
\qquad \text{where } \epsilon_{ij} \sim NID(0, \sigma^2)
$$

Checking assumptions is important

- Normality

- Constant variance

- Independence

- Have we fitted the right model?

- What to do if some of these assumptions are **violated** is beyond the lecture due to limit time

### 2.9.1 Residual

**Residual** for observation $j$ in treatment $i$
$$
e_{ij} = y_{ij} - \hat{y}_{ij}
$$

- where $\hat{y}_{ij}$ is an estimate of observation $j$ in treatment $i$:
  $$
  \hat{y}_{ij} = \hat{\mu} + \hat{\tau_i} 
  = \hat{y}_{\cdot \cdot} + (\hat{y}_{i \cdot} - \hat{y}_{\cdot \cdot})
  = \hat{y}_{i \cdot}
  $$

If the model is adequate, the residuals should be **structureless**; that is, they should contain no obvious patterns

### 2.9.2 The normality assumption test

**Normal probability plot** is a subjective visual examination of the data

- Computer software generates the residuals

- **Residual plots** are very useful

- **Normal probability plot** of residuals

### 2.9.3 Outliers

A very common defect that often shows up on normal probability plot is one residual that is **much larger** than any of the others, called an **outlier**

A **rough** check method: **standardize residuals:** $d_{ij} = e_{ij} / \sqrt{MS_E}$

- If $\epsilon_{ij} \sim \mathcal{N}(0, \sigma^2)$, then $d_{ij} \sim \mathcal{N}(0, 1)$

- $\Pr \left(-1 \le d_{ij} \le 1\right) \approx 0.68$ and $\Pr \left(-2 \le d_{ij} \le 2\right) \approx 0.95$, thus, an outlier if $|d_{ij}| \geq 3$ or $4$.

### 2.9.4 Plot of residuals in time sequence

- Plotting residuals in **time order** of data collection is helpful in detecting **correlations** between the residuals.

- A tendency to have runs of positive and negative residuals indicates positive correction.

- This would imply that the **independence assumption** on the errors has been violated

### 2.9.5 Plot of residuals versus fitted value

- If the model is correct and the assumptions are satisfied, **the residuals should be structureless**; in particular, they should be unrelated to any other variables including the predicted response.

- A simple check is to plot the residuals versus the **fitted value** $\hat{y}_{ij}$

- A defect that occasionally shows up on this plot is **nonconstant variance**.

## 2.10 Determining Sample Size

### 2.10.1 Issues in Designed Experiments

- Answer depends on lots of things

  - including what type of experiment is being contemplated, how it will be conducted, resources, and desired **sensitivity**.

- Sensitivity refers to the **difference in means** that the experimenter wishes to detect.

- Generally, **increasing** the number of **replications increases the sensitivity** or it makes it easier to detect small differences in means.

### 2.10.2 Type I error and Type II error

We can choose the sample size to detect a specific difference in means and achieve desired values of **type I and type II errors**

- **Type I error** $\alpha$ : reject $H_0$ when it is true

- **Type II error** $\beta$ : fail to reject $H_0$ when it is false

- Power $= 1-\alpha$

Calculation of type II error

$$
\begin{split}
\beta &= 1 - \Pr \left( \text{Reject} H_0 \ | \ H_0 \text{ is false} \right) \\
&= 1 - \Pr \left( F_0 > F_{\alpha, a-1, N-a} \ | \ H_0 \text{ is false} \right) \\
&= \Pr \left( F_0 \leq F_{\alpha, a-1, N-a} \ | \ H_0 \text{ is false} \right)
\end{split}
$$

It can be shown that, if $H_0$ is false, the statistic $F_0=MS_{\text{Treatments}}/MS_E$ is distributed as a **noncentral** $F$ random variable with $a-1$ and $N-a$ degrees of freedoms and the **noncentrality** parameter $\delta$, where $\delta=\sum \limits_{i=1}^{n}\tau_i^2 / \sigma^2$

- If $\delta=0$, the noncentral $F$ distribution becomes the usual (central) $F$ distribution.

### 2.10.3 Operating characteristic (QC) curves

These curves plot the probability of type II error $\beta$ against a parameter $\Phi$, where

$$
\Phi^2=\frac{n \sum \limits_{i=1}^a \tau_i^2}{a \, \sigma^2}=\frac{\delta}{a},
\qquad \text{where }
\delta=\frac{\sum \limits_{i=1}^n \tau_i^2}{\sigma^2}
$$

Note: that $\Phi^2$ is related to the noncentrality $\delta$.

- The **OC curves** for the fixed effects model

- A very common way to use these charts is to define a difference in two means $D$ of interest, then the minimum value of $\Phi^2$ is:
$$
\Phi^2 = \frac{n \, D^2}{2 \, a \, \sigma^2}
$$

- Typically work in term of the ratio of $D/\sigma$ and try values of $n$ until the **desired power** is achieved.

- Most statistics software packages can perform power and sample size calculations.

# 3. Random Effects Models

Previous slides have considered **fixed** factors

- A specific set of factor levels is chosen for the experiment

- Inference confined to those levels

- Often **quantitative** factors are fixed (why?)

When factor levels are chosen at random from a larger population of potential levels, the factor is **random**

- Inference is about the **entire population of levels**

- Industrial applications include measurement system studies

## 3.1 Formulation 

- The usual single factor ANOVA model is

$$
Y_{i j} = \mu + \tau_i + \epsilon_{i j},  \qquad
\left\{\begin{array}{l}
i=1,2, \ldots, a \\
i=1,2, \ldots, n
\end{array}\right.
$$

- Both the error term and the treatment effects are random variables:
$$
\epsilon \sim NID(0, \sigma^2), \qquad 
\tau_i \sim NID(0, \sigma_\tau^2)
$$

- Variance components:
$$
V(y_{ij}) = \sigma^2 + \sigma_{\tau}^2
$$

## 3.2 Relevant hypotheses

- In the fixed effects model we test **equality of treatment means**

- This is no longer appropriate because the treatments are randomly selected

  - The individual ones we happen to have are **not of specific interest**

  - we are interested in the **population** of treatments

- The appropriate hypotheses are:

$$
\begin{aligned}
H_0 & : \sigma_{\tau}^2 = 0 \\
H_1 & : \sigma_{\tau}^2 > 0
\end{aligned}
$$

### 3.3 Test hypotheses

- The standard ANOVA partition of the total sum of squares still works; leads to usual ANOVA display.
$$
SS_T = SS_{\text{Treatments}} + SS_E
$$

- Form of the hypothesis test depends on the expected mean squares
$$
\mathbb{E}\left[ MS_E \right] = \sigma^2,
\quad \text{and} \quad
\mathbb{E}\left[ MS_{\text{Treatments}} \right] = \sigma^2 + n \sigma_{\tau}^2
$$

- Proof:

- Therefore, the appropriate test statistic is:
$$
F_0 = \frac{MS_{\text{Treatments}}}{MS_E}
$$

### 3.4 Variance component estimation
$$
\begin{array}{cl}
\hat{\sigma}^2 &= M S_E \\
\hat{\sigma}_\tau^2 &= \dfrac{M S_{\text{Treatments }}-M S_E}{n}
\end{array}
$$

For **unequal sample size**, replace $n$ in the above equation by
$$
n_0 = \frac{1}{a-1} \left(
  \sum_{i=1}^a n_i-\frac{\sum \limits_{i=1}^a n_i^2}{\sum 
  \limits_{i=1}^a n_i}
\right)
$$

# 4 Blocking Principle

Design of Engineering Experiments - The Blocking Principle

## 4.1 The Blocking Principle

- Randomized complete block design (RCBD)

  - RCBD aims to remove the variability between coupons (block) from the experiment error.

  - RCBD: the experimenter test each tip once on each of four coupons. The word "complete" indicates that each block (coupon) contains all the treatments (tips).

- Summary for the hardness testing example

  - To conduct this experiment as a RCBD, assign all 4 tips (treatments) to each coupon (block)

  - Each coupon is called a **"block"**; that is, it is a more homogenous experimental unit on which to test the tips (treatments)

  - Variability **between** blocks can be large, variability **within** a block should be relatively small

  - In general, a **block** is a specific **level of the nuisance factor**

  - A **complete** replicate of the basic experiment is conducted in each block. The word "complete" indicates that each block (coupon) contains all the treatments (tips).

  - A block represents a **restriction on randomization**

  - All runs **within** a block **are randomized**

**Blocking** is a technique for dealing with **nuisance** factors

- A **nuisance** factor is a factor that probably has some effect on the response, but it is of no interest to the experimenter. However, the variability it transmits to the response needs to be minimized

- **Typical nuisance factors** include batches of raw material, operators, pieces of test equipment, time (shifts, days, etc.), different experimental units

- **Many** industrial experiments involve blocking (or should)

- Failure to block is a common flaw in designing an experiment (consequences?)

- If the nuisance variable is **known** and **controllable**, we use blocking

- If the nuisance factor is **known** and **uncontrollable**, sometimes we can use the **analysis of covariance** to remove the effect of the nuisance factor from the analysis

- If the nuisance factor is **unknown** and **uncontrollable** (a **"lurking" variable**), we hope that **randomization** balances out its impact across the experiment

- Sometimes several sources of variability are **combined** in a block, so the block becomes an aggregate variable

## 4.2 Extension of the ANOVA to the RCBD

### 4.2.1 Formulation

- **Assumption**: there are $a$ treatments (factor levels) and $b$ blocks

- A **statistical linear model** (fixed effects model) for the RCBD is:
$$
y_{i j}=\mu+\tau_i+\beta_j+\varepsilon_{i j}, \qquad
\left\{\begin{array}{l}
i=1,2, \ldots, a \\
j=1,2, \ldots, b
\end{array}\right.
$$

- The relevant (fixed effects) hypotheses are:
$$
H_0: \mu_1=\mu_2=\cdots=\mu_a
$$
where
$$
\text { where } \mu_i = \frac{1}{b} \sum_{j=1}^b\left(\mu+\tau_i+\beta_j\right)=\mu+\tau_i
$$

### 4.2.2 ANOVA partitioning of total variability

$$
\begin{split}
\sum_{i=1}^a \sum_{j=1}^b \left( y_{i j} - \bar{y}_{\cdot \cdot} \right)^2 &= 
  \sum_{i=1}^a \sum_{j=1}^b \left[ 
    \left( \bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot} \right)
    + \left( \bar{y}_j-\bar{y}_{\cdot \cdot} \right)
    + \left( y_{i j}-\bar{y}_{i \cdot}-\bar{y}_{\cdot j}+\bar{y}_{\cdot \cdot} \right)
  \right]^2 \\
  &= \underbrace{b \sum_{i=1}^a \left( \bar{y}_{i \cdot} -\bar{y}_{\cdot \cdot} \right)^2}_{SS_{\text {Treatments }}}
    + \underbrace{a \sum_{j=1}^b\left(\bar{y}_{\cdot j}-\bar{y}_{\cdot \cdot} \right)^2}_{SS_{\text {Blocks}}}
    + \underbrace{\sum_{i=1}^a \sum_{j=1}^b \left( y_{i j}-\bar{y}_{i \cdot}-\bar{y}_{\cdot j} + \bar{y}_{\cdot \cdot} \right)^2}_{SS_E} \\
\end{split}
$$

That is:
$$
SS_T = SS_{\text {Treatments }} + SS_{\text {Blocks }} + SS_E
$$

where:
$$
\begin{aligned}
\bar{y}_{i \cdot} &= \frac{\sum \limits_{j=1}^b y_{i j}}{b}, \quad i=1,2, \cdots, a\\
\bar{y}_{\cdot j} &= \frac{\sum \limits_{i=1} y_{i j}}{a}, \quad j=1,2, \cdots, b \\
\bar{y}_{i j} &= \frac{\sum \limits_{i=1}^a \sum \limits_{j=1}^b y_{i j}}{a b}
\end{aligned}
$$
Proof:


### 4.2.3 Degrees of freedom

$$
ab-1 = (a-1) + (b-1) + [(a-1)(b-1)]
$$

Therefore, ratios of sums of squares to their degrees of freedom result in mean squares and the ratio of the mean square for treatments to the error mean square is an $F$ statistic that can be used to test the hypothesis of equal treatment means.
$$
F_0 = \frac{MS_{\text {Treatments }}}{MS_E} \sim f(a-1, (a-1)(b-1))
$$

Three important mean squares and their expected values
$$
\begin{aligned}
MS_{\text {Treatments }} &=\frac{S S_{\text {Treatments }}}{a-1} &
  \mathbb{E} \left[M S_{\text {Treatments }}\right] &= \sigma^2+\frac{b}{a-1} \sum_{i=1}^a \tau_i^2 \\
M S_{\text {Blocks }} &= \frac{S S_{\text {Blocks }}}{b-1} &
  \mathbb{E} \left[M S_{\text {Blocks }}\right] &= \sigma^2+\frac{a}{b-1} \sum_{j=1}^b \beta_j^2 \\
MS_E &=\frac{S S_E}{(a-1)(b-1)} &
  \mathbb{E} \left[M S_E\right] &= \sigma^2
\end{aligned}
$$

**Proof:**

### 4.2 ANOVA display for the RCBD
$$
\begin{array}{l|l|l|l}
\hline 
\text{Source of Var.} & \text{Sum of Squares} & \text{DoF} & \text{Mean Square} & F_0 \\
\hline 
\text{Treatments} & SS_{\text{Treatments}} & a-1 & MS_{\text {Treatments }} = \dfrac{S S_{\text {Treatments }}}{a-1} & \dfrac{MS_{\text {Treatments }}}{MS_E} \\
\hline
\text{Blocks} & SS_{\text{Blocks}} & b-1 & MS_{\text {Blocks}} = \dfrac{SS_{\text {Blocks}}}{b-1} \\
\hline
\text{Error} & SS_E \text{ (By subtraction)} & (a-1)(b-1) & MS_{\text {Error}} = \dfrac{S S_{\text {Error}}}{(a-1)(b-1)} \\
\hline
\text{Total} & SS_T & ab-1 \\
\hline
\end{array}
$$

where
$$
\begin{aligned}
SS_T &=\sum_{i=1}^a \sum_{j=1}^b y_{i j}^2-\frac{y_{\cdot \cdot}^2}{N} \\
SS_{\text {Treatments }} &=\frac{1}{b} \sum_{i=1}^a y_{i \cdot}^2-\frac{y_{\cdot \cdot}^2}{N} \\
SS_{\text {Blocks }} &=\frac{1}{a} \sum_{j=1}^b y_{\cdot j}^2-\frac{y_{\cdot \cdot}^2}{N}
\end{aligned}
$$
the **error sum of squares** is obtained by subtraction as
$$
S S_E=S S_T-S S_{\text {Treatments }}-S S_{\text {Blocks }}
$$

## 4.3 Multiple Comparisons

When the null hypothesis $H_0:\mu_1=\mu_2=\cdots=\mu_a$ is reject in the
ANOVA, we need methods that can identify which means are different, called **multiple comparisons methods**.

- Fisher's least significant difference (LSD) method

- **Hypothesis:** $H_0:\mu_i=\mu_j; \ H_1: \mu_i\neq \mu_j$

- **Test statistic:** $t_0=\dfrac{\bar{y}_{i \cdot}-\bar{y}_{j \cdot}}{\sqrt{2 \, MS_E/b}}$

- $LSD=t_{\alpha / 2, \ (a-1)(b-1)} \sqrt{2 \, MS_E/b}$

- If $|\bar{y}_{i \cdot} - \bar{y}_{j \cdot}| > LSD$, the pair of means $\mu_i$ and $\mu_j$ would be declared significant different.

# References

Montgomery, Douglas C., and George C. Runger. "Chapter 13 Design and Analysis of Single-Factor Experiments: The Analysis of Variance", *Applied Statistics and Probability for Engineers 7th Edition*, Wiley, 2014.