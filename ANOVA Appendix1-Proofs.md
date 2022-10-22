# ANalysis Of VAriance (ANOVA) Appendix 1: Proofs

# 1 Noncentral Chi-square distribution

1 **[Noncentral Chi-square distribution](https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution)** with $k$ degrees of freedom

If $X_i \overset{i.i.d}{\sim} \mathcal{N}(\mu_i, \sigma_i^2), \forall i=1,2,\cdots,k$, then
$$
Y = \sum_{i=1}^{k}\frac{X_i^2}{\sigma^2_i} \sim {\chi'}^2_k(\delta)
$$
where $k$ is the degree of freedom; $\delta$ is **non-centrality parameter** and 
$$
\delta = \sum_{i=1}^{k}\frac{\mu_i^2}{\sigma^2_i} = \mathbb{E}[Y] - k
$$
and $\mathbb{E}[Y]$ is the expectation of the random variable $Y$
$$
\mathbb{E}[Y] = k + \sum_{i=1}^{k}\frac{\mu_i^2}{\sigma^2_i}
$$

2 [**Noncentral $F$ distribution**](https://en.wikipedia.org/wiki/Noncentral_F-distribution) with $(v_1,v_2)$ degrees of freedom and noncentral parameter

If $X_1$ follows a noncentral Chi-square distribution with the noncentral parameter $\delta$, and $X_2$ follows a Chi-square distribution. i.e., $X_1 \sim {\chi'}^2_{k_1}(\delta)$ and $X_1 \sim {\chi'}^2_{k_2}$, then the following random variable $F$ follows a **noncentral $F$ distribution**:
$$
F = \frac{X_1/k_1}{X_2/k_2} \sim{F'}_{k_1, k_2}(\delta)
$$

3 Type II error in ANOVA

$$
\beta = \mathrm{Pr} \left(F_0 \leq F_{\alpha, a-1, N-a} \mid H_0 \text{ is false } \right)
$$
where
$$
F_0 = \frac{ MS_{\text {Treatments }}}{MS_E}
= \frac{ \left(SS_{\text {Treatments }} / \sigma^2 \right) / (a-1)}{ \left(SS_E / \sigma^2 \right) / (N-a)}
\quad \sim \quad 
{F'}_{a-1,N-a} \left(\frac{n}{\sigma^2} \sum_{i=1}^a \tau_i^2 \right)
$$

4 Relation between $\Phi^2$ and $\delta$
$$
\delta = \Phi^2 \cdot a
$$

# 2 Proofs

## 2.1 Fixed effect model

#### Proof of $SS_T = SS_{\text{Treatments}} + SS_E$

$$
\begin{aligned}
SS_T &= 
\sum_{i=1}^a \sum_{j=1}^n ( y_{i j}-\bar{y}_{i \cdot}+\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot} )^2 
\\
&= \sum_{i=1}^a \sum_{j=1}^n \left[ 
  ( \bar{y}_{i \cdot} - \bar{y}_{\cdot \cdot} )^2 
  + 2 (\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot} ) (y_{i j}-\bar{y}_{i \cdot} ) 
  + (y_{i j}-\bar{y}_{i \cdot})^2 
  \right]
\\
&= \underbrace{\sum_{i=1}^a \sum_{j=1}^n ( \bar{y}_{i \cdot} - \bar{y}_{\cdot \cdot} )^2 }_{SS_{\text{Treatments}}}
  + 2 \sum_{i=1}^a \sum_{j=1}^n \left[ (\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot}) (y_{i j} - \bar{y}_{i \cdot} ) \right]
  + \underbrace{ \sum_{i=1}^a \sum_{j=1}^n ( y_{i j} - \bar{y}_{i \cdot} )^2 }_{SS_E}
\\
&= SS_{\text {Treatments }} + SS_E 
  + 2 \sum_{i=1}^a \left\{ \left(\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot}\right) \times \left[\sum_{j=1}^n\left(y_{i j}-\bar{y}_{i \cdot}\right)\right] \right\}
\\
&=S S_{\text {Treatments }}+ SS_E \qquad 
  \text{because }\sum_{j=1}^n\left(y_{i j}-\bar{y}_{i \cdot}\right)=0
\end{aligned}
$$

#### Proof: $SS_{\text{Treatments}} / \sigma^2 \sim \chi^2_{a-1}$ if $H_0$ is true ($\tau_1 = \tau_2 = \cdots = \tau_a$ or $\mu_1 = \mu_2 = \cdots = \mu_a$)

Assume $H_0$ is true, we have $y_{ij} \overset{i.i.d.}{\sim} \mathcal{N} \left(\mu, \sigma^2 \right)
$. Thus, we have
$$
\bar{y}_{i \cdot} = \frac{1}{n} \sum_{j=1}^{n} y_{ij}
\quad \overset{i.i.d.}{\sim} \quad
\mathcal{N} \left(\mu, \frac{\sigma^2}{n} \right)
$$
Then let 
$$
z_i = \frac{\bar{y}_{i \cdot} - \mu}{\sqrt{\sigma^2 / n }} \quad \overset{i.i.d.}{\sim} \quad
\mathcal{N}(0, 1)
$$
Consider
$$
\begin{aligned}
SS_{\text {Treatmens }}
&= n \sum_{i=1}^a (\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot} )^2
\\
&= n \sum_{i=1}^a (\bar{y}_{i \cdot} - \mu + \mu - \bar{y}_{\cdot \cdot} )^2
\\
&= n \left[ \sum_{i=1}^a (\bar{y}_{i \cdot}-\mu)^2
  - 2  \sum_{i=1}^a [(\bar{y}_{i \cdot}-\mu) (\bar{y}_{\cdot \cdot} - \mu ) ] 
  + \sum_{i=1}^a (\bar{y}_{\cdot \cdot} - \mu )^2\right] 
\\
&= n \left[ \sum_{i=1}^a (\bar{y}_{i \cdot} - \mu )^2-a (\bar{y}_{\cdot \cdot}-\mu )^2 \right]
\\
&= n \sum_{i=1}^a\left(\bar{y}_{i \cdot}-\mu\right)^2-n a\left(\bar{y}_{\cdot \cdot}-\mu\right)^2 
\end{aligned}
$$

Thus, 
$$
\begin{aligned}
& \frac{n }{\sigma^2} \sum_{i=1}^a (\bar{y}_{i \cdot}-\mu)^2
\\
=& \frac{n}{\sigma^2} \sum_{i=1}^a (\bar{y}_{i \cdot}-\bar{y}_{\cdot \cdot})^2
	+ \frac{n a}{\sigma^2} \left(\bar{y}_{\cdot \cdot}-\mu \right)^2
\\
=& \underbrace{ \frac{n}{\sigma^2} \sum_{i=1}^a \left[
	(\bar{y}_{i \cdot} - \mu) 
	- \frac{1}{a} \sum \limits_{i=1}^a (\bar{y}_{i \cdot} - \mu ) \right]^2}_{SS_{\text {Treatmens }} / \sigma^2}
	+ \frac{n}{a \sigma^2} \left[\sum_{i=1}^a (\bar{y}_{i \cdot}-\mu ) \right]^2
\\
=& \underbrace{ \sum_{i=1}^{a} (z_i-\bar{z})^2}_{SS_{\text {Treatmens }} / \sigma^2}
   + \frac{1}{a} \left(\sum_{i=1}^a z_i \right)^2 = \sum_{i=1}^{a} z_i^2
\\
=& \left[z_1, z_2, \cdots, z_a\right] 
	\times \left[\mathbf{I}_{a \times a} - \frac{1}{a} \mathbf{1}_{a \times a} \right] \times \left[ \begin{array}{l}
z_1 \\ z_2 \\ \vdots \\ z_a \end{array} \right]
	+ \left[z_1, z_2, \cdots, z_a\right] \times \left[\frac{1}{a} \mathbf{1}_{a \times a} \right] \times \left[\begin{array}{l} z_1 \\ z_2 \\ \vdots \\ z_a \end{array} \right]
\end{aligned}
$$
where $\bar{z} = \frac{1}{a} \sum_{i=1}^{a} z_i$ and matrix $\mathbf{I}_{a \times a}$ is a $a \times a$ identity matrix; $\mathbf{1}_{a \times a}$ is a $a \times a$ all-ones matrix. 

Thus $\text{Rank}\left[\mathbf{I}_{a \times a} - \frac{1}{a} \mathbf{1}_{a \times a} \right] = a - 1$ and $\text{Rank}\left[\frac{1}{a} \mathbf{1}_{a \times a} \right]=1$

Thus, 
$$
\frac{SS_{\text{Treatments }}}{\sigma^2} 
= \sum_{i=1}^a \left(z_i - \bar{z} \right)^2 \sim \chi_{a-1}^2
$$
according to the **Cochran's theorem**.
