Introduction
============

This vigenette is provided to run the models described in the paper “A
joint polytomous logistic regression model for the initiation and
cessation of electronic cigarettes and conventional cigarettes with
time-dependent covariates in a longitudinal study”.

Methods
=======

Simulating data
---------------

``` r
source("model_functions.R")
# Parameter specifications 
nT <- 5 # number of time points (including the baseline)
V <- 4 # number of covariates
N <- 1000

simstat <- Initiate(nT=nT, V=V)
sData <- SimulateData(param=simstat, N=N)

X <- sData[["X"]]
Y <- sData[["Y"]]
```

Data overview
-------------

`simstat` is a list containing parameters like *λ*s, *μ*s, *α*s, *η*s,
*β*s, *γ*s, and *δ*s.

``` r
simstat
```

    $lambda
         [,1] [,2] [,3] [,4] [,5]
    [1,]   NA   -2   -2   -2   -2
    [2,]   NA   -3   -3   -3   -3

    $mu
         [,1] [,2] [,3] [,4] [,5]
    [1,]   NA    1    1    1    1
    [2,]   NA    2    2    2    2

    $alpha
         [,1] [,2] [,3] [,4]
    [1,]  0.1  0.2  0.1  0.2
    [2,]  0.1  0.1  0.3  0.2

    $eta
         [,1] [,2] [,3] [,4]
    [1,] -0.1 -0.1 -0.1 -0.1
    [2,] -0.2 -0.2 -0.2 -0.2

    $beta
    [1] 2 2

    $gamma
    [1] -2 -2

    $delta
         [,1] [,2]
    [1,]    0    0
    [2,]    0    0

`X` is a list of `N` matrices of nT × V, each of which contains `V`
covariates over `nT` time points for one observation.

``` r
head(X, n = 2)
```

    [[1]]
               [,1]        [,2]       [,3]         [,4]
    [1,] -0.8858959 -1.22181813  1.0025334 1.081659e+00
    [2,] -1.1892691 -1.33065021  0.1871055 5.740222e-01
    [3,] -0.7904083 -0.66406286  0.2305506 2.347641e-01
    [4,] -0.5913418  0.13635104  0.1849089 7.204941e-03
    [5,] -0.7507995 -0.03038998 -0.1003550 2.156583e-05

    [[2]]
                [,1]      [,2]      [,3]        [,4]
    [1,] -0.17050876 1.4315800 0.5760126  0.01589943
    [2,] -0.06022801 1.0840697 0.3740817 -0.14956214
    [3,] -0.11953409 0.9576033 0.2168042 -0.16422318
    [4,] -0.21518246 0.8032152 0.2975916 -0.48053650
    [5,] -0.64660521 0.3220836 0.2504134 -0.53618094

`Y` is a list of `N` matrices of nT × V, each of which contains two
outcomes (0: No; 1: Yes) over `nT` time points for one observation.

``` r
head(Y, n = 2)
```

    [[1]]
         [,1] [,2]
    [1,]    0    0
    [2,]    0    1
    [3,]    0    0
    [4,]    0    0
    [5,]    0    0

    [[2]]
         [,1] [,2]
    [1,]    0    0
    [2,]    0    0
    [3,]    0    0
    [4,]    0    0
    [5,]    0    0

Running the joint polytomous logistic regression models
-------------------------------------------------------

``` r
res_joint <- Analyze(simstat, X, Y)
```

Results
=======

The model yields the coefficients and standard errors of all parameters
listed by `simstat`.

``` r
# coefficients 
coef_joint <- vec2list(res_joint$par)
# standard errors 
sigma_joint <- vec2list(sqrt(diag(solve(-res_joint$hessian))))
# create 95% CI plots 
```

We displayed the odds ratio plots for important variables *β*s and *γ*s.

``` r
df <- data.frame(type = factor(
                   rep(c("Initiation", "Cessation"), each = 2), 
                   levels = c("Initiation", "Cessation")
                 ),
                 cig = c("E-cigarettes", "C-cigarettes", 
                         "E-cigarettes", "C-cigarettes"),
                 mean = c(coef_joint$beta,
                          coef_joint$gamma), 
                 sd = c(sigma_joint$beta,
                        sigma_joint$gamma))

df$lb <- exp(df$mean + qnorm(0.025)*df$sd)
df$ub <- exp(df$mean + qnorm(0.975)*df$sd)
df$mean <- exp(df$mean)

ggplot(df, aes(x = cig, y = mean)) +
  geom_point(position=position_dodge(width=0.3), size=2) +
  geom_hline(yintercept = 1.0, linetype="dotted", size=1) + 
  scale_y_continuous(trans='log2', 
                     breaks = c(round(df$lb, 1), 
                                round(df$ub, 1), 1)) +
  scale_color_manual(values = c("grey", "black")) +
  guides(color = guide_legend(reverse = TRUE)) +
  geom_errorbar(aes(ymin = lb, ymax = ub),
                position = position_dodge(0.3), width = 0.25) +
  facet_grid(type~.)+
  coord_flip() +
  ylab("Odds ratio") +
  theme_bw()+
  theme(text = element_text(size=10),
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(), 
        axis.title.y = element_blank())
```

![](vignette_files/figure-markdown_github/unnamed-chunk-7-1.png)
