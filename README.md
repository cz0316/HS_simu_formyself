# HS_simu_formyself

This package was built solely for the purpose of simplifying my own code.

## Installation

```{r}
library(devtools)
devtools::install_github('https://github.com/cz0316/HS_simu_formyself',force = T)
library(HEARTSVG)
```


## demo

```{r}

zi.p=0.6
size1=0.5
mu1=1.5
try_1=new.zi3(spots = 5000,se=1000,ns=5000,type='ZINB',se.mu=mu1*3,ns.mu = mu1,
             lambda=0.6,se.size = size1*3,ns.size = size1,se.p = zi.p/3,ns.p=zi.p,
             ptn='curves_2.png',png_dir =  '~/Desktop/pattern_png/')
```






