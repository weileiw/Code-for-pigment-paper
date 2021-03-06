#  Inverse model for inferring particle exchange rates

The code here is to infering particle respiration rate constants along 
with particle aggregation and disaggragation rate constants, by using Bayesian 
inversion method. The data used here are from MedFlux program, and sampled using 
large volume pumps. 

### Box_model.m

**Use with neglogpost.m**

This is the original model, with objective functions weighted using 
data's own standard deviations. Large particle flux divergence is modeled
based on Krist et al paper. Martin Curve exponential(*b*) is optimized.
Different from Krist et al paper *(**function buildPFD_v2.m**)*, 
we use particle disaggregation rate constant as the loss term of 
large particles.

What we found here is the respiration rate constant of both Chla and
phyeopigment are not well constrained because
1) too much dependency on their priors.
2) too large error bars.

> *b*  = 0.9690$^{0.33}_{0.25}$;

> *d1* = 1.1176$^{1.67}_{0.67}$

> *d2* = 1.0000$^{3.58}_{0.78}$

> *d3* = 1.0000$^{3.58}_{0.78}$

> *a*  = 2.7893$^{4.21}_{1.68}$

> *d*  = 64.0167$^{92.86}_{37.89}$

To compare, here are the priors for each parameter 
> *b_prior*  = 0.84;

> *d1_prior* = 1.00;

> *d2_prior* = 1.00;

> *d3_prior* = 1.00;

> *a_prior*  = 3.00;

> *d_prior*  = 150;

Due to these issues, values of *d2* and *d3* are taken from references
(Wang et al.,2017), and only *d1* (POC) remineralization rate constant 
is optimized. 

Here comes the following model
(parameter denotation has conflicit with that in the paper, fix it!!!)

### Box_model_4p.m

**Use with neglogpost_4p.m**

In this model, only four parameters are optimized, that are 
Martin curve exponential (*b*), POC remineralization rate constant
(*d1*), small particle aggregation rate constant (*a*), and large 
particle disaggregation rate constant (*b*).

>*b*  = 0.91$^{+0.08}_{-0.07}$; 

>*d1* = 1.89$^{0.58}_{0.44}$; 

>*a*  = 4.64$^{1.46}_{1.11}$; 

>*d*  = 95.78$^{28.16}_{21.76}$; 

This is the model that is finally used in the paper.


The following two sets of model code do not have detailed
annotations. Please refer to Box_model.m and Box_model_4p.m
for information.
### Box_model_log.m

**Use with neglogpost_log.m buildPFD_v2.m**

Due to the huge concentration ranges, lognormal distribution of concentrations are
assumed in this model. However, the error bars are huge, and there is not significant 
change to model versus observation correlation. This model is not used in the paper. 


### Box_cons_SV.m 

**Use with negelogpost_cons_SV.m; PFD_cond_SV.m**

Since Armstrong et al., (2009) and Xue et al., (2009) have demonstrated that particles 
caught using sediment traps are sinking at a constant speed. The code here is to test 
if a constant sinking speed for large particles can fit the data better. However, the 
model does not have a unique solution. Sinking speed is positively correlated with particle
exchange rates (POC remineralization and disaggregation). This is not hard to understand 
because fast sinking particles with fast exchange rates have the same effect with slow 
sinking particles having slow exchange rates.

Example

Sinking = 100 m/d

> *d1* = 19.1857

> *d2* = 1.0169

> *d3* = 1.0001

> *a*  = 3.0649

> *d*  = 124.0919

Sinking = 200 m/d

> *d1*  = 38.3702

> *d2*  = 1.0823

> *d3*  = 1.0010

> *a*   = 2.8700

> *d*   = 229.9081

As in the variable sinking speed model, *d2* and *d3* are not well constrained by the model.


