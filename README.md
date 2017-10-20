#  Inverse model for inferring particle exchange rates

The code here is to infering particle respiration rate constants along side 
with particle aggregation and disaggragation rate constant, by using Baysian 
statistic method. The data used are from MedFlux program, and sampled using 
large volume pumps. 

### Box_model.m

**Use with neglogpost.m**

This is the original model, with objective functions weighted using 
data's own standard deviations. Large particle divergence is modeled
based on Krist et al paper. Martin Curve exponential is optimized.
Different from Krist et al paper *(**function buildPFD_v2.m**)* is 
here we use particle 
disaggregation rate constant as the loss term of large particles.

What I found here is the respiration rate constant of both Chla and
phyeopigment are not well constrained because
1) too much dependency on their priors.
2) too large error bars.

> *b*  = 0.9690$^{0.33}_{0.25}$;

> *k1* = 1.1176$^{1.67}_{0.67}$

> *k2* = 1.0000$^{3.58}_{0.78}$

> *k3* = 1.0000$^{3.58}_{0.78}$

> *a*  = 2.7893$^{4.21}_{1.68}$

> *d*  = 64.0167$^{92.86}_{37.89}$

To compare, here are the priors for each parameter 
> *b_prior*  = 0.84;

> *k1_prior* = 1.00;

> *k2_prior* = 1.00;

> *k3_prior* = 1.00;

> *a_prior*  = 3.00;

> *d_prior*  = 150;

Due to these issues, I take values for r2 and r3 from references
(Wang et al.,2017), and only optimize r1 (POC) remineralization
rate constant.
Here comes the following model

### Box_model_4p.m

**Use with neglogpost_4p.m**

In this model, only four parameters are optimized, that are 
Martin curve exponential (*b*), POC remineralization rate constant
(*k1*), small particle aggregation rate constant (*a*), and large 
particle disaggregation rate constant (*b*).

>*b*  = 0.91$^{+0.16}_{-0.14}$; 

>*k1* = 1.12$^{0.81}_{0.47}$; 

>*a*  = 2.50$^{1.80}_{1.05}$; 

>*b*  = 53.00$^{36.88}_{21.75}$; 


### Box_model_log.m

**Use with neglogpost_log.m buildPFD_v2.m**

This time I get slightly better model to observation fit. But 
the finial Hessian matrix is close to zero. try to fix it.

The Hessian problem can be fixed by using prior constrain on 
parameters. However, the error bars are huge. 


### Box_cons_SV_log.m 

Model with a constant sinking speed.



