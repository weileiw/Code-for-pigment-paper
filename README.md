#  Inverse model for inferring particle exchange rates

The code here is to infering particle respiration rate constants along side 
with particle aggregation and disaggragation rate constant, by using Baysian 
statistic method. The data used are from MedFlux program, and sampled using 
large volume pumps. 

### Box_model.m
 
 This is the original model, with objective functions weighted using 
 data's own standard deviations. Large particle divergence is modeled
 based on Krist et al paper. Martin Curve exponential is optimized.
 Different from Krist et al paper *(**function buildPFD_v2**)* is 
 here we use particle 
 disaggregation rate constant as the loss term of large particles.

 What I found here is the respiration rate constant of both Chla and
 phyeopigment are not well constrained because
 1) too much dependency on their priors.
 2) too large error bars.

**Due to these issues, I tried to applied normal distributino to
data. Here come Box_model_log**

This time I get slightly better model to observation fit. But 
the finial Hessian matrix is close to zero. try to fix it.

