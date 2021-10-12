# Changes to the original pmwg package
New functionality for the pmwg package/sampler
Changes to 'old' pmwg:
1) Immediately notable: Broke down into more legible (but less elegant) code

2) Parallelized the sampling of initial startpoints

3) Parallelized creating the conditional multivariate normal

4) Started scaling epsilon adaptively to meet a certain set acceptance rate: pstar. 
Taken from: Garthwaite, P. H., Fan, Y., & Sisson, S. A. (2016). Adaptive optimal scaling of Metropolis–Hastings algorithms using the Robbins–Monro process. Communications in Statistics-Theory and Methods, 45(17), 5098-5111.

This scaling was done to set the acceptance rate for every subject. 
This should also fix the not leaving adaptation problem, since every subject is accepting new particles at a certain rate,
therefore they will get new particles with a fixed speed. 
Related added:
- alphastar (quantity related to pstar)
- n0 (determines after how many iterations were going to start adjusting epsilon)
- updated.epsilon function
- Epsilon per subject is stored in the pmwgs object
*Note that it doesn't work well for the sampling stage, since there the conditional generates the most of the particles,
and the conditional isn't scaled by epsilon. 

5) Started keeping track of where the samples were coming from in the new_particles function.
Their origin is then stored in the pmwgs object

6) Implemented c++ code using Rcpp armadillo for the density and the sampling from the multivariate normal 
(dmvnrm_arma_fast & mvrnorm_arma). Test the density, because now it is calculated at a different way that is less numerically stable. The c++ code also contains an alternative: dmvnorm_arma which is the R mvtnorm::dmvnorm way of calculating it, but then
faster (because it is in C). 
dmvnrm_arma taken from: https://gallery.rcpp.org/articles/dmvnorm_arma/         
mvrnorm_arma taken from: https://gallery.rcpp.org/articles/simulate-multivariate-normal/

7) Everything should still work with the old pmwg objects, functions and function arguments to the run_stage function. Minor
changes to pmwg handling:
- automatic epsilon determination is now done at the init function rather than run_stage. This is because I'm storing epsilon
in the pmwgs object, so I also need one for the init step. 
- You can now set the desired acceptance probability pstar (works best for burn-in and adaptation stages). If no pstar is set
it will not scale epsilon at all, and keep it at it's heuristic or predefined value. 
- If epsilon is not set in the run_stage function, it will take the last epsilon value of the pmwgs object supplied to it. So
for burn-in that would be the init epsilon, and for adaptation this would be burn-in epsilon etc.. This makes sure that when
pstar is set for stage A, but not for stage B, the epsilon that resulted from flexibly adjusting in stage A is carried over to
stage B, but then no longer adjusted flexibly. 
- pmwgs object now also stores epsilon values
- pmwgs object now also stores accepted particle origins

8) Added single subject estimation (non-hierarchical). Run with run_stage_single, no adaptation stage because we can't (reasonably) construct a conditional. Still very much built on the hierarchical architecture, so has quite a lot of code redundancy and storing redundancies. 

9) Fixed a mistake that was in the old pmwg in the inverse gamma mixing weights a_half

10) Added diagonal only estimation of the covariance matrix, so variances only. With new Gibbs step and priors. 

Overall speed-up of 2,3,6 is very roughly estimated around 5%
