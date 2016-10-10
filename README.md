#QERMit

__Quick Evaluation of Risk by Mixtures of t__ - a project originating from my MPhil thesis. Bayesian evaluation of Value at Risk and Expected Shortfall for a given volatility model. Modifiaction and extension of the paper by Hoogerheide and van Dijk (2010). 


_Modification:_ changing the underlying approximation algorithm from AdMit (Hoogerheide et al., 2007) to a more flexible and accurate one, MitISEM (Hoogerheide et al., 2012).

_Extensions:_

* to allow for nonlinear non-Gassian state space models (in the spirit of Barra et al., 2016) with the NAIS state sampler (Koopman et al., 2015);
* to allow for long run risk predictions (like one-month-ahead or one-year-ahead), based on partial candidate construction (with the Partial MitISEM algorithm of Hoogerheide et al., 2012).

To speed up computations, I "MEXed" some of the MATLAB codes (i.e. they are written in C).

###References 

Barra, I.,  L. F. Hoogerheide, S. J. Koopman and A. Lucas (2016), "Joint Bayesian Analysis of Parameters and States in Nonlinear non-Gaussian State Space Models", _Journal of Applied Econometrics_, 31, forthcoming.

Hoogerheide, L. F., J. F. Kaashoek, and H. K. van Dijk (2007), "On the Shape of Posterior Densities and Credible Sets in Instrumental Variable Regression Models with Reduced Rank: an Application of Flexible Sampling Methods using Neural Networks", _Journal of Econometrics_, 139, 154-180.

Hoogerheide, L. F., A. Opschoor, and H. K. van Dijk (2012), "A Class of Adaptive Importance Sampling Weighted EM Algorithms for Effcient and Robust Posterior and Predictive Simulation", _Journal of Econometrics_, 171, 101-120.

Hoogerheide, L. F. and H. K. van Dijk (2010), "Bayesian Forecasting of Value at Risk and Expected Shortfall using Adaptive Importance Sampling", _International Journal of Forecasting_, 26, 231-247.

Koopman, S. J., A. Lucas and M. Scharth (2015), "Numerically Accelerated Importance Sampling for Nonlinear Non-Gaussian State Space Models", _Journal of Business and Economic Statistics_, 33, 114-127.