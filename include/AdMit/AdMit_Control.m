%% AdMit Control
cont.Ns = 10000; %1e5; % number of draws used in the simulation
cont.Np = 1000; %1e3; % number of draws used when estimating the probabilities
cont.Hmax = 10;

cont.CV_tol = 0.1;
cont.norm = true;

% cont.IS.opt = false; % use importance sampling (if false: use usual optimization)
cont.IS.opt = true; % <-- ARCH

cont.IS.scale = [1, 0.25, 4]; % <-- ARCH
% cont.IS.scale = [1, 0.8, 1.2]; % scaling coefficients for the covariance matrix in IS sampling
% cont.IS.perc = [0.1, 0.15, 0.3]; % percentages of weigths used to compute the scale matirix in IS sampling
cont.IS.perc = [0.25, 0.30, 0.35]; % <-- ARCH

cont.pnc = 0.1; % probability of the new component
cont.dfnc = 1;  % degrees of freedom of a new component

cont.resampl_on = true;