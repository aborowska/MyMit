clc;
clear all;
close all;

% Univariate GAS model for volatility. See GasVolaUnivMain.pdf for explanation of the program
% Author: Rutger Lit
tic
global GAUSS STUD_T SIGMA LOG_SIGMA INV_FISHER INV_SQRT_FISHER HESS SAND 
GAUSS = 0; STUD_T = 1; SIGMA = 2; LOG_SIGMA = 3; INV_FISHER = 4; INV_SQRT_FISHER = 5; HESS = 6; SAND = 7;

%*************************************************************************************************
%   START OF USER INPUT	
%   Load data										   
%	Specify scaling, distribution, link function, order of gas model,	   
%	choice of standard errors and starting values 					   
%*************************************************************************************************

% Load data
% mdata = xlsread('DJInd19801999.xls');  
% vy = mdata(2:end,5)';

y = csvread('GSPC_ret_tgarch.csv');
vy = 100*y';

dscaling = 1; % Scaling data can improve stability, 1 for no scaling
vy = vy.*dscaling;
% Distribution: GAUSS, STUD_T
idistribution = GAUSS;
idistribution = STUD_T;

% Link function: SIGMA (f_t=sigma^2_t), LOG_SIGMA (f_t=log(sigma^2_t))
ilinkfunction = SIGMA; %LOG_SIGMA;
% Scaling score: INV_FISHER, INV_SQRT_FISHER
iscalingchoice = INV_FISHER;
% Order of GAS model p, q
ip = 1; iq = 1;
% Standard erros: HESS, SAND
istderr = HESS;	

% Starting values (note the dimensions of ip and iq)
domega = 0.01;  
vA = 0.10; % Extend for higher orders of p, use vector vA = [a ; b ; c ; ..]; 
vB = 0.89; % Extend for higher orders of q, use vector vB = [a ; b ; c ; ..]; B_1 + ... + B_q < 1
dmu = 0;
ddf = 8; % Only estimated if idistribution = STUD_T
	
%*************************************************************************************************
%  END OF USER INPUT							      
%*************************************************************************************************

cT = size(vy,2);
vinput = [idistribution; ilinkfunction; iscalingchoice; ip; iq; istderr];
vp0 = [domega; vA; vB; dmu; ddf];

[vp0, aparnames] = StartingValues(vinput, vp0);
options = optimset('TolX', 0.0001, 'Display', 'iter', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
objfun = @(vp)(-LogLikelihoodGasVolaUniv(vp, vinput, vy));
[vp_mle, dloglik] = fminunc(objfun, vp0, options);

fprintf ('Log Likelihood value = %g \r', -dloglik*cT)
[vpplot, vse] = StandardErrors(objfun, vp_mle, cT, vinput, vy);
horzcat(aparnames, num2cell(horzcat(vpplot, vse)))
PlotSeries(vp_mle, vinput, vy, cT);

toc
