clear all
close all
M = 10000;
p_bar = 0.01;
H = 10;

%     kernel = @(a) posterior_t_garch_noS_hyper_mex(a, y, S, GamMat, hyper);
%     [theta1, accept] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
%     fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept, model, algo);
%     theta1 = theta1(BurnIn+1:M+BurnIn,:);
% 
%     h_T = volatility_t_garch_noS_mex(theta1, y, S);
%     [y_H, eps_H] = predict_t_garch_noS(theta1, y_T, h_T, H);

% load the predicted returns (generated given thetas from the MH algorithm)
% (AR usually around 65%)

y_H = csvread('Returns_t_garch_H10_adapted.csv');

ind_real = (sum(imag(y_H),2)==0);
M_real = sum(ind_real); 
fprintf('M_real = %i.\n',M_real)
y_H = y_H(ind_real,:);

PL_H_ind = fn_PL(y_H);
PL_H = sort(PL_H_ind);
VaR_prelim = PL_H(round(p_bar*M_real));
ES_prelim = mean(PL_H(round(1:p_bar*M)));

% RNE for VaR
ind_prelim = double((PL_H_ind <= VaR_prelim));
RNE_prelim = fn_RNE(ind_prelim, 'MH',[],'Q');
% RNE for ES
ind_prelim = PL_H_ind((PL_H_ind <= VaR_prelim));
RNE_ES_prelim = fn_RNE(ind_prelim, 'MH',[],'Q');