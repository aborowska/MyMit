kernel = @(x) posterior_debug(x, y, a, b, true);
[sigma1, accept ] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
sigma1 = sigma1((BurnIn+1):(M+BurnIn));

%% Future disturbances
eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1

% the returns coresponding to disturbances: y = sqrt(sigma)*eps
y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
% preliminary VaR
[PL, ind] = sort(fn_PL(y_H));
VaR_prelim(bb,1) = PL(p_bar*M);  
ES_prelim = mean(PL(1:p_bar*M));    
fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(bb,1), model, algo);
% theoretical VaR (with the FLAT prior) (or, the likelihood)
VaR_true = norminv(p_bar,0,sqrt(H));
fprintf('Theoretical 100*%4.2f%% VaR (under flat prior): %6.4f.\n', p_bar,VaR_true);


clear y_H
% High loss draws = the target of the truncated H-days-ahead return distibution
draw = [sigma1, eps_H];
clear eps_H sigma1
draw = draw(ind,:);
draw_hl = draw(PL<=VaR_prelim(bb,1),:);  



tic

M = 10000;
MM = 2*M;
kernel = @(x) posterior_debug(x, y, a, b, true);
[sigma1, accept ] = Mit_MH(100*M+BurnIn, kernel, mit1, GamMat);
fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
sigma1 = sigma1((BurnIn+1):(100*M+BurnIn));

eps_H = randn(MM,H); % --> if future disturbances drawn from the target then their weights are 1
% ind_mm = (1:MM)' + (ii-1)*MM; 
sigma_new = sigma1(1:MM);
draw_mm = [sigma_new,eps_H];
y_H = bsxfun(@times,eps_H,sqrt(sigma_new)); 
[PL, ind] = sort(fn_PL(y_H)); 
ind_good = ind(1:M);
ind_bad = ind(M+1:MM);
% ind_good_old = ind_good(ind_good<=M);
ind_good_new = ind_good(ind_good>M);
ind_bad_old = ind_bad(ind_bad<=M);
draw_mm(ind_bad_old,:) = draw_mm(ind_good_new,:);
    
    
for ii = 3:100
    ind_mm = (1:M)' + (ii-1)*M; 
    sigma_new = sigma1(ind_mm,:);
%     eps_H = randn(M,H); 
    draw_mm(M+1:MM,:) = [sigma_new,randn(M,H)];

    y_H = bsxfun(@times,draw_mm(:,2:H+1),sqrt(draw_mm(:,1))); 
    [PL, ind] = sort(fn_PL(y_H)); 
    ind_good = ind(1:M);
    ind_bad = ind(M+1:MM);
    % ind_good_old = ind_good(ind_good<=M);
    ind_good_new = ind_good(ind_good>M);
    ind_bad_old = ind_bad(ind_bad<=M);
    fprintf('\n Length ind_good_new: %d\n',length(ind_good_new))
    draw_mm(ind_bad_old,:) = draw_mm(ind_good_new,:);
end

y_H = bsxfun(@times,draw_mm(:,2:H+1),sqrt(draw_mm(:,1))); 
[PL, ind] = sort(fn_PL(y_H)); 
VaR_est = PL(M);
toc


%% BigDraw function
M = 10000;
BurnIn = 1000;
% WN model
kernel = @(x) posterior_debug(x, y, a, b, true);
% y_H = y_predict(draw_mm); 
y_predict = @(drawu) bsxfun(@times,draw(:,2:end),sqrt(draw(:,1))); 
   

% arch model
kernel = @(xx) posterior_arch(xx, data, S, true);
% y_H = y_predict(draw_mm); 
y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  


[draw_hl, VaR_prelim, y_H, PL] = BigDraw(M, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);

