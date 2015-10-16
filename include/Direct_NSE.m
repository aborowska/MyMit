%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRECT NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort profit/loss values of the candidate draws ascending
%     PL_opt = fn_PL(y_opt);
%     [PL_opt,ind] = sort(PL_opt);
%     % sort the corresponding weights
%     w_opt = w_opt(ind,:);
%     % normalize weights
%     w_opt = w_opt/sum(w_opt);
%     cum_w = cumsum(w_opt);
%     ind_var = min(find(cum_w > p_bar))-1; 
%     VaR_IS = PL_opt(ind_var);
%     ES_IS = sum((w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL_opt(1:ind_var));

%% Check PL function
% y_opt = draw_opt(:,2).*sqrt(h_opt);
% w_opt =  fn_ISwgts(lnk, lnd, false);
% 
% % compute profit/loss at y_opt
% PL_opt = fn_PL(y_opt);
% % compute VaR_IS

% % compute p_hat at VaR_IS
% dens = struct('y',y_opt,'w',w_opt,'c',VaR_IS);
% p_hat = fn_PL(dens, 2);

%% IS NSE & RNE
	
%     g(x) = I{PL(x)<=c} c=hat_VaR
%     g = (fn_PL(y_opt)<=VaR_IS);
G = @(a,b) (fn_PL(a)<=b);
    
%     NSE_prob_IS - NSE of the probability estimator
%     p_hat_IS - the density of PL at VaR_IS  
%     NSE_VaR_IS = NSE_prob_IS/p_hat_IS

%     RNE = (naive variance)/(NSE^2)
%     naive variance, direct sampling ==> only the nominator of NSE_VaR_IS
%     changes, p_hat_IS stays the same ==> RNE for VaR is RNE for prob
[NSE_prob_IS, RNE_VaR_IS] = fn_NSE(y_opt, w_opt, G, VaR_IS) ;
dens_IS = struct('y',y_opt,'w',w_opt,'c',VaR_IS);%,'scale',0.10);
p_hat_VAR_IS = fn_PL(dens_IS, 2);% the density of PL at VaR_IS
NSE_VaR_IS = NSE_prob_IS/p_hat_VAR_IS;

fprintf('(MitISEM) VaR NSE IS: %6.4f. \n',NSE_VaR_IS);
 
% sqrt(p_bar*(1-p_bar)/M)
% naive = sqrt(22*(0.02^2))

%% ES NSE estimation
% First: estimate the density of the ES estimator
%% Step 1: 
% Construct a grid of VaR values
grid_pn = 100;
grid_l = VaR_IS - 4*NSE_VaR_IS;
grid_u = VaR_IS + 4*NSE_VaR_IS;
VaR_grid = linspace(grid_l, grid_u, grid_pn);
VaR_grid = VaR_grid';
ES_grid = zeros(grid_pn,1);
NSE_ES_grid = zeros(grid_pn,1);
p_hat_ES_grid = zeros(grid_pn,1);
p_hat_VaR_grid = zeros(grid_pn,1);
%% Step 2: 
% for each VaR value from the grid evaluate the NSE of the ES estimator
% given the VaR value
% and
% evaluate the asymptotically valid normal density p_hat_ES_grid of the ES
% estimator on a grid
PL_opt = fn_PL(y_opt);
[PL_opt, ind] = sort(PL_opt); 
w = w_opt(ind,:);
w = w/sum(w);
cum_w = cumsum(w);

G = @(x) x;

for ii = 1:grid_pn
    VaR_curr = VaR_grid(ii,1);
    ind_curr = min(find(PL_opt > VaR_curr))-1;
    w_relevant = w(1:ind_curr)/sum(w(1:ind_curr));
    ES_grid(ii,1) = sum(w_relevant.*PL_opt(1:ind_curr));
    NSE_ES_grid(ii,1) = fn_NSE(y_opt(1:ind_curr), w_relevant, G);
    dens_curr = struct('y',y_opt(1:ind_curr),'w',w_relevant,'c',ES_grid(ii,1));%,'scale',10);
    p_hat_ES_grid(ii,1) = fn_PL(dens_curr, 2);% the density of PL at ES_grid 
    dens_IS = struct('y',y_opt,'w',w_opt,'c',VaR_curr);%,'scale',1);
    p_hat_VaR_grid(ii,1) = fn_PL(dens_IS, 2); % weights for ES estimator's density, i.e. estimated density of the VaR estimator
end
%% Step 3: 
% estimate the ES estimator's density p_hat_ES_IS as the weighted
% average of the densities p_hat_ES_grid (given VaR) with weights from the
% estimated density of the VaR estimator p_hat_VaR_grid
w_ES = p_hat_VaR_grid/sum(p_hat_VaR_grid);
p_hat_ES_IS = sum(w_ES.*p_hat_ES_grid);
% NSE_ES_IS is now estimated as the standard deviation of the estimated
% density p_hat_ES_IS
NSE_ES_IS = sqrt(sum(w_ES.*(p_hat_ES_grid-p_hat_ES_IS).^2));
NSE_ES_IS2 = sum(w_ES.*NSE_ES_grid(ii,1));

fprintf('(MitISEM) ES NSE IS: %6.4f. \n',NSE_ES_IS);





























%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRECT NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort profit/loss values of the candidate draws ascending
%     PL_opt = fn_PL(y_opt);
%     [PL_opt,ind] = sort(PL_opt);
%     % sort the corresponding weights
%     w_opt = w_opt(ind,:);
%     % normalize weights
%     w_opt = w_opt/sum(w_opt);
%     cum_w = cumsum(w_opt);
%     ind_var = min(find(cum_w > p_bar))-1; 
%     VaR_IS = PL_opt(ind_var);
%     ES_IS = sum((w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL_opt(1:ind_var));

%% Check PL function
% y_opt = draw_opt(:,2).*sqrt(h_opt);
% w_opt =  fn_ISwgts(lnk, lnd, false);
% 
% % compute profit/loss at y_opt
% PL_opt = fn_PL(y_opt);
% % compute VaR_IS

% % compute p_hat at VaR_IS
% dens = struct('y',y_opt,'w',w_opt,'c',VaR_IS);
% p_hat = fn_PL(dens, 2);

%% IS NSE & RNE
	
%     g(x) = I{PL(x)<=c} c=hat_VaR
%     g = (fn_PL(y_opt)<=VaR_IS);
G = @(a,b) (fn_PL(a)<=b);
    
%     NSE_prob_IS - NSE of the probability estimator
%     p_hat_IS - the density of PL at VaR_IS  
%     NSE_VaR_IS = NSE_prob_IS/p_hat_IS

%     RNE = (naive variance)/(NSE^2)
%     naive variance, direct sampling ==> only the nominator of NSE_VaR_IS
%     changes, p_hat_IS stays the same ==> RNE for VaR is RNE for prob
[NSE_prob_IS, RNE_VaR_IS] = fn_NSE(y_opt, w_opt, G, VaR_IS) ;
dens_IS = struct('y',y_opt,'w',w_opt,'c',VaR_IS);%,'scale',0.10);
p_hat_VAR_IS = fn_PL(dens_IS, 2);% the density of PL at VaR_IS
NSE_VaR_IS = NSE_prob_IS/p_hat_VAR_IS;

fprintf('(AdMit) VaR NSE IS: %6.4f. \n',NSE_VaR_IS);
 
% sqrt(p_bar*(1-p_bar)/M)
% naive = sqrt(22*(0.02^2))

%% ES NSE estimation
% First: estimate the density of the ES estimator
%% Step 1: 
% Construct a grid of VaR values
grid_pn = 100;
grid_l = VaR_IS - 4*NSE_VaR_IS;
grid_u = VaR_IS + 4*NSE_VaR_IS;
VaR_grid = linspace(grid_l, grid_u, grid_pn);
VaR_grid = VaR_grid';
ES_grid = zeros(grid_pn,1);
NSE_ES_grid = zeros(grid_pn,1);
p_hat_ES_grid = zeros(grid_pn,1);
p_hat_VaR_grid = zeros(grid_pn,1);
%% Step 2: 
% for each VaR value from the grid evaluate the NSE of the ES estimator
% given the VaR value
% and
% evaluate the asymptotically valid normal density p_hat_ES_grid of the ES
% estimator on a grid
PL_opt = fn_PL(y_opt);
[PL_opt, ind] = sort(PL_opt); 
w = w_opt(ind,:);
w = w/sum(w);
cum_w = cumsum(w);

G = @(x) x;

for ii = 1:grid_pn
    VaR_curr = VaR_grid(ii,1);
    ind_curr = min(find(PL_opt > VaR_curr))-1;
    w_relevant = w(1:ind_curr)/sum(w(1:ind_curr));
    ES_grid(ii,1) = sum(w_relevant.*PL_opt(1:ind_curr));
    NSE_ES_grid(ii,1) = fn_NSE(y_opt(1:ind_curr), w_relevant, G);
    dens_curr = struct('y',y_opt(1:ind_curr),'w',w_relevant,'c',ES_grid(ii,1));%,'scale',10);
    p_hat_ES_grid(ii,1) = fn_PL(dens_curr, 2);% the density of PL at ES_grid 
    dens_IS = struct('y',y_opt,'w',w_opt,'c',VaR_curr);%,'scale',1);
    p_hat_VaR_grid(ii,1) = fn_PL(dens_IS, 2); % weights for ES estimator's density, i.e. estimated density of the VaR estimator
end
%% Step 3: 
% estimate the ES estimator's density p_hat_ES_IS as the weighted
% average of the densities p_hat_ES_grid (given VaR) with weights from the
% estimated density of the VaR estimator p_hat_VaR_grid
w_ES = p_hat_VaR_grid/sum(p_hat_VaR_grid);
p_hat_ES_IS = sum(w_ES.*p_hat_ES_grid);
% NSE_ES_IS is now estimated as the standard deviation of the estimated
% density p_hat_ES_IS
NSE_ES_IS = sqrt(sum(w_ES.*(p_hat_ES_grid-p_hat_ES_IS).^2));
NSE_ES_IS2 = sum(w_ES.*NSE_ES_grid(ii,1));

fprintf('(AdMit) ES NSE IS: %6.4f. \n',NSE_ES_IS);