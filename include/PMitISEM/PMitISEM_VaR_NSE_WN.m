N_sim = 20;
VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);
M = 10000;
BurnIn = 1000;

% Preliminary
for sim = 1:N_sim   
    fprintf('\nVaR prelim iter: %d\n',sim)
    kernel = @(x) posterior_debug(x, y, a, b, true);
    [sigma1, accept ] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
    sigma1 = sigma1(BurnIn+1:M+BurnIn);

    %% Future disturbances
    eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1
    draw = [sigma1, eps_H];

    % the returns coresponding to disturbances: y = sqrt(sigma)*eps
    y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
    % preliminary VaR
    [PL, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL(p_bar*M);  
    ES_prelim(sim,1) = mean(PL(1:p_bar*M));    
end

% IS
for sim = 1:N_sim   
    fprintf('\nVaR IS iter: %d\n',bb)
     
    sigma1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = randn(M/2,H); 
    draw1 = [sigma1, eps1];

 
%     kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); 
    draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  

%     y_pmit = bsxfun(@times,draw_pmit(:,2:d),sqrt(draw_pmit(:,1)));  % future returns 
%     [PL_pmit, ~] = sort(fn_PL(y_pmit));
% 
%     y_opt = [y_H(1:M/2); y_pmit];
%     PL_opt =  sort(fn_PL(y_opt));
    draw_opt = [draw1; draw_pmit];

    kernel = @(x) posterior_debug(x, y, a, b, true);
    lnk_opt = kernel(draw_opt(:,1));

    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights
    [s1, s2] = fn_partition_ends(partition, d, 1);
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,s1:s2), mit1, true, GamMat));
%     exp_lnd2 = dmvgt(draw_opt(:,s1:s2),pmit(1), true, GamMat);
%     for ss = 2:S
%         [s1,s2] = fn_partition_ends(partition, d, ss);
%         exp_lnd2 = exp_lnd2 + dmvgt(draw_opt(:,s1:s2),pmit(ss), true, GamMat);
%     end
    exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);

    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    y_opt = bsxfun(@times,draw_opt(:,2:d),sqrt(draw_opt(:,1)));  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);   
end



if plot_on
%     VaR_IS_pmit_comp = VaR_IS; 
%     load('figures/PMitISEM/WN1_0.01_VaR_3results20.mat')
    figure(6) 
    set(gcf,'units','normalized','outerposition',[0 0 0.75 0.75]);   
    set(gcf,'defaulttextinterpreter','latex');
%     boxplot([VaR_prelim,VaR_IS,VaR_IS_pmit,VaR_IS_pmit_comp],'labels',{'VaR prelim','VaR Toy PMit','VaR 1comp PMit','VaR multicomp PMit'})
    boxplot([VaR_prelim,VaR_IS],'labels',{'VaR prelim','VaR full pmit'})
    plotTickLatex2D;
    name = ['figures/PMitISEM/',model,'_H', num2str(H),'_VaR_box_B',num2str(B),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end

if save_on
    name = ['figures/PMitISEM/',model,'_', num2str(p_bar),'_VaR_4results',num2str(B),'.mat'];
    save(name,'VaR_prelim','VaR_IS','VaR_IS_pmit','VaR_IS_pmit_comp')
end

if plot_on
    y_pmit0 = [y(end)*ones(M,1),y_pmit];

    ret_pmit = cumsum(y_pmit0,2);
    ind_red = (ret_pmit(1:500,1+H) <= mean(VaR_prelim));

    figure(99)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'defaulttextinterpreter','latex');

    hold on
    plot(0:H,ret_pmit(~ind_red,:)','k')
    plot(0:H,ret_pmit(ind_red,:)','r')
    plot(0:H,mean(VaR_prelim)*ones(1,1+H),'m','LineWidth',2) 
    hold off
    xlabel('Forecast horizon') % x-axis label
    ylabel('Cumulative return') % y-axis label

    [tick_sort,ind_tick] = sort([mean(VaR_prelim), get(gca, 'YTick')]);
    % set(gca, 'YTick', sort([VaR_prelim, get(gca, 'YTick')])); 
    new_label = get(gca, 'YTickLabel');
    new_label = ['VaR';new_label];
    new_label = new_label(ind_tick,:);
    set(gca, 'YTick', tick_sort); 
    set(gca,'YTickLabel',new_label)
    plotTickLatex2D;
end


if plot_on
   
    figure(777)
    set(gcf,'units','normalized','outerposition',[0 0 0.75 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    
    hold on
%     axis(limit)
%     set(gca,'XTick',0:2000:10000)
%     set(gca,'YTick',-10:2:10)   
    
    plot(PL_opt,'b')
    pos =  max(find(PL_opt<=mean(VaR_prelim)));
    scatter(pos, mean(VaR_prelim),'MarkerEdgeColor','green','MarkerFaceColor','green')    
    pos =  max(find(PL_opt<= mean(VaR_IS)));
    scatter(pos,  mean(VaR_IS),'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold off
    plotTickLatex2D; 

    name = ['figures/PMitISEM/',model,'_H', num2str(H),'_predict_B',num2str(B),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end
