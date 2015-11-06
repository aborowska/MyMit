cont2.mit.Hmax = 1;
[mit2_1, summary2_1] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
cont2.mit.Hmax = 2;
[mit2_2, summary2_2] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
cont2.mit.Hmax = 3;
[mit2_3, summary2_3] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
cont2.mit.Hmax = 10;
[mit2_10, summary2_10] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(33)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,2,1)
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.85:0.01:1.15;
    yy = -5:0.01:5;
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;
    Mit2 = MitISEM_plot(mit2_1, 3, xx, yy, GamMat);
    title('(White noise logreturns) Mit2, H=1.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    plotTickLatex2D;
    
        
    subplot(2,2,2)
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.85:0.01:1.15;
    yy = -5:0.01:5;
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;
    Mit2 = MitISEM_plot(mit2_2, 3, xx, yy, GamMat);
    title('(White noise logreturns) Mit2, H=2.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    plotTickLatex2D;

    
    subplot(2,2,3)
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.85:0.01:1.15;
    yy = -5:0.01:5;
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;
    Mit2 = MitISEM_plot(mit2_3, 3, xx, yy, GamMat);
    title('(White noise logreturns) Mit2, H=3.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    plotTickLatex2D;
    
    
    subplot(2,2,4)
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.85:0.01:1.15;
    yy = -5:0.01:5;
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;
    Mit2 = MitISEM_plot(mit2_10, 3, xx, yy, GamMat);
    title('(White noise logreturns) Mit2, H=4 (10 max).')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    plotTickLatex2D;
    

%%%%%%%%%%%%%%%%%%%%

mean_hl_1 = sum(mit2_1.mu.*repmat(mit2_1.p',1,2),1)
PL_mu_hl_1 = fn_PL(sqrt(mean_hl_1(1,1))*mean_hl_1(1,2))

mean_hl_2 = sum(mit2_2.mu.*repmat(mit2_2.p',1,2),1)
PL_mu_hl_2 = fn_PL(sqrt(mean_hl_2(1,1))*mean_hl_2(1,2))

mean_hl_3 = sum(mit2_3.mu.*repmat(mit2_3.p',1,2),1)
PL_mu_hl_3 = fn_PL(sqrt(mean_hl_3(1,1))*mean_hl_3(1,2))

mean_hl_10 = sum(mit2_10.mu.*repmat(mit2_10.p',1,2),1)
PL_mu_hl_10 = fn_PL(sqrt(mean_hl_10(1,1))*mean_hl_10(1,2))














mean_hl_2 = sum(mit2.mu.*repmat(mit2.p',1,2),1)
PL_mu_hl_2 = fn_PL(sqrt(mean_hl_2(1,1))*mean_hl_2(1,2))
