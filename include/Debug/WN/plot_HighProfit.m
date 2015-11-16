  if plot_on2
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    if strcmp(model,'WN_hp')
        kernel = @(x) posterior_debug_hp(x, y, a, b, VaR_prelim, false);
        xx = 0.9:0.01:1.1;  n = length(xx);
        yy = -4:0.01:4;       m = length(yy);
    else
        kernel = @(x) posterior_arch_hp(x, data, S, VaR_prelim, false);
        xx = 0:0.01:0.5;    n = length(xx);
        yy = -5:0.01:5;     m = length(yy);
    end

    [XX,YY] = meshgrid(xx,yy);
    XX1 = reshape(XX,n*m,1);
    YY1 = reshape(YY,n*m,1);

    Val1 = kernel([XX1,YY1]);
    Val = reshape(Val1,m,n);
    P = surf(XX,YY,Val);
    set(P,'LineStyle','none')
    if strcmp(model,'WN_hp')
        title('(WN) The high profit density $$p_{2}(\sigma^2,\varepsilon_{T+1})$$.')
        xlabel('$$\sigma^2$$')
    else
        title('(arch) The high profit density $$p_{2}(\alpha,\varepsilon_{T+1})$$.')
        xlabel('$$\alpha$$')
    end
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');

    if strcmp(model,'WN_hp')
        xx = 0.9:0.01:1.1;
        yy = -4:0.01:4;
        Mit3 = MitISEM_plot(mit3, 3, xx, yy, GamMat);
        title('(WN) Approximation to the high profit density $$q_{2,Mit}(\sigma^2,\varepsilon_{T+1})$$.')
        xlabel('$$\sigma^2$$')
    else
        xx = 0:0.01:0.5;    n = length(xx);
        yy = -5:0.01:5;     m = length(yy);    
        Mit3 = MitISEM_plot(mit3, 3, xx, yy, GamMat);     
        title('(arch) Approximation to the high profit density $$q_{2,Mit}(\alpha,\varepsilon_{T+1})$$.')
        xlabel('$$\alpha$$')        
    end
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp   
  end