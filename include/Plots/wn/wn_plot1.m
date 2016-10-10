if plot_on
%     http://www.cs.tut.fi/kurssit/SGN-84007/slides/Lecture4.pdf
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    kernel = @(x) posterior_debug(x, y, a, b, false);

    xx = 0.9:0.01:1.1;  n = length(xx);
    yy = -4:0.01:4;       m = length(yy);
    %     xx = 0.95:0.01:1.05;
    %     yy = -4:0.01:-2;

    [XX,YY] = meshgrid(xx,yy);
    XX1 = reshape(XX,n*m,1);
    YY1 = reshape(YY,n*m,1);

    Val1 = kernel(XX1).*normpdf(YY1);
    Val = reshape(Val1,m,n);
    P = surf(XX,YY,Val);
    set(P,'LineStyle','none')

    title('(WN) The whole density $$p_{1}(\sigma^2,\varepsilon_{T+1})$$.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp 

    %%%

    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim, false);

    xx = 0.9:0.01:1.1;  n = length(xx);
    yy = -4:0.01:4;       m = length(yy);
    %     xx = 0.95:0.01:1.05;
    %     yy = -4:0.01:-2;

    [XX,YY] = meshgrid(xx,yy);
    XX1 = reshape(XX,n*m,1);
    YY1 = reshape(YY,n*m,1);

    Val1 = kernel([XX1,YY1]);
    Val = reshape(Val1,m,n);
    P = surf(XX,YY,Val);
    set(P,'LineStyle','none')

    title('(WN) The high loss density $$p_{2}(\sigma^2,\varepsilon_{T+1})$$.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp 

    %%%
    figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.9:0.01:1.1;
    yy = -4:0.01:4;
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;
    Mit2 = MitISEM_plot(mit2, 3, xx, yy, GamMat);
    title('(WN) Approximation to the high loss density $$q_{2,Mit}(\sigma^2,\varepsilon_{T+1})$$.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp 
end