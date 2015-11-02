if plot_on
    figure(100)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.9:0.01:1.1;
    yy = [];
    Mit1 = MitISEM_plot(mit1, 2, xx, yy, GamMat);
    xlabel('$$\sigma^2$$')
    title('(WN) Approximation to the posterior density $$p(\sigma^2)$$')
end