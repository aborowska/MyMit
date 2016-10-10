if (plot_on && (sim == N_sim))
    figure(6)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.9:0.01:1.1; xx = xx';
    yy = -4:0.01:4;
    n = length(xx);
    mit_ep = normpdf(yy);
    Mit1_aug = diag(Mit1)*repmat(mit_ep,n,1); Mit1_aug = Mit1_aug';
    [X1,X2] = meshgrid(xx,yy);
    Mit_opt = 0.5*Mit1_aug + 0.5*Mit2;
    surf(X1,X2,Mit_opt,'EdgeColor','interp'); %colormap(bone); colormap(hot)
    title('(WN) Approximation to the optimal posterior density')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
    shading interp 
end