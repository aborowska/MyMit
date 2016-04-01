log_chi2_one = @(aa) -0.5*log(2*pi) + 0.5*(aa - exp(aa));
log_chi2_one_init = @(aa) 0.5*log(2*pi) - 0.5*(aa - exp(aa));

mit_chi.mu = 0;
mit_chi.Sigma = 1;
mit_chi.df = 1;
mit_chi.p = 1;

MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

[mit1_chi, summary1_chi] = MitISEM_new(log_chi2_one_init, log_chi2_one, -1, cont, GamMat);


cont.mit.CV_tol = 0.01;
[mit2_chi, summary2_chi] = MitISEM_new(log_chi2_one_init, log_chi2_one, -1, cont, GamMat);


figure(1)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
xx = -5:0.01:5;
yy = -5:0.01:5;
Mit2 = MitISEM_plot(mit2_chi, 2, xx, yy, GamMat);
hold on
plot(xx,exp(log_chi2_one(xx)),'r')
plot(xx,omori_plot,'g');
hold off
if v_new
    set(gca,'TickLabelInterpreter','latex')
else
    plotTickLatex2D;
end













Omori = [0.00609	1.92677	0.112650
0.04775	1.34744	0.177880
0.13057	0.73504	0.267680
0.20674	0.02266	0.406110
0.22715	-0.85173	0.626990
0.18842	-1.97278	0.985830
0.12047	-3.46788	1.574690
0.05591	-5.55246	2.544980
0.01575	-8.68384	4.165910
0.00115	-14.65000	7.333420]

Omori_p = Omori(:,1);
Omori_m = Omori(:,2);
Omori_v2 = Omori(:,3);

fn_omori = @(aa) (Omori_p')*(exp(((aa-Omori_m).^2)./(2*Omori_v2))./sqrt(2*pi*Omori_v2));

omori_plot = 0*xx;
for ii = 1:length(xx)
    omori_plot(ii) = fn_omori(xx(ii));
end
