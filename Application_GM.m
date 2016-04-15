clear all
addpath(genpath('include/'));
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x=(0:0.00001:50)'+0.00001;
GamMat = gamma(x);

A = 1; B = 0; C1 = 3; C2 = 3; 
L = true;
% mu0 = [3, 4];
% mu_init = [3, 1/3];
mu_init = [0, 0.1];
kernel_init = @(a) -GelmanMeng(a,A,B,C1,C2,L);
kernel = @(a) GelmanMeng(a,A,B,C1,C2,L);

%% AdMit
AdMit_Control 
cont.resampl_on = false;
cont.IS.opt = false;
cont.IS.scale = [1, 0.8, 1.2];
cont.IS.perc = [0.1, 0.15, 0.3]; 
 
[mit_admit, summary_admit] = AdMit_old(kernel_init, kernel, mu_init, cont, GamMat);
save('results/GM_admit.mat','mit_admit', 'summary_admit');

%% MitISEM
MitISEM_Control 

[mit_mitisem, summary_mitisem] = MitISEM(kernel_init, kernel, mu_init, cont, GamMat);
save('results/GM_mitisem.mat','mit_mitisem', 'summary_mitisem');


%% Plot    
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');

x = -1:0.05:6;
n = length(x);
[X1,X2] = meshgrid(x,x);
V1 = reshape(X1,n*n,1); V2 = reshape(X2,n*n,1);
V = [V1,V2];

GM_fix = @(a) GelmanMeng(a,A,B,C1,C2,false);
GM = arrayfun(@(ii) GM_fix(V(ii,:)), 1:n*n, 'un', 0);
GM = reshape(GM,n,n);
GM = cell2mat(GM);

for ii = 1:2
    if ii == 1
        mit = mit_admit;
    else
        mit = mit_mitisem;
    end
    [H,d]  = size(mit.mu);
    Mit = zeros(n,n);
    for h=1:H
        p_h = mit.p(h);
        mu_h = mit.mu(h,:);
        Sigma_h = reshape(mit.Sigma(h,:),d,d);
        df_h = mit.df(h);

        MitX_h = @(a) p_h*dmvt(a,mu_h,Sigma_h,df_h,GamMat);
        MitX_h = arrayfun(@(ii) MitX_h(V(ii,:)), 1:n*n, 'un', 0);
        MitX_h = reshape(MitX_h,n,n);
        MitX_h = cell2mat(MitX_h);
        Mit = Mit + MitX_h;
    end
    subplot(2,2,ii)
    surf(x,x,Mit);
    campos([40.1365   -3.6617    1.5353])
%     contour(x,x,Mit)
    if ii == 1
        title('AdMit approximation')
        Mit_admit = Mit;
        H_admit = H;
    else
        title('MitISEM approximation')
        Mit_mitisem = Mit;
        H_mitisem = H;
    end
    set(gca,'ZTickLabel',[])
end
subplot(2,2,3)
surf(x,x,GM)
campos(1.0e+03 *[0.0319   -0.0024    1.2423])
% contour(x,x,GM)
title('Target Gelman-Meng')

subplot(2,2,4)
hold on
plot(1:H_admit, summary_admit.CV, 'b-o');
plot(1:H_mitisem, summary_mitisem.conv.CV, 'k-o');
hold off
title('CoV');
legend('AdMit', 'MitISEM')


name = 'GM_results.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%% Time 
% time_draw = sum(summary.time.draw);
% time_max = sum(summary.time.max);
% time_opt = sum(summary.time.opt);
% time_list = time_draw + time_max + time_opt;
% time_all = summary.time.all;