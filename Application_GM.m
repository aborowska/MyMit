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
% AdMit_Control 
% cont.resampl_on = false;
% cont.IS.opt = false;
% cont.IS.scale = [1, 0.8, 1.2];
% cont.IS.perc = [0.1, 0.15, 0.3]; 
%  
% [mit_admit, summary_admit] = AdMit_old(kernel_init, kernel, mu_init, cont, GamMat);
% save('results/GM/GM_admit.mat','mit_admit', 'summary_admit');

%% MitISEM
cont = MitISEM_Control 

[mit_mitisem, CV_mitisem] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
save('results/GM/GM_mitisem.mat','mit_mitisem', 'CV_mitisem');
MIT{3} = mit_mitisem;

load('results/GM/GM_H1.mat','mit_new','CV_new');
MIT{1} = mit_new;
CV1 = CV_new;
load('results/GM/GM_H2.mat','mit_new','CV_new');
MIT{2} = mit_new;
CV2 = CV_new;


%% Plot    
ff = figure(1);
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
subplot(2,2,1)
GM_surf = surf(x,x,GM);
set(GM_surf,'LineStyle','none')
% xlabel('$$x$$')
% ylabel('$$y$$')
campos(1.0e+03 *[0.0319   -0.0024    1.2423])
% contour(x,x,GM)
title('Target Gelman-Meng','fontsize',12)
% saveas(gcf,['results/GM/GM_target.fig'])    

for ii = 1:3
    mit = MIT{ii};
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
    subplot(2,2,ii+1)
%     figure(ii+1)
%     set(gcf,'defaulttextinterpreter','latex');
    Mit_surf = surf(x,x,Mit);
    set(Mit_surf,'LineStyle','none')
%     campos(1.0e+03 *[0.0319   -0.0024    1.2423])
    campos([40.1365   -3.6617    1.5353])

    title(['MitISEM approximation, H = ',num2str(ii)],'fontsize',12)
    set(gca,'ZTickLabel',[])
%     saveas(gcf,['results/GM/GM_H',num2str(ii),'.fig'])    
end

name = 'results/GM/GM_comb.eps';
set(gcf,'PaperPositionMode','auto');
print_fail = 1;
while print_fail 
    try                   
        print(gcf,name,'-depsc','-r0')
        print_fail = 0;
    catch
        print_fail = 1;
    end
end

%% Combine
figure(5)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');

figure(1)
mit = open(['results/GM/GM_target.fig']);  
mit_ax = gca;

figure(5)
sp = subplot(2,2,1);
sp_pos = get(sp,'position');
P = copyobj(mit_ax,sp);
set(P,'position',sp_pos)
for ii = 1:3
%     sp = subplot(2,2,ii+1);
    figure(ii+1)
    mit = open(['results/GM/GM_H',num2str(ii),'.fig']);  
%     copyobj(get(mit,'children'),sp)
end    



mit0 = open(['results/GM/GM_target.fig']);  
mit_ax0 = gca;
mit1 = open(['results/GM/GM_H',num2str(1),'.fig']);  
mit_ax1 = gca;
mit2 = open(['results/GM/GM_H',num2str(2),'.fig']);  
mit_ax2 = gca;
mit3 = open(['results/GM/GM_H',num2str(3),'.fig']);  
mit_ax3 = gca;


figure(5);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');

sp0 = subplot(2,2,1);
copyobj(allchild(mit_ax0),sp0);

sp1 = subplot(2,2,2);
copyobj(allchild(mit_ax1),sp1);

sp2 = subplot(2,2,3);
copyobj(allchild(mit_ax2),sp2);

sp3 = subplot(2,2,4);
copyobj(allchild(mit_ax3),sp3);