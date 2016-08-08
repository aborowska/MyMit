close all
clear all
model = 'arch';
% Data
data = csvread('GSPC_ret.csv');
data = 100*data;
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data); % data variance for the variance targeting
         
horizons = [10, 20, 40, 100, 250];
p_bar = 0.01;
N_sim = 20;
M = 10000;
VaRs = zeros(length(horizons),1);
ii = 0;
for H = horizons;
    ii = ii+1;
    name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'draw_hl','VaR_est');
    Draw_hl{ii} = draw_hl;
    VaRs(ii,1) = VaR_est;
end

alpha = zeros(M,length(horizons));
for ii = 1:length(horizons)
    alpha(:,ii) = Draw_hl{1,ii}(:,1);
end
hist(alpha)

eps10 = zeros(M,length(horizons));
for ii = 1:length(horizons)
    eps10(:,ii) = Draw_hl{1,ii}(:,11);
end
hist(eps10)


draw_hl_10 = Draw_hl{1,5}(:,1:11);
draw_hl_20 = Draw_hl{1,5}(:,1:21);
draw_hl_40 = Draw_hl{1,5}(:,1:41);
draw_hl_100 = Draw_hl{1,5}(:,1:101);
draw_hl_250 = Draw_hl{1,5}(:,1:251);

y_10 = predict_arch(draw_hl_10(:,1), y_T, S, 10, draw_hl_10(:,2:10+1));  
y_20 = predict_arch(draw_hl_20(:,1), y_T, S, 20, draw_hl_20(:,2:20+1));  
y_40 = predict_arch(draw_hl_40(:,1), y_T, S, 40, draw_hl_40(:,2:40+1));  
y_100 = predict_arch(draw_hl_100(:,1), y_T, S, 100, draw_hl_100(:,2:100+1));  
y_250 = predict_arch(draw_hl_250(:,1), y_T, S, 250, draw_hl_250(:,2:250+1));  

PL10 = sort(fn_PL(y_10));
PL20 = sort(fn_PL(y_20));
PL40 = sort(fn_PL(y_40));
PL100 = sort(fn_PL(y_100));
PL250 = sort(fn_PL(y_250));

plot(PL10);
hold on
plot(PL10*0+VaRs(1,1),'r')
hold off


plot(PL100);
hold on
plot(PL100*0+VaRs(4,1),'r')
hold off



plot(PL250);
hold on
plot(PL250*0+VaRs(5,1),'r')
hold off



bar([-VaRs,sqrt(horizons'/10)*(-VaRs(1,1))])
