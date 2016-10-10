clear all
cd ../
addpath('include/');
addpath('results/arch');


load('arch_mitisem_1.mat', 'VaR_IS', 'PL_opt', 'sim') 
VaR_mitisem = VaR_IS(sim,1);
PL_mitisem = PL_opt;

load('arch_admit_1_small.mat', 'VaR_IS', 'PL_opt', 'sim') 
VaR_admit = VaR_IS(sim,1);
PL_admit = PL_opt;

load('arch_direct_1.mat', 'VaR_direct', 'PL_direct', 'sim') 
VaR_direct = VaR_direct(sim,1);
 
clear sim VaR_IS PL_opt
N = size(PL_admit,1);
xx = 1:1:N;
limit = [1, N,-15, 15];
%%
for ii = 1:3
    figure(ii)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    switch ii
        case 1
            axis(limit)
            plot(xx,sort(PL_mitisem))
            pos = max(find(sort(PL_mitisem)<= VaR_mitisem));
            scatter(pos, VaR_mitisem,'MarkerFaceColor','red')
            name = 'figures/arch_predict.png';
        case 2
            axis(limit)
            plot(xx,sort(PL_admit))
            pos = max(find(sort(PL_admit)<= VaR_admit));
            scatter(pos, VaR_admit,'MarkerFaceColor','red')
            name = 'figures/arch_predict_admit.png';
        case 3
            axis(limit)
            plot(xx,PL_direct)
            pos = max(find(PL_direct<= VaR_direct));
            scatter(pos, VaR_mitisem,'MarkerFaceColor','red')
            name = 'figures/arch_predict_direct.png';
    end
    hold off
%     subplot(3,1,1)
%     subplot(3,1,2)
%     subplot(3,1,3)
    title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
    set(gca,'TickLabelInterpreter','latex')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end