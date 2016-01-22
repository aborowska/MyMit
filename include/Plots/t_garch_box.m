v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

load(['results/',model,'_Direct.mat'], 'VaR_direct', 'ES_direct');
load(['results/',model,'_IS_AdMit.mat'], 'VaR_IS', 'ES_IS');
VaR_IS_admit = VaR_IS;
ES_IS_admit = ES_IS;
load(['results/',model,'_IS_MitISEM.mat'], 'VaR_IS', 'ES_IS');

%%%%%%%%%%%%%%%%%%%%%%%%

figure(590)
if v_new
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
else
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
end
set(gcf,'defaulttextinterpreter','latex');
boxplot([VaR_direct, VaR_IS_admit, VaR_IS],'labels',{'VaR Direct','VaR AdMit','VaR MitISEM'})        
%     title(['100*', num2str(p_bar),'\% VaR estimates: Direct, AdMit and MitISEM (',strrep(model,'_','\_'),', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])  
if v_new
    set(gca,'TickLabelInterpreter','latex')
else
    plotTickLatex2D;
end
if print_on
    name = ['figures/(',model,')', num2str(p_bar),'_VaR_3box_',num2str(M),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end

%%%%%%%%%%%%%%%%%%%%%%%%
figure(690)
if v_new
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
else
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
end
set(gcf,'defaulttextinterpreter','latex');
boxplot([ES_direct, ES_IS_admit, ES_IS],'labels',{'ES Direct','ES AdMit','ES MitISEM'})        
%     title(['100*', num2str(p_bar),'\% ES estimates: Direct, AdMit and MitISEM (',strrep(model,'_','\_'),', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])  
if v_new
    set(gca,'TickLabelInterpreter','latex')
else
    plotTickLatex2D;
end
if print_on
    name = ['figures/(',model,')', num2str(p_bar),'_ES_3box_',num2str(M),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end