%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

data = csvread('GSPC_ret.csv');
data = 100*data;

% QERMit ARCH 
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data);

mu_init = 0.03;
M = 10000;
N_sim = 20;

plot_on = false;
print_on = false;

if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(1998,2000+104/365,T);
    plot(xx, data)
    set(gca,'XTickLabel',{'1998',' ','1999',' ','2000',' '})
	title('S$$\&$$P 500 log-returns $$y$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/arch_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

AdMit_Control
cont.IS.opt = true;
cont.IS.scale = [1, 0.25, 4]; 
cont.IS.perc = [0.25, 0.30, 0.35];

cont.IS.opt = true;
cont.Hmax = 10;
cont.dfnc = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.: 
L = true;
kernel_init = @(a) - posterior_arch(a, data, S, L);
kernel = @(a) posterior_arch(a, data, S, L);

[mit1, summary1] = AdMit(kernel_init, kernel, mu_init, cont, GamMat);

% for QERMit 2.: draws from the joint
[draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1, kernel);
draw1 = [draw1,randn(M,1)];


%% Fig. 5
if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
	xx = 0:0.01:0.5;
	yy = -5:0.01:5;
	Mit1 = MitISEM_plot(mit1, 2, xx, yy, GamMat);
	xlabel('$$\alpha$$')
	title('(AdMit) Approximation to the posterior density $$p(\alpha)$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/q1_mit_admit.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1b.:
% generate set of draws of alpha using independence MH with
% candidate from AdMit; 
% then simulate returns from normal with variance based on the draw of alpha 

[alpha, accept] = Mit_MH(M+1000, kernel, mit1, GamMat);
fprintf('(AdMit) MH acceptance rate: %4.2f. \n',accept);
alpha = alpha(1001:M+1000);

f_stdev = @(aa) sqrt(S+(y_T^2-S)*aa);
stdev = f_stdev(alpha);
eps_T1 = randn(M,1);
y_T1 = stdev.*eps_T1;

% get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
p_bar = 0.01; % p_bar = 1-alpha, 100alpha% VaR
[PL_T1, ind] = sort(fn_PL(y_T1));
alpha_hl = alpha(ind,:);
alpha_hl = alpha_hl(1:p_bar*M,:);
eps_hl = eps_T1(ind,:);
eps_hl = eps_hl(1:p_bar*M,:);

VaR_prelim = PL_T1(p_bar*M);

% PL_hp= sort(fn_PL(y_hp));
% VaR_prelim_hp = PL_hp(p_bar*M);
fprintf('(AdMit) Preliminary VAR estimate: %6.4f. \n',VaR_prelim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% high loss
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)
c = 100*log(1 + VaR_prelim/100);
eps_bord = c./stdev;
high_loss_subspace = [alpha, eps_bord];

%% Fig. 6.1
if plot_on
	figure(2)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
	xx = min(draw1(:,1)):0.01:max(draw1(:,1));
	yy = c./f_stdev(xx);
	hold on
	xlabel('$$\alpha$$')
	ylabel('$$\varepsilon_{T+1}$$')
	% scatter(high_loss_subspace(:,1),high_loss_subspace(:,2),10,'r.')
	plot(xx, yy, 'r')
	% h = area(xx,yy);
	% h.FaceColor = [1,0.4,0.3];
	% h.EdgeColor = 'red';
	scatter(draw1(:,1),draw1(:,2),10,'k.')
	hold off
	title('(AdMit) Joint density $$p(\alpha,\varepsilon_{T+1}|y)$$ and hight loss subspace')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')

	if print_on
		name = 'figures/high_loss_admit.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1c.: 
% High loss density approximation
% L = true;
kernel_init = @(a) - posterior_arch_hl(a, data, S, VaR_prelim, L);
kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, L);

mu_hl = [0.15,-3];

cont2 = cont; % to write in the output later on
cont2.Ns = 1e4;
cont2.Hmax = 10;
cont2.CV_tol = 0.01;
cont2.dfnc = 5;
cont2.IS.scale = [  0.35  ]; 
cont2.IS.perc = [0.3 0.25   ]; 
[mit2, summary2] = AdMit(kernel_init, kernel, mu_hl, cont2, GamMat);

%% Plot hl
if plot_on
	figure(3)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
	xx = 0:0.01:0.5;
	yy = -5:0.01:5;
	Mit2 = MitISEM_plot(mit2, 3, xx, yy, GamMat);
	title('(AdMit) Approximation to the high loss density $$q_{2,Mit}(\alpha,\varepsilon_{T+1})$$.')
	xlabel('$$\alpha$$')
	ylabel('$$\varepsilon_{T+1}$$','Position',[0.5,5.1,-10.6])
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	campos([-3,-30,65]);
	if print_on
		name = 'figures/q2_mit_admit.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 2: 
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estimate VaR and ES for theta and y (or alpha in eps)

plot_on = false;
print_on = false;
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sim = 1:N_sim
    kernel = @(a) posterior_arch(a, data, S, L);
    [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M/2, mit1, kernel);
    draw1 = [draw1,randn(M/2,1)];
    
    kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, L);
	[draw2, lnk2, ind_red2] = fn_rmvgt_robust(M/2, mit2, kernel);
	draw_opt = [draw1; draw2];

    if (sim == N_sim)
        plot_on = true;
        print_on = true;
    end

	%% Figures
	if plot_on
		figure(4)
		set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
		set(gcf,'defaulttextinterpreter','latex');
		hold on
		scatter(draw_opt(:,1),draw_opt(:,2),10,'k.')
		xlabel('$$\alpha$$')
		ylabel('$$\varepsilon_{T+1}$$')
		hold off
		title('(AdMit) Draws from the optimal importance density $$q_{opt}(\alpha,\varepsilon_{T+1}|y)$$.')
%        plotTickLatex2D;
		set(gca,'TickLabelInterpreter','latex')
		if print_on
			name = 'figures/q_opt_admit.png';
			fig = gcf;
			fig.PaperPositionMode = 'auto';
			print(name,'-dpng','-r0')
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		figure(5)
		set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
		set(gcf,'defaulttextinterpreter','latex');
		xx = 0:0.01:0.5; xx = xx';
		yy = -5:0.01:5;
		n = length(xx);
		mit_ep = normpdf(yy);
		Mit1 = diag(Mit1)*repmat(mit_ep,n,1); Mit1 = Mit1';
		[X1,X2] = meshgrid(xx,yy);
		MP = surf(X1,X2,Mit1);
		set(MP,'LineStyle','none')
		xlabel('$$\alpha$$')
		ylabel('$$\varepsilon_{T+1}$$','Position',[0.57,6.3,-3.6])
		title('(AdMit) Joint density $$p(\alpha,\varepsilon_{T+1}|y)$$.')
%         plotTickLatex2D;
		set(gca,'TickLabelInterpreter','latex')
		campos([-3,-32,20]);

		if print_on
			name = 'figures/joint_admit.png';
			fig = gcf;
			fig.PaperPositionMode = 'auto';
			print(name,'-dpng','-r0')
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%% 

		figure(6)
		set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
		set(gcf,'defaulttextinterpreter','latex');
		Mit_opt = 0.5*Mit1 + 0.5*Mit2;
		surf(X1,X2,Mit_opt,'EdgeColor','interp'); %colormap(bone); colormap(hot)
		title('(AdMit) Approximation to the optimal posterior density')
%       plotTickLatex2D;
		set(gca,'TickLabelInterpreter','latex')
		campos([-3,-32,40]);
		xlabel('$$\alpha$$')
		ylabel('$$\varepsilon_{T+1}$$','Position',[0.13,1.8,-2.4])

		if print_on
			name = 'figures/q_opt_mit_admit.png';
			fig = gcf;
			fig.PaperPositionMode = 'auto';
			print(name,'-dpng','-r0')
		end    
		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		figure(7)
		set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
		set(gcf,'defaulttextinterpreter','latex');
		subplot(2,1,1)
		plot(draw_opt(:,1))
		title('$$\alpha_1$$')
		subplot(2,1,2)
		plot(draw_opt(:,2))
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% IS weights
	L = true;
	kernel = @(a) posterior_arch_whole(a, data, S, L);
	lnk = kernel(draw_opt);
	exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1),mit1,false, GamMat);
	exp_lnd2 = 0.5*dmvgt(draw_opt,mit2,false, GamMat);
	exp_lnd = exp_lnd1 + exp_lnd2;
	lnd = log(exp_lnd);
	w_opt =  fn_ISwgts(lnk, lnd, false);

	%% VaR and ES IS estimates
	h_opt = S*(1-draw_opt(:,1)) + draw_opt(:,1)*(y_T^2);
	y_opt = draw_opt(:,2).*sqrt(h_opt);
	dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
	IS_estim = fn_PL(dens, 1);
	VaR_IS(sim,1) = IS_estim(1,1);
	ES_IS(sim,1) = IS_estim(1,2);
    PL_opt= fn_PL(y_opt);

    if plot_on    
        figure(8)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gcf,'defaulttextinterpreter','latex');
        hold on
        plot(sort(PL_opt))
        pos =  max(find(sort(PL_opt)<=VaR_IS(sim,1)));
        scatter(pos, VaR_IS(sim,1),'MarkerFaceColor','red')
        hold off
    %     subplot(3,1,1)
    %     subplot(3,1,2)
    %     subplot(3,1,3)
        title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
        set(gca,'TickLabelInterpreter','latex')

        if print_on
            name = 'figures/arch_predict_admit.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
	fprintf('AdMit IS VAR estimate: %6.4f. \n',VaR_IS(sim,1));
	fprintf('AdMit IS ES estimate: %6.4f. \n',ES_IS(sim,1));
end

mean_VaR_IS = mean(VaR_IS);
mean_ES_IS = mean(ES_IS);
    
NSE_VaR_IS = std(VaR_IS);
NSE_ES_IS = std(ES_IS);
    
fprintf('AdMit IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);
fprintf('AdMit IS VAR (mean) estimate: %6.4f. \n',mean_ES_IS);
    
fprintf('AdMit IS NSE VaR estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('AdMit IS NSE ES estimate: %6.4f. \n',NSE_ES_IS);

cont = cont2;
model = 'arch';
s = 'admit';
hp = 1;
gen_out;


clear GamMat x fig
save(['results/arch_admit_',num2str(hp),'.mat']);

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);