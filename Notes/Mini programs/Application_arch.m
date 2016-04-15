clear all
addpath('include/');
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

data = csvread('GSPC_ret.csv');
data = 100*data;
% QERMit ARCH 
I_arch = find(data<=-5.5);
data = data(1:576,1);
S = var(data);

mu_init = 0.03;

L = true;
kernel_init = @(a) - posterior_arch(a, data, S, L);
kernel = @(a) posterior_arch(a, data, S, L);

%% AdMit
% AdMit_Control
% 
% [mit_admit, summary_admit] = AdMit(kernel_init, kernel, mu_init, cont);
% save('results/arch_admit.mat','mit_admit', 'summary_admit');

%% MitISEM 
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

MitISEM_Control

[mit_mitisem, summary_mitisem] = MitISEM(kernel_init, kernel, mu_init, cont);
save('results/arch_mitisem.mat','mit_mitisem', 'summary_mitisem');
