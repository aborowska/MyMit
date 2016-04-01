mu1 = [0,3]';
mu2 = [3,0]';
mu3 = [-3,0]';

sigma1 = [2,0.5;0.5,0.5];
sigma2 = [1,0;0,0.1];
sigma3 = [2,-0.5;-0.5,0.5];

% pi = 1/3;
n = 100;
n_n = 50;
s_n = 10;

ind = randsample((1:3),n,true);
n1 = sum(ind(ind==1));
n2 = sum(ind(ind==2));
n3 = n-n1-n2;

smpl1 = mvnrnd(mu1,sigma1,n1);
smpl2 = mvnrnd(mu2,sigma2,n2);
smpl3 = mvnrnd(mu3,sigma3,n3);

noise = s_n*2*(-0.5 + rand(n_n,2));

obs = [smpl1;smpl2;smpl3;noise];
hold all
scatter(smpl1(:,1),smpl1(:,2),'+');
scatter(smpl2(:,1),smpl2(:,2),'s');
scatter(smpl3(:,1),smpl3(:,2),'o');
scatter(noise(:,1),noise(:,2),'.');
hold off



df = 2;

m = 20;
mu_hat = zeros(1,m+1);
mu_hat(1,1) = mean(obs);
sigma_hat = zeros(1,m+1);
sigma_hat(1,1) = sqrt(var(obs));
for ii = 1:m
%       aux = ((obs - mu_hat(1,ii)).^(2))./(2*sigma_hat(1,ii).^(2));
    aux = ((obs - mu_hat(1,ii)).^(2))./(sigma_hat(1,ii).^(2));
    w = (df+1)./(df+aux); % mean of posterior Gamma
    mu_hat(ii+1) = sum(obs.*w)/sum(w);
    sigma_hat(ii+1) = sqrt(sum(w.*(obs-mu_hat(ii+1)).^(2))/n);
end
