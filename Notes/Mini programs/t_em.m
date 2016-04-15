clear all;
n = 100;
mu = 1;
sigma = 0.8;
obs = mu + sigma*randn(n,1);
n_out = 10;
mu_out = 10;
out = randsample(n,n_out);
out = sort(out);
obs(out) = mu_out + 2*sigma*rand(n_out,1);

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

subplot(2,2,1)
plot(obs)
subplot(2,2,3)
plot(w)
subplot(2,2,2)
plot(mu_hat)
subplot(2,2,4)
plot(sigma_hat)
