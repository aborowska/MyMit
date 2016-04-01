N = 500;
alpha= 0.91;
sigma = 1.0;
beta= 0.5;
V = randn(N,1);
W = randn(N,1);
X = zeros(N,1);
Y = zeros(N,1);
X(1,1) = ((sigma^2)/(1-alpha^2))*randn;

for ii = 2:N
    X(ii,1) = alpha*X(ii-1,1) + sigma*V(ii,1);
    Y(ii,1) = beta*exp(X(ii,1)/2)*W(ii,1);
end

hold on
plot(X, 'k')
plot(Y,'r')
hold off

M = 1000;