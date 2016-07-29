clear;
name = 'consumption1.data';

% Preparation of variables
S = load(name);
[Start_row, Start_col] = find(S==1953);
[End_row, End_col] = find(S==1996);
S = S(Start_row(1,1)-2:End_row(4,1),4:5);
Y = S(:,1);
C = S(:,2);
y = log(Y);
c = log(C);
dy = diff(y);
dc = diff(c);
ldy = dy(1:length(dy)-1,1);
dy = dy(2:length(dy),1);
dc = dc(2:length(dc),1);
n = length(dc);
iota = ones(n,1);

% OLS estiamtion
X = [iota, dy, ldy];

X = ones(n,1);
dc = double(ind_direct);
[n,k] = size(X);

beta = (X'*X)^(-1)*X'*dc;

sigma2 = (dc-X*beta)'*(dc-X*beta)/(n-k);
D = inv(X'*X);
d = diag(D);
se = sqrt(sigma2*d);

t = beta./se;

% Test for serial correlation
res = dc - X*beta;
lres = [0 ; res(1:length(res)-1)]; 

X_test = [X, lres];
[n_test,k_test] = size(X_test);

beta_test = (X_test'*X_test)^(-1)*X_test'*dc;

sigma2_test = (dc-X_test*beta_test)'*(dc-X_test*beta_test)/(n_test-k_test);
D_test = inv(X_test'*X_test);
d_test = diag(D_test);
se_test = sqrt(sigma2_test*d_test);

t_test = beta_test./se_test;
% Two-tailed t-test
p_test = 2*(1-tcdf(abs(t_test),n_test-k_test));


% HAC Newey-White estimator
% HAC covariance matrices p = 1,...,8
 Gamma = zeros(k,k,9);

for p = 0:8
    for t = p+1:n
        Gamma(:,:,p+1) = Gamma(:,:,p+1) + res(t)*res(t-p)*X(t,:)'*X(t-p,:);
    end
end
Gamma = Gamma/n;

Sigma = zeros(k,k,8);
Var_beta = zeros(k,k,8);
HAC_se = zeros(8,1);

for p = 1:8
    Sigma(:,:,p) = Gamma(:,:,1);
    for j = 1:p
        Sigma(:,:,p) = Sigma(:,:,p) + (1-j/(p+1))*(Gamma(:,:,j+1) + Gamma(:,:,j+1)');
    end
    Var_beta(:,:,p) = n*(X'*X*Sigma(:,:,p)^(-1)*X'*X)^(-1);
    HAC_se(p,:) = sqrt(diag(Var_beta(:,:,p)));
end



[N, d] = size(theta);
x = ones(N,1); % one regressand --> nvar = 1!
k = 1;

xpxi = 1/N; % xpxi = inv(x'*x);
% mean
beta = mean(theta,1); % beta    = xpxi*(x'*y);
% yhat    = x*results.beta;
% demeaned
resid = theta - repmat(beta,N,1); % resid   = y - results.yhat;
% demeaned'*demeaned
sigu = sum(resid.^2,1); %sigu = results.resid'*results.resid;
sige    = sigu/(N-1);
% perform Newey-West correction
% emat = [];
% for i=1:nvar;
% emat = [emat
%         results.resid'];
% end;
% emat = resid';
% demeaned' * ones'
hhat = resid';%hhat=emat.*x';

% G = zeros(1,d); %G=zeros(nvar,nvar);

w = zeros(lag,1);
% a=0;
G = sum(resid.^2,1);
for a = 1:nlag
% while a~=nlag+1;
%     ga = zeros(1,d);
%     ga = zeros(nvar,nvar);
    w(a,1) = (lag+1-a)/(lag+1);
%     w(nlag+1+a,1)=(nlag+1-a)/(nlag+1);
    za = sum(resid((a+1):N,:).*resid(1:N-a,:),1);
%     za=hhat(:,(a+1):nobs)*hhat(:,1:nobs-a)';
%     if a==0;
%         ga=ga+za;
%         G = za; 
%     else
        G = G + w(a,1)*(2*za);
%         ga=ga+za+za';
%     end;
%     G = G +w(1+a,1)*ga;
%     a=a+1;
end; % end of while

V_NW = G/(N^2) %V=xpxi*G*xpxi;