mu_hl = [0.2,-2.2];

[mu, Sigma] = fn_initopt(kernel_init, mu_hl);
mu
Sigma

xx = 0:0.01:0.5;
n=length(xx);
yy = -5:0.01:0; 
m=length(yy);
[X1,X2] = meshgrid(xx,yy);
V1 = reshape(X1,n*m,1); V2 = reshape(X2,n*m,1);
V = [V1,V2];
kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, true);

hl = arrayfun(@(ii) kernel(V(ii,:)), 1:n*m, 'un', 0);
hlr = reshape(hl,m,n);
hlr = cell2mat(hlr);
hlr=hlr-max(max(hlr));
hlre = exp(hlr);     
surf(xx,yy,hlre)

[max_hlre,ind] = max(hlre);
[max_hlre2,ind2] = max(max_hlre);

[max_hlre_v,ind_v] = max(hlre');
[max_hlre2_v,ind2_v] = max(max_hlre_v);



