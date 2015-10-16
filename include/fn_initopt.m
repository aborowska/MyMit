function [mu, Sigma, val] = fn_initopt(kernel, mu0) % , options)
    d = length(mu0);
    [x,val,~,~,~,hessian] = fminunc(kernel,mu0); %,options);
%     options = optimset('Display','iter');
%     x = fminsearch(kernel_init,mu_init, options);
%     x = fminsearch(kernel_init,mu_hl), options);
%     [x,~,~,~,~,hessian] = fminunc(kernel,mu0,options)
    mu = x;
	Sigma = inv(hessian);
    Sigma = reshape(Sigma,1,d^2);
    
end