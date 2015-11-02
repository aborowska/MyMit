function d = posterior_arch_whole2(theta, data, S, L)
    alpha = theta(:,1);
    eps = theta(:,2);
    
    T = length(data);
    [N,~] = size(theta);
    ind = 2:T;
    
   
    d = -Inf*ones(N,1);
    h = zeros(T,1); h(1,1) = S;
    omega = S*(1-alpha); % variance targeting constraint

    for ii = 1:N
        pdf = zeros(T,1);
        h(ind) =  omega(ii,1) + alpha(ii,1)*(data(ind-1,1)).^2;
        pdf(ind) = -0.5*(log(2*pi) + log(h(ind)) + (data(ind).^2)./h(ind));
        d(ii,1) = sum(pdf) + prior(ii,2)  - 0.5*(log(2*pi) + (eps(ii,1)).^2); 
    end

end
