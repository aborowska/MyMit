function draw = rmvt_old(mu,Sigma,df,n)
%  location/scale multivariate student t distribution 
    d = length(mu);
    draw = zeros(n,d);
    mt = @(a) rmvt_one(a,mu,Sigma,df);
    mt = arrayfun(@(ii) mt(draw(ii,:)), (1:n)', 'un', 0);
    draw = cell2mat(mt);
%     for ii = 1:n
%         draw(ii,:) = rmvt_one(mu,Sigma,df);
%     end
end

function draw = rmvt_one(draw,mu,Sigma,df)
%draws one varaible from multivariate t distribution with location mu,
% scale Sigma and df number of degrees of freedom
    d = length(mu);
    Z = randn(d,1);
    B = chol(Sigma); B = B';
    Z = B*Z; Z = Z';
    R2 = chi2rnd(df);
    draw = mu + Z/sqrt(df/R2);
end