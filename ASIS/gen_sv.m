function [y,h] = gen_sv(T, par_SV, parametrisation)

    c = par_SV(1,1);
    phi = par_SV(1,2);
    sigma2 = par_SV(1,3);
    sigma = sqrt(sigma2);
    
    y = zeros(T,1);
    h = zeros(T+1,1);
    
    if strcmp(parametrisation,'C') % centred parametrisation
        h(1,1) = c + sqrt(1/(1-phi^2))*randn;
    
        for ii = 2:T+1
            h(ii,1) = c + phi*(h(ii-1) - c) + sigma*randn;        
            y(ii-1,1) = exp(h(ii,1)/2)*randn; 
        end
    
    else % noncentred parametrisation
            
    end
end