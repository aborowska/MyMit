function mit_opt = fn_qermit_opt(mit1, mit2)
% combine approximation to the whole density and to the high loss density:
% rescale probabilities and concatenate the rest

% mit_opt - optimal 50-50% importance density
% mit1 - approximation to the whole density  
% mit2 - approximation to the high loss density
    mit_opt.p = [0.5*mit1.p, 0.5*mit2.p];
    mit_opt.mu = [mit1.mu; mit2.mu];
    mit_opt.Sigma = [mit1.Sigma; mit2.Sigma];
    mit_opt.df = [mit1.df, mit2.df];
end