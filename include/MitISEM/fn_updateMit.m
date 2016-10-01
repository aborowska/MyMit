function mit_m = fn_updateMit(mit, mit_nc, p)
% combine old mixture and new component: rescale probabilities and
% concatenate the rest

% mit - old mixture with H components
% mit_nc - new mit density 
% mit_m - update mixture including all components

    mit_m.mu = [mit.mu; mit_nc.mu];
    mit_m.Sigma = [mit.Sigma; mit_nc.Sigma];
    mit_m.df = [mit.df, mit_nc.df];
    if nargin > 2
        mit_m.p = p;
    else
        p_nc = mit_nc.p; % probability of the new component
        mit_m.p = [(1-p_nc)*mit.p, p_nc];
    end
end