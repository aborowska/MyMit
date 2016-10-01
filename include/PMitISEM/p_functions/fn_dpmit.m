function lnd = fn_dpmit(theta, pmit, partition, fn_const_X, L, GamMat)
% density evaluation of theta on the Partial MitISEM structure
% if L == true than in logs
% theta is partitioned according to partition

    d = size(theta,2);
    S = length(partition);
    
    for ss = 1:S
        [s1,s2] = fn_partition_ends(partition, d, ss);
%         lnd = lnd + dmvgt(theta(:,s1:s2), pmit(ss), true, GamMat);
        if (ss == 1)
            lnd = dmvgt(theta(:,s1:s2), pmit(ss), L, GamMat);
        else
            X = fn_const_X(theta(:,1:s1-1));
            lnd = lnd + p_dmvgt(theta(:,s1:s2), pmit(ss), L, GamMat, X);
        end
    end   

    if ~L
        lnd = exp(lnd);
    end
end