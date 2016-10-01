function [lnd, input_X] = fn_dpmit3(input_X, pmit, partition, fn_const_X, L, GamMat)
% density evaluation of theta on the Partial MitISEM structure
% in L == true then in logs
% theta is partitioned according to partition
    theta = input_X.theta;
    d = size(theta,2);
    SS = length(partition);
    
    for ss = 1:SS
        [s1,s2] = fn_partition_ends(partition, d, ss);
%         lnd = lnd + dmvgt(theta(:,s1:s2), pmit(ss), true, GamMat);
        if (ss == 1)
            lnd = dmvgt(theta(:,s1:s2), pmit(ss), L, GamMat);
        else
            input_X.theta = theta(:,1:s1-1);                         
            [X, input_X] = fn_const_X(input_X);
            lnd = lnd + p_dmvgt(theta(:,s1:s2), pmit(ss), L, GamMat, X);
        end
    end   
    
    if (nargout == 2)
        input_X.theta = theta;                         
        [~, input_X] = fn_const_X(input_X);
    end
    
    if ~L
        lnd = exp(lnd);
    end
end