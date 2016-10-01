function X = WN_const_X(theta)
    [N, h] = size(theta);
    if (h > 1) % not the first forecast
        y = theta(:,2:h);
        y = bsxfun(@times, y, sqrt(theta(:,1)));
        X = [ones(N,1), fn_PL(y)];
    else
        X = ones(N,1);
    end
end

