
function Xmax = fn_Xmax(x)
    Xmax = 2.0*max(x,[],2);
    ord = floor(log(Xmax)./log(10));
    next = Xmax./10^(ord);
    next = ceil(next*2)./2;
    Xmax = next.*10.^(ord);
end