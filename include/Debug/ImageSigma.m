function Sigma = ImageSigma(mit,h)
    d2 = size(mit.Sigma(h,:),2);
    d = sqrt(d2);
    Sigma = reshape(mit.Sigma(h,:),d,d);
    surf(Sigma);
end