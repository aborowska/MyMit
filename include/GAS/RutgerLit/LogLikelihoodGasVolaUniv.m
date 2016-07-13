function [dloglik, vf, vscore, vscaledsc, vfunc] = LogLikelihoodGasVolaUniv(vp, vinput, vy)
	global GAUSS STUD_T SIGMA SAND;
	idistribution = vinput(1);
	ilinkfunction = vinput(2);
	iscalingchoice = vinput(3);
	ip = vinput(4);
	iq = vinput(5);
	istderr = vinput(6);
	cT = size(vy, 2);
	dloglik = 0;
	vf = zeros(1, cT); vscore = zeros(1, cT); vscaledsc = zeros(1, cT);
	[domega, vA, vB, dmu, ddf] = TransformPar(vp, ilinkfunction, idistribution, ip, iq);
	imaxpq = max(ip, iq);
	vf(1:imaxpq) = domega/(1 - sum(vB)); % Unconditional mean of factor f
	[vf, vscore, vscaledsc] = Init(vf, vscore, vscaledsc, dmu, imaxpq, ddf, ilinkfunction, vy, idistribution, iscalingchoice); 
	vfunc = zeros(cT-imaxpq, 1);
	for t = imaxpq+1:cT				
		vf(t) = domega + vscaledsc(t-ip:t-1)*vA + vf(t-iq:t-1)*vB;
		if(ilinkfunction == SIGMA)
			dsigma2_t = vf(t);
		else
			dsigma2_t = exp(vf(t));
		end
		[dscore, dinvfisher, dllik] = SystemMatrices(vy(t)-dmu, dsigma2_t, ddf, idistribution, ilinkfunction);		
		if(istderr == SAND) 
			vfunc(t-imaxpq) = dllik; 
		end
		vscore(t) = dscore;
		dloglik = dloglik + dllik;
		dscale = Scaling(dinvfisher, iscalingchoice);
		vscaledsc(t) = dscale*dscore;
	end
	if(idistribution == GAUSS)
		dloglik = double((-log(2*pi)*(cT-imaxpq)/2 + dloglik)/cT);
	elseif(idistribution == STUD_T)
		dloglik = double(dloglik/cT);
	end
	if(isfinite(dloglik) == 0 || isreal(dloglik) == 0 || sum(vB) > 1)
		dloglik = -100;
	end    
end

    
function [domega, vA, vB, dmu, ddf] = TransformPar(vP, ilinkfunction, idistribution, ip, iq)
	global GAUSS STUD_T SIGMA;
	if ilinkfunction == SIGMA
		domega = exp(vP(1));
	else
		domega = vP(1);
	end
	vA = flipud(vP(2:ip+1));
	vB = flipud(vP(1+ip+1:1+ip+iq));
	if idistribution == GAUSS
		dmu = vP(1+ip+iq+1);
		ddf = 0;	
	elseif idistribution == STUD_T	
		dmu = vP(1+ip+iq+1);
		ddf = vP(1+ip+iq+2);
	end
end
    

function [vf, vscore, vscaledsc] = Init(vf, vscore, vscaledsc, dmu, imaxpq, ddf, ilinkfunction, vy, idistribution, iscalingchoice)
	global SIGMA;
	for p = 1:imaxpq
		if ilinkfunction == SIGMA
			dsigma2_t = vf(p);
		else
			dsigma2_t = exp(vf(p));
	end 
		[dscore, dinvfisher, dllik] = SystemMatrices(vy(p)-dmu, dsigma2_t, ddf, idistribution, ilinkfunction);
		vscore(p) = dscore;
		dscale = Scaling(dinvfisher, iscalingchoice);
		vscaledsc(p) = dscale*dscore;
	end
end


function [dscore, dinvfisher, dllik] = SystemMatrices(dy, dsigma2_t, ddf, idistribution, ilinkfunction)
	global GAUSS STUD_T LOG_SIGMA;
	if idistribution == GAUSS
		dscore = (dy^2-dsigma2_t)/(2*dsigma2_t^2);
		dinvfisher = 2*dsigma2_t^2;
		dllik = -log(dsigma2_t)/2 - dy^2/(2*dsigma2_t);
	elseif idistribution == STUD_T
		dw = (1 + (dy^2/((ddf-2)*dsigma2_t)))^-1*((ddf+1)/(ddf-2));	
		dscore = dw*dy^2/(2*dsigma2_t^2) - 1/(2*dsigma2_t);
		dinvfisher = (2*dsigma2_t^2*(ddf+3))/ddf;
		dllik = log(gamma((ddf+1)/2)) - log(gamma(ddf/2)) - log((ddf-2)*pi)/2 - log(dsigma2_t)/2 - ((ddf+1)/2)*log(1+(dy^2/((ddf-2)*dsigma2_t)));
	else 
		error('Specify correct distribution type');
	end    
	if ilinkfunction == LOG_SIGMA 
		dscore = dscore*dsigma2_t; % Apply chain rule
		dinvfisher = dinvfisher/dsigma2_t^2;
	end
end


function dscale = Scaling(dinvfisher, iscalingchoice)
	global INV_FISHER INV_SQRT_FISHER;
	if iscalingchoice == INV_FISHER
		dscale = dinvfisher;
	elseif iscalingchoice == INV_SQRT_FISHER
		dscale = sqrt(dinvfisher);
	else 
		error('Specify correct scaling type');
	end
end
