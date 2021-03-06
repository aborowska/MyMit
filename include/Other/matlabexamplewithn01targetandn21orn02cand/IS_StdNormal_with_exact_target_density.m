
% We are interested in (dP * 100)% percentile:

dP  = 0.05;


% Target distribution: N( 0, 1 ):

dTargetMean     = 0;
dTargetStdev    = 1;


% Candidate distribution: N( -2, 1 ):

dCandidateMean  = -2;
dCandidateStdev = 1;


% Number of draws:

iN = 100000;


% Obtain candidate draws, compute candidate en target density and IS weights:

vCandidateDraws = normrnd( dCandidateMean, dCandidateStdev, iN, 1 );
vCandidatePdf   = normpdf(vCandidateDraws, dCandidateMean, dCandidateStdev );
vTargetPdf      = normpdf(vCandidateDraws, dTargetMean, dTargetStdev );
vISweights      = (1/iN) * vTargetPdf ./ vCandidatePdf;


% Sort draws (and corresponding weights) using "sortrows":

mCandidateDraws_ISweights       = [ vCandidateDraws  vISweights];
mSortedCandidateDraws_ISweights = sortrows( mCandidateDraws_ISweights, 1 );


% Compute cumulative sum of IS weights and

vCumulativeISweights = cumsum( mSortedCandidateDraws_ISweights(:,2) );

iNumberOfCumulativeWeightsSmallerThanP = sum( vCumulativeISweights < dP );

dISestimateOfPercentile = mSortedCandidateDraws_ISweights( iNumberOfCumulativeWeightsSmallerThanP, 1 )

    