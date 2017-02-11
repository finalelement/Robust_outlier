function JS = calc_JSD(P,Q)
% JSD measures distance betwween two probability distributions. In this
% case, the probability distribution is the ODF defined by discrete
% sampling on the sphere. 
%           JSD(P,Q) = DKL(P,M)+DKL(Q,M)/2
% where M(i)=P(i)+Q(i)/2
% where P(i) and Q(i) are distributions of the discrete ODF reqconstructed
% from HARDI data, i is the index of sampling of ODF (i=1...724) and
% DKL(X,Y) is the kullback leibler divergence
%           DKL(P,Q) = sumi(P(i)log(P(i)/Q(i))
% JSD is symettric

% first need to make sure they are each probability distributions (sum to
% 1)
P=P/sum(P);
Q=Q/sum(Q);

% calc M
M = (P+Q)/2; % same size as P and Q

% calc DKL(P,M)
DKLPM = sum(P.*log(P./M));
% calc DKL(Q,M)
DKLQM = sum(Q.*log(Q./M));
% calc JS
JS = (DKLPM+DKLQM)/2;









