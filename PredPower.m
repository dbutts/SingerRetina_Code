function [pp,r2] = PredPower( data, rt, dt )
%
% Usage: [pp,r2] = PredPower( data, rt, dt )
%
% data is repeated spike trails
% rt is predicted [model] probability per bin
% dt is resolution of pred, default to data already in units of bins

% Leave out how many initial bins
BUFF = 30;
% Leave out how many repeats
REP_OMIT = 2;

Nreps = repeat_detect(data);
NT = length(rt);
T = NT*dt;

Trials.D = zeros( Nreps-REP_OMIT, NT-BUFF );
for n = 1:(Nreps-REP_OMIT)
  spks = clip_repeats(data,n+REP_OMIT);
  h = histc(spks,((1:NT)-1)*dt);
  Trials.D(n,:) = h((BUFF+1):NT);
end
Trials.SR = 1;
Pred.SR = 1;

Pred.D = rt((BUFF+1):NT);
Pred.D = Pred.D(:)';
[r2,pp] = assessPerformance( Pred, Trials, [],2 );
