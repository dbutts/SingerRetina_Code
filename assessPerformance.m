function [VarExplained,Beta,VarSignal,VarNoise,VarVarSignal,VarExplainedCorr] ...
  = assessPerformance(Pred,Trials,P,compVar,Sigma)
% assessPerformance computes the fraction of explained variance and the beta error
% for a prediction Pred and the Repetitions of the experiment in Trials.
% 
% The Beta error makes only little assumptions (see Sahani & Linden 2003a):
% Assumed: additive noise, possibly with time or level dependent distribution

%global U; if isempty(U) U = units; end

if ~exist('compVar','var') compVar = 0; end
if isfield(P,'Tselect') & ~isempty(P.Tselect) 
  Pred = HF_selectRange(Pred,P.Tselect); 
  Trials = HF_selectRange(Trials,P.Tselect);
end
SR = Pred.SR; Pred = Pred.D; Trials = Trials.D; 
Beta = NaN; VarSignal = NaN; VarNoise = NaN; VarVarSignal = NaN;

NTrials = size(Trials,1); %if NTrials>1 fprintf([n2s(NTrials),' Trials.\n']); end
NSteps = length(Pred);
Resp = mean(Trials,1); % compute average Rate

% Compute the Coefficient of Determination / Fraction of explained Variance 
% on the level of the average, i.e.
% - compares the prediction with the average trial
% - if Pred comes from a linear regression estimate it should also hold that: 
%    VarExplained = VarPred/VarResp (see Coefficient of Determination)

VarError = var(Resp - Pred);
VarResp = var(Resp);
VarExplained = 1-VarError/VarResp;

% BETA is first defined in Machens2004, but had been introduced in Sahani2003
% The two versions are not exactly equivalent, since Machens uses the 
% in Machens:       B = ( <P(r^n)> - <P(e^n)> ) / P(S)   
% whereas in Sahani B = ( P(<r^n>) - P(<e^n>) ) / P(S)
% The terms in the numerator are the same in expectation 
% and Sahani's version should converge faster
% Also: Sahani's Version can be computed more efficiently, 
% since the trial errors do not have to be computed.
%
% Compute the Fraction of explained Variance for individual trials (beta-error)
% - takes unreliability/noise into account
if NTrials > 1
  VarErrorTrial = zeros(NTrials,1); mPred = mean(Pred);
  mTrials = sum(Trials,2)/NSteps; 
  if Trials.*Trials==Trials % Bernoulli Case
    for i=1:NTrials
      VarErrorTrial(i) = sum((Trials(i,:) - Pred - mTrials(i) + mPred).^2)./(NSteps-1);
      %VarRespTrial(i)  = sum((Trials(i,:) - mTrials(i)).^2)./(NSteps-1);
    end
    VarRespTrial = mTrials.*(1-mTrials)*(NSteps/(NSteps-1)); % Using the upper equation instead of the this one gives the unbiased variance estimator
  else
    VarRespTrial = zeros(NTrials,1);
    for i=1:NTrials % Efficient computation of variances
      VarErrorTrial(i) = sum((Trials(i,:) - Pred - mTrials(i) + mPred).^2)./(NSteps-1);
      VarRespTrial(i)  = sum((Trials(i,:) - mTrials(i)).^2)./(NSteps-1);
    end
  end
  
  % VarErrorInd captures the average Error in the prediction of single trials
  VarErrorInd = mean(VarErrorTrial);
  
  % VarRespInd captures the average Variance within a trial
  VarRespInd = mean(VarRespTrial);
  
  % SIGNAL POWER ESTIMATE
  % i.e. Noise, which cannot be captured given always the same stimulus
  VarSignal = 1/(NTrials-1) * (NTrials*VarResp - VarRespInd); % Sahani2003a
  % VarResidual captures trial-to-trial variation, 
  %VarResidual = NTrials/(NTrials-1) * (VarRespInd - VarResp); % Machens2004
  %VarSignal = VarRespInd - VarResidual;
  VarNoise = VarRespInd - VarSignal;

  % BETA
  % Enumerator of Beta: Explained Variance / Power
  % Denominator of Beta: Explainable Variance / Power
  Beta(1,1) = (VarResp - VarError)/VarSignal; % Sahani2003a
  %Beta(1,2) = (VarRespInd - VarErrorInd)/VarSignal; % Machens2004
  
  
  VarExplainedCorr=[];
  
  %{
  % HAEFNER MEASURE: Corrected Explained Variance
  NSigma = NSteps*(NTrials-1);
  Sigma2 = 1/(NTrials*NSteps*(NTrials-1))*sum(sum((Trials - repmat(Resp,NTrials,1)).^2));
  LDD = (NSigma-2)/NSigma*sum((Resp-mean(Resp)).^2)/Sigma2 - (NTrials-1);
  LDM = (NSigma-2)/NSigma*sum((Resp-Pred).^2)/Sigma2 - (NTrials-P.NParams);
  VarExplainedCorr = 1 - LDD/LDM;
  %}
  
  % VARIANCE OF SIGNAL POWER (Sahani2003a)
  if exist('Sigma','var') compVar = 3; end
  if compVar>0
    switch compVar
      case 1; % Compute the reduced variance over only 5 ms
        Noise = Trials - repmat(Resp,NTrials,1);
        Noise = Noise - repmat(mean(Noise,2),1,NSteps);
        Sigma = sparse(NSteps,NSteps);
        MaxCorrDist = round(5*U.ms*SR);
        for i=1:NSteps
          iStop = min([i+MaxCorrDist,NSteps]);
          Sigma(i,i:iStop) = Noise(:,i:iStop)'*Noise(:,i);
        end
        Sigma = (Sigma + Sigma')/NTrials;
        Sigma([1:NSteps+1:NSteps^2]) = Sigma([1:NSteps+1:NSteps^2])/2;
        %Sigma = cov(Noise); %Sigma = diag(Resp.*(1-Resp));
      case 2; % Compute the full variance
        Noise = Trials - repmat(Resp,NTrials,1);
        Noise = Noise - repmat(mean(Noise,2),1,NSteps);
        Sigma = cov(Noise);
      case 3; % Use prescribed Variancematrix Sigma
    end
    
    mSigma1 = sum(Sigma,1)/NSteps; mSigma = sum(mSigma1)/NSteps;
    mResp = mean(Resp);
    %trSigma = trace(Sigma*Sigma);  
    trSigma = sum(Sigma(:).^2); % = Sigma(:)'*Sigma(:); (very memory intensive)
    VarVarSignal = ...
      4/NTrials*(1/NSteps^2 * Resp*Sigma*Resp' - 2/NSteps * mResp*mSigma1*Resp' + mResp*mSigma*mResp) ...
      + 2/(NTrials*(NTrials-1))*(1/NSteps^2*trSigma - 2/NSteps*mSigma1*mSigma1' + mSigma.^2);
 
  end
end % Only one Trial, no Variance computable

