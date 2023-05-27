function [logLikelihoodSumSum,logLikelihoodMeanSum,bic,aic] = bciMultiCondition(parameters,parameterNames,data,selection)
% This function calculates the negloglikelihood of a general Bayesian
% Causal Inference Model that allow some parameters to vary across
% conditions. Here, we assume that only Pcommon, sig1(sigA), sig2(sigV),
% muP could change.
% INPUT:
%     params - parameter values in an array
%     parameterNames - Names of paramters. The order needs to match the params.
%     data - a struct, including N, space and trials. Trials should include
%     specific condition indice to allow data selection.
%     selection - the decision strategy
% OUTPUT:
%     neglikelihood
%     bic
%     aic



%% Get parameters
if iscolumn(parameters) parameters = parameters';end


Pnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'p[0-9]*')));
sigAnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigA[0-9]*')));
sigVnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*')));
sigPnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigP[0-9]*')));
muPnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'muP[0-9]*')));


%% initialize something
logLikelihoodSumSum = 0;
logLikelihoodMeanSum = 0;


%% fit across various condition levels

dataHeading=data.trials.Properties.VariableNames;
conditionP=[];conditionA=[];conditionV=[];conditionMu=[];
if length(Pnames)>1
    conditionP = find(cellfun(@isempty,regexp(dataHeading,'conditionP'))==0);
end
if length(sigAnames)>1
    conditionA = find(cellfun(@isempty,regexp(dataHeading,'conditionA'))==0);
end
if length(sigVnames)>1
    conditionV = find(cellfun(@isempty,regexp(dataHeading,'conditionV'))==0);
end
if length(muPnames)>1
    conditionMu = find(cellfun(@isempty,regexp(dataHeading,'conditionMu'))==0);
end

conditions = sortrows(unique(data.trials(:,[conditionP,conditionA,conditionV,conditionMu])));

for iCond = 1:size(conditions,1)
    
    % choose a subset of data
    actData = data;
    actData.trials = actData.trials(ismember(actData.trials(:,[conditionP,conditionA,conditionV,conditionMu]),conditions(iCond,:),'rows'),:);


    % ------------ choose corresponding parameters - P, sigA, sigV, sigP, mu
    if ~isempty(conditionP)
        actP = parameters(ismember(parameterNames,Pnames{conditions(iCond,:).conditionP}));
    else
        actP=parameters(ismember(parameterNames,Pnames));
    end
    if ~isempty(conditionA)
        actSigA = parameters(ismember(parameterNames,sigAnames{conditions(iCond,:).conditionA}));
    else
        actSigA = parameters(ismember(parameterNames,sigAnames));
    end
    if ~isempty(conditionV)
        actSigV = parameters(ismember(parameterNames,sigVnames{conditions(iCond,:).conditionV}));
    else
        actSigV = parameters(ismember(parameterNames,sigVnames));
    end
    actSigP = parameters(ismember(parameterNames,sigPnames));
    if ~isempty(conditionMu)
        actMu = parameters(ismember(parameterNames,muPnames{conditions(iCond,:).conditionMu}));
    else
        actMu=parameters(ismember(parameterNames,muPnames));
    end
    
    actParameters = [actP,actSigA, actSigV,actSigP,actMu]; % [pcommon, sig1, sig2, sigP (,muP)]
    % ------------

    [logLikelihoodSum,~,~] = bciBasicContinuous(actParameters,actData,selection);
    logLikelihoodSumSum = logLikelihoodSumSum + logLikelihoodSum; % negative

end

%% output



k = numel(parameterNames);
bic = -logLikelihoodSum-0.5*k*log(data.Nsample);
aic = 2*k + 2*logLikelihoodSum;


end
