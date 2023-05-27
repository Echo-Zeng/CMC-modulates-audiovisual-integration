%% Setting

nmin = 15; % extract simulations with the 15 min negloglikelihood
decisionLabels = {'ms','ma','pm'};

path = pwd;

folder = fullfile(pwd,'refineResult');
if ~exist(folder)
    mkdir(folder);
end


%% Refine Setting

% Options for bads
options = bads('defaults');
options.TolFun = '1e-12';
options.Display = 'final';
options.MaxIter = 4000;

% set upper and lower bounds on parameters
parameterNames = {'p1','p2','p3','p4','sigA','sigV','sigP'}; 
LB = [0 0 0 0 0.01 0.01, 0.01];
UB = [1 1 1 1  30   30    40];

for subid = 2:29

    for decision=1:3
    load([path '\gridSearchResult\' num2str(subid) '_Simulations_gridsearch_acrossconditions_' decisionLabels{decision} '.mat'])

    % Creating anonymous function for input to fminsearch
    fun = @(param)  bciMultiCondition(param,parameterNames,bcidata,decisionLabels{decision});

    %% Finding the n minima of the log-likelihood values and the corresponding parameter combination
    % sort gridsearch log-likelihoods from the 1st best (global minimum)
    logLike_sort = sort(logLike_allSum);

    % get id number and parameter combination for the n best log-likelihoods
    for i = 1:nmin
        idx = find(logLike_allSum==logLike_sort(i));
        if length(idx)>1
            idx_temp=randi(length(idx)); %randomly select one
            idx=idx(idx_temp);
        end
        bestParam = paramCombinations(idx,1:end);
        % save in subj-specific structure
        subj_results(i).idx = idx;
        subj_results(i).bestParam = bestParam;
        subj_results(i).strategy = decisionLabels{decision};
    end
    % get n best log-likelihoods and save in subj-specific structure
    bestlogLike(:,1) = logLike_sort(1:nmin);
    for i=1:nmin
        subj_results(i).bestlogLike = bestlogLike(i);
    end

    %% Starting Refining: time and details
    cStart = clock;

    % Performing fminsearch for each parameter combination
    parfor i = 1:nmin
        fprintf(['Performing bads optimization ' num2str(i) ' of subject ' num2str(subid) '... \n']);
        %         [fmin_parameters,errorval] = fminsearchbnd(fun,subj_results(i).bestParam,LB,UB,opts);

        [bads_parameters,errorval] = bads(fun,subj_results(i).bestParam,LB,UB,[],[],[],options);

        subj_results(i).bads_parameters = bads_parameters;
        subj_results(i).errorval = errorval;
    end


    % Finishing timer and printing elapsed time
    fprintf('badssearch elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));


    %% Saving
    save([path '\refineResult\' num2str(subid) '_Refined_' decisionLabels{decision} '.mat'],"bcidata","subj_results")
    end
end

