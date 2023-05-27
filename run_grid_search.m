
%% Preparations

load("data_total_exp1.mat")

folder = 'gridSearchResult';
path = fullfile(pwd,folder);
if ~exist(path)
    mkdir(path);
end


%% Grid search parameter Setting
% grid resolution
n = 5;
p_common1=linspace(0.1,0.9,n);
p_common2=linspace(0.1,0.9,n);
p_common3=linspace(0.1,0.9,n);
p_common4=linspace(0.1,0.9,n);
sigP = linspace(0.1,40,n);
sigV = linspace(0.1,30,n);
sigA = linspace(0.1,30,n);

parameterNames = {'p1','p2','p3','p4','sigA','sigV','sigP'};
gridVectors = {p_common1,p_common2,p_common3,p_common4,sigA,sigV,sigP,};
nParameters = numel(gridVectors);
% Full factorial expansion of the specified parameter vectors
coords = cell(1,nParameters);
[coords{:}] = ndgrid(gridVectors{:});
coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
% Matrix of all possible parameter combinations: nRows = number of
% combinations, columns correspond to parameters
paramCombinations = cat(2,coords{:});

decisionLabels = {'ms','ma','pm'};

%% Performing Grid Search

for subid = 2:29
    bcidata.trials = data{subid};
    bcidata.trials = bcidata.trials(bcidata.trials.block>0,:);

    % Preprocssing
    dataHeading=bcidata.trials.Properties.VariableNames;
    dataHeading{find(cellfun(@isempty,regexp(dataHeading,'loc_a'))==0)}='stim1';
    dataHeading{find(cellfun(@isempty,regexp(dataHeading,'loc_v'))==0)}='stim2';
    dataHeading{find(cellfun(@isempty,regexp(dataHeading,'response'))==0)}='re1';
    %     dataHeading{find(cellfun(@isempty,regexp(dataHeading,'reV'))==0)}='re2';
    bcidata.trials.Properties.VariableNames=dataHeading;
    bcidata.space=linspace(-18.5,18.5,50);
    bcidata.N_trials=[35,nan,nan] % Each kind of auditory stimulus appeared 35 times at baseline condition and 80 times in mental imagery and physical stimulation conditions
    bcidata.N = 10000;
    bcidata.Nsample = size(bcidata.trials,1)
    clear dataHeading

    % decision loop
    for decision = 1:length(decisionLabels)

        % Pre-allocating array for data collection
        logLike_allSum = NaN(size(paramCombinations,1),1);
        logLike_allMean = logLike_allSum;

        % Starting the timer and printing details
        cStart = clock;
        fprintf('Performing gridsearch... \n');

        %% Performing grid search on the specified parameter space
        parfor in = 1:size(paramCombinations,1)
            [logLike_allSum(in),logLike_allMean(in),~] = bciMultiCondition(paramCombinations(in,:),parameterNames,bcidata,decisionLabels{decision});
        end

        % Printing elapsed time
        cTemp = clock;
        fprintf('Gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
            datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));

        % Saving subject specific bci simulations
        fprintf('\n\nSaving data...\n\n');
        save([path '\' num2str(subid) '_Simulations_gridsearch_acrossconditions_' decisionLabels{decision} '.mat'], ...
            'subid', 'logLike_allSum', 'paramCombinations', 'parameterNames', 'bcidata');

    end % end of decision loop , save grid search results

end
