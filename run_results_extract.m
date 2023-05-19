
path=pwd;

decisionLabels = {'ms','ma','pm'};
nmin=15;
Labels={'p1','p2','p3','p4','sigA','sigV','sigP','negloglike','id','Ntrials','strategy'}; % Label和refine文件的label对齐

refineList=[];
summaryList=[];

for subid=2:29
    for decision=1:3
    load([path '\refineResult\' num2str(subid) '_Refined_' decisionLabels{decision} '.mat']);
    id=bcidata.trials.ID(1);
    
    refineResult=struct2table(subj_results);
    refineResult=refineResult(:,5:end);
    refineResult.id=id*ones(nmin,1);
    refineResult.N=size(bcidata.trials,1)*ones(nmin,1);
    refineResult.strategy=repmat(decisionLabels{decision},nmin,1);
    refineResult=splitvars(refineResult,"bads_parameters");
    refineResult.Properties.VariableNames=Labels;

    negLogNull = bciNullModel(bcidata);
    refineResult.nullLog=negLogNull*ones(nmin,1);

    refineList=[refineList;refineResult];
    
    idx=find(refineResult.negloglike==min(refineResult.negloglike));
    summaryList=[summaryList;refineResult(idx,:)];

    clear bcidata subj_results
    end
end

k=7; % 模型中的自由参数个数
summaryList.bic = -summaryList.negloglike-0.5*k*log(summaryList.Ntrials);
summaryList.aic = 2*k + 2*summaryList.negloglike;

writetable(summaryList,"bciSummary.xlsx");