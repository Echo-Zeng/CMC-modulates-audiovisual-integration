function [negLoglikeSum,negLoglikeMean] = bciNullModel(data)

% This function calculates negloglikelihood for a null model which assumes the estimates x orignates from a uniform distribution

rng('shuffle')

s1_min = min(data.space);
s1_max = max(data.space);
s1_interval = s1_max - s1_min;

patch = -1000;
stim_pairs = [data.trials.stim1,data.trials.stim2]; 
stim_pairs(isnan(stim_pairs)) = patch;
stim_pairs = unique(stim_pairs,'rows');
stim_pairs(stim_pairs==patch) = nan;
stim1 = stim_pairs(:,1);
stim2 = stim_pairs(:,2);
N_stim_pairs = height(stim_pairs);
N_stim1 = length(unique(stim1));
N_stim2 = length(unique(stim2));
N_bi = N_stim_pairs - N_stim1 - N_stim2;

s1_hat = s1_interval * rand(N_stim1,data.N) + s1_min;

h = hist(s1_hat, data.space);
freq_pred1 = bsxfun(@rdivide,h,sum(h)); % size:N_space * N_bi, each column refers to a function

% negloglikelihood for the auditory estimate for the unisensory condition
% sub = sortrows(data.trials(isnan(data.trials.stim2),["stim1","re1"]));
% s = spline(data.space,freq_pred1',sub.re1);
% idx = sortrows(repmat(eye(N_stim1),data.N_trials(1),1),'descend')';
% ss = s(logical(idx));
% ss(ss<1e-5) = 1e-5;
% spline1_uni = sum(sum(log(ss),'omitnan'),'omitnan');

% negloglikelihood for the auditory estimate for the bimodal condition
sub = sortrows(data.trials(~isnan(data.trials.stim1) & ~isnan(data.trials.stim2),["stim1","stim2","re1"]));
s = spline(data.space,freq_pred1',sub.re1);
idx = gen_idx(sub);
ss = s(logical(idx));
ss(ss<1e-5) = 1e-5;
spline1_biSum = sum(sum(log(ss),'omitnan'),"omitnan");

negLoglikeSum = -(spline1_biSum);
end

%%
function idx = gen_idx(tbl)


if sum(find(strcmp(tbl.Properties.VariableNames,'cond')))==0
    tbl.stim_idx = tbl.stim1*100 + tbl.stim2;
    ft = tabulate(tbl.stim_idx);
    ft(:,end+1) = cumsum(ft(:,2)); % end point
    ft(2:end,end+1) = ft(1:end-1,end)+1;
    ft(1,end) = 1;
    idx_mat = cell2mat(arrayfun(@(x) [repmat(x,ft(x,2),1),(ft(x,end):ft(x,end-1))',ones(ft(x,2),1)],1:size(ft,1),'UniformOutput',false)');
    idx = sparse(idx_mat(:,1),idx_mat(:,2),idx_mat(:,3));
else
    cond_list = unique(tbl.cond);
    N_cond = length(cond_list);
    idx_list = cell(N_cond,1);
    size_list = nan(N_cond,2);
    for i = 1:length(cond_list)
        tbl_sub = tbl(tbl.cond==cond_list(i),:);
        idx_list{i} = gen_idx(removevars(tbl_sub,"cond"));
        size_list(i,:) = size(idx_list{i});
    end

    idx = cell(N_cond,1);
    for i = 1:length(cond_list)
        temp = zeros(size_list(i,1),sum(size_list(:,2)));
        temp(:,sum(size_list(1:i-1,2))+1:sum(size_list(1:i,2))) = logical(idx_list{i});
        idx{i} = temp;
    end
    idx = cell2mat(idx);
end
return
end
