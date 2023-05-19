function [negloglikeSum,negloglikeMean,sim] = bciBasicContinuous(params,data,selection)
% This function calculates the negloglikelihood of the basic Bayesian
% Causal Inference Model
% INPUT:
%   params - [pcommon, sig1, sig2, sigP (,muP)]
%   data - a struct, including simulation times N, location range space,
%   and experiment data trials. Trials should contain at least four
%   columns: stim1, stim2, re1, re2
%   selection - the decision strategy, -1 means model selection, 0 means
%   model averaging, 1 means probability matching
% OUTPUT:
%   negloglikelihood
%   sim



rng('shuffle')

%Basic parameters
p_common = params(1);
sig1 = params(2); var1 = sig1^2;
sig2 = params(3); var2 = sig2^2;
sigP = params(4); varP = sigP^2;
if length(params) < 5
    muP = 0;
else
    muP = params(5);
end

patch = -1000;
stim_pairs = [data.trials.stim1,data.trials.stim2];
stim_pairs(isnan(stim_pairs)) = patch;
stim_pairs = unique(stim_pairs,'rows');
stim_pairs(stim_pairs==patch) = nan;
stim1 = stim_pairs(:,1);
stim2 = stim_pairs(:,2);
N_stim_pairs = size(stim_pairs,1);
N_stim1 = sum(isnan(stim2));
N_stim2 = sum(isnan(stim1));
N_bi = N_stim_pairs - N_stim1 - N_stim2;

%Generation of different conditions/stim locations
var12_hat = 1/(1/var1 + 1/var2 + 1/varP);
var1_hat = 1/(1/var1 + 1/varP);
var2_hat = 1/(1/var2 + 1/varP);

var_common = var1 * var2 + var1 * varP + var2 * varP;
var1_indep = var1 + varP;
var2_indep = var2 + varP;

x1 = stim1 + sig1*randn(N_stim_pairs,data.N);
x2 = stim2 + sig2*randn(N_stim_pairs,data.N);


% Bimodal estimates of perceived location of causes
biInds = (~isnan(stim1) & ~isnan(stim2));
s_hat_common = ((x1(biInds,:)/var1) + (x2(biInds,:)/var2) + repmat(muP,sum(biInds),data.N)/varP) * var12_hat;
s1_hat_indep = ((x1(biInds,:)/var1) + repmat(muP,sum(biInds),data.N)/varP) * var1_hat;
s2_hat_indep = ((x2(biInds,:)/var2) + repmat(muP,sum(biInds),data.N)/varP) * var2_hat;

% Unimodal Estimates - organize them according to their locations
uniInds1 = isnan(stim2);
s1_hat_uni = ((x1(uniInds1,:)/var1) + repmat(muP/varP,sum(uniInds1),data.N)) * var1_hat;
uniInds2 = isnan(stim1);
s2_hat_uni = ((x2(uniInds2,:)/var2) + repmat(muP/varP,sum(uniInds2),data.N)) * var2_hat;

%Distance between stim locations
quad_common = (x1(biInds,:)-x2(biInds,:)).^2 * varP + (x1(biInds,:)-muP).^2 * var2 + (x2(biInds,:)-muP).^2 * var1;
quad1_indep = (x1(biInds,:)-muP).^2;
quad2_indep = (x2(biInds,:)-muP).^2;

%Likelihood calculations and
likelihood_common = exp(-quad_common/(2*var_common))/(2*pi*sqrt(var_common));
likelihood1_indep = exp(-quad1_indep/(2*var1_indep))/sqrt(2*pi*var1_indep);
likelihood2_indep = exp(-quad2_indep/(2*var2_indep))/sqrt(2*pi*var2_indep);
likelihood_indep =  likelihood1_indep .* likelihood2_indep;
post_common = likelihood_common * p_common;
post_indep = likelihood_indep * (1-p_common);
pC = post_common./(post_common + post_indep);

%Decision Strategies for bimodal trials
%Model Selection
if selection=='ms'
    s1_hat_bi = (pC>0.5).*s_hat_common + (pC<=0.5).*s1_hat_indep;
    s2_hat_bi = (pC>0.5).*s_hat_common + (pC<=0.5).*s2_hat_indep;

elseif selection=='ma'
    s1_hat_bi = (pC).*s_hat_common + (1-pC).*s1_hat_indep;
    s2_hat_bi = (pC).*s_hat_common + (1-pC).*s2_hat_indep;
    %Model Averaging
elseif selection=='pm' %Model Matching
    match = 1 - rand(sum(biInds),data.N);
    s1_hat_bi = (pC>match).*s_hat_common + (pC<=match).*s1_hat_indep;
    s2_hat_bi = (pC>match).*s_hat_common + (pC<=match).*s2_hat_indep;
end


%Prediction of location estimates
uniInds1_sim = (repmat((find(uniInds1)*data.N)-data.N,1,data.N) + repmat((1:data.N),N_stim1,1))';
uniInds2_sim = (repmat((find(uniInds2)*data.N)-data.N,1,data.N) + repmat((1:data.N),sum(uniInds2),1))';
biInds_sim = (repmat((find(biInds)*data.N)-data.N,1,data.N) + repmat((1:data.N),sum(biInds),1))';
uniInds1_sim = uniInds1_sim(:);
uniInds2_sim = uniInds2_sim(:);
biInds_sim = biInds_sim(:);


sim = nan(N_stim_pairs*data.N,4);
sim(uniInds1_sim,1) = sortrows(repmat((1:N_stim1)',data.N,1));%cat(1,ones(data.N,1),2*ones(data.N,1),3*ones(data.N,1),4*ones(data.N,1),5*ones(data.N,1));
sim(uniInds1_sim,3) = reshape(s1_hat_uni',[],1);%s1_hat_uni(:);
sim(uniInds2_sim,2) = sortrows(repmat((1:N_stim2)',data.N,1));%cat(1,ones(data.N,1),2*ones(data.N,1),3*ones(data.N,1),4*ones(data.N,1),5*ones(data.N,1));
sim(uniInds2_sim,4) = reshape(s2_hat_uni',[],1);%s2_hat_uni(:);
sim(biInds_sim,1) = sortrows(repmat((1:N_bi)',data.N,1));%cat(1,ones(data.N*5,1),2*ones(data.N*5,1),3*ones(data.N*5,1),4*ones(data.N*5,1),5*ones(data.N*5,1));
sim(biInds_sim,3) = reshape(s1_hat_bi',[],1);%s1_hat_bi(:);
sim(biInds_sim,2) = sortrows(repmat((1:N_bi)',data.N,1));%repmat(cat(1,ones(data.N,1),2*ones(data.N,1),3*ones(data.N,1),4*ones(data.N,1),5*ones(data.N,1)),5,1);
sim(biInds_sim,4) = reshape(s2_hat_bi',[],1);%s2_hat_bi(:);

% %% bimodal prediction
h = hist(s1_hat_bi', data.space);
freq_pred1_bi = bsxfun(@rdivide,h,sum(h)); % size:N_space * N_bi, each column refers to a function


% %% unimodal prediction
% h = hist(s1_hat_uni', data.space);
% freq_pred1_uni = bsxfun(@rdivide,h,sum(h));



% sub = sortrows(data.trials(isnan(data.trials.stim1),["stim2","re2"]));
% s = spline(data.space,freq_pred2_uni',sub.re2);
% idx = sortrows(repmat(eye(N_stim1),data.N_trials(2),1),'descend')';
% ss = s(logical(idx));
% ss(ss<1e-5) = 1e-5;
% spline2_uni = sum(sum(log(ss),'omitnan'),'omitnan');


sub = sortrows(data.trials(~isnan(data.trials.stim1) & ~isnan(data.trials.stim2),["stim1","stim2","re1"]));
s = spline(data.space,freq_pred1_bi',sub.re1);
idx = gen_idx(sub);
ss = s(logical(idx));
ss(ss<1e-5) = 1e-5;
spline1_biSum = sum(sum(log(ss),'omitnan'),"omitnan");
spline1_biMean = log(mean(mean(ss,'omitnan'),"omitnan"));

% negloglike = -(spline1_bi+spline1_uni);
negloglikeSum = -(spline1_biSum);
negloglikeMean = -(spline1_biMean);
if p_common < 0 || p_common > 1 || sig1 < 0 || sig2 < 0 || sigP < 0
    negloglikeSum = Inf;
end
end
% --------------------------------------------------------------------

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