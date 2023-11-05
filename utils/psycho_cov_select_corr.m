function [trialinfo,explVar] = psycho_cov_select_corr(trialinfo,covIdx,varargin)

% fields = trialinfo(:,4:end-2).Properties.VariableNames;
% tmp = trialinfo{:,4:end-2};
% fields = fieldnames(rmfield(imageFeatures,'name'));
% tmp = struct2array(rmfield(imageFeatures,'name'));
if ~isempty(varargin) && strcmp(varargin{1},'visual check')
    visualCheck = varargin{2};
else
    visualCheck = 1;
end
    
%%
fields = trialinfo(:,covIdx).Properties.VariableNames;
tmp = trialinfo{:,covIdx};

to_del = [];
for i = 1:length(tmp)
% for i = 1:size(tmp,2)
    if any(isnan(tmp(i,:)))
        to_del = [to_del i];
    end
end
tmp(to_del,:) = NaN;
covariates = [];
idx = cell2mat(cellfun(@(x) ~contains(x,'primer'),fields,'UniformOutput',false));
% for i = 2:2:length(fields)
for i = find(idx)
    tmp2 = tmp(:,i);
    covariates = [covariates tmp2(~isnan(tmp2(:,1)))];
%     covariates = [covariates tmp(:,i)];
end

if visualCheck
    % [c1,p1] = corr(covariates,covariates,'Tail','left')
    r = corr(covariates(~isnan(covariates(:,1)),:));
    % r(5,7) = -0.4617;
    % r(7,5) = r(5,7);
    % r(6,7) = 0.4158;
    % r(7,6) = r(6,7);
    % r(abs(r)<0.5) = 0;
    figure();imagesc(r,[-1,1]);colorbar
    cmap = interp1([-1;0;1],[1 1 0;0 0 1;1 1 0],linspace(-1,1)');
    colormap(cmap)
    % labels = strrep(fields(1:2:end),'_','-');
    % labels = strrep(fields,'_','-');
    % for i = 1:length(labels)
    %     labels{i}(1:7) = [];
    % end

    % labels = {'Entropy', 'Contrast', 'Corr.', 'Energy', 'Homog.', 'Compact.', 'Ratio',...
    %     '#Spect.\newlineClust.', 'Freq.\newlineEnergy', 'MaxFreq', 'Max.\newlineFreq.\newlineDist'};
    labels = {'#Phon.', 'Lexical.\newlineFreq.', 'Movie.\newlineFreq.', 'AoA', 'Visual.\newlineComplex.', 'Familiar.',...
        'Imageab.'};

    xticks(1:length(labels))
    xticklabels(labels)
    yticklabels(labels)
    % title('Covariates correlation')
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;

    toMerge = input('Which covariates to merge? (e.g. {[1 2]; [5 6 7]})');
else
    toMerge = {[2 3];[5 6]};
end

explVar = [];
newRegressTarget = NaN(size(trialinfo,1),size(toMerge,1));
for j = 1:size(toMerge,1)
%     cov = covariates(:,toMerge{j});
%     [coef, score, ~,~,explained,~] = pca(zscore(cov(~isnan(cov(:,1)),:)));
    [coef, score, ~,~,explained,~] = pca(zscore(covariates(:,toMerge{j})));
    newRegressTarget(~isnan(tmp(:,1)),j) = score(:,1);
    explVar(j) = explained(1);
    compWeights{j} = coef;
end
% covariates(:,unique([toMerge{:}])) = [];

covariates = [];
idx = cell2mat(cellfun(@(x) contains(x,'primer'),fields,'UniformOutput',false));
% for i = 2:2:length(fields)
for i = find(idx)
    tmp2 = tmp(:,i);
    covariates = [covariates tmp2(~isnan(tmp2(:,1)))];
end
newRegressPrimer= NaN(size(trialinfo,1),size(toMerge,1));
for j = 1:size(toMerge,1)
    [coef, score, ~,~,explained,~] = pca(zscore(covariates(:,toMerge{j})));
    newRegressPrimer(~isnan(tmp(:,1)),j) = score(:,1);
    explVar(end+1) = explained(1);
    compWeights{j} = coef;
end

trialinfo = table2struct(trialinfo);
% trialinfo(to_del) = [];
trialinfo = rmfield(trialinfo,{'primer_lexical_freq','target_lexical_freq',...
    'primer_visual_freq','target_visual_freq','primer_visual_complexity',...
    'target_visual_complexity','primer_familiarity','target_familiarity'});
trialinfo = [struct2table(trialinfo)...
    table(newRegressPrimer(:,1), 'VariableNames', {'primer_psycho_freq'}),...
    table(newRegressTarget(:,1), 'VariableNames', {'target_psycho_freq'}),...
    table(newRegressPrimer(:,2), 'VariableNames', {'primer_familiarity'}),...
    table(newRegressTarget(:,2), 'VariableNames', {'target_familiarity'})];

% trialinfo = rmfield(trialinfo,{'name','entropyValue','contrastValue','energyValue','homogeneityValue','numberOfCluster','maxClusterFreq','maxClusterDist'});
% trialinfo = table2struct([struct2table(trialinfo) table(contrastPCA, 'VariableNames', {'contrastPCA'}) table(freqPCA, 'VariableNames', {'freqPCA'})]);

