function [r] = var_selection(config)
% This function creates a trialinfo table with all the information needed
% about the studied confounders
% /!\/!\/!\
% This function has to be modified if you want to study other confounders
% /!\/!\/!\

%% Define the raw structure containing information about the covariates in the field "trialinfo"
dinfo = dir(fullfile(config.BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    dinfo = dir(fullfile(config.BIDS_FOLDER,subj_name,'eeg',['*',config.ext]));
    EEG_FILE = fullfile(config.BIDS_FOLDER, subj_name, 'eeg', dinfo.name);
    
    cfg = [];
    cfg.dataset = EEG_FILE;
    raw_eeg = ft_preprocessing(cfg);
    
    cfg                     = [];
    cfg.dataset             = EEG_FILE;
    cfg.root                = config.BIDS_FOLDER;
    cfg.trialdef.eventtype  = config.trialdef_eventtype;
    cfg.trialfun            = config.trial_function; % define the extraction of covariates information (cf. utils/SenSem_trialfun_trial.m)
    cfg.cov_description     = config.cov_description;
    cfg.trialdef.eventvalue = config.trialdef_eventvalue;
    cfg.trialdef.prestim    = config.trialdef_prestim;
    cfg.trialdef.poststim   = config.trialdef_poststim; %Note: the minimal time between the target appearance and the next primer is 2s
    cfg = ft_definetrial(cfg);
    eeg = ft_redefinetrial(cfg,raw_eeg);
    
    if config.save_choice
        save(fullfile(config.PATH_TO_DERIV, subj_name, 'eeg', sprintf('%s_task-%s_raw.mat',subj_name,config.task_name)))
    end
end

%% Extract the covariates
imageFeatures = extract_cov(config.PATH_TO_DERIV,config.PATH_TO_IMAGES,config.task_name); % this process takes few minutes
if config.save_choice
    save(fullfile(config.PATH_TO_DERIV,'image_features.mat'),'imageFeatures')
end

visualSimilarity = extract_visual_sim(config.PATH_TO_SIMILARITY, config.PATH_TO_ITEMS);

tmp = struct2array(rmfield(imageFeatures,'name'));
correl = corrcoef(tmp);
correl(correl<0.55 & correl>-0.55) = 0;
figure;imagesc(correl)

%% compute PCA for homogeneity-contrast related features
tmp = [zscore(imageFeatures.entropyValue) zscore(imageFeatures.contrastValue) zscore(imageFeatures.energyValue)...
    zscore(imageFeatures.homogeneityValue)];
[coef, score, ~,~,explained,~] = pca(tmp);
contrastPCA = score(:,1);
fprintf('contrast PCA coefficients = %d with corresponding explained variances = %d\n', coef, explained)

%% compute PCA for frequency related features
tmp = [zscore(imageFeatures.numberOfCluster) zscore(imageFeatures.maxClusterFreq) zscore(imageFeatures.maxClusterDist)];
[coef, score, ~,~,explained,~] = pca(tmp);
freqPCA = score(:,1);
fprintf('frequency PCA coefficients = %d with corresponding explained variances = %d\n', coef, explained)

%% Create the trialinfo structure
imageFeaturesReduced = rmfield(imageFeatures,{'name','entropyValue','contrastValue','energyValue','homogeneityValue','numberOfCluster','maxClusterFreq','maxClusterDist'});
imageFeaturesReduced.contrastPCA = contrastPCA;
imageFeaturesReduced.freqPCA = freqPCA;
f = fieldnames(imageFeaturesReduced);

dinfo = dir(fullfile(config.BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    trialinfo = load(fullfile(config.BIDS_FOLDER,'derivatives',subj_name, 'eeg','trialinfo_psycho.mat'));
    trialinfo = trialinfo.(cell2mat(fieldnames(trialinfo)));
    for j = 1:size(trialinfo,1)
        idxTarget = cell2mat(cellfun(@(x) strcmp(x(1:end-4),trialinfo(j).target{1}),cellstr(char(imageFeatures.name)), 'UniformOutput',false));
        idxPrimer= cell2mat(cellfun(@(x) strcmp(x(1:end-4),trialinfo(j).primer{1}),cellstr(char(imageFeatures.name)), 'UniformOutput',false));
        for k = 1:length(f)
            try
                trialinfo(j).([f{k} '_primer']) = imageFeaturesReduced.(f{k})(idxPrimer);
                trialinfo(j).([f{k} '_target']) = imageFeaturesReduced.(f{k})(idxTarget);
            catch
                trialinfo(j).([f{k} '_primer']) = NaN;
                trialinfo(j).([f{k} '_target']) = NaN;
            end
            if isempty(trialinfo(j).([f{k} '_primer']))
                trialinfo(j).([f{k} '_primer']) = NaN;
            end
            if isempty(trialinfo(j).([f{k} '_target']))
                trialinfo(j).([f{k} '_target']) = NaN;
            end
                
        end
        idx = cell2mat(cellfun(@(x,y) strcmp(x,trialinfo(j).target{1}) & strcmp(y,trialinfo(j).primer{1}),...
            cellstr(char(visualSimilarity.target)),cellstr(char(visualSimilarity.primer)), 'UniformOutput',false));
        try
            trialinfo(j).visualSimilarity = visualSimilarity.similarity{idx};
        catch
            trialinfo(j).visualSimilarity = NaN;
        end
    end
    
    trialinfo = struct2table(trialinfo);
    if config.save_choice
        trialinfo_filename = 'trialinfo_psycho_image.mat';
        save(fullfile(config.BIDS_FOLDER,'derivatives',subj_name, 'eeg',trialinfo_filename),'trialinfo');
    end
end

%% Total correlation analysis
covIdx = 4:26;
fields = trialinfo(:,covIdx).Properties.VariableNames;
% tmp = struct();
% for varName = fields
%     tmp.(varName{1}) = trialinfo.(varName{1});
% end
tmp = table2array(trialinfo(:,covIdx));

to_del = [];
for i = 1:length(tmp)
    if any(isnan(tmp(i,:)))
        to_del = [to_del i];
    end
end
tmp(to_del,:) = [];
covariates = [];
idx = cell2mat(cellfun(@(x) ~contains(x,'primer'),fields,'UniformOutput',false));
% for i = 2:2:length(fields)
for i = find(idx)
    tmp2 = tmp(:,i);
    covariates = [covariates tmp2(~isnan(tmp2(:,1)))];
end
r = corr(covariates);
figure();imagesc(r,[-1,1]);colorbar
cmap = interp1([-1;0;1],[1 1 0;0 0 1;1 1 0],linspace(-1,1)');
colormap(cmap)
% labels = {'Entropy', 'Contrast', 'Corr.', 'Energy', 'Homog.', 'Compact.', 'Ratio',...
%     '#Spect.\newlineClust.', 'Freq.\newlineEnergy', 'MaxFreq', 'Max.\newlineFreq.\newlineDist'};
labels = {'#Phon.', 'AoA', 'Imageab.','Psycho.\newlineFreq.', 'Familiar.',...
    'Corr.', 'Compact.', 'Ratio', 'Freq.\newlineEn.', 'Contrast', 'Image\newlineFreq.','Visual\newlineSim.'};

xticks(1:length(labels))
xticklabels(labels)
yticks(1:length(labels))
yticklabels(labels)
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;

%% Model with image features only
dinfo = dir(fullfile(config.BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    trialinfo = load(fullfile(BIDS_FOLDER,'derivatives',subj_name, 'eeg',trialinfo_filename));
    trialinfo = trialinfo.(cell2mat(fieldnames(trialinfo)));
    f = trialinfo.Properties.VariableNames;
    for k = 4:15 % it has to be adapted following your model
        trialinfo = removevars(trialinfo,f(k));
    end
    if config.save_choice
        trialinfo_filename = 'trialinfo_image.mat';
%         trialinfo_filename = 'new_image_trialinfo.mat';
        save(fullfile(config.BIDS_FOLDER,'derivatives',subj_name, 'eeg',trialinfo_filename),'trialinfo');
    end
end

%% Show covariates difference between categories
% trialinfo.condition(ismember(trialinfo.condition,[1 3 5])) = 1;
% trialinfo.condition(ismember(trialinfo.condition,[2 4 6])) = 2;

x1 = find(trialinfo.condition==1);
x2 = find(trialinfo.condition==2);
g = [zeros(length(x1), 1); ones(length(x2), 1); 2;...
     3*ones(length(x1), 1); 4*ones(length(x2), 1); 5;...
     6*ones(length(x1), 1); 7*ones(length(x2), 1); 8;...
     9*ones(length(x1), 1); 10*ones(length(x2), 1); 11;...
     12*ones(length(x1), 1); 13*ones(length(x2), 1); 14];

fields = {'target_phoneme_nb','target_AoA','target_imageability','target_psycho_freq','target_familiarity'};
boxMat = [];
for f = fields
    f = f{1};
    meanInfo = nanmean(trialinfo.(f));
    stdInfo = nanstd(trialinfo.(f));
    zInfo = (trialinfo.(f) - meanInfo) ./ stdInfo;
%     zInfo = zscore(trialinfo.(f));
    boxMat = [boxMat;zInfo(trialinfo.condition==1);zInfo(trialinfo.condition==2);nan];
end
figure;
boxplot(boxMat,g)
title('confounder pairwise comparison natural vs manufactured')
ylabel('zscore')

% Initialize the new cell array
newFields = cell(1, length(fields) * 3);
% Fill the new cell array
for i = 1:length(fields)
    newFields{(i-1)*3 + 1} = fields{i};
    newFields{(i-1)*3 + 2} = fields{i};
    newFields{(i-1)*3 + 3} = "";
end
xticklabels(newFields); % Apply the newFields labels
% Rotate the x-tick labels by 90 degrees
xtickangle(90);
