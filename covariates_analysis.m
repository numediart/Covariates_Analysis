%% Covariates Analysis

%% (DATA SPECIFIC) Set paths and variable names
% required data (set to true for the first use of ft2limo)
model_name = 'model_psycho_and_image_new';

SOURCE_ANALYSIS             = false;

task_name                   = 'semanticPriming';

% paths
PATH_TO_SOURCE              = 'D:\__EEG-data\BIDS_source';
PATH_TO_FIELDTRIP           = 'D:\_MSoC\fieldtrip\fieldtrip';
PATH_TO_LIMO                = 'C:\Users\luca-\OneDrive - UMONS\_PhD\_Matlab\3-LIMO\limo_tools';
PATH_TO_FT2LIMO             = 'C:\Users\luca-\OneDrive - UMONS\_PhD\_Matlab\3-LIMO\LIMO-for-FieldTrip';
PATH_TO_ROOT                = 'D:\__EEG-data\BIDS_files';
PATH_TO_COV_ANALYSIS        = 'C:\Users\luca-\OneDrive - UMONS\_PhD\_Matlab\3-LIMO\Covariates_Analysis';

% output folder (derivatives)
PATH_TO_DERIV               = fullfile(PATH_TO_ROOT, 'derivatives');

% add toolboxes to path
ft_defaults
addpath(PATH_TO_LIMO)
addpath(genpath(fullfile(PATH_TO_LIMO,'external')))
addpath(genpath(fullfile(PATH_TO_LIMO,'limo_cluster_functions')))
addpath(genpath(PATH_TO_FT2LIMO))
addpath(genpath(PATH_TO_COV_ANALYSIS))

BIDS_FOLDER = PATH_TO_ROOT;

cd(PATH_TO_ROOT)


%% Extract the covariates
PATH_TO_IMAGES = 'D:\_ARC-images\ARC t?che Final\ARC t?che Final\300x300';
% [psychoFeatures,imageFeatures] = extract_cov(PATH_TO_IMAGES);
imageFeatures = extract_cov(PATH_TO_IMAGES);

PATH_TO_SIMILARIITY = 'C:\Users\luca-\OneDrive - UMONS\_PhD\ARC\visual_similarity.csv';
PATH_TO_ITEMS = 'C:\Users\luca-\OneDrive - UMONS\_PhD\ARC\cible_amorce.csv';
visualSimilarity = extract_visual_sim(PATH_TO_SIMILARIITY, PATH_TO_ITEMS);

tmp = struct2array(rmfield(imageFeatures,'name'));
% correl = corrcoef([freqPCA tmp]);
correl = corrcoef(tmp);
% correl(correl<0.55 & correl>-0.55) = 0;
figure;imagesc(correl)

%% compute PCA for homogeneity-contrast related features
tmp = [zscore(imageFeatures.entropyValue) zscore(imageFeatures.contrastValue) zscore(imageFeatures.energyValue)...
    zscore(imageFeatures.homogeneityValue)];
[~, score, ~] = pca(tmp);
contrastPCA = score(:,1);

% compute PCA for frequency related features
tmp = [zscore(imageFeatures.numberOfCluster) zscore(imageFeatures.maxClusterFreq) zscore(imageFeatures.maxClusterDist)];
[~, score, ~] = pca(tmp);
freqPCA = score(:,1);


imageFeaturesReduced = rmfield(imageFeatures,{'name','entropyValue','contrastValue','energyValue','homogeneityValue','numberOfCluster','maxClusterFreq','maxClusterDist'});
imageFeaturesReduced = table2struct([struct2table(imageFeaturesReduced) table(contrastPCA, 'VariableNames', {'contrastPCA'}) table(freqPCA, 'VariableNames', {'freqPCA'})]);
% f = imageFeaturesReduced.Properties.VariableNames;
f = fieldnames(imageFeaturesReduced);

%%
dinfo = dir(fullfile(BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    load(fullfile(BIDS_FOLDER,'derivatives',subj_name, 'eeg','nat_man_trialinfo.mat'));
    trialinfo = table2struct(trialinfo);
    
    for j = 1:size(trialinfo,1)
%         names = cellfun(@(x) x(1:end-4),cellstr(char(imageFeatures.name)), 'UniformOutput',false);

        idxTarget = cell2mat(cellfun(@(x) strcmp(x(1:end-4),trialinfo(j).target{1}),cellstr(char(imageFeatures.name)), 'UniformOutput',false));
        idxPrimer= cell2mat(cellfun(@(x) strcmp(x(1:end-4),trialinfo(j).primer{1}),cellstr(char(imageFeatures.name)), 'UniformOutput',false));
%         idx = find(cellfun(@any, cellfun(@(x) strcmp(x(1:end-4),trialinfo(j).target{1}),imageFeatures.name, 'UniformOutput',false)));
        for k = 1:length(f)
            try
                trialinfo(j).([f{k} '_primer']) = imageFeaturesReduced(idxPrimer).(f{k});
                trialinfo(j).([f{k} '_target']) = imageFeaturesReduced(idxTarget).(f{k});
            catch
                trialinfo(j).([f{k} '_primer']) = NaN;
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
    
    save(fullfile(BIDS_FOLDER,'derivatives',subj_name, 'eeg','new_psycho_and_image_trialinfo.mat'),'trialinfo');
end

%% Model with image features only
dinfo = dir(fullfile(BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    load(fullfile(BIDS_FOLDER,'derivatives',subj_name, 'eeg','new_psycho_and_image_trialinfo.mat'));
    f = trialinfo.Properties.VariableNames;
    for k = 4:15
        trialinfo = removevars(trialinfo,f(k));
    end
    
    save(fullfile(BIDS_FOLDER,'derivatives',subj_name, 'eeg','new_image_trialinfo.mat'),'trialinfo');
end


%% Create the complete model

contrast.mat = [1 -1 0];
regress_cat = {1:2 ,1;
                0  ,0};

% select the desired regressors
% my_trialinfo = 'new_image_trialinfo.mat';
% selected_regressors = 4:28; %selection from trialinfo.Properties.VariableNames

my_trialinfo = 'new_psycho_and_image_trialinfo.mat';
selected_regressors = [4,5,8:11,14,15:28]; %selection from trialinfo.Properties.VariableNames
trial_start = -200; %starting time of the trial in ms
trial_end = 500; %ending time of the trial in ms


% model is a structure that specifiy information to build a model
model = create_model(PATH_TO_DERIV,PATH_TO_SOURCE,SOURCE_ANALYSIS,task_name,my_trialinfo,trial_start,trial_end,selected_regressors,regress_cat);
if isempty(model.cont_files{1})
    model.cont_files = {};
end
save(fullfile(BIDS_FOLDER,'derivatives',[model_name,'.mat']),'model')


%% Remove previous GLM folder

dinfo = dir(fullfile(BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    del_folder = fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg','GLM_OLS_Time_Channels');
    [root,name,ext] = fileparts(del_folder);
    cd(root)
    if exist(del_folder,'dir')
        rmdir(name,'s')
    end
end

%% Run limo_batch on complete model

cd(PATH_TO_ROOT)
option = 'both'; % or 'model specification', 'contrast only' or 'both'
[LIMO_files, procstatus] = limo_batch(option,model,contrast);

%% Rename GLM folder

dinfo = dir(fullfile(BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    del_folder = fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg','GLM_OLS_Time_Channels');
    [root,name,ext] = fileparts(del_folder);
    cd(root)
    if exist(del_folder,'dir')
        new_name = 'new_psycho_image_GLM_OLS_Time_Channels';
        if exist(fullfile(root,new_name),'dir')
            rmdir(new_name,'s')
        end
        movefile(name,new_name);
        cd(fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg',new_name))
        load('LIMO.mat')
        LIMO.dir = pwd;
        save('LIMO.mat','LIMO')
    end
end

%%
% limo_review(LIMO)
tmp = LIMO.design.X;
tmp(tmp(:,1:2)==0) = -1;
figure;imagesc(tmp);colormap(gray);caxis([-1,1]);

% xticks(1:3);
% xticklabels({'Manufactured\newlineitem', 'Natural\newlineitem', 'Error'});

xticks(1:11);
xticklabels({'Manufactured\newlineitem', 'Natural\newlineitem', '#Phoneme\newlinePrimer','#Phoneme\newlineTarget', ...
    'Movie\newlineFrequency\newlinePrimer', 'Movie\newlineFrequency\newlineTarget', 'AoA\newlinePrimer','AoA\newlineTarget',...
    'Familiarity\newlinePrimer','Familiarity\newlineTarget', 'Error'});

% xticks(1:16);
% xticklabels({'Manufactured\newlineitem', 'Natural\newlineitem', 'Correlation\newlinePrimer', 'Correlation\newlineTarget',...
%     'Compactness\newlinePrimer', 'Compactness\newlineTarget', 'Ratio\newlinePrimer','Ratio\newlineTarget',...
%     'FreqEnergy\newlinePrimer', 'FreqEnergy\newlineTarget', 'Contrast\newlinePrimer','Contrast\newlineTarget',...
%     'frequency\newlinePrimer', 'frequency\newlineTarget', 'Visual\newlineSimilarity', 'Error'});

% xticks(1:24);
% xticklabels({'Man.\newlineitem', 'Nat.\newlineitem', '#Phon.\newlinePrimer','#Phon.\newlineTarget', ...
%     'MovFreq.\newlinePrimer', 'MovFreq.\newlineTarget', 'AoA\newlinePrimer','AoA\newlineTarget',...
%     'Fam.\newlinePrimer','Fam.\newlineTarget',...
%     'Corr.\newlinePrimer', 'Corr.\newlineTarget',...
%     'Compact.\newlinePrimer', 'Compact.\newlineTarget', 'Ratio\newlinePrimer','Ratio\newlineTarget',...
%     'FreqEn.\newlinePrimer', 'FreqEn.\newlineTarget', 'Contr.\newlinePrimer','Contr.\newlineTarget',...
%     'freq.\newlinePrimer', 'freq.\newlineTarget', 'Visual\newlineSim.', 'Error'});

% labels = {'#Phoneme', 'LexicalFreq', 'MovieFreq', 'AoA', 'VisualComplexity', 'Familiarity',...
%     'Imageability'};

ylabel('trials / subjects')


%% Create the naive model
my_name = 'psycho_image';

load(fullfile(BIDS_FOLDER,'derivatives',[model_name,'.mat']))
numberOfCov = size(model.cont_files{1},2);

model_complete = model;
cd(PATH_TO_ROOT)
dinfo = dir(fullfile(PATH_TO_DERIV,'sub-*'));
subj = {dinfo.name};
rng('shuffle')
option = 'both'; % or 'model specification', 'contrast only' or 'both'

tmp_betas = cell(length(subj),1);
tmp_con = cell(length(subj),1);
betas_mean = cell(length(subj),1);
con_mean = cell(length(subj),1);

% repeat xxx times to get mean naive
for j=1:10
    i = 1;
    for subj_name = subj
        subj_name = subj_name{1};
        load(fullfile(PATH_TO_DERIV,subj_name, 'eeg', [my_name '_GLM_OLS_Time_Channels'], 'LIMO.mat'))
        idx = ~isnan(model_complete.cont_files{i});
        model_tmp = model_complete.cont_files{i};
        delete(fullfile(BIDS_FOLDER,'limo_batch_report','GLM_OLS_Time_Channels',['subject' num2str(i)],'PIPE.lock'))
        naive_path = fullfile(PATH_TO_DERIV,subj_name, 'eeg', [my_name '_naive_GLM_OLS_Time_Channels']);
        glm_path = fullfile(PATH_TO_DERIV,subj_name, 'eeg', 'GLM_OLS_Time_Channels');
        if ~exist(glm_path,'dir')
            mkdir(glm_path)
        end
        LIMO.dir = glm_path;
        if isfield(LIMO,'contrast')
            LIMO = rmfield(LIMO,'contrast');
            save(fullfile(glm_path,'LIMO.mat'),'LIMO')
        end
        null_mat = mvnrnd(zeros(1,numberOfCov),cov(LIMO.design.X(:,3:numberOfCov+2)),size(LIMO.design.X,1));
        model_tmp(sum(idx,2)==numberOfCov,:) = null_mat;
        model.cont_files{i} = model_tmp;
        
        i = i + 1;
    end
    
    [LIMO_files, procstatus] = limo_batch(option,model,contrast);
    
    i=1;
    for subj_name = subj
        subj_name = subj_name{1};
        tmp = load(LIMO_files.mat{i});
        load(fullfile(tmp.LIMO.dir,'Betas.mat'));
        load(fullfile(tmp.LIMO.dir,'con_1.mat'));
        if j==1
            betas_mean{i} = Betas;
            con_mean{i} = con;
        else
            betas_mean{i} = ((j-1)*betas_mean{i} + Betas)./j;
            con_mean{i} = ((j-1)*con_mean{i} + con)./j;
        end
%         tmp.LIMO.data.Cont = model_tmp;
%         tmp.LIMO.design.X = 
        LIMO = tmp.LIMO;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')

        i = i + 1;
    end	
    save(fullfile(PATH_TO_DERIV,[my_name '_betas_mean.mat']), 'betas_mean')
    save(fullfile(PATH_TO_DERIV,[my_name '_con_mean.mat']), 'con_mean')
end

i=1;
for subj_name = subj
    subj_name = subj_name{1};
%     load(fullfile(PATH_TO_DERIV,subj_name, 'eeg', 'GLM_OLS_Time_Channels', 'LIMO.mat'))

%     Betas = betas_mean{i};
%     save(fullfile(LIMO.dir,'Betas.mat'),'Betas')
%     con = con_mean{i};
%     save(fullfile(LIMO.dir,'con_1.mat'),'con')
    
    Betas = betas_mean{i};
    save(fullfile(PATH_TO_DERIV,subj_name, 'eeg', 'GLM_OLS_Time_Channels','Betas.mat'),'Betas')
    con = con_mean{i};
    save(fullfile(PATH_TO_DERIV,subj_name, 'eeg', 'GLM_OLS_Time_Channels','con_1.mat'),'con')    
    
    i = i + 1;
end

save(fullfile(PATH_TO_DERIV, [model_name,'_naive.mat']), 'model')


%% Rename naive GLM folder

dinfo = dir(fullfile(BIDS_FOLDER,'sub-*'));
subj = {dinfo.name};
for subj_name = subj
    subj_name = subj_name{1};
    del_folder = fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg','GLM_OLS_Time_Channels');
%     del_folder = fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg','new_nat_man_naive_GLM_OLS_Time_Channels');
    [root,name,ext] = fileparts(del_folder);
    cd(root)
    if exist(del_folder,'dir')
        new_name = [my_name '_naive_GLM_OLS_Time_Channels'];
%         new_name = 'GLM_OLS_Time_Channels';
        if exist(fullfile(root,new_name),'dir')
            rmdir(new_name,'s')
        end
        movefile(name,new_name);
%         copyfile(name,new_name);
        cd(fullfile(BIDS_FOLDER,'derivatives',subj_name,'eeg',new_name))
        if exist('LIMO.mat','file')
            load('LIMO.mat')
            LIMO.dir = pwd;
            save('LIMO.mat','LIMO')
        end
    end
end

%% Find cluster of significant contrast within complete model

% modify Beta_files_GLM_OLS_Time_Channels.txt files to access complete and
% naive Betas

cd(PATH_TO_ROOT)
expected_chanlocs = limo_avg_expected_chanlocs(PATH_TO_DERIV, model.defaults);

my_name = 'nat_man_simple';
my_param = 'con_1';

%     my_con = 'con_1';
LIMOfiles = fullfile(PATH_TO_ROOT,sprintf('%s_files_GLM_OLS_Time_Channels_%s.txt',my_param,my_name));
if ~exist(fullfile(PATH_TO_ROOT,[my_name '_' my_param]),'dir')
    mkdir(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
end
cd(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
LIMOPath = limo_random_select('one sample t-test',expected_chanlocs,'LIMOfiles',... 
    LIMOfiles,'analysis_type','Full scalp analysis',...
    'type','Channels','nboot',100,'tfce',1,'skip design check','yes');

% limo_results %find regions of significant contrast through clustering algo with p=0.05
% pause()
% mask_complete = mask;

cd(fullfile(PATH_TO_ROOT,'nat_man_simple_con_1'))
limo_results %find regions of significant contrast in categorical data only model
pause()
mask_simple = mask;

%% Find cluster of significant beta within both models (where do the covariates have a significant influence on the model)
numberOfCov = size(model.cont_files{1},2);
mask_beta = {};

my_name = 'psycho_image_naive';
my_param = 'Beta';

LIMOfiles = fullfile(PATH_TO_ROOT,sprintf('%s_files_GLM_OLS_Time_Channels_%s.txt',my_param,my_name));
if ~exist(fullfile(PATH_TO_ROOT,[my_name '_' my_param]),'dir')
    mkdir(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
end
cd(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
LIMOPath = limo_random_select('one sample t-test',expected_chanlocs,'LIMOfiles',... 
    LIMOfiles,'analysis_type','Full scalp analysis',...
    'parameters',{3:numberOfCov+2},'type','Channels','nboot',100,'tfce',1,'skip design check','yes');

% find regions of significant effect through clustering algo with p=0.05 for each covariate
j=1;
for i = numberOfCov+2:-1:3
    if j==1
        disp('complete model')
    else
        disp('naive model')
    end
    disp(i)
    mask = zeros(size(mask));
    limo_results
    pause()
    mask_beta{i,j} = mask;
end

mask_beta_naive = mask_beta;
save(fullfile(PATH_TO_ROOT,[my_name '_' my_param],'mask_beta_naive.mat'),'mask_beta_naive')

%%

my_name = 'image_naive';
my_param = 'Beta';

LIMOfiles = fullfile(PATH_TO_ROOT,sprintf('%s_files_GLM_OLS_Time_Channels_%s.txt',my_param,my_name));
if ~exist(fullfile(PATH_TO_ROOT,[my_name '_' my_param]),'dir')
    mkdir(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
end
cd(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
LIMOPath = limo_random_select('one sample t-test',expected_chanlocs,'LIMOfiles',... 
    LIMOfiles,'analysis_type','Full scalp analysis',...
    'parameters',{3:numberOfCov+2},'type','Channels','nboot',100,'tfce',1,'skip design check','yes');


% find regions of significant effect through clustering algo with p=0.05 for each covariate
j=2;
for i = numberOfCov+2:-1:3
    if j==1
        disp('complete model')
    else
        disp('naive model')
    end
    disp(i)
    mask = zeros(size(mask));
%     limo_results
%     pause()
    mask_beta{i,j} = mask;
end


%% Plot unthresholded contrast maps from each model
close all

SE = strel('rectangle',[6 15]);
new_mask = imdilate(mask_simple,SE);
SE = strel('rectangle',[4 13]);
new_mask2 = imerode(new_mask,SE);
[B_simple,~] = bwboundaries(new_mask2);

figure;
p=1;
% modelList = {'nat_man_simple_con_1','new_nat_man_con_1','image_con_1','psycho_image_con_1'};
modelList = {'nat_man_simple_con_1','new_nat_man_con_1','image_con_1','psycho_image_con_1'};
titleList = {'control model', 'psycho model', 'image model', 'psycho-image model'};
for m = modelList
    tmpModel = m{1};
    load(fullfile(BIDS_FOLDER,tmpModel,'one_sample_ttest_parameter_1.mat'))
    unthresh_map = one_sample(:,:,1);
    subplot(2,2,p)
    imagesc(unthresh_map);
    caxis([-1.5,1.5])
    title(titleList{p}, 'FontSize', 14)
    p=p+1;
    hold on
    colorbar
    
    %categorical effect (simple model)
    for k = 1:length(B_simple)
       boundary = B_simple{k};
       plot(boundary(:,2), boundary(:,1), 'LineWidth', 2,'Color',[0, 1, 1, 0.6])
    end
    
    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
    xlabel('time(ms)', 'FontSize', 14)
    ylabel('electrode', 'FontSize', 14)
%     hold off
end
% suptitle('Contrast unthresholded maps from each model')
%  axes( 'Position', [0, 0.95, 1, 0.05] ) ;
%  set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%  text( 0.5, 0, 'Contrast unthresholded maps from each model', 'FontSize', 18, 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

h = plot(NaN,NaN,'-', 'LineWidth', 2,'Color',[0, 1, 1, 0.6]);
my_legend = legend(h(1), 'control clusters');
% newPosition = [0.47 0.48 0.08 0.05];
% newPosition = [0.47 0.02 0.08 0.05];
newPosition = [0.78 0.01 0.08 0.05];
newUnits = 'normalized';
% my_legend = legend([f(3),f(2),f(1)],'covariate effect','random covariate effect','categorical effect');
set(my_legend,'Position', newPosition,'Units', newUnits, 'FontSize',16);

%% Plot covariates influence on the model
% my_folder = 'new_nat_man_complete_Beta';
% my_folder = 'new_image_Beta';
my_folder = 'psycho_image_naive_Beta';
numberOfCov = size(model.cont_files{1},2);
% fields = {'#Phoneme\newlinePrimer','#Phoneme\newlineTarget', ...
%     'Movie\newlineFrequency\newlinePrimer', 'Movie\newlineFrequency\newlineTarget', 'AoA\newlinePrimer','AoA\newlineTarget',...
%     'Familiarity\newlinePrimer','Familiarity\newlineTarget'};
fields = {'Correlation\newlinePrimer', 'Correlation\newlineTarget',...
    'Compactness\newlinePrimer', 'Compactness\newlineTarget', 'Ratio\newlinePrimer','Ratio\newlineTarget',...
    'FreqEnergy\newlinePrimer', 'FreqEnergy\newlineTarget', 'Contrast\newlinePrimer','Contrast\newlineTarget',...
    'frequency\newlinePrimer', 'frequency\newlineTarget', 'Visual Similarity'};

% mask of categorical effect in complete model
% SE = strel('rectangle',[3 15]);
% new_mask = imdilate(mask_complete,SE);
% SE = strel('rectangle',[2 10]);
% new_mask2 = imerode(new_mask,SE);
% [B_complete,~] = bwboundaries(new_mask2);

% SE = strel('rectangle',[3 15]);
% new_mask = imdilate(mask_simple,SE);
% SE = strel('rectangle',[2 10]);
% new_mask2 = imerode(new_mask,SE);
% [B_simple,~] = bwboundaries(new_mask2);

SE = strel('rectangle',[6 15]);
new_mask = imdilate(mask_simple,SE);
SE = strel('rectangle',[4 13]);
new_mask2 = imerode(new_mask,SE);
[B_simple,~] = bwboundaries(new_mask2);

% figure;
% hold on
% boundary = B_complete{1};
% plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2,'Color',[1, 0, 0, 0.5])
% plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2,'Color',[0, 1, 0, 0.5])
% plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2,'Color',[0, 0, 1, 0.5])
% f=get(gca,'Children');
% my_legend = legend([f(3),f(2),f(1)],'covariate effect','random covariate effect','categorical effect');
% close all

% primer influence
figure;
p=1;
for i = 3:2:numberOfCov+2
    load(fullfile(my_folder,['parameter_' num2str(i)],['one_sample_ttest_parameter_' num2str(i) '.mat']))
    unthresh_map = one_sample(:,:,1);
    subplot(ceil(numberOfCov/4),2,p)
    p=p+1;
    imagesc(unthresh_map); colorbar;
    caxis([-1,1])
    hold on

    %covariate effect (psycho or image only)
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta{i,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6])
    end
    
    %covariate effect (psycho only)
%     SE = strel('rectangle',[3 15]);
%     new_mask = imdilate(mask_beta_all{i,1},SE);
    %covariate effect (image only)
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta_all{i+8,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6])
    end
    
%     %random covariate effect
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta_naive{i,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6])
    end
    
%     %categorical effect (complete model
%     for k = 1:length(B_complete)
%        boundary = B_complete{k};
%        plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6])
%     end
    
    %categorical effect (simple model)
    for k = 1:length(B_simple)
       boundary = B_simple{k};
       plot(boundary(:,2), boundary(:,1), 'LineWidth', 2,'Color',[0, 1, 1, 0.6])
    end
    
    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
%     xlabel('time(ms)')
%     ylabel('electrode')
    xlabel('time(ms)','FontSize',16)
    ylabel('electrode','FontSize',16)
    f = fields{i-2};
    idx = strfind(f,'\newline');
    if isempty(idx)
        idx = length(f)+1;
    end
    title(strrep(f(1:idx(end)-1),'\newline',' '))
%     switch i
%         case 3
%             title('phoneme number')
%         case 5
%             title('movie occurence')
%         case 7
%             title('AoA')
%         case 9
%             title('familiarity')
%     end
end

% get(gca,'Children');
h = zeros(3,1);
h(1) = plot(NaN,NaN,'-', 'LineWidth', 2,'Color',[0, 1, 1, 0.6]);
h(2) = plot(NaN,NaN, '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6]);
h(3) = plot(NaN,NaN, '-.b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6]);
% h(4) = plot(NaN,NaN, 'LineWidth', 2,'Color',[0, 1, 1, 0.9]);
% my_legend = legend(h, 'covariate effect','random covariate effect','categorical effect (complete)','categorical effect (simple)');
% my_legend = legend(h, 'control clusters','psycho clusters','psycho-image clusters');
my_legend = legend(h, 'control clusters','image clusters','psycho-image clusters');
% newPosition = [0.468 0.422 0.1 0.1];
newPosition = [0.46 0.01 0.1 0.08];
newUnits = 'normalized';
% my_legend = legend([f(3),f(2),f(1)],'covariate effect','random covariate effect','categorical effect');
set(my_legend,'Position', newPosition,'Units', newUnits);
% suptitle('Influence of primer covariates on the regression')


% target influence
figure;
p=1;
for i = 4:2:numberOfCov+2
    load(fullfile(my_folder,['parameter_' num2str(i)],['one_sample_ttest_parameter_' num2str(i) '.mat']))
    unthresh_map = one_sample(:,:,1);
    subplot(floor(numberOfCov/4),2,p)
    p=p+1;
    imagesc(unthresh_map); colorbar;
    caxis([-1,1])
    hold on

    %covariate effect (psycho or image only)
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta{i,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6])
    end
    
    %covariate effect (psycho only)
%     SE = strel('rectangle',[3 15]);
%     new_mask = imdilate(mask_beta_all{i,1},SE);
    %covariate effect (image only)
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta_all{i+8,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6])
    end
    
    %categorical effect (simple model)
    for k = 1:length(B_simple)
       boundary = B_simple{k};
       plot(boundary(:,2), boundary(:,1), 'LineWidth', 2,'Color',[0, 1, 1, 0.9])
    end
    
    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
    xlabel('time(ms)')
    ylabel('electrode')
    f = fields{i-2};
    idx = strfind(f,'\newline');
    title(strrep(f(1:idx(end)-1),'\newline',' '))
%     switch i
%         case 4
%             title('phoneme number')
%         case 6
%             title('visual frequency')
%         case 8
%             title('AoA')
%         case 10
%             title('familiarity')
%     end
end

h = zeros(3,1);
h(1) = plot(NaN,NaN,'-', 'LineWidth', 2,'Color',[0, 1, 1, 0.6]);
h(2) = plot(NaN,NaN, '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6]);
h(3) = plot(NaN,NaN, '-.b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6]);
% my_legend = legend(h, 'control clusters','psycho clusters','psycho-image clusters');
my_legend = legend(h, 'control clusters','image clusters','psycho-image clusters');
newPosition = [0.46 0.01 0.1 0.08];
newUnits = 'normalized';
set(my_legend,'Position', newPosition,'Units', newUnits);


% h = zeros(4,1);
% h(1) = plot(NaN,NaN,'-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6]);
% h(2) = plot(NaN,NaN, '--k', 'LineWidth', 1,'Color',[0, 0.4, 0, 0.8]);
% h(3) = plot(NaN,NaN, 'b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6]);
% h(4) = plot(NaN,NaN, 'LineWidth', 2,'Color',[0, 1, 1, 0.9]);
% my_legend = legend(h, 'covariate effect','random covariate effect','categorical effect (complete)','categorical effect (simple)');
% newPosition = [0.465 0.422 0.1 0.1];
% newUnits = 'normalized';
% set(my_legend,'Position', newPosition,'Units', newUnits);
% suptitle('Influence of target covariates on the regression')


%% image features influence
paramNumber = 3:9; %the numbers  corresponding to the parameters we want to study
figure;
p=1;
nRows = 4;
nCols = 2;
for i = paramNumber
    load(fullfile(['parameter_' num2str(i)],['one_sample_ttest_parameter_' num2str(i) '.mat']))
    unthresh_map = one_sample(:,:,1);
    subplot(nRows,nCols,p)
    p=p+1;
    imagesc(unthresh_map);
    caxis([-1,1])
    hold on

    %covariate effect
    SE = strel('rectangle',[3 15]);
    new_mask = imdilate(mask_beta{i,1},SE);
    SE = strel('rectangle',[2 14]);
    new_mask2 = imerode(new_mask,SE);
    [B,~] = bwboundaries(new_mask2);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6])
    end
    
%     %random covariate effect
%     SE = strel('rectangle',[3 15]);
%     new_mask = imdilate(mask_beta{i,2},SE);
%     SE = strel('rectangle',[2 14]);
%     new_mask2 = imerode(new_mask,SE);
%     [B,~] = bwboundaries(new_mask2);
%     for k = 1:length(B)
%        boundary = B{k};
%        plot(boundary(:,2), boundary(:,1), '--k', 'LineWidth', 1,'Color',[0, 0.4, 0, 0.8])
%     end
    
    %categorical effect (complete model)
    for k = 1:length(B_complete)
       boundary = B_complete{k};
       plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6])
    end
    
    %categorical effect (simple model)
    for k = 1:length(B_simple)
       boundary = B_simple{k};
       plot(boundary(:,2), boundary(:,1), 'LineWidth', 2,'Color',[0, 1, 1, 0.9])
    end
    
    
    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
    xlabel('time(ms)')
    ylabel('electrode')
    switch i
        case paramNumber(1)
            title('correlation')
        case paramNumber(2)
            title('compactness')
        case paramNumber(3)
            title('ratio')
        case paramNumber(4)
            title('spectral energy') 
        case paramNumber(5)
            title('contrast')
        case paramNumber(6)
            title('frequency')
        case paramNumber(7)
            title('visual similarity')
    end
end

% h = zeros(3,1);
h = zeros(2,1);
h(1) = plot(NaN,NaN,'-.r', 'LineWidth', 4,'Color',[1, 0, 0, 0.6]);
% h(2) = plot(NaN,NaN, '--k', 'LineWidth', 1,'Color',[0, 0.4, 0, 0.8]);
% h(3) = plot(NaN,NaN, 'b', 'LineWidth', 2,'Color',[0, 0, 1, 0.6]);
h(2) = plot(NaN,NaN, 'b', 'LineWidth', 4,'Color',[0, 0, 1, 0.6]);
h(3) = plot(NaN,NaN, 'b', 'LineWidth', 4,'Color',[0, 1, 1, 1]);
% my_legend = legend(h, 'covariate effect','random covariate effect','categorical effect');
my_legend = legend(h, 'covariate effect','categorical effect (complete)','categorical effect (simple)');
% newPosition = [0.465 0.422 0.1 0.1];
newPosition = [0.7 0.12 0.1 0.1];
newUnits = 'normalized';
set(my_legend,'Position', newPosition,'Units', newUnits,'FontSize',16);
suptitle('Influence of image features on the regression')

%% Plot covariate influence on time/region of interest
% idx = mask_complete==3; %cluster 1
idx = ~ismember(mask_complete,[0,3]); %cluster 2

group_partial_coef = [];
group_R2 = [];
i=1;
for my_path = model.set_files'
    disp(i)
    my_path = char(my_path);
    [root,~,~] = fileparts(my_path);
    cd([root '\new_nat_man_GLM_OLS_Time_Channels'])
    LIMO = load('LIMO.mat');
    LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
    R2 = load('R2.mat');
    R2 = R2.(cell2mat(fieldnames(R2)));
    R2 = R2(:,:,1);
    group_R2(i,1) = mean(mean(R2(idx)));

    limo_semi_partial_coef(LIMO); %output dim: channel*time*(R2,F-value,p-value)
    for j = 1:9
        semi_partial_coef = load([root '\new_nat_man_GLM_OLS_Time_Channels\semi_partial_coef_' num2str(j) '.mat']);
        semi_partial_coef = semi_partial_coef.(cell2mat(fieldnames(semi_partial_coef)));
        semi_partial_coef = semi_partial_coef(:,:,1);
        avg_partial_coef = squeeze(mean(mean(abs(semi_partial_coef(idx)))));
        group_partial_coef(i,j,1) = avg_partial_coef;
    end
    
    cd([root '\new_nat_man_naive_GLM_OLS_Time_Channels'])
    LIMO = load('LIMO.mat');
    LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
    load('R2.mat')
    R2 = R2(:,:,1);
    group_R2(i,2) = mean(mean(R2(idx)));
% 
    limo_semi_partial_coef(LIMO); %output dim: channel*time*(R2,F-value,p-value)
    for j = 1:9
        semi_partial_coef = load([root '\new_nat_man_naive_GLM_OLS_Time_Channels\semi_partial_coef_' num2str(j) '.mat']);
        semi_partial_coef = semi_partial_coef.(cell2mat(fieldnames(semi_partial_coef)));
        semi_partial_coef = semi_partial_coef(:,:,1);
        avg_partial_coef = squeeze(mean(mean(abs(semi_partial_coef(idx)))));
        group_partial_coef(i,j,2) = avg_partial_coef;
    end
    i=i+1;
end

% complete_group_R2 = group_R2;
% complete_partial_coef = group_partial_coef;

%% Plot the results
fontSize = 10;

%R2 simple vs R2 cov
target_group = group_R2;
isout = isoutlier(target_group,'quartiles');
xClean = target_group;
xClean(isout) = NaN;

[est,HDI]=data_plot(xClean,'estimator','trimmed mean'); % test with estimator 
xticks([1,2.25])
xticklabels({'complete', 'naive'})
title(sprintf("Explained variance by model\n(R^2 values)"),'FontSize',fontSize)
ylim([0 0.25])

figure()

target_group = [group_partial_coef(:,:,1), group_partial_coef(:,:,2)];
isout = isoutlier(target_group,'quartiles');
xClean = target_group;
xClean(isout) = NaN;
subplot(2,1,1)
[est2,HDI2]=data_plot(xClean(:,1:9),'estimator','trimmed mean','figure','on');
my_ticks = 1;
for i = 1:8
    my_ticks = [my_ticks my_ticks(end)+1.25];
end
xticks(my_ticks)
xticklabels({'cat', 'cov1', 'cov2', 'cov3', 'cov4', 'cov5', 'cov6', 'cov7', 'cov8'})
title(sprintf("Semi-partial variance\ncomplete model"),'FontSize',fontSize)
ylim([0 0.05])

subplot(2,1,2)
[est3,HDI3]=data_plot(xClean(:,10:end),'estimator','trimmed mean','figure','on');
my_ticks = 1;
for i = 1:8
    my_ticks = [my_ticks my_ticks(end)+1.25];
end
xticks(my_ticks)
xticklabels({'cat', 'cov1', 'cov2', 'cov3', 'cov4', 'cov5', 'cov6', 'cov7', 'cov8'})
title(sprintf("Semi-partial variance\nnaive model"),'FontSize',fontSize)
ylim([0 0.05])

% mean and std of R2 from complete model within categorical clusters
figure;
mean_R2 = [mean(complete_group_R2(:,1)); mean(group_R2(:,1))];
std_R2 = [std(complete_group_R2(:,1)); std(group_R2(:,1))];
errorbar(mean_R2,std_R2,'x')
xticks([1,2])
xlim([0,3])
ylim([0 0.25])
xticklabels({'cluster1','cluster2'})
title('Total explained variance (R^2) of the complete model\newline                    within the categorical clusters\newline')

% mean and std of partial coef from complete model within categorical clusters
figure;
mean_R2 = [mean(complete_partial_coef(:,:,1))'; mean(group_partial_coef(:,:,1))'];
std_R2 = [std(complete_partial_coef(:,:,1))'; std(group_partial_coef(:,:,1))'];
errorbar(mean_R2,std_R2,'x')
xticks(1:18)
xlim([0,19])
ylim([0 0.025])
xticklabels({'cat', 'cov1', 'cov2', 'cov3', 'cov4', 'cov5', 'cov6', 'cov7', 'cov8','cat', 'cov1', 'cov2', 'cov3', 'cov4', 'cov5', 'cov6', 'cov7', 'cov8'})
title('Explained variance by covariate (partial coef) of the complete model\newline                    within the categorical clusters\newline')

