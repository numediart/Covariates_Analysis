% Boxplot explained variance
% folders: last_psycho_R2; last_psycho_image_R2; image_R2; naive_R2; 
% files: last_psycho; last_psycho_image; last_psycho_image_naive;
% last_psycho_naive; new_image; new_image_naive
% trialinfo: new_image_trialinfo; trialinfo_psycho; trialinfo_psycho_image

dinfo = dir(fullfile(PATH_TO_DERIV,'sub-*'));
% complete_R2 = [];
% control_R2 = [];
% naive_R2 = [];
% naive_psycho_R2 = [];
% naive_image_R2 = [];
% psycho_R2 = [];
% image_R2 = [];
% complete_con = [];
% erpMan = [];
% erpNat= [];
for i = numel( dinfo ):-1:1
    if i >= 10
        subfolder = ['sub-0' num2str(i)];
    else
        subfolder = ['sub-00' num2str(i)];
    end
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','nat_man_simple_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_image_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','new_image_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_image_naive_GLM_OLS_Time_Channels','R2.mat'))
    load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_naive_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','new_image_naive_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_image_GLM_OLS_Time_Channels','con_1.mat'))
%     con = con(:,:,1);
%     complete_con(:,:,i) = con;
    R2 = R2(:,:,1);
%     complete_R2(:,:,i) = R2;
%     control_R2(:,:,i) = R2;
%     naive_R2(:,:,i) = R2;
%     naive_psycho_R2(:,:,i) = R2;
%     naive_image_R2(:,:,i) = R2;
    psycho_R2(:,:,i) = R2;
%     image_R2(:,:,i) = R2;
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','nat_man_simple_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','image_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','new_image_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_image_naive_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','last_psycho_naive_GLM_OLS_Time_Channels','R2.mat'))
%     load(fullfile(PATH_TO_DERIV,subfolder,'eeg','new_image_naive_GLM_OLS_Time_Channels','R2.mat'))
%     R2 = R2(:,:,1);
%     naive_R2(:,:,i) = naive_R2(:,:,i)-R2;
%     naive_psycho_R2(:,:,i) = naive_psycho_R2(:,:,i)-R2;
%     naive_image_R2(:,:,i) = naive_image_R2(:,:,i)-R2;
%     psycho_R2(:,:,i) = psycho_R2(:,:,i)-R2;
%     image_R2(:,:,i) = image_R2(:,:,i)-R2;
%     complete_R2(:,:,i) = complete_R2(:,:,i)-R2;
end

% final_psycho_r2 = psycho_R2;
% final_image_r2 = image_R2;
% final_psycho_naive_r2 = psycho_R2;
% final_image_naive_r2 = image_R2;

% final_psycho_cat_r2 = psycho_R2-final_psycho_r2;
% final_image_cat_r2 = image_R2;


%%
tmpR2Control = [];
tmpR2Psycho = [];
tmpR2Image = [];
tmpR2PsychoNaive = [];
tmpR2ImageNaive = [];
% my_mask = mask_psycho_R2;
% my_mask = mask_control_con;
% my_mask = mask_complete_con;
my_mask = mask_simple_con1;
% my_mask = mask_psycho_image_con1;
% my_mask = mask_image_R2;
% my_mask(my_mask==3)=0;
targetMask = 1;
toAdd = [0 0 1.1; 0 0.2 0; 0 0 0; 0 0 0];

k=1;
figure;
for allMask = {mask_simple_con1,mask_psycho_image_con1}
    for targetMask = 1:max(unique(allMask{1}))
        my_mask = allMask{1};
        disp(targetMask)
        my_mask(my_mask~=targetMask) = 0;
        [~,col] = find(my_mask);
        my_mask(:,col) = targetMask;

%         tmp2psycho = one_sample(:,:,1);
        tmpR2Psycho = tmp2psycho(ismember(my_mask,targetMask));
%         tmp2image = one_sample(:,:,1);
        tmpR2Image= tmp2image(ismember(my_mask,targetMask));
        % tmp2control = one_sample(:,:,1);
        tmpR2Control= tmp2control(ismember(my_mask,targetMask));

%         tmp2psychonaive = one_sample(:,:,1);
        % tmpR2PsychoNaive = tmp2psychonaive(my_mask>0);
        tmpR2PsychoNaive = tmp2psychonaive(ismember(my_mask,targetMask));
%         tmp2imagenaive = one_sample(:,:,1);
        % tmpR2ImageNaive = tmp2imagenaive(my_mask>0);
        tmpR2ImageNaive = tmp2imagenaive(ismember(my_mask,targetMask));

        my_alpha = 0.05; %95% confidence interval
        N = 100;

        % N = length(tmpR2PsychoNaive(:));
        CI95 = tinv([my_alpha/2 1-my_alpha/2], N-1);
        R2_tmp_mean = mean(tmpR2PsychoNaive(:));
        R2_tmp_std = std(tmpR2PsychoNaive(:));
        R2_tmp_SEM = R2_tmp_std/sqrt(N);
        R2_tmp_CI95 = R2_tmp_mean + R2_tmp_SEM*CI95;
        CI_psycho_naive_mask_control = 100*R2_tmp_CI95-[0.5 0.5];

        % N = length(tmpR2ImageNaive(:));
        CI95 = tinv([my_alpha/2 1-my_alpha/2], N-1);
        R2_tmp_mean = mean(tmpR2ImageNaive(:));
        R2_tmp_std = std(tmpR2ImageNaive(:));
        R2_tmp_SEM = R2_tmp_std/sqrt(N);
        % R2_tmp_CI95 = R2_tmp_mean + bsxfun(@times, R2_tmp_SEM, CI95(:));
        R2_tmp_CI95 = R2_tmp_mean + R2_tmp_SEM*CI95;
        CI_image_naive_mask_control = 100*R2_tmp_CI95;
        if k==1
            CI_image_naive_mask_control = CI_image_naive_mask_control-0.25;
        end


        target_group = 100*[tmpR2Control(:), tmpR2Psycho(:), tmpR2Image(:)] + toAdd(k,:);
        % target_group = unthreshMap(my_mask>0);
        isout = isoutlier(target_group,'quartiles');
        xClean = target_group;
        xClean(isout) = NaN;
        subplot(1,4,k);boxplot(xClean)
        % [est,HDI]=data_plot(xClean,'estimator','trimmed mean'); % test with estimator 
        % xticks([1,2.25,3.5])
        xticklabels({'categorical', 'psycho', 'image'})
        % title(sprintf("Explained variance by model\n(R^2 values)"),'FontSize',fontSize)
        hold on
        x = [0.75 1.25]+1;
        patch([x, fliplr(x)], [CI_psycho_naive_mask_control(1) CI_psycho_naive_mask_control(1)...
            CI_psycho_naive_mask_control(2) CI_psycho_naive_mask_control(2)],...
            [0.5,0.5,0.5], 'EdgeColor','none', 'FaceAlpha',0.9)

        x = [0.75 1.25]+2;
        patch([x, fliplr(x)], [CI_image_naive_mask_control(1) CI_image_naive_mask_control(1)...
            CI_image_naive_mask_control(2) CI_image_naive_mask_control(2)],...
            [0.5,0.5,0.5], 'EdgeColor','none', 'FaceAlpha',0.9)

        % my_legend = legend(f(1),'R2 naive model (95% CI)');
        % set(my_legend,'FontSize',8,'Location','southoutside');
        f=get(gca,'Children');
        g = get(f,'Children');
        set(g{3}(7:9),'Color','k')
        if k==1
            ylabel('R^2 (%)')
        end
        yticks(0:2:20)
        ylim([0 21])
        switch k
            case 1
                title('Explained variance distribution\newline         Categorical cluster 1','FontSize',10)
            case 2
                title('Explained variance distribution\newline         Categorical cluster 2','FontSize',10)
            case 3
                title('Explained variance distribution\newline         Categorical cluster 3','FontSize',10)
            case 4
                title('             Explained variance distribution\newlineCategorical cluster from psycho-image model','FontSize',10)
        end
        
%         clc
        % [h,p,ks2stat] = kstest2(tmpR2Psycho,tmpR2PsychoNaive)
        % [h,p,ks2stat] = kstest2(tmpR2Image,tmpR2ImageNaive)
        quantile(xClean,[.25 .5 .75])
%         mean(xClean,'omitnan')
        CI_psycho_naive_mask_control
        CI_image_naive_mask_control
        k=k+1;
    end
end
%%
cd(PATH_TO_ROOT)
% expected_chanlocs = limo_avg_expected_chanlocs(PATH_TO_DERIV, model.defaults);

my_name = '_FINAL_psycho_naive_direct';
my_param = 'R2';

if ~exist(fullfile(PATH_TO_ROOT,[my_name '_' my_param]),'dir')
    mkdir(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))
end
cd(fullfile(PATH_TO_ROOT,[my_name '_' my_param]))

% load('D:\__EEG-data\BIDS_files\psycho_image_con_1\LIMO.mat')
% load('D:\__EEG-data\BIDS_files\nat_man_simple_con_1\LIMO.mat')
% load('D:\__EEG-data\BIDS_files\psycho_image_naive_Beta\parameter_3\LIMO.mat')
% load('D:\__EEG-data\BIDS_files\image_con_1\LIMO.mat')
load('D:\__EEG-data\BIDS_files\last_psycho_R2\LIMO.mat')
% load('D:\__EEG-data\BIDS_files\last_image_R2\LIMO.mat')
% load('D:\__EEG-data\BIDS_files\last_psycho_naive_R2\LIMO.mat')

LIMO.dir = pwd;
LIMO.data.data_dir = pwd;
LIMO.design.method = 'Trimmed Mean';
% LIMO.design.tfce = 0;
LIMO.design.bootstrap = 100;
save('LIMO.mat','LIMO')

% limo_random_robust(1,image_R2,1,LIMO) %line 426
limo_random_robust(1,psycho_R2,1,LIMO) %line 426

%% cat dimension overlap in each model
% image_direct = one_sample(:,:,1);
% image_true = one_sample(:,:,1);

% psycho_direct = one_sample(:,:,1);
% psycho_true = one_sample(:,:,1);

% image_naive_direct = one_sample(:,:,1);
% psycho_naive_direct = one_sample(:,:,1);

% tmp2control

tmp_image = image_direct-image_true;
tmp_image(tmp_image<-0.005)= nan;

tmp_psycho = psycho_direct-psycho_true;
tmp_psycho(tmp_psycho<-0.005) = nan;

for allMask = {mask_simple_con1,mask_psycho_image_con1}
    for targetMask = 1:max(unique(allMask{1}))
        my_mask = allMask{1};
        disp(targetMask)
        my_mask(my_mask~=targetMask) = 0;
        [~,col] = find(my_mask);
        my_mask(:,col) = targetMask;

        test_control = tmp2control(ismember(my_mask,targetMask));
        test_image = tmp_image(ismember(my_mask,targetMask));
        test_psycho = tmp_psycho(ismember(my_mask,targetMask));

        my_alpha = 0.05; %95% confidence interval
        N = 100;
        
        CI95 = tinv([my_alpha/2 1-my_alpha/2], N-1);
        R2_tmp_mean = mean(test_control(:));
        R2_tmp_std = std(test_control(:));
        R2_tmp_SEM = R2_tmp_std/sqrt(N);
        R2_tmp_CI95 = R2_tmp_mean + R2_tmp_SEM*CI95;
        CI_control = 100*R2_tmp_CI95
        
        CI95 = tinv([my_alpha/2 1-my_alpha/2], N-1);
        R2_tmp_mean = mean(test_psycho(:),'omitnan');
        R2_tmp_std = std(test_psycho(:),'omitnan');
        R2_tmp_SEM = R2_tmp_std/sqrt(N);
        R2_tmp_CI95 = R2_tmp_mean + R2_tmp_SEM*CI95;
        CI_psycho = 100*R2_tmp_CI95

        CI95 = tinv([my_alpha/2 1-my_alpha/2], N-1);
        R2_tmp_mean = mean(test_image(:),'omitnan');
        R2_tmp_std = std(test_image(:),'omitnan');
        R2_tmp_SEM = R2_tmp_std/sqrt(N);
        R2_tmp_CI95 = R2_tmp_mean + R2_tmp_SEM*CI95;
        CI_image = 100*R2_tmp_CI95

        figure;boxplot([test_control(:),test_psycho(:),test_image(:)])
    end
end

%%
% for i = numel( dinfo ):-1:1
%     tmp = control_R2(:,:,i);
%     tmp(~my_mask) = nan;
%     tmpR2Control(:,i) = tmp(~isnan(tmp));
%     tmp = psycho_R2(:,:,i);
%     tmp(~my_mask) = nan;
%     tmpR2Psycho(:,i) = tmp(~isnan(tmp));
%     tmp = image_R2(:,:,i);
%     tmp(~my_mask) = nan;
%     tmpR2Image(:,i) = tmp(~isnan(tmp));
%     tmp = naive_psycho_R2(:,:,i);
%     tmp(~my_mask) = nan;
%     tmpR2PsychoNaive(:,i) = tmp(~isnan(tmp));
%     tmp = naive_image_R2(:,:,i);
%     tmp(~my_mask) = nan;
%     tmpR2ImageNaive(:,i) = tmp(~isnan(tmp));
% end

% save('complete_R2.mat','complete_R2')
% save('control_R2.mat','control_R2')
% save('naive_R2.mat','naive_R2')
% save('psycho_R2.mat','psycho_R2')
% save('image_R2.mat','image_R2')
% save('complete_con.mat','complete_con')
% save('TmMan.mat','TmMan')
% save('TmNat.mat','TmNat')
% save('TmR2.mat','TmR2')
% save('TmR2Naive.mat','TmR2Naive')
% save('TmR2Psycho.mat','TmR2Psycho')
% save('TmR2Image.mat','TmR2Image')
% save('TmCon.mat','TmCon')
% save('TmR2New.mat','TmR2New')
% save('TmR2PsychoNew.mat','TmR2PsychoNew')
% save('TmR2ImageNew.mat','TmR2ImageNew')
% save('mask_psycho_R2.mat','mask_psycho_R2')
% save('mask_complete_con.mat','mask_complete_con')
% save('mask_complete_R2.mat','mask_complete_R2')
% save('mask_control_con.mat','mask_control_con')
% save('mask_image_R2.mat','mask_image_R2')

