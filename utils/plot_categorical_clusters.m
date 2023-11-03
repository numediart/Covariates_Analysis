% plot categorical clusters
% categorical model: nat_man_simple_con_1
% naive model: last_psycho_image_naive_con_1
% complete model: last_psycho_image_con_1

% Find the masks
% limo_results
% mask_psycho_image_con1 = mask;
% mask_simple_con1 = mask;
% mask_simple_con1(ismember(tmp,[1,2])) = 3;
% mask_psycho_image_naive_con1 = mask;
% save('D:\__EEG-data\BIDS_files\masks\mask_psycho_image_con1.mat','mask_psycho_image_con1')
% save('D:\__EEG-data\BIDS_files\masks\mask_simple_con1.mat','mask_simple_con1')
% save('D:\__EEG-data\BIDS_files\masks\mask_psycho_image_naive_con1.mat','mask_psycho_image_naive_con1')

myMask = {mask_simple_con1,mask_psycho_image_con1,mask_psycho_image_naive_con1};
myTitle = {'Categorical Model','Psycho-Image Model','Naive Psycho-Image Model'};
figure;
cmap_pos = hot(64);
cmap_neg = flipud(cmap_pos);
cmap = [cmap_neg; cmap_pos];

i=1;
for myFile = {'nat_man_simple_con_1','last_psycho_image_con_1','last_psycho_image_naive_con_1'}
%     load('one_sample_ttest_parameter_1.mat')
    load(['D:\__EEG-data\BIDS_files\' myFile{1} '\one_sample_ttest_parameter_1.mat'])
    unthreshMap = squeeze(one_sample(:,:,1));
    mask = myMask{i};
    unthreshMap(~mask) = nan;
%     figure;imagesc(unthreshMap);colorbar;
    subplot(1,3,i)
    imagesc(unthreshMap*100,'AlphaData',mask);%colorbar;
    colormap(cmap)
%     colormap hot
%     set(gca,'color',1*[.93 .93 .93]);
    set(gca,'color',1*[.9 .9 .9]);
%     caxis([0 100*max(unthreshMap(:))])
    caxis([-180 180])
    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
    xlabel('time(ms)', 'FontSize', 14)
    ylabel('channels', 'FontSize', 14)

    title(myTitle{i},'FontSize',14)
    i = i+1;
end