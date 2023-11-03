% plot explained variance with ROIs (patch)
% unthreshMap = squeeze(one_sample(:,:,1));
% figure;imagesc(unthreshMap)
% 
% simpleMap = unthreshMap; %last_simple_R2
% 
% imageMap = unthreshMap; %new_image_R2
% imageMap(imageMap<0) = 0;
% 
% psychoMap = unthreshMap; %last_last_psycho_R2
% psychoMap(psychoMap<0.01) = 0;

% % R2 psycho-image model
% test = unthreshMap*100;
% test(unthreshMap<0) = 0;
% figure;imagesc(test)
% load('D:\__EEG-data\BIDS_files\masks\mask_r2_complete.mat')
% hold on
% for i = 1:length(unique(mask_r2_complete))-1
%     [~,col] = find(mask_r2_complete == i);
%     patch_start = min(col);
%     patch_end = max(col);
%     rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[1 0 0 0.4], 'EdgeColor', [1, 0, 0, 0])
% end
% title('Explained Variance of Psycho-Image Model','FontSize',14)

figure;
k=1;
myMask = {0,mask_psycho_R2,mask_image_R2};
myTitle = {'Explained Variance of Categorical Model','Explained Variance of Psycho Model','Explained Variance of Image Model'};
for myMap = {simpleMap, psychoMap, imageMap}
    
    subplot(1,3,k);imagesc(myMap{1}*100)
    caxis([0 4])
    title(myTitle{k},'FontSize',11)

    % significant regions

    hold on    
    for i = 1:length(unique(mask_simple_con1))-1
        [~,col] = find(mask_simple_con1 == i);
        patch_start = min(col);
        patch_end = max(col);
        rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[0 1 0 0.4], 'EdgeColor', [1, 0, 0, 0])
    end

%     mask_psycho_image_con1(:,300:end) = 0;
    for i = 1:length(unique(mask_psycho_image_con1))-1
        [~,col] = find(mask_psycho_image_con1 == i);
        patch_start = min(col);
        patch_end = max(col);
        rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[1 1 0 0.6], 'EdgeColor', [1, 0, 0, 0])
    end
    
    if k>1
        for i = 1:length(unique(myMask{k}))-1
            [~,col] = find(myMask{k} == i);
            if i==1 && k==2
                patch_start = min(col)+25;
            else
                patch_start = min(col);
            end
            patch_end = max(col);
            rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[1 0 0 0.5], 'EdgeColor', [1, 0, 0, 0])
        end
    end

    x = linspace(-200,500,8);
    xticks(1:50:351)
    xticklabels(num2str(x'))
    xlabel('time(ms)', 'FontSize', 10)
    ylabel('channels', 'FontSize', 10)

    k=k+1;
end