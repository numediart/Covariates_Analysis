% plot categorical clusters with ROIs (patch)

myMask = mask_psycho_image_con1;
myTitle = 'Categorical Contrast from Psycho-Image Model';
figure;
cmap_pos = parula(64);
cmap_neg = flipud(cmap_pos);
cmap = [cmap_neg; cmap_pos];

myFile = 'last_psycho_image_con_1';
load(['D:\__EEG-data\BIDS_files\' myFile '\one_sample_ttest_parameter_1.mat'])
unthreshMap = squeeze(one_sample(:,:,1));
imagesc(unthreshMap);colorbar;
colormap(cmap)
caxis([-2 2])
x = linspace(-200,500,8);
xticks(1:50:351)
xticklabels(num2str(x'))
xlabel('time(ms)', 'FontSize', 14)
ylabel('channels', 'FontSize', 14)
title(myTitle,'FontSize',14)

hold on
for i = 1:length(unique(mask_r2_complete))-1
    [~,col] = find(mask_r2_complete == i);
    patch_start = min(col);
    patch_end = max(col);
    rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[1 0 0 0.4], 'EdgeColor', [1, 0, 0, 0])
end

mask_psycho_image_con1(:,300:end) = 0;
for i = 1:length(unique(mask_psycho_image_con1))-1
    [~,col] = find(mask_psycho_image_con1 == i);
    patch_start = min(col);
    patch_end = max(col);
    rectangle('Position',[patch_start 0 patch_end-patch_start 65],'FaceColor',[0 1 0 0.4], 'EdgeColor', [1, 0, 0, 0])
end
