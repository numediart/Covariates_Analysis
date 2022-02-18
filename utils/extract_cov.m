function [imageFeatures] = extract_cov(PATH_TO_IMAGES)
%% extract_cov extracts psycholinguistic and image features from the database


%% Covariates on image features
% PATH_TO_IMAGES = 'D:\_ARC-images\ARC tâche Final\ARC tâche Final\300x300';
dinfo = dir(fullfile(PATH_TO_IMAGES,'*.jpg'));
img_name = {dinfo.name};
img_name=img_name(~startsWith(img_name,{'Forme','T','masque'}));
n = numel(img_name);

% define features
entropyValue = NaN(n,1);
contrastValue = NaN(n,1);
correlationValue = NaN(n,1);
energyValue = NaN(n,1);
homogeneityValue = NaN(n,1);
compactness = NaN(n,1);
ratio = NaN(n,1);
numberOfCluster = NaN(n,1);
maxClusterEnergy = NaN(n,1);
maxClusterFreq = NaN(n,1);
maxClusterDist = NaN(n,1);

% im = find(strcmp(img_name,'biberon.jpg'));
for im = 1:n
% Read image
    img = imread(fullfile(PATH_TO_IMAGES,img_name{im}));
    Igray = rgb2gray(img);
    
    % Entropy
    entropyValue(im) = entropy(Igray);
    
    % Other properties
    imProps = graycoprops(Igray,'all');
    %contrast
    contrastValue(im) = imProps.Contrast;
    %correlation
    correlationValue(im) = imProps.Correlation;
    %energy
    energyValue(im) = imProps.Energy;
    %homogeneity
    homogeneityValue(im) = imProps.Homogeneity;
%     continue
    
    % compactness

    [counts,x] = imhist(Igray,256);  
    buffer = 3;
    prevcount = Inf;
    localmin = Inf;
    i = max(x);
    while buffer > 0
        count = counts(i);
        if count > prevcount
            buffer = buffer - 1;
            if localmin > prevcount
                localmin = prevcount;
                thresh = x(i) + 1;
            end
        elseif count < localmin
            buffer = 3;
        end
        prevcount = count;
        i = i-1;
    end
    
    thresh = max(thresh, prctile(x,90));
    bw = Igray < thresh;
    A = bwarea(bw);
    B = bwboundaries(bw);
    [~,idx] = max(cellfun('size',B,1));
    boundary = B{idx};
    P = sum(diag(pdist2(boundary,circshift(boundary,1))));
%     figure('Position', [10 10 1500 800]);
%     plot(boundary(:,2), -boundary(:,1), '-.r', 'LineWidth', 2,'Color',[1, 0, 0, 0.6])
    compactness(im) = 4*pi*A/P^2;


% ration length/width

    [row,col] = find(bw);
    vert = max(row)-min(row) > max(col)-min(col);
    switch(vert)
        case 1
            y1 = min(row);
            x1 = round(mean(find(bw(y1,:))));
            y2 = max(row);
            x2 = round(mean(find(bw(y2,:))));
        case 0
            x1 = min(col);
            y1 = round(mean(find(bw(:,x1))));
            x2 = max(col);
            y2 = round(mean(find(bw(:,x2))));
    end
    
    if x1==x2
        barAngle = 0;
    elseif y1==y2
        barAngle = -90;
    else
        x = 1:1e-3:300;
        y = interp1([x1,x2],[y1,y2],x);
        x(y>300 | y<1 | isnan(y)) = [];
        y(y>300 | y<1 | isnan(y)) = [];
        x = round(x);
        y = round(y);
        orient_mat = zeros(300);
        idx = sub2ind(size(orient_mat), y, x);
        orient_mat(idx) = 1;
    %     figure;imagesc(orient_mat)

        % Perform the Hough transform
        [H, theta, ~] = hough(orient_mat);
        % Find the peak pt in the Hough transform
        peak = houghpeaks(H);
        % Find the angle of the bars
        barAngle = theta(peak(2));
    end
    
    img_rot = imrotate(bw,barAngle,'bilinear','crop');

    % crop
    [row,col] = find(img_rot);
    rect = [min(col) min(row) max(col)-min(col) max(row)-min(row)];
    img_crop = imcrop(img_rot,rect);
%     figure(),imshow(img_crop);
%     pause()
%     close all

% length/width ratio

    [~,idx] = max(size(img_crop));
    ratio(im) = size(img_crop,idx)/size(img_crop,mod(idx,2)+1);
    
    
% frequencial features

    invertedIgray = 255-double(Igray);
    bwIgray=invertedIgray;
    bwIgray(bw) = max(invertedIgray(:));
    
%     figure;imagesc(test2)
    F = abs(fft2(invertedIgray));
    F2 = abs(fft2(bwIgray));
    Fdiff = F-F2;
    Fdiff = max(Fdiff,0);
    FdiffReduced = Fdiff(1:round(size(Fdiff,1)/2),1:round(size(Fdiff,2)/2)); %keep first quadrant
    
    cumulEnergy = 0;
    peakFreq = [];
    i = 1;
    Ftmp = FdiffReduced;
    totalEnergy = sum(Ftmp(:));
    while cumulEnergy < 0.1*totalEnergy
        [val,peakFreq(i)] = max(Ftmp(:));
        cumulEnergy = cumulEnergy + val;
        Ftmp(peakFreq(i)) = 0;
        i = i + 1;
    end
    
    Fdiff10percent = zeros(size(FdiffReduced));
    Fdiff10percent(peakFreq) = 1;
    
    SE = strel('square',ceil(length(Fdiff10percent)/10));
    Fdiff10percentClosed = logical(imclose(Fdiff10percent,SE));

    clc;
    [B,~,numberOfCluster(im),~] = bwboundaries(Fdiff10percentClosed,'noholes');
        
    centroid = cell2mat(struct2cell(regionprops(Fdiff10percentClosed,'Centroid'))');
    in = cell(1,numberOfCluster(im));
    energyByCluster = NaN(numberOfCluster(im),1);
    for i = 1:numberOfCluster(im)
        [row,col] = find(~isnan(Fdiff10percentClosed));
        tmpIn = inpolygon(col,row,B{i}(:,2),B{i}(:,1));
        in{i} = find(tmpIn);
        energyByCluster(i) = sum(FdiffReduced(in{i}));
    end
    [maxClusterEnergy(im),idx] = max(energyByCluster);
    maxClusterFreq(im) = sqrt(centroid(idx,1)^2+centroid(idx,2)^2)/length(FdiffReduced);

    centroidDist = pdist2(centroid,centroid);
    maxClusterDist(im) = max(centroidDist(:));

end

imageFeatures = struct('name',{img_name},'entropyValue',entropyValue,'contrastValue',contrastValue,'correlationValue',correlationValue,...
    'energyValue',energyValue,'homogeneityValue',homogeneityValue,'compactness',compactness,'ratio',ratio,...
    'numberOfCluster',numberOfCluster,'maxClusterEnergy',maxClusterEnergy,'maxClusterFreq',maxClusterFreq,...
    'maxClusterDist',maxClusterDist);

% imageFeatures = struct('name',img_name,'entropyValue',num2cell(entropyValue),'contrastValue',num2cell(contrastValue),...
%     'correlationValue',num2cell(correlationValue),'energyValue',num2cell(energyValue),'homogeneityValue',num2cell(homogeneityValue),...
%     'compactness',num2cell(compactness),'ratio',num2cell(ratio),'numberOfCluster',num2cell(numberOfCluster),...
%     'maxClusterEnergy',num2cell(maxClusterEnergy),'maxClusterFreq',num2cell(maxClusterFreq),'maxClusterDist',num2cell(maxClusterDist));
