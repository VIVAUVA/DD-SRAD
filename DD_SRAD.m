% Scott Acton
% for fun!
% August 2015
% Nazia Tabassum
% July 2018

clear;
clc;

%% read in the images from a dataset (3D stack)

% only use with RSS distance
load('SyntheticSquare.mat');
TheStack = SyntheticSquare + 100; % so that output doesn't have NaN
% load('realData.mat');

startNum=0;
[m,n,slices] = size(TheStack);

I = TheStack(:,:,1);
Image=double((I));
[height, width]=size(Image);
Image=Image*255/max(max(Image));
Ilog=log(Image+1);
Ilog=(Ilog-min(min(Ilog)))*255/max(max(Ilog-min(min(Ilog))));

%% homogeneous region detection

[PixelList,idx] = homogeneousRegionDetection(TheStack);

%% set flags for which distance you'd like to use

% 0 = RSS
% 1 = KS
% 2 = Bhattacharyya
distanceFlag = 1;

% 0 = RSS unweighted
% 1 = RSS weighted
weightingFlag = 0; 

%% iterations

close;

% can iterate more or less depending on how smooth you want your result
for round=1:60,
    
    % printing iteration number
    disp('Iteration');
    disp(round);
    
    % deciding whether to use weighted approach or not
    if weightingFlag == 1 && distanceFlag == 0;
        [outimage, CE, CS]=sradVectorWeights(TheStack,PixelList,idx);
    else
        [outimage, CE, CS]=sradVector(TheStack,PixelList,idx,distanceFlag);
    end
    
    TheStack=outimage;
    C=(CE/max(max(CE))+CS/max(max(CS)))/2;

    outimage=outimage(:,:,end);
    
    outimage255=(outimage-min(min(outimage)))*255/(max(max(outimage))-min(min(outimage)));

    C=1-C;
    g=graythresh(C);
    Cmap=C>g;
    OIlog=log(outimage255+1);
    OIlog=(OIlog-min(min(OIlog)))*255/(max(max(OIlog))-min(min(OIlog)));
    C255=(C-min(min(C)))*255/(max(max(C))-min(min(C)));
    
    BW = edge(OIlog,'Canny',.25,1);
    figure(startNum+1),
    hold on;
    subplot(1,3,1),
    imagesc(Ilog),colormap(gray),axis image, axis off, axis tight; 
    title('Original Image');
    subplot(1,3,2),
    imagesc(OIlog),colormap(gray),axis image, axis off, axis tight; 
    title('Denoised Result');
    subplot(1,3,3),
    imagesc(BW),colormap(gray),axis image, axis off, axis tight; 
    title('Edgemap');
    drawnow
    
end

hold off;
drawnow

disp('Finished');