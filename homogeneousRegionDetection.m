% Nazia Tabassum 
% 20 July 2018
% Homogeneous Region Detection

function [ PixelList, idx ] = homogeneousRegionDetection( TheStack )

%% take the median

med = median(TheStack,3);

%% edge detection on the median

BW = edge(med,'canny');
% [m,n] = size(BW);

%% filtering

trial = BW;
window_size = 3;
binary_white = @(x) sum(x(:)) == 0; % if any edge pixels, result will
% be 0
homog = nlfilter(trial, [window_size, window_size], binary_white);

%% displaying homogeneous region

CC = bwconncomp(homog);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);

for i = 1:length(numPixels);
    if i ~= idx
        homog(CC.PixelIdxList{i}) = 0;
    end
end

% stats = regionprops('table',homog,'Area','BoundingBox','PixelIdxList');

[B,L] = bwboundaries(homog,'noholes');
imshow(trial);

hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

%% getting the list of pixels in the homogeneous region

PixelList = CC.PixelIdxList;

end