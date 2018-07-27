function [ ks_Sdrive, ks_Edrive, ks_Wdrive, ks_Ndrive ] = ksMap( time_series )
% ks_map Generates edge map based on Kolmogorov-Smirnov distance
%   Takes a stack of of images and generates edge map based on distance
%   between the time series of each pixel
%   Distance taken in the NS and EW directions

[m,n,t] = size(time_series); % dimensions of image stack

% after making 3D matrix of time series, I now want to calculate a ks test 
% map in NS & EW directions
ks_Sdrive = zeros(m-2,n-2); % smaller size than one image
ks_Edrive = zeros(m-2,n-2);
ks_Wdrive = zeros(m-2,n-2);
ks_Ndrive = zeros(m-2,n-2);

[m,n] = size(ks_Sdrive); % dimensions of ks map

for i = 1:m-1
    for j = 1:n-1
        pixel_cur = squeeze(time_series(i,j,:)); % current pixel
        pixel_S = squeeze(time_series(i+1,j,:)); % pixel to the south
        pixel_E = squeeze(time_series(i,j+1,:)); % pixel to the east
        
        
        % comparing pixels below and to the right
        % placing result of test as current pixel value in our ks map
        [h,p,ks_Sdrive(i,j)] = kstest2(pixel_cur,pixel_S);
        [h,p,ks_Edrive(i,j)] = kstest2(pixel_cur,pixel_E);
        
    end
end

for i = 2:m
    for j = 2:n
        pixel_cur = squeeze(time_series(i,j,:)); % current pixel
        pixel_W = squeeze(time_series(i,j-1,:)); % pixel to the west
        pixel_N = squeeze(time_series(i-1,j,:)); % pixel to the north
        
        % comparing pixels above and to the left
        [h,p,ks_Wdrive(i,j)] = kstest2(pixel_cur,pixel_W);
        [h,p,ks_Ndrive(i,j)] = kstest2(pixel_cur,pixel_N);
    end
end

end