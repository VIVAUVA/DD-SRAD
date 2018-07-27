function [ bhatta_Sdrive, bhatta_Edrive ] = bhattacharyyaMap( time_series )
%bhattacharyya_map Calculates Bhattacharyya distance between each pixel in
%the time series
%   Takes a stack of of images and generates edge map based on distance
%   coefficient between the time series of each pixel
%   Distance taken in the NS and EW directions

[m,n,t] = size(time_series); % dimensions of image stack

% after making 3D matrix of time series, I now want to calculate a bhattacharyya 
% map in NS & EW directions
bhatta_Sdrive = zeros(m,n); % same size as one image
bhatta_Edrive = zeros(m,n);

for i = 1:m-1
    for j = 1:n-1
        pixel_cur = squeeze(time_series(i,j,:)); % current pixel
        pixel_S = squeeze(time_series(i+1,j,:)); % pixel to the south
        pixel_E = squeeze(time_series(i,j+1,:)); % pixel to the east
        
        % comparing pixels below and to the right
        % placing result of test as current pixel value in our ks map
        bhatta_Sdrive(i,j) = bhattacharyya(pixel_cur,pixel_S);
        bhatta_Edrive(i,j) = bhattacharyya(pixel_cur,pixel_E);
    end
end

end
