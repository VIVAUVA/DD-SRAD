function [ DI_Wv, DI_Ev, DI_Nv, DI_Sv ] = RSS( Image )
% RSS_weighted This is a function that calculates the Root Sum of Squares
% distance between a stack of images

[height,width,Slices]=size(Image); % Images is stack of images

for s = 1:Slices
    
    DI_Wv(:,:,s) = Image(2:height-1,1:width-2,s)-Image(2:height-1,2:width-1,s);
    DI_Ev(:,:,s) = Image(2:height-1,3:width,s)-Image(2:height-1,2:width-1,s);
    DI_Nv(:,:,s) = Image(1:height-2,2:width-1,s)-Image(2:height-1,2:width-1,s);
    DI_Sv(:,:,s) = Image(3:height,2:width-1,s)-Image(2:height-1,2:width-1,s);
   
end

end
