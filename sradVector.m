%*****************************************************************************************
%  This function performes speckle filtering using SRAD, Speckle Reducing Anisotropic Diffusion
%
%  With modified the morrior boundary conditions
%  Y. Yu
%  June 28,2002.
% Modified
% Vector version
% Scott Acton
% August 2015
% Added distances
% Nazia Tabassum
% July 2018
%****************************************************************************************

function [Image, East, South]=sradVector(Image,PixelList,idx,distanceFlag)

[height,width,Slices]=size(Image); % Image is stack of images

mn0=mean2(Image(:,:,1));

I = Image(:,:,1);
variance = var(I(PixelList{idx}));
average = mean(I(PixelList{idx}));
constant = 1;
Cu2 = constant*variance/(average^2);

if distanceFlag == 0
    [ DI_Wv, DI_Ev, DI_Nv, DI_Sv ] = RSS( Image );
elseif distanceFlag == 1
    % calculate ks maps for stack
    [ks_Sdrive, ks_Edrive, ks_Wdrive, ks_Ndrive] = ksMap(Image);
else
    % bhatta
    [bhatta_Sdrive, bhatta_Edrive] = bhattacharyyaMap(Image);
end

for s = 1:Slices
    
    notanumber=isnan(Image(:,:,s));
    Img = Image(:,:,s);
    Img(notanumber)=0;
    Image(:,:,s) = Img;
    
    Image1=zeros(height,width);
    Image_RD=zeros(height-2,width-2);
    
    Image_RD(1:height-2,1:width-2)=Image(2:height-1,2:width-1,s);
    
    % Compute four first-order differences
    DI_W=Image(2:height-1,1:width-2,s)-Image(2:height-1,2:width-1,s);
    
    DI_E=Image(2:height-1,3:width,s)-Image(2:height-1,2:width-1,s);
    
    DI_N=Image(1:height-2,2:width-1,s)-Image(2:height-1,2:width-1,s);
    
    DI_S=Image(3:height,2:width-1,s)-Image(2:height-1,2:width-1,s);
    
    if distanceFlag == 0
        
        DI_Wv_med = median(DI_Wv,3);
        DI_Ev_med = median(DI_Ev,3);
        DI_Nv_med = median(DI_Nv,3);
        DI_Sv_med = median(DI_Sv,3);
        
        rootsumsquaresW = sqrt(sum(DI_Wv_med.*DI_Wv_med,3));
        rootsumsquaresE = sqrt(sum(DI_Ev_med.*DI_Ev_med,3));
        rootsumsquaresN = sqrt(sum(DI_Nv_med.*DI_Nv_med,3));
        rootsumsquaresS = sqrt(sum(DI_Sv_med.*DI_Sv_med,3));
        
        % Compute diffusion coefficient
        Laplace=(rootsumsquaresW+rootsumsquaresE+rootsumsquaresN+rootsumsquaresS)./Image_RD;
        magrad=(rootsumsquaresN.^2+rootsumsquaresS.^2+rootsumsquaresW.^2+rootsumsquaresE.^2)./Image_RD.^2;
        
    elseif distanceFlag == 1
        
        Laplace=(ks_Wdrive+ks_Edrive+ks_Ndrive+ks_Sdrive)./Image_RD;
        magrad=(ks_Ndrive.^2+ks_Sdrive.^2+ks_Wdrive.^2+ks_Edrive.^2)./Image_RD.^2;
    else
        
        bhatta_Sdrive = bhatta_Sdrive(1:height-2,1:width-2);
        bhatta_Edrive = bhatta_Edrive(1:height-2,1:width-2);
        Laplace = (bhatta_Sdrive + bhatta_Edrive)./Image_RD;
        magrad = (bhatta_Sdrive.^2 + bhatta_Edrive.^2)./Image_RD.^2;
    end
  
    Ci2=abs(0.25*magrad-(Laplace.^2)/17)./(1+0.25*Laplace).^2;

    D=(Ci2-Cu2)./(Cu2+1);

    D=D/max(max(D));
    ci=1./(1+D./Cu2);
    
    ciS=zeros(size(ci));
    ciE=ciS;
    ciS(1:height-3,1:width-2)=ci(2:height-2,1:width-2);
    ciS(height-2,1:width-2)=ciS(height-3,1:width-2);
    ciE(1:height-2,1:width-3)=ci(1:height-2,2:width-2);
    ciE(1:height-2,width-2)=ciE(1:height-2,width-3);
    
    minc=min(min(ciS));
    maxc=max(max(ciS));
    maxc=maxc-minc;
    ciS=(ciS-minc)/maxc;
    
    minc=min(min(ciE));
    maxc=max(max(ciE));
    maxc=maxc-minc;
    ciE=(ciE-minc)/maxc;
    
    minc=min(min(ci));
    maxc=max(max(ci));
    maxc=maxc-minc;
    ci=(ci-minc)/maxc;
    
    if s == 1,
        East = ciE;
        South = ciS;
    end
    
    %Ye old diffusion equation
    Image_RD1=Image_RD+.05*(1/4)*(ci.*(DI_N+DI_W)+ciS.*DI_S+ciE.*DI_E);
    
    Pos=find(Image_RD1<=0|ci==0);
    Image_RD1(Pos)=Image_RD(Pos);
    
    % Treat the boundaries of the image
    Image1(2:height-1,2:width-1)=Image_RD1(1:height-2,1:width-2);
    Image1(1,:)=Image1(2,:);
    Image1(height,:)=Image1(height-1,:);
    Image1(:,1)=Image1(:,2);
    Image1(:,width)=Image1(:,width-1);
    
    Image(:,:,s)=Image1;
    
    % mn=mean2(Image(:,:,s));
    
    notanumber=isnan(Image(:,:,s));
    Img = Image(:,:,s);
    Img(notanumber)=0;
    Image(:,:,s) = Img;
    
    % Image(:,:,s)=Image(:,:,s)/mn*mn0; % normalizing
    Image(:,:,s)=Image(:,:,s)/mn0; % normalizing
    
end

return;