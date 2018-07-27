%*****************************************************************************************
%  This function performes speckle filtering using SRAD, Speckle Reducing Anisotropic Diffusion
%
%  With modified mirror boundary conditions
%  Y. Yu
%  June 28,2002.
% Modified
% Vector version
% Scott Acton
% August 2015
% Modified weighted version
% Nazia Tabassum
% July 2018
%****************************************************************************************

function [Image, East, South]=sradVectorWeights(Image,PixelList,idx)

[height,width,Slices]=size(Image); % Image is stack of images

mn0=mean2(Image(:,:,1));

I = Image(:,:,1);
variance = var(I(PixelList{idx}));
average = mean(I(PixelList{idx}));
constant = 1;
Cu2 = constant*variance/(average^2);

g = gausswin(2*Slices); % gaussian used for weighting images

[ DI_Wv, DI_Ev, DI_Nv, DI_Sv ] = RSS( Image );

for s = 1:Slices
    
    notanumber=isnan(Image(:,:,s));
    Img = Image(:,:,s);
    Img(notanumber)=0;
    Image(:,:,s) = Img;

    Image1=zeros(height,width);
    Image_RD=zeros(height-2,width-2);
    % Image_RD1=zeros(height-2,width-2);
    % D=zeros(height-2,width-2);
    
    % Compute four first-order differences
    DI_W=Image(2:height-1,1:width-2,s)-Image(2:height-1,2:width-1,s);
    
    DI_E=Image(2:height-1,3:width,s)-Image(2:height-1,2:width-1,s);
    
    DI_N=Image(1:height-2,2:width-1,s)-Image(2:height-1,2:width-1,s);
    
    DI_S=Image(3:height,2:width-1,s)-Image(2:height-1,2:width-1,s);
    
    [m,n,l] = size(DI_Wv);
    
    DI_Wv_avg = zeros(m,n);
    DI_Ev_avg = zeros(m,n);
    DI_Nv_avg = zeros(m,n);
    DI_Sv_avg = zeros(m,n);
    
    % weighting
    start = (Slices+1)-(s-1);
    stop = start + (Slices-1);
    g_subsec = g(start:stop);
    
    for i = 1:Slices
        DI_Wv_avg = DI_Wv_avg + DI_Wv(:,:,i)*g_subsec(i);
        DI_Ev_avg = DI_Ev_avg + DI_Ev(:,:,i)*g_subsec(i);
        DI_Nv_avg = DI_Nv_avg + DI_Nv(:,:,i)*g_subsec(i);
        DI_Sv_avg = DI_Sv_avg + DI_Sv(:,:,i)*g_subsec(i);
    end
    
    DI_Wv_avg = DI_Wv_avg/sum(g_subsec);
    
    rootsumsquaresW = sqrt(sum(DI_Wv_avg.*DI_Wv_avg,3));
    rootsumsquaresE = sqrt(sum(DI_Ev_avg.*DI_Ev_avg,3));
    rootsumsquaresN = sqrt(sum(DI_Nv_avg.*DI_Nv_avg,3));
    rootsumsquaresS = sqrt(sum(DI_Sv_avg.*DI_Sv_avg,3));

    Image_RD(1:height-2,1:width-2) = Image(2:height-1,2:width-1,s);
    
    % Compute diffusion coefficient
    Laplace=(rootsumsquaresW+rootsumsquaresE+rootsumsquaresN+rootsumsquaresS)./Image_RD;
    magrad=(rootsumsquaresN.^2+rootsumsquaresS.^2+rootsumsquaresW.^2+rootsumsquaresE.^2)./Image_RD.^2;
  
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

    I = Image(:,:,s);
    variance = var(I(PixelList{idx}));
    average = mean(I(PixelList{idx}));
    Cu2 = variance/(average^2);
    
    Image(:,:,s)=Image1;
    
    notanumber=isnan(Image(:,:,s));
    Img = Image(:,:,s);
    Img(notanumber)=0;
    Image(:,:,s) = Img;
    
    mn=mean2(Image(:,:,s));
    
    Image(:,:,s)=Image(:,:,s)/mn*mn0;
    
end

return;