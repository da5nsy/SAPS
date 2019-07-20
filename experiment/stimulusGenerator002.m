%% LAB

stimsize=2188;
stim=zeros(stimsize,stimsize,3);

LABstim=stim;
LABstim(:,:,1)=60; %L
for i=1:stimsize %a and b
    LABstim(i,:,2)=linspace(-50,50,stimsize);
    LABstim(:,i,3)=linspace(-50,50,stimsize);
end
clear i stim 
 
%% LAB --> XYZ
%XYZstim=lab2xyz(LABstim,'WhitePoint','d65');

XYZstim=LABstim;

Xn=357.882;
Yn=389.386;
Zn=432.084;

for i=1:stimsize 
    for j=1:stimsize

        L = LABstim(i,j,1);
        a = LABstim(i,j,2);
        b = LABstim(i,j,3);

        Y_Yn = (L+16) / 116;
        X_Xn = (a/500) + Y_Yn; 
        Z_Zn = Y_Yn -(b/200);

        X=Xn*(X_Xn^3);
        Y=Yn*(Y_Yn^3);
        Z=Zn*(Z_Zn^3);
        
        XYZstim(i,j,:)=[X,Y,Z];
    end
end
% XYZstim_X=XYZstim(:,:,1);
% XYZstim_Y=XYZstim(:,:,2);
% XYZstim_Z=XYZstim(:,:,3);

clear X Y Z L a b Y_Yn X_Xn Z_Zn Xn Yn Zn i j 
%imshow(xyz2rgb(XYZstim))
%% XYZ to R'G'B'

M = [127.86,155.33,79.68;
    75.50,259.18,58.32;
    12.17,46.44,377.15];

M_in=inv(M);

RGBstim=XYZstim;

%[XYZstim(1,1,1);XYZstim(1,1,2);XYZstim(1,1,3)]

for i=1:stimsize
    for j=1:stimsize
        RGBstim(i,j,:)=M_in*[XYZstim(i,j,1);XYZstim(i,j,2);XYZstim(i,j,3)];
    end
end

RGBstim_R=RGBstim(:,:,1);
RGBstim_G=RGBstim(:,:,2);
RGBstim_B=RGBstim(:,:,3);

% RGBstim255=uint8(RGBstim.*255);
%imshow(RGBstim255)

% RGBstim255_R=RGBstim255(:,:,1);
% RGBstim255_G=RGBstim255(:,:,2);
% RGBstim255_B=RGBstim255(:,:,3);

clear M M_in i j

%% R'G'B' --> RGB (Linearization

% Generate LUT
red=[0.81	0.98	1.27	2.08	3.25	5.09	7.54	10.67	14.36	18.81	23.65	29.32	35.15	41.54	48.32	56.06	64.44	75.56];
red=red./max(red);
green=[0.98	1.38	2.79	5.52	9.99	16.44	25.33	36.59	50	65.11	81.72	100.13	120	141.64	163.34	189.38	222.15	258.75];
green=green./max(green);
blue=[1.1	1.09	1.39	1.99	3.01	4.35	6.17	8.68	11.67	15.19	18.87	22.87	27.6	32.41	37.62	42.86	49.39	58.29];
blue=blue./max(blue);

x = 0:1/17:1;
yr = red;
yg = green;
yb = blue;
xx = 0:1/9999:1;
yyr = spline(x,yr,xx);
yyg = spline(x,yg,xx);
yyb = spline(x,yb,xx);

LUT=zeros(4,length(xx)); %allocate memory
LUT(1,:)=xx; %linear space

for i=1:length(xx)
    [~, index] = min(abs(yyr-(i/length(xx)))); 
   LUT(2,i)= index/length(xx);
end

for i=1:length(xx)
    [~, index] = min(abs(yyg-(i/length(xx))));
   LUT(3,i)= index/length(xx);
end

for i=1:length(xx)
    [~, index] = min(abs(yyb-(i/length(xx))));
   LUT(4,i)= index/length(xx);
end

LUT(2,1:150)=0;
LUT(3,1:100)=0;
LUT(4,1:200)=0;

RGBstim2=RGBstim;

for i=1:stimsize
    for j=1:stimsize
        
        [~, index] = min(abs(RGBstim(i,j,1)-LUT(1,:)));
        RGBstim2(i,j,1)=LUT(2,index);
        
        [~, index] = min(abs(RGBstim(i,j,2)-LUT(1,:)));
        RGBstim2(i,j,2)=LUT(3,index);
        
        [~, index] = min(abs(RGBstim(i,j,3)-LUT(1,:)));
        RGBstim2(i,j,3)=LUT(4,index);
    end
    disp(i)
end
%%
% RGBstim255=RGBstim.*255;
% imshow(RGBstim255)

RGBstim255=uint8(RGBstim2.*255);
figure,imshow(RGBstim255)

imwrite(RGBstim255,'zazzle_60_50_8bit.tif','compression','none')

%% output string 

out=zeros(stimsize*stimsize*3+10,1);

out=RGBstim(:);
csvwrite('STIM.csv',out)

%%
subplot(2,2,1)
imshow(lab2rgb(LABstim))
title('LAB')

subplot(2,2,2)
imshow(xyz2rgb(XYZstim))
title('XYZ')

subplot(2,2,3)
imshow(RGBstim255)
title('RGB')

subplot(2,2,4)
plot(1:10,1:10)
title('Fourth subplot')