% Analyse Spectral Measurements of Tablet Under Various Lighting Conditions

% The aim here is to understand how much influence the amient surroundings
% have upon the measured (and thus viewed) output of the screen. 

clear, clc, close all


%% Pre run commands:
rootdir = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\Tablet_Characterization';

titles={'WW','CW','MH','WW2'};

for i=1:4
    filename{i}=fullfile(rootdir,titles{i},'White');
end
clear i j

% %% Produce thumbnails for spectra
% %Uncomment saveas for this to actually work
% %Untested save function
% 
% for j=1:length(filename)
%     cd(filename{j});
%     files= dir('*.mat');
%     
%     for file = 1:length(files)
%         load(files(file).name)
%         fig=figure; plot(spc(:,1),spc(:,2));
%         title(files(file).name);
%         
%         iname = fullfile(rootdir,sprintf('%s.tif',files(file).name(1:end-4)));
%         %saveas(fig,strcat(files(file).name(1:end-4),'.tif'))
%         %close
%         
%     end
% end

%% Load all data

PV=0:15:255;
data = zeros(101,length(PV),4); % 4 lighting conditions


for j=1:length(filename)    
    for i = 0:15:255
        load(fullfile(filename{j},sprintf('%03d', i)));
        data(:,i/15+1,j)=spc(:,2);
    end
end

% %Plot data
% figure,
% for j=1:length(filename)
%     subplot(2,2,j), hold on
%     for i = 0:15:255
%         plot(spc(:,1),data(:,i/15+1,j))
%         axis([380,780,0,0.0015]) %Focusing on low levels
%         title(titles{j})
%         %legend('show')
%         %suptitle('Raw data')
%     end
% end

lambda=spc(:,1);

%% Dark Field Correction

%Load all DF (dark field) measurements
for j=1:length(filename)
    cd(filename{j})
    d=dir('DF*mat');
    for i=1:length(d)
        load(fullfile(filename{j},d(i).name));
        DF(:,j,i)=spc(:,2);
    end
end, clear d

%Calculate Average Dark Fields
DF(DF==0)=nan; 
    %Some measurements have more DF measurements than others. 
    %Since this data is saved in a matrix the ones with less 
    %measurements act as though there is an additional zero measurement.
    %including zeros will be bias down the averages
    %so change 0s to NaNs, then use nanmean
DFmean=zeros(101,4);

for i=2:4 %#1 doesn't have DF measurements
    DFmean(:,i)=nanmean(DF(:,i,:),3);
end
DFmean(:,1)=nanmean(DF(:,4,:),3); %steal #4 for #1
%figure, hold on; plot(DFmean); %plot DFCs

data_dfc=data;
for j=1:length(filename) 
    for i=1:size(data,2)
        data_dfc(:,i,j)=data(:,i,j)-DFmean(:,j);        
    end
end

%Plot DFC data
figure,
for j=1:length(filename)
    subplot(2,2,j), hold on
    for i = 0:15:255
        plot(spc(:,1),data_dfc(:,i/15+1,j))
    end
    axis([380,780,0,0.0015])
    title(titles{j})
    legend({'0','15','30','45','60','75','90','105','120','135','150','165','180','195','210','225','240','255'})
end
%suptitle('DFC data')

clear i j

% %% TEMP, kill DFC
% 
% data_dfc=data;

% %% Investigating the difference between WW1 and WW2
% figure, hold on
% 
% for j=1:18
%     for i=[1,4]
%         if i==1
%             plot(data_dfc(:,j,i),'r')
%         elseif i==4
%             plot(data_dfc(:,j,i),'b')
%         end
%     end
% end

%% Calculate Chromaticities

ciefile = fullfile('C:','Users','cege-user','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

cie2= xlsread(ciefile,'1931 col observer','A6:D86');
CIEwavelength=cie2(:,1);
xbar=cie2(:,2);
ybar=cie2(:,3);
zbar=cie2(:,4);

%Calculate XYZ then convert to u'v'

data_dfc_int=zeros(length(cie2),18,4);
XYZ=zeros(3,18,4); xy=zeros(2,18,4); uv=zeros(2,18,4);

for i=1:4
    for j=1:18
        data_dfc_int(:,j,i)=interp1...
            (lambda,...
            data_dfc(:,j,i),...
            CIEwavelength,'spline');
        
        XYZ(1,j,i)=xbar'*data_dfc_int(:,j,i);
        XYZ(2,j,i)=ybar'*data_dfc_int(:,j,i);
        XYZ(3,j,i)=zbar'*data_dfc_int(:,j,i);
        
        xy(1,j,i)=XYZ(1,j,i)/sum(XYZ(:,j,i));
        xy(2,j,i)=XYZ(2,j,i)/sum(XYZ(:,j,i));
        
        uv(1,j,i)=4*XYZ(1,j,i)/(XYZ(1,j,i)+15*XYZ(2,j,i)+3*XYZ(3,j,i));
        uv(2,j,i)=9*XYZ(2,j,i)/(XYZ(1,j,i)+15*XYZ(2,j,i)+3*XYZ(3,j,i));
        
    end
end

%save('C:\Users\cege-user\Dropbox\Documents\MATLAB\DFC_Y.mat')

%% Plot chromaticities
figure, hold on

% ubar=4.*xbar ./ (xbar + 15.*ybar + 3.*zbar);
% vbar=9.*ybar ./ (xbar + 15.*ybar + 3.*zbar);
% plot(ubar,vbar)

% %2D scatter, size scaled
% for i=1:4
%     scatter(uv(1,:,i),uv(2,:,i),[1:18]*4,'filled')
% end

% %3D scatter
% for i=1:4
%     scatter3(uv(1,:,i),uv(2,:,i),XYZ(2,:,i),'filled')
% end
% 
% xlabel('u'''),ylabel('v'''),axis equal

% %2D with annotations
% annot={'0','15','30','45','60','75','90','105','120','135','150','165','180','195','210','225','240','255'};
% for i=1:4
%     scatter(uv(1,:,i),uv(2,:,i),'filled')
%     maxTxtVl=60;
%     text(uv(1,1:(maxTxtVl/15+1),i)+0.001,uv(2,1:(maxTxtVl/15+1),i),annot(1:(maxTxtVl/15+1)));
% end
% xlabel('u'''),ylabel('v'''),zlabel('Y'),axis equal

%2D, just above 60
annot={'0','15','30','45','60','75','90','105','120','135','150','165','180','195','210','225','240','255'};
for i=1:4
    scatter(uv(1,6:end,i),uv(2,6:end,i),'filled')
    minTxtVl=210;
    maxTxtVl=255;
    text(uv(1,(minTxtVl/15+1):(maxTxtVl/15+1),i)+0.001,uv(2,(minTxtVl/15+1):(maxTxtVl/15+1),i),annot((minTxtVl/15+1):(maxTxtVl/15+1)));
end
xlabel('u'''),ylabel('v'''),zlabel('Y'),axis equal

legend(titles,'Location','Best')

%% u'v' distance from average

uv_top_mean=mean(uv(:,end-5:end,:),2);
scatter(uv_top_mean(1,:,:),uv_top_mean(2,:,:),'k','filled')

uv_top_mean2=mean(mean(uv(:,end-5:end,:),2),3);
scatter(uv_top_mean2(1,:,:),uv_top_mean2(2,:,:),'r','filled')

distance_from_mean=sqrt((uv(1,:,:)-uv_top_mean(1)).^2+((uv(2,:,:)-uv_top_mean(2)).^2));

figure, hold on
for i=1:4
    plot(PV,distance_from_mean(:,:,i))
end
    
xlabel('PV')
ylabel('Delta u''v''')


%% Histograms for each channel
%clear

I=imread('C:\Users\cege-user\Dropbox\Documents\PsychoPy\zazzle_luv_60_50_8bit_BM.tif');
R=imhist(I(:,:,1));
G=imhist(I(:,:,2));
B=imhist(I(:,:,3));
figure, plot(R,'r')
hold on, plot(G,'g')
plot(B,'b'), legend(' Red channel','Green channel','Blue channel');

%% Plot gamut

clear, clc, close all 

rootdir = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\Tablet_Characterization';

gamut_data = zeros(101,3); %101 spectral bins, 3=RGB
chan={'Red','Green','Blue','White'};

for i=1:length(chan)
    gamut_filename{i}=fullfile(rootdir,'WW',chan{i});
end

lambda=[380:4:780]'; %lambda=spc(:,1);

% Create and save spectral plots for quick identification
for i=1:length(chan)
    cd(gamut_filename{i});
    clear files
    files= dir('*.mat');
    
    for file = 1:length(files)
        load(files(file).name)
        %         % Create and save spectral plots for quick identification
        %         fig=figure; plot(spc(:,1),spc(:,2));
        %         title(files(file).name);
        %
        %         iname = fullfile(rootdir,sprintf('%s.tif',files(file).name(1:end-4)));
        %         %saveas(fig,strcat(files(file).name(1:end-4),'.tif'))
        %         %drawnow
        gamut_data(:,i,file)=spc(:,2);
    end
    %     %clear spc
end

gamut_data=gamut_data(:,:,1:18);

% figure, hold on
% for i=1:size(gamut_data,3)
%     plot(lambda,gamut_data(:,1,i),'r')
%     plot(lambda,gamut_data(:,2,i),'g')
%     plot(lambda,gamut_data(:,3,i),'b')
%     %plot(lambda,gamut_data(:,4,i),'k')
%     ylim([0, max(gamut_data(:,3,i))])
%     title(files(i).name)
%     drawnow
%     pause(0.02)
%     %if i ~= size(gamut_data,3); cla; end
% end

% Load CIE 2deg obs

ciefile = fullfile('C:','Users','cege-user','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');
cie2= xlsread(ciefile,'1931 col observer','A6:D86');
CIEwavelength=cie2(:,1);
xbar=cie2(:,2);
ybar=cie2(:,3);
zbar=cie2(:,4);

% Calculate XYZ then convert to u'v'

gamut_data_int=zeros(length(cie2),length(chan),18);
gamut_XYZ=zeros(3,length(chan),18); 
gamut_xy=zeros(2,length(chan),18); 
gamut_uv=zeros(2,length(chan),18);

for i=1:length(chan)
    for d=1:18 %d for drive
    gamut_data_int(:,i,d)=interp1...
        (lambda,...
        gamut_data(:,i,d),...
        CIEwavelength,'spline');
    
    gamut_XYZ(1,i,d)=xbar'*gamut_data_int(:,i,d);
    gamut_XYZ(2,i,d)=ybar'*gamut_data_int(:,i,d);
    gamut_XYZ(3,i,d)=zbar'*gamut_data_int(:,i,d);
    
    gamut_xy(1,i,d)=gamut_XYZ(1,i,d)/sum(gamut_XYZ(:,i,d));
    gamut_xy(2,i,d)=gamut_XYZ(2,i,d)/sum(gamut_XYZ(:,i,d));
    
    gamut_uv(1,i,d)=4*gamut_XYZ(1,i,d)/(gamut_XYZ(1,i,d)+15*gamut_XYZ(2,i,d)+3*gamut_XYZ(3,i,d));
    gamut_uv(2,i,d)=9*gamut_XYZ(2,i,d)/(gamut_XYZ(1,i,d)+15*gamut_XYZ(2,i,d)+3*gamut_XYZ(3,i,d));
    end    
    
end

% Plot
figure('Position',[[100,100], [500,400]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'defaultAxesFontName', 'Courier',...
    'Renderer','Painters',...
    'color','white'); 
hold on

ubar=4.*xbar ./ (xbar + 15.*ybar + 3.*zbar);
vbar=9.*ybar ./ (xbar + 15.*ybar + 3.*zbar);

sRGB_dcs = XYZToSRGBPrimary([xbar,ybar,zbar]'); %sRGB display colours
sRGB_dcs(sRGB_dcs>1) = 1;
sRGB_dcs(sRGB_dcs<0) = 0;
for i=1:3
    for j=1:size(sRGB_dcs,2)-1
        t(i,j) = (sRGB_dcs(i,j)+sRGB_dcs(i,j+1))/2;
    end    
end
sRGB_dcs = t;

for i = 1:size(xbar,1)-1
    plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i));
end
plt(1) = plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i),'DisplayName','Spectral Locus');

for d=1:18
    scatter3(gamut_uv(1,:,d),gamut_uv(2,:,d),gamut_XYZ(2,:,d),[],[1,0,0;0,1,0;0,0,1;0,0,0],'.')
    % Plot lines between gamut points:
%     if d==18
%         plot3([gamut_uv(1,1:3,d),gamut_uv(1,1,d)],[gamut_uv(2,1:3,d),gamut_uv(2,1,d)],...
%             [gamut_XYZ(2,1:3,d),gamut_XYZ(2,1,d)],'k','DisplayName','Device Gamut')
%     end
end

plt(2) = plot([gamut_uv(1,1:3,18),gamut_uv(1,1,18)],[gamut_uv(2,1:3,18),gamut_uv(2,1,18)],'k--','DisplayName','Device Gamut');

xlabel('u'''),ylabel('v'''),zlabel('Y')
axis equal
xlim([0 0.65])
ylim([0 0.65])

% Add sRGB
sRGB_R_xy = [0.6400,0.3300];
sRGB_G_xy = [0.3000,0.6000];
sRGB_B_xy = [0.1500,0.0600];

sRGB_R_uv = xyTouv(sRGB_R_xy'); %Uses PsychToolbox function
sRGB_G_uv = xyTouv(sRGB_G_xy');
sRGB_B_uv = xyTouv(sRGB_B_xy');

% scatter(sRGB_R_uv(1),sRGB_R_uv(2),'r*');
% scatter(sRGB_G_uv(1),sRGB_G_uv(2),'g*');
% scatter(sRGB_B_uv(1),sRGB_B_uv(2),'b*');

plt(3) = plot([sRGB_R_uv(1),sRGB_G_uv(1)],[sRGB_R_uv(2),sRGB_G_uv(2)],'k:','DisplayName','sRGB Gamut');
plot([sRGB_G_uv(1),sRGB_B_uv(1)],[sRGB_G_uv(2),sRGB_B_uv(2)],'k:')
plot([sRGB_B_uv(1),sRGB_R_uv(1)],[sRGB_B_uv(2),sRGB_R_uv(2)],'k:')
%save('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\SAPS_TabletGamut.mat')

legend(plt,'Location','southeast')

cleanTicks
%save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\gamut.pdf')


%% Plot 'representative gamut' (the gamut of the points actually displayed)
%Requires previous section to have been run prior 
close all
for kspec=1%:30
    
    %clearvars -except kspec rootdir chan gamut_filename files file gamut_data lambda ciefile cie2 CIEwavelength xbar ybar zbar
    clc, 
    %close all
        
    % 1. Load the screenshots
    % 2. Flatten the screenshot RGBs into one long 3 channel vector (n,3)
    %   Take the (768,1366,3,30) and make it (n,3)
    % 3. Convert the RGBs via lookup to spectra
    % 4. to xy and then u'v'
    
    cd(fullfile(rootdir,'Screenshots of PAMELA SAPS','Second Attempt'))
    files= dir('*.png'); 
        
    if ~exist('kspec','var')
        % Old version for bringing in multiple images at once:
        % 1. Load the screenshots
        for i=1:30
            images(:,:,:,i)=imread(files(i).name);
            %figure,imshow(a(:,:,:,i))
        end
        
        %Split by channel
        aR2(:,:,:)=images(:,:,1,:);
        aG2(:,:,:)=images(:,:,2,:);
        aB2(:,:,:)=images(:,:,3,:);
    else
        image=imread(files(kspec).name);
        aR2=image(:,:,1);
        aG2=image(:,:,2);
        aB2=image(:,:,3);
        
    end
    
    % 2. Flatten the screenshot RGBs into one long 3 channel vector (n,3)
    i2=[aR2(:),aG2(:),aB2(:)]; 
    clear iR2 iG2 iB2 i %tidy up
    
    % New code for creating a 4D matrix with a cube for each image, so that
    % I can see how many images each chromaticity is in.
    cube=zeros(256,256,256,30);
    for i=1:length(i2)
        cube(i2(i,1)+1,i2(i,2)+1,i2(i,3)+1,ceil(i/1049088))=...
        cube(i2(i,1)+1,i2(i,2)+1,i2(i,3)+1,ceil(i/1049088))+1;
    end
    
    %This code is really really slow. For now I'll save and load cubennzsum
    %but I really should try to find a better way to do it.
    %load('cubennzsum.mat')
    cubennzsum= zeros(256,256,256);
    for i=1:256
        for j=1:256
            for k=1:256
                cubennzsum(i,j,k)=nnz(cube(i,j,k,:));
                
            end
        end
        i
    end
    
    %create a cube and add each time a specific triplet is seen to that area
    %(like a 3D histogram)
    cube=zeros(256,256,256);
    for i=1:length(i2)
        cube(i2(i,1)+1,i2(i,2)+1,i2(i,3)+1)=...
        cube(i2(i,1)+1,i2(i,2)+1,i2(i,3)+1)+1;
    end
    
    %Calculate spectrum for each pixel value
    gamut_data_int=zeros(101,4,256);
    for i=1:4
        for j=1:101
            gamut_data_int(j,i,:)=interp1(0:15:255,squeeze(gamut_data(j,i,:)),0:255);
        end
    end
    
    %Check the interp1 worked as planned
    % figure, hold on
    % for i=1:size(gamut_data_int,3)
    %     plot(lambda,gamut_data_int(:,1,i),'r')
    %     plot(lambda,gamut_data_int(:,2,i),'g')
    %     plot(lambda,gamut_data_int(:,3,i),'b')
    %     plot(lambda,gamut_data_int(:,4,i),'k')
    %     ylim([0, max(gamut_data_int(:,3,i))])
    %     drawnow
    %     %if i ~= size(gamut_data,3); cla; end
    % end
    
    xbar2_101=interp1(CIEwavelength,xbar,380:4:780,'spline');
    ybar2_101=interp1(CIEwavelength,ybar,380:4:780,'spline');
    zbar2_101=interp1(CIEwavelength,zbar,380:4:780,'spline');
    
    % Calculate for each channel (RGB) the XYZ tristimulus values at
    % each pixel value (0:255)
    
    % I think that the following general procedure for calculating colorimetry
    % should be checked against 'Measuring Colour' and 'Computational colour
    % science using MATLAB'. For now, this will do.
    
    XYZ_LUT=zeros(4,3,256);
    xy_LUT=zeros(4,2,256);
    uv_LUT=zeros(4,2,256);
    
    for i=1:4
        for j=1:256
            XYZ_LUT(i,1,j)=xbar2_101*gamut_data_int(:,i,j);
            XYZ_LUT(i,2,j)=ybar2_101*gamut_data_int(:,i,j);
            XYZ_LUT(i,3,j)=zbar2_101*gamut_data_int(:,i,j);
            
            xy_LUT(i,1,j)=XYZ_LUT(i,1,j)/sum(XYZ_LUT(i,:,j));
            xy_LUT(i,2,j)=XYZ_LUT(i,2,j)/sum(XYZ_LUT(i,:,j));
            
            uv_LUT(i,1,j)=4*XYZ_LUT(i,1,j)/(XYZ_LUT(i,1,j)+15*XYZ_LUT(i,2,j)+3*XYZ_LUT(i,3,j));
            uv_LUT(i,2,j)=9*XYZ_LUT(i,2,j)/(XYZ_LUT(i,1,j)+15*XYZ_LUT(i,2,j)+3*XYZ_LUT(i,3,j));
        end
    end
    
    %Check (note that these are chromaticities based on interpolated spectra,
    %which I think explains why there are some funky bits where it seems to
    %change direction for a small period of time: this really is just one
    %different recorded value that'd being interpolated to make it look like
    %many more points.)
    
    figure, hold on
    scatter3(uv_LUT(1,1,1:15:256),uv_LUT(1,2,1:15:256),XYZ_LUT(1,2,1:15:256),'r*');
    scatter3(uv_LUT(2,1,1:15:256),uv_LUT(2,2,1:15:256),XYZ_LUT(2,2,1:15:256),'g*');
    scatter3(uv_LUT(3,1,1:15:256),uv_LUT(3,2,1:15:256),XYZ_LUT(3,2,1:15:256),'b*');
    scatter3(uv_LUT(4,1,1:15:256),uv_LUT(4,2,1:15:256),XYZ_LUT(4,2,1:15:256),'k*');
    %xlim([0.1 0.4]);
    %ylim([0.25,0.6]);
    
    % %Normalise cube (frequency of occurance of pixel value: n/30).
    %cube2=cubennzsum/max(max(max(max(cubennzsum))));
    cube2=cubennzsum;
    
    %This is based on the total number of pixels: n/numel(images(:,:,:,1))
    % %Normalise cube (frequency of occurance of pixel value.
%     cube2=cube/max(max(max(cube)));
    
%     occ_NonNormalised=cube(:);
%     [ans1,ans2]=max(occ==1);
%     occ_NonNormalised(ans2-100:ans2+100);
    
    cube_XYZ_LUT=zeros(256,256,256,3);
    cube_u_LUT=zeros(256,256,256);
    cube_v_LUT=zeros(256,256,256);
    
    for i=1:256
        for j=1:256
            for k=1:256
                if cube2(i,j,k)~=0
                    cube_XYZ_LUT(i,j,k,:)=XYZ_LUT(1,:,i)+XYZ_LUT(2,:,j)+XYZ_LUT(3,:,k);
                    cube_u_LUT(i,j,k)=4*cube_XYZ_LUT(i,j,k,1)/(cube_XYZ_LUT(i,j,k,1)+15*cube_XYZ_LUT(i,j,k,2)+3*cube_XYZ_LUT(i,j,k,3));
                    cube_v_LUT(i,j,k)=9*cube_XYZ_LUT(i,j,k,2)/(cube_XYZ_LUT(i,j,k,1)+15*cube_XYZ_LUT(i,j,k,2)+3*cube_XYZ_LUT(i,j,k,3));
                end
            end
        end
    end
    
    %figure, hold on
    
    ubar=4.*xbar ./ (xbar + 15.*ybar + 3.*zbar);
    vbar=9.*ybar ./ (xbar + 15.*ybar + 3.*zbar);
    plot(ubar,vbar,'k')
  
    xlabel('u'''),ylabel('v'''),zlabel('Y')
    axis equal
    
    u=cube_u_LUT(:);
    v=cube_v_LUT(:);
    occ=cube2(:);
    %occ=cube2(:); %occurance
    Y=squeeze(cube_XYZ_LUT(:,:,:,2)); Y=Y(:);
    
    %scatter(u,v,'k.');
    %scatter(u(1:25:end),v(1:25:end),20,occ(1:25:end),'filled')
    %scatter(0.1915,0.4717,'k*')
    scatter3(u(1:50:end),v(1:50:end),Y(1:50:end),20,occ(1:50:end),'filled')
    %scatter3(u(1:end),v(1:end),Y(1:end),20,occ(1:end),'filled')
    
    % xlim([0.1482 0.2403])
    % ylim([0.4217 0.5007])
    % zlim([0.046 0.051])
end
colorbar

%% %Plot single stimulus

scatter(u(1:end),v(1:end),'k','filled')
xlabel('u''')
ylabel('v''')
%axis([0.14 0.24, 0.41,0.51])
axis equal
grid on
xlim([0.14 0.25]),ylim([0.41 0.52])

%% Focus on the first stimulus
xlim([0.1482 0.2403])
ylim([0.4217 0.5007])

%% Save/load data and plot fresh
% Quickfire way to avoid all of the code above

%save('SAPS_SelectableColoursGamut.mat','u','v','Y','occ')
clc, clear, close all

% figure('Position',[[100,100], [500,309]],...
%     'defaultLineLineWidth',2,...
%     'defaultAxesFontSize',12,...
%     'defaultAxesFontName', 'Courier',...
%     'Renderer','Painters',...
%     'color','white'); 
% hold on

 load('SAPS_SelectableColoursGamut.mat')
% s=scatter3(u(1:50:end),v(1:50:end),Y(1:50:end),20,occ(1:50:end),'filled');
% view(2)
% axis('equal')
% daspect([1,1,0.1]) %This isn't exactly what I want to do but it'll suffice
% xlim([0.14 0.25]),ylim([0.41 0.52]),zlim([0.046 0.06])
% colorbar
% xlabel('u'''),ylabel('v'''),zlabel('Y') 
% 
% cleanTicks
%save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\practical_gamut.pdf')


% with gamut:

load('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\SAPS_TabletGamut.mat')

figure('Position',[[100,100], [500,400]],...
    'defaultLineLineWidth',2,...
    'defaultAxesFontSize',12,...
    'defaultAxesFontName', 'Courier',...
    'Renderer','Painters',...
    'color','white'); 
hold on

ubar=4.*xbar ./ (xbar + 15.*ybar + 3.*zbar);
vbar=9.*ybar ./ (xbar + 15.*ybar + 3.*zbar);

sRGB_dcs = XYZToSRGBPrimary([xbar,ybar,zbar]'); %sRGB display colours
sRGB_dcs(sRGB_dcs>1) = 1;
sRGB_dcs(sRGB_dcs<0) = 0;
for i=1:3
    for j=1:size(sRGB_dcs,2)-1
        t(i,j) = (sRGB_dcs(i,j)+sRGB_dcs(i,j+1))/2;
    end    
end
sRGB_dcs = t;

for i = 1:size(xbar,1)-1
    plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i));
end

% for d=1:18
%     scatter3(gamut_uv(1,:,d),gamut_uv(2,:,d),gamut_XYZ(2,:,d),[],[1,0,0;0,1,0;0,0,1;0,0,0],'.')
%     % Plot lines between gamut points:
%     if d==18
%         plot3([gamut_uv(1,1:3,d),gamut_uv(1,1,d)],[gamut_uv(2,1:3,d),gamut_uv(2,1,d)],...
%             [gamut_XYZ(2,1:3,d),gamut_XYZ(2,1,d)],'k:')
%     end
% end
%s=scatter3(u(u~=0),v(u~=0),Y(u~=0),10,occ(u~=0),'filled');
s=scatter(u(u~=0),v(u~=0),2,occ(u~=0),'filled');
plot([gamut_uv(1,1:3,18),gamut_uv(1,1,18)],[gamut_uv(2,1:3,18),gamut_uv(2,1,18)],'k--','DisplayName','Device Gamut');

xlabel('u'''),ylabel('v'''),zlabel('Y')
axis equal
xlim([0 0.65])
ylim([0 0.65])

cleanTicks
colorbar
save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\practical_gamut.pdf')
