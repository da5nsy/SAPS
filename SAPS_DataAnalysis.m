clc, clear, close all

%General purpose script assimilating previous advances in different files
%Many pieces cut, see prev versions for missing sections

% To do list

% Add Grant data (different run length, should just be able to change n
% Replace dist function in BM plotting
% Replace full paths with relative ones using rootdir and something like
% 'fig_output_folder'

set(groot,'defaultfigureposition',[100 100 500 400]); 
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesFontName', 'Courier');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vectors
set(groot,'defaultfigurecolor','white')

saveFigs = 0;

%% Load Data

%-% Tablet Data
% (Gone back to having 2 xlsread calls, despite being slower, because having
% both raw and numerical data is useful, but they require different excel
% addresses to specify data)

location='BM';  %changing this changes EVERYTHING
%'PAMELA' or 'GRANT' or 'BM' or '200s' or
%'PAMELA_20180205' or 'basement_rgby_test' or
%'basement_greyCard_test'

%Warning: As far as I recall, for the 'GRANT' and '200s' data, a different
%stimulus, defined in CIELAB was used, and no accounting for this is made
%in the following analysis. If required, contact DG for access to previous
%versions (pre-git) of this code.

rootdir = fullfile('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data',location);
cd(rootdir)

try
    load(sprintf('TabletData_%s',location))
catch
    disp('Reading from excel. Wait 42 seconds (PAMELA), and 171 seconds for BM')
    tic
    TBfilename = [rootdir,'\',location,'.xlsx'];
    [~,sheets] = xlsfinfo(TBfilename);
    if strcmp(location,'200s')
        TabletData = zeros(200,8,numel(sheets));
    else
        TabletData = zeros(40,8,numel(sheets));
    end
    
    for i=1:numel(sheets) %65 for BM, 20 for PAMELA
        disp(i) %to provide a progress update in command window, should take ~2mins
        if strcmp(location,'200s')
            TabletData(:,:,i)=xlsread(TBfilename, i,'A2:H201');
        else
            TabletData(:,:,i)=xlsread(TBfilename, i,'A2:H41');
        end
        
        date=str2num(sheets{i}(6:end));
        [t_stamp(i).y,t_stamp(i).m,t_stamp(i).d...
            ,t_stamp(i).h,t_stamp(i).mi,t_stamp(i).s] ...
            = datevec(datenum([1970 1 1 0 0 date]));
        
        files(i).Time=t_stamp(i);
        files(i).filename=cellstr(sheets(i));
        files(i).data=TabletData(:,:,i) ;
        
        if ~strcmp(location,'200s')  %Ignore extra parameters for 200s for now, build in as needed
            [~,~,raw{i}]=xlsread(TBfilename, i,'A2:H49');
            files(i).extraInfo=raw{1,i}{42,2};
            files(i).checker_switch=raw{1,i}{43,2};
            files(i).repeats=raw{1,i}{44,2};
            files(i).participant=raw{1,i}{45,2};
            files(i).Gallery=raw{1,i}{46,2};
            files(i).Light=raw{1,i}{48,2};
        end
    end
    toc
    
    clear i
    save(sprintf('TabletData_%s',location))
end

if strcmp(location,'BM') %correct for poor data collection
    files(33).participant = 'DG';
    files(14).participant = 'DG';
end

% Set up baseline dummy data
TabletData(:,:,end+1)=TabletData(:,:,end);
TabletData(:,7:8,end)=zeros;

%-% Supplementary variables
n=size(TabletData,1); %how many trials
n2=size(TabletData,3);

files(length(files)+1)=files(length(files));
for i=1:n2 %add date info onto top level for readability (has to come here)
    files(i).date=files(i).Time.d;
end

files(length(files)).date=0; %Set the final dataset (dummy data) date to 0 so that it can be specially treated when plotting
files(length(files)).participant='dummy'; %As above

%% Processing


whiteXYZ=[360.476; 394.674; 416.714];

L_star=65;

% Get white point u and v
uv0 = XYZTouv(whiteXYZ);

for i=1:n2
    
    %Calculate systematic offset (spatial calibration)
    avXdp(i) = median(TabletData(end-9:end,7,i)-TabletData(end-9:end,1,i));
    avYdp(i) = median(TabletData(end-9:end,8,i)-TabletData(end-9:end,2,i));
    
    Xdp(:,i) = TabletData(end-9:end,7,i)-TabletData(end-9:end,1,i);
    Ydp(:,i) = TabletData(end-9:end,8,i)-TabletData(end-9:end,2,i);
    
    avXdp(n2)=0;
    avYdp(n2)=0;
    
    %Adjust data by spatial calibration, and recorded offsets
    TabletData(:,7,i)=TabletData(:,7,i)-avXdp(i)-TabletData(:,1,i);
    TabletData(:,8,i)=TabletData(:,8,i)-avYdp(i)-TabletData(:,2,i);
    
    %Recalculate direction and magnitude (d=4, m=5)
    for j=1:n
        TabletData(j,4,i)=(atan2(TabletData(j,8,i),TabletData(j,7,i)))*...
            (180/pi)+TabletData(j,3,i);
        TabletData(j,5,i)= sqrt(TabletData(j,7,i)^2+TabletData(j,8,i)^2);
    end
    
    %Calculate spatial x and spatial y from mag/dir
    TabletData(:,9,i) = (TabletData(:,5,i)).*cosd(TabletData(:,4,i)); %x
    TabletData(:,10,i) = (TabletData(:,5,i)).*sind(TabletData(:,4,i)); %y
    
    % spatial xy --> u*v*
    % (the figure 21.88 relates the pixel dimensions of the stimulus to the 
    % colorimetric definition of that stimulus)
    TabletData(:,9,i)   =   (TabletData(:,9,i)./21.88); %u* 
    TabletData(:,10,i)  =   (TabletData(:,10,i)./21.88); %v*
    
    % u*v* --> u'v'
    TabletData(:,11,i) = (TabletData(:,9,i) ./ (13.0 * L_star) ) + uv0(1); %u'
    TabletData(:,12,i) = (TabletData(:,10,i)./ (13.0 * L_star) ) + uv0(2); %v'
    
end

%% Time analysis

timeA = 0; %time analysis
if timeA    
    count = 1;
    for i=1:length(raw)
        if ~strcmp(files(i).participant,'dud')
            if ~strcmp(files(i).participant,'dummy')
                if and(strcmp(location,'BM'),i~=40) %a specific dud from the BM set
                    time(:,count) = (cell2mat(raw{1,i}(1:40,6)));
                    count = count+1;
                end
            end
        end
    end
    
    figure,
    %plot(time) %plots everything - messy
    plot(nanmedian(time,2),'k-o') %takes the median, ignoring NaNs (which occur when it times out because someone took too long)
    ylim([0,10])
    xlabel('Trial number')
    ylabel('Median response time (seconds)')
    
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\mediantime.pdf')
    end
    
    figure,
    %hist(sum(time(2:end,:)),50) %seconds
    histogram(sum(time(2:end,:)/60),50,'FaceColor','k') %minutes
    xlabel('Total time taken (minutes)')
    %yticks([0 max(ylim)])
    ylabel('number of observers')
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\timeHistogram.pdf')
    end
end

%% Plot u'v'

ciefile = fullfile('C:','Users','cege-user','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

cie2= xlsread(ciefile,'1931 col observer','A6:D86');

lambdaCie2=cie2(:,1);
xbar2=cie2(:,2);
ybar2=cie2(:,3);
zbar2=cie2(:,4);

figure
hold on
axis equal

ubar=4.*xbar2 ./ (xbar2 + 15.*ybar2 + 3.*zbar2);
vbar=9.*ybar2 ./ (xbar2 + 15.*ybar2 + 3.*zbar2);

sRGB_dcs = XYZToSRGBPrimary([xbar2,ybar2,zbar2]'); %sRGB display colours
sRGB_dcs(sRGB_dcs>1) = 1;
sRGB_dcs(sRGB_dcs<0) = 0;
for i=1:3
    for j=1:size(sRGB_dcs,2)-1
        t(i,j) = (sRGB_dcs(i,j)+sRGB_dcs(i,j+1))/2;
    end    
end
sRGB_dcs = t;

for i = 1:size(xbar2,1)-1
    plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i));
end
plt(1) = plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i),'DisplayName','Spectral Locus'); %Repeat just for legened


%% Pull light measurements

if strcmp(location,'BM')
    
    clear out
    
    %cd('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\BM\GL'); %old measurements
    cd('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\BM\GL 20180427') %new measurements
    mmg=dir('*.mmg');
    mmgl=length(mmg);
    
    out=zeros(256,mmgl); % Pull in data, store to 'out'
    
    %Grab individual datasets and store in out
    for j=1:mmgl
        xFile = xmlread(mmg(j).name);
        for k=0:255
            out(k+1,j+1)=str2double(xFile.getElementsByTagName('row').item(k).getAttribute('value'));
        end
    end
    
    % Grab wavelength data
    for i=0:255
        out(i+1,1)=str2double(xFile.getElementsByTagName('row').item(i).getAttribute('wavelength'));
    end
    
    spd=out(28:end,:);
    
    %for new data, where this is 10 of each datasets
    spd_all=spd; %backup all data before averaging
    spd(:,1)=spd(:,1);
    spd(:,2)=mean(spd(:,2:11),2);
    spd(:,3)=mean(spd(:,12:21),2);
    spd(:,4)=mean(spd(:,22:end),2);
    spd=spd(:,1:4);
    
    plt_spd = 0;
    if plt_spd
        figure, hold on
        for i=[4,3,2]
            plot(spd(:,1),spd(:,i)/max(spd(:,i)))
        end
        
        legend({'Rooms 77/78','Room 25','Great Court'},'Location','northoutside')
        
        xlabel('Wavelength (nm)')
        ylabel('Normalised power')
        axis tight
        ylim([0 1])
        yticks(ylim)
    end
    
    xbar2_GL=interp1(lambdaCie2,xbar2,spd(:,1),'spline');
    ybar2_GL=interp1(lambdaCie2,ybar2,spd(:,1),'spline');
    zbar2_GL=interp1(lambdaCie2,zbar2,spd(:,1),'spline');
    
    GLXYZ=[xbar2_GL,ybar2_GL,zbar2_GL]'*spd(:,2:4);
    GLxy=[GLXYZ(1,:)./sum(GLXYZ);GLXYZ(2,:)./sum(GLXYZ)];
    GLuv=[4*GLxy(1,:)./(-2*GLxy(1,:)+12*GLxy(2,:)+3);9*GLxy(2,:)./(-2*GLxy(1,:)+12*GLxy(2,:)+3)];
    
    %scatter(spd_uv(:,1),spd_uv(:,2),'*') %all
    plt(2) = scatter(GLuv(1,1),GLuv(2,1),'*','DisplayName','Great Court');
    plt(3) = scatter(GLuv(1,2),GLuv(2,2),'*','DisplayName','Room 25');
    plt(4) = scatter(GLuv(1,3),GLuv(2,3),'*','DisplayName','Room 77/78');
    
    %text(spd_uv(:,1)'+([1:size(spd_uv,1)]/300),spd_uv(:,2)'+([1:size(spd_uv,1)]/300),string([1:size(spd_uv,1)]))
    
    axis equal
    xlim([0 0.65])
    ylim([0 0.65])
    cleanTicks
    legend(plt,'Location','southeast')
    xlabel('u'''),ylabel('v''')
    if saveFigs
    save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\BM_chromaticities.pdf')
    end
    
    %     for i=2:mmgl+1
    %         figure
    %         plot(out(:,1),out(:,i))
    %         title(mmg(i-1).name)
    %         drawnow
    %         %pause(0.5)
    %         %saveas(gcf,strcat(mmg(i-1).name(1:end-4),'.tif'))
    %     end
    
    %     for i=2:4
    %         figure
    %         plot(spd(:,1),spd(:,i))
    %     end
    
    %figure, scatter(GLuv(1,:),GLuv(2,:),'k*');
    %figure,
    
    %     scatter(flip(GLuv(1,:)),flip(GLuv(2,:)),'k*');
    %     text(flip(GLuv(1,:)),flip(GLuv(2,:)),{'1','2','3'})
    
    %save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\BM_SPD.pdf')
end

if strcmp(location,'PAMELA')
    try
        load('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\20170331 Spectra\PAMELA_SPD_data.mat')
    catch
        [spd_data,~,~,spd_uv] = read_UPRtek('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\20170331 Spectra',0,0,0);
        save('PAMELA_SPD_data.mat','spd_data','spd_uv')
    end    
    
    %scatter(spd_uv(:,1),spd_uv(:,2),'*') %all
    plt(2) = scatter(spd_uv([1,5:7],1),spd_uv([1,5:7],2),'*','DisplayName','CW');
    plt(3) = scatter(spd_uv([2:4],1),spd_uv([2:4],2),'*','DisplayName','WW');
    plt(4) = scatter(spd_uv([11:18],1),spd_uv([11:18],2),'*','DisplayName','MH');
    
    %text(spd_uv(:,1)'+([1:size(spd_uv,1)]/300),spd_uv(:,2)'+([1:size(spd_uv,1)]/300),string([1:size(spd_uv,1)]))
    
    axis equal
    xlim([0 0.65])
    ylim([0 0.65])
    cleanTicks
    legend(plt,'Location','southeast')
    xlabel('u'''),ylabel('v''')
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_chromaticities.pdf')
    end
        
    % [1,5:7]   == CW
    % [2:4]     == WW
    % [11:18]   == MH   
    
    spd_uv_CW=mean(spd_uv([1,5:7],:));
    spd_uv_WW=mean(spd_uv(2:4,:));
    spd_uv_MH=mean(spd_uv(11:18,:));
    
    % Note that the following pulls data from a specified measurement
    % session, not neccessarily the one specified at the start of this
    % script
    
    plt_ind_spd = 1; %plot individual SPD
    if plt_ind_spd
        try
            load('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\2017 Spectra\PAMELA_SPD_data.mat')
        catch
            [spd_data,peak,lux,spd_uv]=read_UPRtek('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\2017 Spectra',0,0,0);
            save('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\2017 Spectra\PAMELA_SPD_data.mat','spd_data','spd_uv')
        end
        
        figure, hold on
        h= plot(360:760,spd_data(:,[1,3:8]));
        set(h, {'color'}, {'k';'r';'g';'b';[0.9,0.9,0];'k';'k'});
        set(h, {'LineStyle'}, {'-';'-';'-';'-';'-';':';'--'});
        legend({'High output','Red','Green','Blue','Amber','Warm White','Cool White'},'Location','northoutside')
        
        xlabel('Wavelength (nm)')
        ylabel('Normalised power')
        axis tight
        ylim([0 1])
        yticks(ylim)
        %title('Normalised SPDs of PAMELA lighting channels')
        
        if saveFigs
            save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_SPD.pdf')
        end
    end
end

% Abandoned this as it already exists in
% CalculateChromaticitiesOfPamelaLighting.m
%
% if strcmp(location,'PAMELA_20180205')
%     try
%         load('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\20180205 Spectra\PAMELA_ChromData.mat')
%     catch
%         [spectral_data,peak,lux_fromExcel,spd_uv,spd_xy] = read_UPRtek('C:\Users\cege-user\Dropbox\Documents\MATLAB\SAPS\data\PAMELA\20180205 Spectra',0,0,0);
%         S_spectral_data = [360,1,401];
%         save('PAMELA_ChromData.mat','spectral_data','peak','lux_fromExcel','spd_uv','S_spectral_data','spd_xy')
%     end
% end

%% Assess variance

for i=1:n2
    kstd_u(i)=nanstd(TabletData(1:end-10,11,i));
    kstd_v(i)=nanstd(TabletData(1:end-10,12,i));
    Mu_ind_u(i)=nanmean(TabletData(1:end-10,11,i));
    Mu_ind_v(i)=nanmean(TabletData(1:end-10,12,i));
end

kstd_mean=mean([kstd_u;kstd_v]);
kstd_min=min([kstd_u;kstd_v]);
diffs=squeeze(sqrt((TabletData(3,11,:)-TabletData(8,11,:)).^2+(TabletData(3,12,:)-TabletData(8,12,:)).^2))';

if strcmp(location,'PAMELA')
    %calc means for DG data for specific lighting conditions
    Mu_u_CW=mean(Mu_ind_u([1,2,6,8,15,17]));
    Mu_u_WW=mean(Mu_ind_u([3,5,9,10]));
    Mu_u_MH=mean(Mu_ind_u([12,14,19]));
    Mu_v_CW=mean(Mu_ind_v([1,2,6,8,15,17]));
    Mu_v_WW=mean(Mu_ind_v([3,5,9,10]));
    Mu_v_MH=mean(Mu_ind_v([12,14,19]));
end

if strcmp(location,'PAMELA_20180205')
    figure, hold on
    %scatter(1,kstd_mean(1),'k','filled');
    plot([1,9],[kstd_mean(1),kstd_mean(1)],'k--','DisplayName','Real touch baseline data');
    plot([1,9],[kstd_mean(end),kstd_mean(end)],'k:','DisplayName','Computed baseline data');
    sc(1)=scatter(1:9,kstd_mean(2:10),'r','filled','DisplayName','MH1');
    sc(2)=scatter(1:9,kstd_mean(11:19),'g','filled','DisplayName','ML');
    sc(3)=scatter(1:9,kstd_mean(20:28),'b','filled','DisplayName','MH2');
    xlim([0 10])
    xticks(1:9)
    xlabel('Observer')
    ylabel('Mean SD in u'' and v''')
    ylim([0, 0.02])
    yticks([0, 0.02])
    % %Replace numbers with participant identifiers
    %xticklabels({[],files([2:10]).participant})
    legend
    
    % %Plot DBUR
    %     figure, hold on
    %     scatter(1:9,diffs(2:10),'r','filled','DisplayName','MH1');
    %     scatter(1:9,diffs(11:19),'g','filled','DisplayName','ML');
    %     scatter(1:9,diffs(20:28),'b','filled','DisplayName','MH2');
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\exp3excl.pdf')
    end
end


MinOrMean = 'Mean'; %'Min' or 'Mean'
thresh_SD = 0.01;
thresh_DBUR = 0.018;

if strcmp(location,'BM')
    if strcmp(MinOrMean,'Mean')
        kstd=kstd_mean;
    else
        kstd=kstd_min;
    end
    %figure
    %bar([kstd_mean;diffs]')
    
    figure, hold on
    %scatter(kstd_mean,diffs); %plot all data, incl LM data, which is not
    %included below
    
    ylim([0 0.07])
    
    scatter(kstd(5:28),diffs(5:28),'r','filled','DisplayName','Room 77/78');
    scatter(kstd(29:54),diffs(29:54),'g','filled','DisplayName','Room 25');
    scatter(kstd(55:end),diffs(55:end),'b','filled','DisplayName','Great Court');
    
    %scatter(kstd_mean(4),diffs(4),'k','filled','DisplayName','known dud');
    scatter(kstd(66),diffs(66),'y','filled','DisplayName','Baseline Data');
    
    %plot([kstd(end),kstd(end)],[0,0.07],'k:','DisplayName','Logical Threshold')
    
    plot([thresh_SD,thresh_SD],[0,0.07],'k:','DisplayName','SD > 0.01')
    plot([min(xlim),max(xlim)],[thresh_DBUR,thresh_DBUR],'k:','DisplayName','DBUR > 0.018')
    
    f(1)=fill([min(xlim),max(xlim),max(xlim),min(xlim)],[thresh_DBUR,thresh_DBUR,max(ylim),max(ylim)],'k','LineStyle','none','FaceAlpha','0.1','DisplayName','Excluded Data');
    f(2)=fill([thresh_SD,max(xlim),max(xlim),thresh_SD],[min(ylim),min(ylim),max(ylim),max(ylim)],'k','LineStyle','none','FaceAlpha','0.1','DisplayName','Excluded Data');
    
    if strcmp(MinOrMean,'Min')
        xlabel('Min SD')
    else
        xlabel('Mean SD')
    end
    ylabel('DBUR')  
    
    %set( get( get( f(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    %set( get( get( f(2), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    legend('Location','east')
    
    %cleanTicks
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\excl3.pdf')
    end

    plt_outL = 0;
    if plt_outL
        figure, hold on
    end
    clear sdcoll
    for i=1:n2
        if strcmp(files(i).participant, 'Public') %only include public (also excludes dummy)
            if avXdp(i) < 100 % exclude ones which were clearly duds
                if ~exist('sdcoll','var'),sdcoll=[];end %collect sd
                %sdcoll=[sdcoll; nanstd(TabletData(end-9:end,11,i)),nanstd((TabletData(end-9:end,12,i)))];
                sdcoll=[sdcoll; nanstd(Xdp(:,i)),nanstd(Ydp(:,i))];
                if plt_outL
                    if sdcoll(end,1)>15 %annotate outliers
                        text(sdcoll(end,1),sdcoll(end,2),string(i),'FontSize',14)
                    end
                end
            end
        end
    end
    if plt_outL
        scatter(sdcoll(:,1),sdcoll(:,2))
        axis equal
        xlim([0 inf])
        ylim([0 inf])
        xlabel('SD in x axis during checkerboard (pixels)')
        ylabel('SD in y axis during checkerboard (pixels)')
    end
    
    figure, hold on  %scatter excluding outliers
    scatter(sdcoll(sdcoll(:,1)<15,1),sdcoll(sdcoll(:,1)<15,2),'k','filled','MarkerFaceAlpha',.5)
    axis equal
    xlim([0 12])
    ylim([0 12])
    %title('(Excluding 4 outliers)')
    xlabel({'SD in x axis', 'during checkerboard (pixels)'})
    ylabel({'SD in y axis', 'during checkerboard (pixels)'})    
    plot([0,12],[0,12],'k--')

    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\BMtouch.pdf')
    end
end

%%
if strcmp(location,'BM')
    figure, hold on
    for i=1:n2
        if strcmp(files(i).Gallery,'Africa')
            if strcmp(files(i).participant,'LM')
                h1=scatter(kstd_u(i),kstd_v(i),'r','filled','DisplayName','LM');
            elseif strcmp(files(i).participant,'DG')
                h2=scatter(kstd_u(i),kstd_v(i),'g','filled','DisplayName','DG');
            elseif strcmp(files(i).participant,'Public')                
                h3=scatter(kstd_u(i),kstd_v(i),'b','filled','DisplayName','Public');
            end
        end
    end
    legend([h1 h2 h3])
    axis equal
    xlim([0 0.015])
    ylim([0 0.015])
    xlabel('SD in u'' dimension')
    ylabel('SD in v'' dimension')
end

cleanTicks
if saveFigs
    save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\naive.pdf')
end

% hold on
% scatter(kstd([10,18,54,48]),diffs([10,18,54,48]),'k', 'filled')

% clc %print out the ones which are thrown up as anomolies above
% Xdp(:,[10,18,54,48])
% Ydp(:,[10,18,54,48])

% figure, hold on 
% for i=1:n2
%     scatter(Xdp(:,i),Ydp(:,i),'filled')
%     text(Xdp(end,i),Ydp(end,i),string(i))
% end
% legend

% %% Analyse repeat values
%
% fig=figure('Position', [50, 50, 800, 800]); hold on;
% for i = 1:size(xbar2,1)-1
%     plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i));
% end
% axis equal;
%
% load('SAPS_SelectableColoursGamut.mat')
% s=scatter(u(1:50:end),v(1:50:end),20,occ(1:50:end));
% view(2)
% axis('equal')
% xlim([0.14 0.25]),ylim([0.41 0.52]) %close to selectable gamut boundary
% colorbar
% xlabel('u'''),ylabel('v'''),zlabel('Y')
%
% repeats=[3,8];
% kstd_cutoff=0.012;
% kstd_ind=kstd_mean<0.012;
%
% for i=1:n2
%     plot(TabletData(repeats,11,i),TabletData(repeats,12,i),'k')
%     scatter(TabletData(repeats(2),11,i),TabletData(repeats(2),12,i),kstd_mean(i)*10000,'k','filled')
%     %text(TabletData(repeats(2),11,i)+0.003,TabletData(repeats(2),12,i),files(i).participant)
%     if kstd_ind(i)
%         scatter(TabletData(repeats(2),11,i),TabletData(repeats(2),12,i),kstd_mean(i)*10000,'r','filled')
%     end
% end
%
% % For PAMELA data
% %for i=2:10
% %for i=11:19
% %for i=20:28
%
% % % Add Baseline data (essentially zero apart from the difference from
% % % input variability)
% % plot(TabletData(repeats,11,1),TabletData(repeats,12,1),'r')
% % scatter(TabletData(repeats(2),11,1),TabletData(repeats(2),12,1),'r')

%% Plot

figure, hold on

for i = 1:size(xbar2,1)-1
    plot([ubar(i),ubar(i+1)],[vbar(i),vbar(i+1)],'Color',sRGB_dcs(:,i));
end

%scatter(GLuv(1,:),GLuv(2,:),'k*');

load('SAPS_SelectableColoursGamut.mat')
s=scatter(u(1:50:end),v(1:50:end),20,occ(1:50:end),'.');
view(2)
axis equal
xlim([0.14 0.25]),ylim([0.41 0.52]) %close to selectable gamut boundary
colorbar
xlabel('u'''),ylabel('v'''),zlabel('Y')

% from http://stackoverflow.com/questions/3417028/ellipse-around-the-data-in-matlab

G = 1*ones(n-10,1);
for i=2:n2
    G = [G ; i*ones(n-10,1)];
end

X=reshape(permute([ TabletData(1:end-10,11,:), ...
    TabletData(1:end-10,12,:)],[1 3 2])...
    ,[(n-10)*n2,2]);

%PlotZeros = 1;
if exist('PlotZeros','var')
    clear G X
    G = 1*ones(n,1);
    for i=2:n2
        G = [G ; i*ones(n,1)];
    end
    
    X=reshape(permute([ TabletData(1:end,11,:), ...
        TabletData(1:end,12,:)],[1 3 2])...
        ,[(n)*n2,2]);
end

%gscatter(X(:,1), X(:,2), G)
%gscatter(X(1:n-10,1), X(1:n-10,2), G(1:n-10))

% for i=1:n2
%     kstd_u(i)=nanstd(TabletData(1:end-10,11,i));
%     kstd_v(i)=nanstd(TabletData(1:end-10,12,i));
% end

for k=1:n2
    %# indices of points in this group
    idx = ( G == k );
    %     if kstd(k) > .013
    %         continue
    %     end
    
    %# substract mean
    Mu = nanmean( X(idx,:) );
    X0 = bsxfun(@minus, X(idx,:), Mu);
    
    STD = 1;
    conf = 2*normcdf(STD)-1;
    scale = chi2inv(conf,2);
    Cov = cov(X0,'omitrows') * scale;
    [V, D] = eig(Cov);
    
    %# eigen decomposition [sorted by eigen values]
    [D, order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);
    
    %# Make elipse
    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %# project circle back to orig space
    
    e2=zeros(100,1);
    for i=1:100
        e2(i)=sqrt(e(1,i)+e(1,i));
    end
    
    %# plot
    if strcmp(location,'200s')
        p1=plot(e(1,:), e(2,:),'Color','k');
        if strcmp(files(k).participant,'dummy')
            p1=plot(e(1,:), e(2,:),'Color','r');
        end
        
        %         %Mod to pick out short adapt vs long adap
        %     elseif strcmp(location,'PAMELA')
        %         P='DG'; %specify participant
        %         if strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
        %             p1{i}=plot(e(1,:), e(2,:),'DisplayName',files(k).Light);
        %         end
        
    elseif strcmp(location,'PAMELA')
        OM = 1; %Original or mod? O=0, M=1;
        %Mod plots means as squares, with low alpha so that overlapping
        %points are visible
        
        P='DG'; %specify participant
        %Original
        if OM
            if strcmp(files(k).Light(1:2),'CW') && strcmp(files(k).participant,P)
                p(1)=plot(e(1,:), e(2,:), 'Color','g','DisplayName','CW');
            elseif strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
                p(2)=plot(e(1,:), e(2,:), 'Color','b','DisplayName','WW');
            elseif strcmp(files(k).Light(1:2),'MH') && strcmp(files(k).participant,P)
                p(3)=plot(e(1,:), e(2,:), 'Color','r','DisplayName','MH');
            elseif strcmp(files(k).participant,'dummy')
                p(4)=plot(e(1,:), e(2,:),':','Color','k','DisplayName','BL');
            end
        end
        %Mod
        if ~OM
            if strcmp(files(k).Light(1:2),'CW') && strcmp(files(k).participant,P)
                p(1)=scatter(Mu(1),Mu(2),'g','filled','DisplayName','CW','MarkerFaceAlpha',.5);
                if ~exist('CW_catch','var'), CW_catch=[]; end
                CW_catch = [CW_catch; Mu];
            elseif strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
                p(2)=scatter(Mu(1),Mu(2),'b','filled','DisplayName','WW','MarkerFaceAlpha',.5);
            elseif strcmp(files(k).Light(1:2),'MH') && strcmp(files(k).participant,P)
                p(3)=scatter(Mu(1),Mu(2),'r','filled','DisplayName','MH','MarkerFaceAlpha',.5);
%             elseif strcmp(files(k).participant,'dummy')
%                 p4=scatter(Mu(1),Mu(2),'k','filled','DisplayName','BL','MarkerFaceAlpha',.5);
            end
        end
        %text(Mu(1),Mu(2),string(k))
        
    elseif strcmp(location,'BM') && kstd(k) < thresh_SD && diffs(k) < thresh_DBUR
        %elseif strcmp(location,'BM') %without exclusions
        %if strcmp(files(k).participant,'Public') %Exclude LM and DG
        if k==55 %from time-stamp this appears to be from before the data collection started, and so I assume it is not 'real' data
            continue
        elseif files(k).date==10||files(k).date==11 %77/78
            
            %scatter mean
            p(1) = scatter(Mu(1),Mu(2),'rs','filled','DisplayName','77/78');
            
            %plot elipse
            %scatter(e(1,:), e(2,:))
            %plot(e(1,:), e(2,:), 'Color','r');
            
            %plot line
            %                 e2=dist(e);
            %                 [~,I] = max(e2(:));
            %                 [I_row, I_col] = ind2sub(size(e2),I);
            %                 plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'r:')
            
        elseif files(k).date==12||files(k).date==13 %25 Africa
            p(2) = scatter(Mu(1),Mu(2),'kv','filled','DisplayName','25');
            %plot(e(1,:), e(2,:), 'Color','g');
            %                 e2=dist(e);
            %                 [~,I] = max(e2(:));
            %                 [I_row, I_col] = ind2sub(size(e2),I);
            %                 plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'k--')
        elseif files(k).date==14 %Great Court
            p(3) = scatter(Mu(1),Mu(2),'bo','filled','DisplayName','GC');
            %plot(e(1,:), e(2,:), 'Color','b');
            %                 e2=dist(e);
            %                 [~,I] = max(e2(:));
            %                 [I_row, I_col] = ind2sub(size(e2),I);
            %                 plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'b-.')
        elseif strcmp(files(k).participant,'dummy') %dummy data
            %scatter(Mu(1),Mu(2),'gv','filled');
            %scatter(X(idx,1),X(idx,2),'g','filled')
        end
        %             if k==55 %testing specific values
        %                 scatter(X(idx,1),X(idx,2),'g','filled')
        %             end
        %             text(Mu(1),Mu(2),num2str(k))
        %end
%     elseif strcmp(location,'PAMELA_20180205')
%         end
    elseif strcmp(location,'PAMELA_20180205')
        plt_all = 0;
        if plt_all
            %plots all, excluding SP, BG, NPG
            if ~(strcmp(files(k).participant,'SP') ||...
                    strcmp(files(k).participant,'BG') ||...
                    strcmp(files(k).participant,'NPG'))
                if k==1
                    %p1{k}=plot(e(1,:), e(2,:),'k');
                    p(1)=scatter(Mu(1),Mu(2),'k','filled','MarkerFaceAlpha',.5,'DisplayName','BL');
                elseif (1<k) && (k<11)
                    %p1{k}=plot(e(1,:), e(2,:),'r');
                    p(2)=scatter(Mu(1),Mu(2),'r','filled','MarkerFaceAlpha',.5,'DisplayName','MH1');
                elseif (10<k) && (k<20)
                    %p1{k}=plot(e(1,:), e(2,:),'g');
                    p(3)=scatter(Mu(1),Mu(2),'g','filled','MarkerFaceAlpha',.5,'DisplayName','ML');
                elseif (19<k)
                    %p1{k}=plot(e(1,:), e(2,:),'b');
                    p(4)=scatter(Mu(1),Mu(2),'b','filled','MarkerFaceAlpha',.5,'DisplayName','MH2');
                end
            end
        else
            %specify a single participant (+dummy)
            if strcmp(files(k).participant,'dummy') ||...
                    strcmp(files(k).participant,'PK')
                
                % plots all, excluding SP, BG, NPG
                %         if ~(strcmp(files(k).participant,'SP') ||...
                %                 strcmp(files(k).participant,'BG') ||...
                %                 strcmp(files(k).participant,'NPG'))
                
                if k==29 || k==1
                    p(1)=plot(e(1,:), e(2,:),'k','DisplayName','BL');
                    scatter(X(idx,1),X(idx,2),'k','filled')
                    %comet(X(idx,1),X(idx,2))                    
                elseif (1<k) & (k<11)
                    p(2)=plot(e(1,:), e(2,:),'r','DisplayName','MH1');
                    scatter(X(idx,1),X(idx,2),'r','filled')
                    %comet(X(idx,1),X(idx,2))
                elseif (10<k) & (k<20)
                    p(3)=plot(e(1,:), e(2,:),'g','DisplayName','ML');
                    scatter(X(idx,1),X(idx,2),'g','filled')
                    %comet(X(idx,1),X(idx,2))
                elseif (19<k<29)
                    p(4)=plot(e(1,:), e(2,:),'b','DisplayName','MH1');
                    scatter(X(idx,1),X(idx,2),'b','filled')
                    %comet(X(idx,1),X(idx,2))
                end
                %title(files(k).participant)
                %saveas(fig,strcat('bg',files(k).participant,'.tif'))
            end
        end
    elseif strcmp(location,'basement_rgby_test')
        if k==1
            p1{k}=plot(e(1,:), e(2,:),'r');
            scatter(X(idx,1),X(idx,2),...
                'MarkerFaceColor','r','MarkerEdgeColor','k')
        elseif k==2
            p1{k}=plot(e(1,:), e(2,:),'g');
            scatter(X(idx,1),X(idx,2),...
                'MarkerFaceColor','g','MarkerEdgeColor','k')
        elseif k==3
            p1{k}=plot(e(1,:), e(2,:),'b');
            scatter(X(idx,1),X(idx,2),...
                'MarkerFaceColor','b','MarkerEdgeColor','k')
        elseif k==4
            p1{k}=plot(e(1,:), e(2,:),'Color',[1,.9,0]);
            scatter(X(idx,1),X(idx,2),...
                'MarkerFaceColor',[1,.9,0],'MarkerEdgeColor','k')
        elseif k==5
            p1{k}=plot(e(1,:), e(2,:),'Color',[0.5,0.5,0.5]);
            scatter(X(idx,1),X(idx,2),...
                'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','k')
        end
    elseif strcmp(location,'basement_greyCard_test')
        
        if k==4 %LM hole
            p1{k}=plot(e(1,:), e(2,:),'r');
            scatter(X(idx,1),X(idx,2),...
                'r','DisplayName','LM - hole method')
        elseif k==5 %LM default
            p1{k}=plot(e(1,:), e(2,:),'r');
            scatter(X(idx,1),X(idx,2),'s',...
                'MarkerFaceColor','r','DisplayName','LM - default method')
        elseif k==6 %DG hole
            p1{k}=plot(e(1,:), e(2,:),'b');
            scatter(X(idx,1),X(idx,2),...
                'b','DisplayName','DG - hole method')
        elseif k==7 %DG default
            p1{k}=plot(e(1,:), e(2,:),'b');
            scatter(X(idx,1),X(idx,2),'s',...
                'MarkerFaceColor','b','DisplayName','DG - default method')
            %         elseif k<4 % baseline data and data from pre-runs
            %             p1{k}=plot(e(1,:), e(2,:),'k');
            %             scatter(X(idx,1),X(idx,2),...
            %                 'k','DisplayName','Test data')
        end
        
    end
    
    xlabel('u'''),ylabel('v''')
end

% Here goes all the stuff that you don't need to be applied for each run of
% the data, like graph formatting
if strcmp(location,'PAMELA_20180205')
    %legend([h1 h2 h3 h4],{'Baseline','MH1','ML','MH2'})
    legend(p)
    %     xlim([0.15 0.25])
    %     ylim([0.43 0.51])
    xlim([0.17 0.22])
    ylim([0.45 0.50])
    cleanTicks
    if saveFigs
        if plt_all
            save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_20180205_results.pdf')
        else
            save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_20180205_Individual2.pdf')
        end
    end
end

if strcmp(location,'BM')
    legend(p,'Location','northwest')
    colorbar off
    xlim([0.16, 0.28])
    ylim([0.46, 0.56])
    
    scatter(GLuv(1,3),GLuv(2,3),'r*','DisplayName','77/78');
    scatter(GLuv(1,2),GLuv(2,2),'k*','DisplayName','25');
    scatter(GLuv(1,1),GLuv(2,1),'b*','DisplayName','GC');
    
    cleanTicks
    
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\exp2.pdf')
    end
    
end
    

if strcmp(location,'basement_rgby_test')
    axis equal
    xlim([0.14 0.24])
    ylim([0.41 0.52])
    cleanTicks
    if saveFigs
        save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\basement_rgby_test.pdf')
    end
end

if strcmp(location,'PAMELA')
    legend(p,'Location','southeast')
%     xlim([0.175 0.23])
%     ylim([0.455 0.5])
   % xticks(0.15:0.01:0.25)
   % yticks(0.43:0.01:0.51)
    
    %Plot light chromaticities
    if exist('spd_uv_CW','var')
        scatter(spd_uv_CW(1),spd_uv_CW(2),'*g','DisplayName','CW')        
        scatter(spd_uv_WW(1),spd_uv_WW(2),'*b','DisplayName','WW')        
        scatter(spd_uv_MH(1),spd_uv_MH(2),'*r','DisplayName','MH')
        xlim([0.175 0.255])
        ylim([0.455 0.53])
    end
    
    stats=0;
    if stats
        scatter(Mu_u_CW,Mu_v_CW,'g^','DisplayName','CW')
        scatter(Mu_u_WW,Mu_v_WW,'b^','DisplayName','WW')
        scatter(Mu_u_MH,Mu_v_MH,'r^','DisplayName','MH')
        
    end
   % offset
   disp(Mu_u_WW-Mu_u_CW)
   disp(Mu_v_WW-Mu_v_CW)
   
   if ~OM
       xlim([0.18 0.26])
       ylim([0.46 0.535])
   else
       xlim([0.17 0.22])
       ylim([0.45 0.50])
   end
   cleanTicks
   
   if saveFigs
       if ~OM
           save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_DG_results.pdf')
       else
           save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\tablet\PAMELA_DG_results_ellipses.pdf')
       end
   end
       
end

%         if strcmp(files(k).Light(1:2),'CW') && strcmp(files(k).participant,P)
%             p1=plot(e(1,:), e(2,:), 'Color','g','DisplayName','CW');
%         elseif strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
%             p2=plot(e(1,:), e(2,:), 'Color','b','DisplayName','WW');
%         elseif strcmp(files(k).Light(1:2),'MH') && strcmp(files(k).participant,P)
%             p3=plot(e(1,:), e(2,:), 'Color','r','DisplayName','MH');
%             

%close

% %% Assessing intra-observer variation
% % Nabbed script from 'AptialCAtestAnyalysis003_WIP.m'
% 
% clc, close all
% 
% subsample_median=   zeros(n-10,2,n2-1); %n2-1 would exclude dummy
% subsample_SD=       zeros(n-10,2,n2-1); %n2-1 would exclude dummy
% 
% for j=1:n2-1                                                         %For all the datasets (including dummy)
%     for i=1:n-10                                      %Using 'i' as a value from 1 to the 190 (the max number of trials)
%         clear subsample
%         subsample=TabletData(1:i,11,j);
%         subsample_median(i,1,j)=nanmedian(subsample);
%         subsample_SD(i,1,j)=nanstd(subsample);                           %calculate the standard deviation of that sample
%     end
% end
% 
% for j=1:n2-1
%     for i=1:n-10
%         clear subsample
%         subsample=TabletData(1:i,12,j);
%         subsample_median(i,2,j)=nanmedian(subsample);
%         subsample_SD(i,2,j)=nanstd(subsample);
%     end
% end
% 
% figure, hold on
% plot(10:n-10,mean(subsample_SD(10:end,:,1),2),':r','LineWidth',2);
% plot(10:n-10,mean(subsample_SD(10:end,:,2),2),'--b','LineWidth',2);
% 
% %axis([10 n-1 3.2*10^-3 5.5*10^-3])
% xlabel('# of data points in subsample')
% ylabel({'Standard Deviation (SD)','Median Across All Datasets'})
% ax = gca; ax.XTick = [10:20:190];
% legend({'u''','v'''})
% 
% % %% Single run demo
% %
% % figure,
% % plot(1:190,u_prime(:,10))
% % axis([1 190 -inf inf])
% % xlabel('Trial #','FontSize',20)
% % ylabel('u''','FontSize',20)
% %
