clc, clear, close all
% Requires Psychtoolbox and read_UPRtek

%Add some text here that says what this does

%% Load lighting measurements

cd('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\20180205 Spectra')

try
    load('PAMELA_ChromData')
catch
    [spectral_data,peak,lux_fromExcel,spd_uv,spd_xy]=read_UPRtek('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\20180205 Spectra',[],[],[]);
    S_spectral_data=[360,1,401];
    save('PAMELA_ChromData')
end

for i=1:size(spectral_data,2)
    spectral_data_UnNormalised(:,i) = spectral_data(:,i).*peak(i);
end

% figure; plot(SToWls(S_spectral_data),spectral_data)
% figure; plot(SToWls(S_spectral_data),spectral_data_UnNormalised)

%comment: mel_low = 0; 
%comment: mel_high = 1
index=[ones(1,6),zeros(1,7),ones(1,6)];

% %Plot each SPD
% figure; plot(SToWls(S_spectral_data),spectral_data(:,index==0))
% figure; plot(SToWls(S_spectral_data),spectral_data(:,index==1))

ML.mean = mean(spectral_data_UnNormalised(:,index==0),2);
MH.mean = mean(spectral_data_UnNormalised(:,index==1),2);

ML.int81=SplineSpd(S_spectral_data, ML.mean,[380:5:780]', 1);
MH.int81=SplineSpd(S_spectral_data, MH.mean,[380:5:780]', 1);
ML.int441=SplineSpd(S_spectral_data, ML.mean,[390:830]', 1);
MH.int441=SplineSpd(S_spectral_data, MH.mean,[390:830]', 1);
ML.int271=SplineSpd(S_spectral_data, ML.mean,[390:660]', 1);
MH.int271=SplineSpd(S_spectral_data, MH.mean,[390:660]', 1);

% % Test of switch from old interp to new interp
% ML.int81a=SplineSpd(S_spectral_data, ML.mean, [380,5,81], 1);
% figure, hold on
% plot(380:5:780,ML.int81a,'r') %SplineSpd
% ML.int81b=interp1(SToWls(S_spectral_data),ML.mean,380:5:780,'spline',0);
% plot(380:5:780,ML.int81b,'b') %interp1
% plot(SToWls(S_spectral_data),ML.mean,'g') %original data

% ML.int81=interp1(SToWls(S_spectral_data),ML.mean,380:5:780,'spline',0);
% MH.int81=interp1(SToWls(S_spectral_data),MH.mean,380:5:780,'spline',0);
% ML.int441=interp1(SToWls(S_spectral_data),ML.mean,390:830,'spline',0);
% MH.int441=interp1(SToWls(S_spectral_data),MH.mean,390:830,'spline',0);
% ML.int271=interp1(SToWls(S_spectral_data),ML.mean,390:660,'spline',0);
% MH.int271=interp1(SToWls(S_spectral_data),MH.mean,390:660,'spline',0);

% figure; hold on
% plot(SToWls(S_spectral_data),ML.mean);
% plot(380:5:780,ML.int81);
% plot(390:830,ML.int441,'--');

%% Plot SPDs

norm = 1; %normalise SPDs to max=1

figure; hold on
if norm
    plot(SToWls(S_spectral_data),ML.mean/max(ML.mean),'r','DisplayName','ML');
    plot(SToWls(S_spectral_data),MH.mean/max(MH.mean),'b','DisplayName','MH');
    xlabel('Wavelength (nm)')
    ylabel('Normalised Power')
else
    plot(SToWls(S_spectral_data),ML.mean,'r','DisplayName','ML');
    plot(SToWls(S_spectral_data),MH.mean,'b','DisplayName','MH');
    xlabel('Wavelength (nm)')
    %ylabel('?')!!!!!!!!
end
legend

%% Calculate chromaticities

load T_xyz1931
ML.XYZ_1931 = T_xyz1931*ML.int81;
MH.XYZ_1931 = T_xyz1931*MH.int81;
ML.xy_1931 = [ML.XYZ_1931(1)/sum(ML.XYZ_1931),ML.XYZ_1931(2)/sum(ML.XYZ_1931)];
MH.xy_1931 = [MH.XYZ_1931(1)/sum(MH.XYZ_1931),MH.XYZ_1931(2)/sum(MH.XYZ_1931)];
ML.lum = ML.XYZ_1931(2)*683/1000; %not sure where the factor of 1000 comes from
MH.lum = MH.XYZ_1931(2)*683/1000; %not sure where the factor of 1000 comes from

ML.lum_fromExcel = mean(lux_fromExcel(index==0));
MH.lum_fromExcel = mean(lux_fromExcel(index==1));

% %How representative are the averages of the actual measurements over time?
%figure, plot(lux_fromExcel)

load T_xyz1964
ML.XYZ_1964 = T_xyz1964*ML.int81;
MH.XYZ_1964 = T_xyz1964*MH.int81;
ML.xy_1964 = [ML.XYZ_1964(1)/sum(ML.XYZ_1964),ML.XYZ_1964(2)/sum(ML.XYZ_1964)];
MH.xy_1964 = [MH.XYZ_1964(1)/sum(MH.XYZ_1964),MH.XYZ_1964(2)/sum(MH.XYZ_1964)];

load T_cones_ss2
ML.LMS_ss2 = T_cones_ss2*ML.int441;
MH.LMS_ss2 = T_cones_ss2*MH.int441;

load T_cones_ss10
ML.LMS_ss10 = T_cones_ss10*ML.int441;
MH.LMS_ss10 = T_cones_ss10*MH.int441;

load T_melanopsin
ML.mel = T_melanopsin*ML.int271;
MH.mel = T_melanopsin*MH.int271;

% %-% Trying to work out the scaling for mel vs photopic lux
% 
% % Are they normalised differently?
% figure, hold on
% plot(SToWls(S_xyz1931),T_xyz1931(2,:))
% plot(SToWls(S_melanopsin),T_melanopsin)
% 
% figure, hold on
% plot(SToWls(S_xyz1931),MH.int81)
% plot(SToWls(S_melanopsin),MH.int271)
% 
% 
% %-%


%% Plot CIE 1931 xy chromaticity diagram

% I think the data that ships with Psychtoolbox disagrees with the CIE data
% and creates some funky bits at the end of the spectrum because of
% rounding error bits
plotx_31=T_xyz1931(1,1:73)./sum(T_xyz1931(:,1:73));    %plotx(end+1)=plotx(1);
ploty_31=T_xyz1931(2,1:73)./sum(T_xyz1931(:,1:73));    %ploty(end+1)=ploty(1);

figure, hold on
plot(plotx_31,ploty_31,'k');
scatter(plotx_31,ploty_31,'k');
axis equal

scatter(ML.xy_1931(1),ML.xy_1931(2),'k','filled')
text(ML.xy_1931(1),ML.xy_1931(2),'ML\_1931')

scatter(MH.xy_1931(1),MH.xy_1931(2),'b','filled')
text(MH.xy_1931(1),MH.xy_1931(2),'MH\_1931')

%% 1964

plotx_64=T_xyz1964(1,1:73)./sum(T_xyz1964(:,1:73));    %plotx(end+1)=plotx(1);
ploty_64=T_xyz1964(2,1:73)./sum(T_xyz1964(:,1:73));    %ploty(end+1)=ploty(1);

figure, hold on
plot(plotx_64,ploty_64,'k');
scatter(plotx_64,ploty_64,'k');
axis equal

scatter(ML.xy_1964(1),ML.xy_1964(2),'k','filled')
text(ML.xy_1964(1),ML.xy_1964(2),'ML\_1964')

scatter(MH.xy_1964(1),MH.xy_1964(2),'b','filled')
text(MH.xy_1964(1),MH.xy_1964(2),'MH\_1964')

%% Version which calculates values for each individual SPD
% Are the values recorded for lux in excel well matched by those calculated
% here?
% A: Yes, only out by <1.3

% Load lighting measurements
[spectral_data,peak,lux_fromExcel]=read_UPRtek('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\20180205 Spectra',[],[]);
S_spectral_data=[360,1,401];

for i=1:size(spectral_data,2)
    spectral_data_UnNormalised(:,i) = spectral_data(:,i).*peak(i);
    
    M(i).int81=SplineSpd(S_spectral_data, spectral_data_UnNormalised(:,i),[380:5:780]', 1);
    M(i).int441=SplineSpd(S_spectral_data, spectral_data_UnNormalised(:,i),[390:830]', 1);
    M(i).int271=SplineSpd(S_spectral_data, spectral_data_UnNormalised(:,i),[390:660]', 1);
end


% Calculate chromaticities

load T_xyz1931
load T_xyz1964
load T_cones_ss2
load T_cones_ss10
load T_melanopsin

for i=1:size(spectral_data,2)
    M(i).XYZ_1931 = T_xyz1931*M(i).int81;
    M(i).xy_1931 = [M(i).XYZ_1931(1)/sum(M(i).XYZ_1931),M(i).XYZ_1931(2)/sum(M(i).XYZ_1931)];
    M(i).lum = M(i).XYZ_1931(2)*683;
    
    M(i).XYZ_1964 = T_xyz1964*M(i).int81;
    M(i).xy_1964 = [M(i).XYZ_1964(1)/sum(M(i).XYZ_1964),M(i).XYZ_1964(2)/sum(M(i).XYZ_1964)];
    
    M(i).LMS_ss2 = T_cones_ss2*M(i).int441;    
    M(i).LMS_ss10 = T_cones_ss10*M(i).int441;
    
    M(i).mel = T_melanopsin*M(i).int271;
    
    M(i).lux=lux_fromExcel(i);    
%    M(i).lux-(M(i).lum/1000)
end

%% Compare spreadsheet chromaticity values with computed values
uv_mean_ML=mean(spd_uv(index==0,:));
uv_mean_MH=mean(spd_uv(index==1,:));

xy_mean_ML=mean(spd_xy(index==0,:));
xy_mean_MH=mean(spd_xy(index==1,:));

figure, hold on
plot(plotx_31, ploty_31)

scatter(ML.xy_1931(1),ML.xy_1931(2),'k','filled')
text(ML.xy_1931(1),ML.xy_1931(2),'ML\_1931')

scatter(xy_mean_ML(1),xy_mean_ML(2),'r','filled')
text(xy_mean_ML(1),xy_mean_ML(2),'ML\_1931')

scatter(MH.xy_1931(1),MH.xy_1931(2),'b','filled')
text(MH.xy_1931(1),MH.xy_1931(2),'MH\_1931')

scatter(xy_mean_MH(1),xy_mean_MH(2),'r','filled')
text(xy_mean_MH(1),xy_mean_MH(2),'ML\_1931')

axis equal

%% Assessing how much light the tablet emits compared to the ambient lighting

clc, clear, close all

% I measured 3 repeats of 3 different stimuli. 
% Then 3 with the display turned off. 
% Then 3 with the display on, but lights off.
% Then 3 with display and lights off.

try
    load('PAMELA_Tablet_Emission')
catch
    [spectral_data,peak,lux_fromExcel,spd_uv,spd_xy]=read_UPRtek('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\20180207 Spectra',1,1,1);
    S_spectral_data=[360,1,401];
    save('PAMELA_Tablet_Emission')
end

spectral_data_UnNormalised = spectral_data.*peak;
%figure, plot(spectral_data)
figure, hold on

%Two sets of plot calls to make the legend work correctly.
%There's probably a neater way to do this.
p1 = plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,1),'r');
p2 = plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,10),'g'); %Looks like I accidentally only took 2 of these rather than 3
p3 = plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,12),'b');
p4 = plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,15),'k');

plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,1:9),'r');
plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,10:11),'g'); %Looks like I accidentally only took 2 of these rather than 3
plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,12:14),'b');
plot(SToWls(S_spectral_data),spectral_data_UnNormalised(:,15:17),'k');

legend([p1 p2 p3 p4], {'Light, with stimuli','Light, no stimuli','No light, with stimuli','No light, no stimuli'})

xlabel('Wavelength (nm)')
ylabel('Spectral irradiance (W·m?2·nm?1)')

