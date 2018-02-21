clc, clear, close all
% Requires Psychtoolbox and read_UPRtek

%% Load lighting measurements
[spectral_data,peak]=read_UPRtek('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\20180205 Spectra',[]);
S_spectral_data=[360,1,401];

for i=1:size(spectral_data,2)
    spectral_data_UnNormalised(:,i) = spectral_data(:,i).*peak(i);
end

% figure; plot(SToWls(S_spectral_data),spectral_data)
% figure; plot(SToWls(S_spectral_data),spectral_data_UnNormalised)

%mel_low = 0; 
%mel_high = 1
index=[ones(1,6),zeros(1,7),ones(1,6)];

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

% figure; plot(SToWls(S_spectral_data),ML.mean);title('ML');
% figure; plot(SToWls(S_spectral_data),MH.mean);title('MH');

%% Calculate chromaticities

load T_xyz1931
ML.XYZ_1931 = T_xyz1931*ML.int81;
MH.XYZ_1931 = T_xyz1931*MH.int81;
ML.xy_1931 = [ML.XYZ_1931(1)/sum(ML.XYZ_1931),ML.XYZ_1931(2)/sum(ML.XYZ_1931)];
MH.xy_1931 = [MH.XYZ_1931(1)/sum(MH.XYZ_1931),MH.XYZ_1931(2)/sum(MH.XYZ_1931)];
ML.lum = ML.XYZ_1931(2)*683;
MH.lum = MH.XYZ_1931(2)*683;

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

%% Plot CIE 1931 xy chromaticity diagram

clear plotx ploty

% I think the data that ships with Psychtoolbox disagrees with the CIE data
% and creates some funky bits at the end of the spectrum because of
% rounding error bits
plotx=T_xyz1931(1,1:73)./sum(T_xyz1931(:,1:73));    %plotx(end+1)=plotx(1);
ploty=T_xyz1931(2,1:73)./sum(T_xyz1931(:,1:73));    %ploty(end+1)=ploty(1);

figure, hold on
plot(plotx,ploty,'k');
scatter(plotx,ploty,'k');
axis equal

scatter(ML.xy_1931(1),ML.xy_1931(2),'k','filled')
text(ML.xy_1931(1),ML.xy_1931(2),'ML\_1931')

scatter(MH.xy_1931(1),MH.xy_1931(2),'b','filled')
text(MH.xy_1931(1),MH.xy_1931(2),'MH\_1931')

%% 1964

clear plotx ploty
plotx=T_xyz1964(1,1:73)./sum(T_xyz1964(:,1:73));    %plotx(end+1)=plotx(1);
ploty=T_xyz1964(2,1:73)./sum(T_xyz1964(:,1:73));    %ploty(end+1)=ploty(1);

figure, hold on
plot(plotx,ploty,'k');
scatter(plotx,ploty,'k');
axis equal

scatter(ML.xy_1964(1),ML.xy_1964(2),'k','filled')
text(ML.xy_1964(1),ML.xy_1964(2),'ML\_1964')

scatter(MH.xy_1964(1),MH.xy_1964(2),'b','filled')
text(MH.xy_1964(1),MH.xy_1964(2),'MH\_1964')
