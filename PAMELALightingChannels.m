% Code for plotting the different lights at PAMELA

clc, clear, close all

%Specify folder, load data
rootdir=('C:\Users\cege-user\Dropbox\UCL\Data\Tablet\PAMELA\2017 Spectra');
[data,peak,lux,spd_uv]=read_UPRtek(rootdir,0,0,0);

%% Plot data

h= plot(360:760,data(:,[1,3:8]));
set(h, {'color'}, {'k';'r';'g';'b';[0.9,0.9,0];'k';'k'});
set(h, {'LineStyle'}, {'-';'-';'-';'-';'-';':';'--'});
legend({'High output white','Red','Green','Blue','Amber','Warm White','Cool White'})

xlabel('Wavelength (nm)')
ylabel('Normalised SPD')
title('Normalised SPDs of PAMELA lighting channels')