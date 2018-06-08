%Untested
%Writing to sub out from SAPS_DataAnalysis.m but not got time
%Ideally should be similar to read_UPRtek

% SEAHA has a GL Optis device, which I used to take measurements at the
% British Museum and elsewhere

function [spd] = read_GL(folder,plt,sv_plt)

cd(folder)
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

for i=0:255
    out(i+1,1)=str2double(xFile.getElementsByTagName('row').item(i).getAttribute('wavelength'));
end

spd=out(28:end,:);

if plt==1
    %nabbed from read_UPRtek
    %     if i>1 %plot previous in black, to make new ones show up in red
    %         plot();
    %     end
    %     plot();
    %     if i==length(xlsx) %overplot the last one in black
    %         plot(360:760,data(:,i)*peak(i),'k');
    %     end
    
end

%nabbed from read_UPRtek
% if sv_plt
%     drawnow
%     saveas(gcf,strcat(filename(1:end-5),'.tif'))
% end

end