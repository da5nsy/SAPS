% Function for reading data from a UPRtek MK350
% https://www.uprtek.com/en/index.html

% Note, before use, all files were opened in excel and resaved as '.xlsx'
% files due to unusual '.XLS' files (probably a quicker way to do this
% could be found).

% I use this for reading the spectral data from PAMELA, since they have one
% of these devices.

function [data,peak,lux,spd_uv] = read_UPRtek(folder,plt,sv_plt,norm)

cd(folder)
xlsx=dir('*.xlsx');

if plt==1
    figure, hold on
end

for i=1:length(xlsx)
    %figure,
    filename = xlsx(i).name;
    xlRange = 'B10:B410';
    
    data(:,i) = xlsread(filename,xlRange);
    
    [~,txt,~] = xlsread(filename,'A3:A3');
    words = regexp(txt,' ','split');
    spd_uv(i,1) = str2double(words{1,1}{1,3});
    clear txt words
    
    [~,txt,~] = xlsread(filename,'A4:A4');
    words = regexp(txt,' ','split');
    spd_uv(i,2) = str2double(words{1,1}{1,3});
    clear txt words
    
    [~,txt,~] = xlsread(filename,'A7:A7');
    words = regexp(txt,' ','split');
    peak(i) = str2double(words{1,1}{1,4});
    clear txt words
    
    [~,txt,~] = xlsread(filename,'A7:A7');
    words = regexp(txt,' ','split');
    peak(i) = str2double(words{1,1}{1,4});
    clear txt words
    
    [~,txt,~] = xlsread(filename,'A9:A9');
    words = regexp(txt,' ','split');
    lux(i) = str2double(words{1,1}{1,3});
    clear txt words
    
    if plt==1
        if norm            
            if i>1 %plot previous in black, to make new ones show up in red
                plot(360:760,data(:,i-1),'k');
            end
            plot(360:760,data(:,i),'r');
            if i==length(xlsx) %overplot the last one in black
                plot(360:760,data(:,i),'k');
            end
        else
            if i>1 %plot previous in black, to make new ones show up in red
                plot(360:760,data(:,i-1)*peak(i-1),'k');
            end
            plot(360:760,data(:,i)*peak(i),'r');
            if i==length(xlsx) %overplot the last one in black
                plot(360:760,data(:,i)*peak(i),'k');
            end            
        end        
    end
    %plot(360:760,av_data*peak,'b');
    title(xlsx(i).name,'Interpreter', 'none')
    
    if sv_plt
        drawnow
        saveas(gcf,strcat(filename(1:end-5),'.tif'))
    end
end



%%
%av_data=mean(data(:,4:end),2);


end
