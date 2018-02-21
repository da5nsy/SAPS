% Function for reading data from a UPRtek MK350
% https://www.uprtek.com/en/index.html

% Note, before use, all files were opened in excel and resaved as '.xlsx'
% files due to unusual '.XLS' files (probably a quicker way to do this
% could be found).

function [data,peak] = read_UPRtek(folder,plt)

cd(folder)
xlsx=dir('*.xlsx');

if plt==1
    figure, hold on
end

for i=1:length(xlsx)
    filename = xlsx(i).name;
    xlRange = 'B10:B410';
    
    data(:,i) = xlsread(filename,xlRange);
    [~,txt,~] = xlsread(filename,'A7:A7');
    words = regexp(txt,' ','split');
    peak(i) = str2double(words{1,1}{1,4});
    clear txt words
    
    if plt==1
        plot(360:760,data(:,i)*peak,'k');
    end
    %plot(360:760,av_data*peak,'b');
    %title(xlsx(i).name)
    
    drawnow
end
end

%%
%av_data=mean(data(:,4:end),2);