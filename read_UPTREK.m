function [data,peak] = read_UPTREK(folder,plt)

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