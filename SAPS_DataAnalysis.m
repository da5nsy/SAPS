clc, clear, close all

%General purpose script assimilating previous advances in different files
%Many pieces cut, see prev versions for missing sections

% To do list

% Add Grant data (different run length, should just be able to change n)
% Clean up data, find a way to threshold and recognise duds
% Threshold BM data

%% Load Data

%-% Tablet Data
% (Gone back to having 2 xlsread calls, despite being slower, because having
% both raw and numerical data is useful, but they require different excel
% addresses to specify data)

location='PAMELA_20180205';  %changing this changes EVERYTHING
                    %'PAMELA' or 'GRANT' or 'BM' or '200s' or
                    %'PAMELA_20180205' or 'basement_rgby_test'

rootdir = fullfile('C:','Users','cege-user','Dropbox','UCL','Data','Tablet',location);
cd(rootdir)

try
    load(sprintf('TabletData_%s',location))
catch
    disp('Reading from excel. Wait ~42.3 seconds (PAMELA), and 171.3 seconds for BM')
    tic
    TBfilename = fullfile(rootdir,location);
    [~,sheets] = xlsfinfo(TBfilename);
    if strcmp(location,'200s')
        TabletData = zeros(200,8,numel(sheets));
    else
        TabletData = zeros(40,8,numel(sheets));
    end
    
    for i=1:numel(sheets) %65 for BM, 20 for PAMELA
        i %to provide a progress update in command window, should take ~2mins
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
    
    %-% Supplementary variables
    n=size(TabletData,1); %how many trials
    
    clear i
    save(sprintf('TabletData_%s',location))
end

%% Processing

whiteXYZ=[360.476; 394.674; 416.714];

L_star=65;

% Get white point u and v
uv0 = XYZTouv(whiteXYZ);

for i=1:numel(sheets)
    
    %Calculate systematic offset (spatial calibration)
    avXdp(i) = median(TabletData(end-10:end,7,i)-TabletData(end-10:end,1,i));
    avYdp(i) = median(TabletData(end-10:end,8,i)-TabletData(end-10:end,2,i));
    
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
    TabletData(:,9,i)   =   (TabletData(:,9,i)./21.88); %u*
    TabletData(:,10,i)  =   (TabletData(:,10,i)./21.88); %v*
    
    % u*v* --> u'v'
    TabletData(:,11,i) = (TabletData(:,9,i) ./ (13.0 * L_star) ) + uv0(1); %u'
    TabletData(:,12,i) = (TabletData(:,10,i)./ (13.0 * L_star) ) + uv0(2); %v'
    
end

%% Plot u'v'

ciefile = fullfile('C:','Users','cege-user','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

cie2= xlsread(ciefile,'1931 col observer','A6:D86');

lambdaCie2=cie2(:,1);
xbar2=cie2(:,2);
ybar2=cie2(:,3);
zbar2=cie2(:,4);

ciedata2_10001=interp1(lambdaCie2,cie2(:,2:4),380:.04:780);
xbar2_10001=interp1(lambdaCie2,xbar2,380:.04:780);
ybar2_10001=interp1(lambdaCie2,ybar2,380:.04:780);
zbar2_10001=interp1(lambdaCie2,zbar2,380:.04:780);
ubar=4*xbar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
vbar=9*ybar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);

figure, hold on
locus=scatter(ubar,vbar,100,'filled');
LCol=xyz2rgb(ciedata2_10001);
LCol(LCol<0)=0;LCol(LCol>1)=1;
locus.CData=LCol;
axis square;

%% Pull light measurements

if strcmp(location,'BM')
    
    clear out
    
    rootdir = ('C:\Users\ucesars\Dropbox\UCL\Data\Tablet\BM\GL');
    %rootdir = uigetdir; %select folder where .mmg (GL Optis) files stored
    cd(rootdir)
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
    
    % for i=2:mmgl+1
    %     figure
    %     plot(out(:,1),out(:,i))
    %     title(i-1)
    %     drawnow
    %     pause(0.5)
    % end
    
    spd=out(28:end,:);
    
    xbar2_GL=interp1(lambdaCie2,xbar2,spd(:,1),'spline');
    ybar2_GL=interp1(lambdaCie2,ybar2,spd(:,1),'spline');
    zbar2_GL=interp1(lambdaCie2,zbar2,spd(:,1),'spline');
    
    GLXYZ=[xbar2_GL,ybar2_GL,zbar2_GL]'*spd(:,2:4);
    GLxy=[GLXYZ(1,:)./sum(GLXYZ);GLXYZ(2,:)./sum(GLXYZ)];
    GLuv=[4*GLxy(1,:)./(-2*GLxy(1,:)+12*GLxy(2,:)+3);9*GLxy(2,:)./(-2*GLxy(1,:)+12*GLxy(2,:)+3)];
    scatter(GLuv(1,:),GLuv(2,:),'k*');
end

%% Plot
%figure, hold on
%axis square

fig=figure('Position', [50, 50, 800, 800]); hold on
locus=scatter(ubar,vbar,100,'filled');
LCol=xyz2rgb(ciedata2_10001);
LCol(LCol<0)=0;LCol(LCol>1)=1;
locus.CData=LCol;
axis square;

load('SAPS_SelectableColoursGamut.mat')
s=scatter(u(1:50:end),v(1:50:end),20,occ(1:50:end));
view(2)
axis('equal')
xlim([0.14 0.25]),ylim([0.41 0.52]) %close to selectable gamut boundary
colorbar
xlabel('u'''),ylabel('v'''),zlabel('Y')

% from http://stackoverflow.com/questions/3417028/ellipse-around-the-data-in-matlab

G = 1*ones(n-10,1);
for i=2:numel(sheets)
    G = [G ; i*ones(n-10,1)];
end

X=reshape(permute([ TabletData(1:end-10,11,:), ...
    TabletData(1:end-10,12,:)],[1 3 2])...
    ,[(n-10)*numel(sheets),2]);

%PlotZeros = 1;
if exist('PlotZeros','var')
    clear G X
    G = 1*ones(n,1);
    for i=2:numel(sheets)
        G = [G ; i*ones(n,1)];
    end
    
    X=reshape(permute([ TabletData(1:end,11,:), ...
        TabletData(1:end,12,:)],[1 3 2])...
        ,[(n)*numel(sheets),2]);
end

%gscatter(X(:,1), X(:,2), G)
%gscatter(X(1:n-10,1), X(1:n-10,2), G(1:n-10))

% for i=1:numel(sheets)
%     kstd_u(i)=nanstd(TabletData(1:end-10,11,i));
%     kstd_v(i)=nanstd(TabletData(1:end-10,12,i));
% end

for k=1:numel(sheets)
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
    
        %Mod to pick out short adapt vs long adap
    elseif strcmp(location,'PAMELA')
        P='DG'; %specify participant
        if strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
            p1{i}=plot(e(1,:), e(2,:),'DisplayName',files(k).Light);
        end
        
%     elseif strcmp(location,'PAMELA')
%         P='DG'; %specify participant
%         if strcmp(files(k).Light(1:2),'CW') && strcmp(files(k).participant,P)
%             p1=plot(e(1,:), e(2,:), 'Color','b','DisplayName','CW');
%         elseif strcmp(files(k).Light(1:2),'WW') && strcmp(files(k).participant,P)
%             p2=plot(e(1,:), e(2,:), 'Color','r','DisplayName','WW');
%         elseif strcmp(files(k).Light(1:2),'MH') && strcmp(files(k).participant,P)
%             p3=plot(e(1,:), e(2,:), 'Color','g','DisplayName','MH');
%         end
        
    elseif strcmp(location,'BM')
        %if strcmp(files(k).participant,'Public') %Exclude LM and DG
            if files(k).Time.d==10||files(k).Time.d==11 %
                
                %scatter mean
                scatter(Mu(1),Mu(2),'rs','filled');
                
                %plot elipse
                %scatter(e(1,:), e(2,:))
                %plot(e(1,:), e(2,:), 'Color','r');
                
                %plot line
                e2=dist(e);
                [~,I] = max(e2(:));
                [I_row, I_col] = ind2sub(size(e2),I);
                plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'r:')
                
            elseif files(k).Time.d==12||files(k).Time.d==13 %Africa
                scatter(Mu(1),Mu(2),'kv','filled');
                %plot(e(1,:), e(2,:), 'Color','g');
                e2=dist(e);
                [~,I] = max(e2(:));
                [I_row, I_col] = ind2sub(size(e2),I);
                plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'k--')
            elseif files(k).Time.d==14 %Great Court
                scatter(Mu(1),Mu(2),'bo','filled');
                %plot(e(1,:), e(2,:), 'Color','b');
                e2=dist(e);
                [~,I] = max(e2(:));
                [I_row, I_col] = ind2sub(size(e2),I);
                plot([e(1,I_row),e(1,I_col)],[e(2,I_row),e(2,I_col)],'b-.')
            end
            %end
    elseif strcmp(location,'PAMELA_20180205')
        if strcmp(files(k).participant,'test, corners, colour & bw') ||...
                strcmp(files(k).participant,'KW')
            %specify a single participant (also plots test case)
            if k==1
                p1{k}=plot(e(1,:), e(2,:),'r'); 
                scatter(X(idx,1),X(idx,2),'r','filled')
                %comet(X(idx,1),X(idx,2))
                
            elseif (1<k) && (k<11)
                p1{k}=plot(e(1,:), e(2,:),'g'); 
                scatter(X(idx,1),X(idx,2),'g','filled')
                %comet(X(idx,1),X(idx,2))
            elseif (10<k) && (k<20)
                p1{k}=plot(e(1,:), e(2,:),'b'); 
                scatter(X(idx,1),X(idx,2),'b','filled')
                %comet(X(idx,1),X(idx,2))
            elseif (19<k)
                p1{k}=plot(e(1,:), e(2,:),'k');
                scatter(X(idx,1),X(idx,2),'k','filled') 
                %comet(X(idx,1),X(idx,2))
            end
            %title(files(k).participant) 
            %saveas(fig,strcat('bg',files(k).participant,'.tif'))
        end
%     elseif strcmp(location,'PAMELA_20180205')
%         if k==1
%            p1{k}=plot(e(1,:), e(2,:),'r');
%         elseif (1<k) && (k<11)
%             p1{k}=plot(e(1,:), e(2,:),'g');
%         elseif (10<k) && (k<20)
%             p1{k}=plot(e(1,:), e(2,:),'b');
%         elseif (19<k) 
%             p1{k}=plot(e(1,:), e(2,:),'k');
%         end
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
        
    end
    
    xlabel('u'''),ylabel('v''')
end

% if strcmp(location,'PAMELA')
%     legend([p1 p2 p3])
% end
%legend 

%close

%% Plot standard deviation for each observer across runs

for i=1:numel(sheets)
    kstd_u(i)=nanstd(TabletData(1:end-10,11,i));
    kstd_v(i)=nanstd(TabletData(1:end-10,12,i));
end

kstd_mean=mean([kstd_u;kstd_v]);

figure, hold on
scatter(0,kstd_mean(1),'r','filled');
plot([0,9],[kstd_mean(1),kstd_mean(1)],'r--');
scatter(1:9,kstd_mean(2:10),'g','filled');
scatter(1:9,kstd_mean(11:19),'b','filled');
scatter(1:9,kstd_mean(20:28),'k','filled');

xlabel('Observer')
%Replace numbers with participant identifiers
%xticklabels({files([2:10]).participant}) 
ylabel('Mean SD in u'' and v''')

%% Analyse repeat values


fig=figure('Position', [50, 50, 800, 800]); hold on; 
locus=scatter(ubar,vbar,100,'filled');
LCol=xyz2rgb(ciedata2_10001);
LCol(LCol<0)=0;LCol(LCol>1)=1;
locus.CData=LCol;
axis square;

load('SAPS_SelectableColoursGamut.mat')
s=scatter(u(1:50:end),v(1:50:end),20,occ(1:50:end));
view(2)
axis('equal')
xlim([0.14 0.25]),ylim([0.41 0.52]) %close to selectable gamut boundary
colorbar
xlabel('u'''),ylabel('v'''),zlabel('Y')


repeats=[3,8];

% % Test I've got the right ones
%TabletData(repeats,1:3,:);

for i=length(sheets)
%for i=2:10
%for i=11:19
%for i=20:28
    plot(TabletData(repeats,11,i),TabletData(repeats,12,i),'k')
    scatter(TabletData(repeats(2),11,i),TabletData(repeats(2),12,i),kstd_mean(i)*10000,'k','filled')
    text(TabletData(repeats(2),11,i)+0.003,TabletData(repeats(2),12,i),files(i).participant)
end

% % Add Baseline data (essentially zero apart from the difference from
% % input variability)
% plot(TabletData(repeats,11,1),TabletData(repeats,12,1),'r')
% scatter(TabletData(repeats(2),11,1),TabletData(repeats(2),12,1),'r')


