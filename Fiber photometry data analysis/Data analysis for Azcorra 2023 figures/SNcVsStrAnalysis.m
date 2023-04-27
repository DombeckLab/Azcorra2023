function SNcVsStrAnalysis(data6,deconvolve)

%% defaults, check vars
if ~exist('deconvolve','var') || deconvolve ~= 1
    deconvolve = 0;
    v5 = 'no';
else
    v5 = 'yes';
end
disp(['    Deconvolve? ' v5]);

%% load data
filePath = matlab.desktop.editor.getActiveFilename;
i = strfind(filePath, '\');
filePath = filePath(1:i(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder');


%% switch channels for recordings where str is in chG and snc in chR (few)
data6(data6.dup,:) = []; %remove previous dups
flip = false(size(data6,1),1);
flip(data6.exStr(:,1) == 1 & data6.exStr(:,2) == 0) = true;
flip2 = find(flip);
data6 = flipRG(data6,flip2);
data6.flip(flip2) = true;
% sort
data6 = sortrows(data6,1);


%% get exclusion criteria
good = double(data6.exBothCh & data6.exSNcStr & data6.exRG405 & sum(data6.Bad405All,2) == 0);

% export table
export = data6(logical(good),[1 2 3 4 5 6 7 9 10 11 12 16]);
export = [cell(size(export,1),1), export];
export.Properties.VariableNames{1} = 'Idx';
for i = 1:size(export,1)
    export.Idx{i} = export.data{i}.Properties.RowNames{1};
end
export.data = [];
exportFig4 = export;
save([dataProcessingFolder '\exportFig4_p1.mat'], 'exportFig4');

%% asign numbers for different subtypes
subtypes = {'VGlut' 'Calb' 'Raldh' 'Anxa' 'Dat'};
load colors.mat C
colors = [C.red{1}; C.orange{1}; C.green{1}; C.aqua{1}; C.grey{1}];
s = length(subtypes);
for i = 1:s
    sIdx = strcmp(data6.Exp,subtypes{i});
    good(good > 0 & sIdx) = i;
end


%% get crosscorr snc vs str
% Fig. 7E and Fig. S10C
win = 100;
ccData = nan(size(data6,1),win*2+1);
ccData405 = nan(size(data6,1),win*2+1);
figure
set(gcf, 'Position',  [1, 42, 320, 954])
for i = 1:s
    temp = crossCorr2SNcStr(data6, good == i, deconvolve, win,[s,2,i*2-1]);
    temp2 = cell2mat(temp(:,1));
    ccData(good == i,:) = temp2(good == i,:);
    temp2 = cell2mat(temp(:,2));
    ccData405(good == i,:) = temp2(good == i,:);
    subplot(s,2,i*2)
    caxis([-0.1,0.7])
end

%% Distribution of cross correlation peaks by subtype - histograms
% Fig. 7F and Fig. S10D
peaks = max(ccData,[],2);

bins = linspace(-0.1,1,10);
figure
set(gcf, 'Position',  [800 780 680 200])
for i = 1:s
    subplot(1,s,i)
    histogram(peaks(good == i),bins,'FaceColor',colors(i,:),'EdgeColor','k');
    view([90 -90])
    xlabel('Peak cross corr')
    ylabel('# recordings')
    
    xlim([-0.1 1]);
    yl = ylim;
    ylim([0 yl(2)+1])
    
    if i ~= 5 %Dat
    	p = ranksum(peaks(good == i),peaks(good == 5));
        title([subtypes{i} ' - p=' num2str(p)])
        disp(p)
    else
        title(subtypes{i})
    end
    
    disp([subtypes{i} ' - mean = ' num2str(round(mean(peaks(good==i)),2))])
end
peaksD = peaks(good == 5);
save([dataProcessingFolder '\PeaksSNcStrDat.mat'], 'peaksD');


%% Split Aldh to get only dStr vs vSNc in relative numbers
% Fig. S10E
threshStr = 2.25;
goodStrR = good == 3 & cell2mat(data6.depthR) <= threshStr;
goodStrD = good == 5 & cell2mat(data6.depthR) <= threshStr;
goodStr = {goodStrD;goodStrR};

colors = [C.grey{1}; C.green{1}];
tit = {'D-dStr' 'R-dStr'};
lineSt = {'-' '--' '-.' ':'};
xlabels = {'Normalized depth' 'Relative depth (mm)'};

figure
set(gcf, 'Position',  [854   336   393   636])
for j = 1:2
    idx = find(goodStr{j});
    recDay = cell(size(idx));
    for i = 1:length(idx)
        e = idx(i);
        try
            recDay{i} = [data6.Mouse{e}, data6.Date{e}{1}];
        catch
            recDay{i} = [data6.Mouse{e}, data6.Date{e}];
        end
        idx(i,3) = data6.depthG{e};
    end
    [~,~,xx] = unique(recDay);
    idx(:,2) = xx;

    for i = 1:length(idx)
        e = idx(i);
        g = data6.data{e}.chGreen{1};
        r = data6.data{e}.chRed{1};
        bad = isnan(r)|isnan(g);
        r(bad) = [];g(bad)=[];
        g = fastsmooth(g,10);
        r = fastsmooth(r,10);
        c2 = crosscorr(r,g,100);
        peak = max(c2);
        idx(i,4) = peak;
    end
    
    for i = 1:max(idx(:,2))
        % stretch depth
        idx2 = find(idx(:,2) == i);
        idx(idx2,5) = linspace(0,1,length(idx2));
        % relative depth
        idx(idx2,6) = idx(idx2,3) - mean([max(idx(idx2,3)),min(idx(idx2,3))]);
    end
    
    for ss = 1:2
        subplotMN(2,2,ss,j)
        hold on
        for i = 1:max(idx(:,2))
            tempY = idx(idx(:,2) == i,4);
            tempX = idx(idx(:,2) == i,4+ss);
            if length(tempY) > 1
                plot(tempX,tempY,[lineSt{rem(i-1,4)+1} 'o']','LineWidth',1.5,'Color',colors(j,:),'MarkerFaceColor',colors(j,:))
            end
        end
        title(tit{j})       
        if ss == 1
            xlim([-0.2 1.2])
            xticks([0 1])
            xticklabels({'d-SNc' 'vSNc'})
        else
            xlim([-0.55 0.55])
        end
        
        ylabel('peak cross corr')
        ylim([-0.1 0.8])
        xlabel(xlabels{ss})
    end
end


%% Distribution of cross correlation peaks by subtype - histograms
% only ALDH dorsal str vs ventral snc
% Only in bioRxiv paper
threshSNcRel = 0.0999;%threshSNcRel = 0.1999;

goodX = good;
goodX(idx(idx(:,6) >= threshSNcRel,1)) = 31;

peaks = max(ccData,[],2);

bins = linspace(-0.1,1,10);
figure
set(gcf, 'Position',  [1480 780 130 200])
histogram(peaks(goodX == 31),bins,'FaceColor',C.green{1},'EdgeColor','k');
view([90 -90])
xlabel('Peak cross corr')
ylabel('# recordings')

xlim([-0.1 1]);
yl = ylim;
ylim([0 yl(2)+1])

p = ranksum(peaks(goodX == 31),peaks(good == 5));
title([subtypes{3} ' - p=' num2str(p)])
disp(p)

disp(['   Raldh-split - mean = ' num2str(round(mean(peaks(goodX == 31)),2))])


%% for aldh only dorsalStr-ventralSNc - mean cross corr and heatmap
% Only in bioRxiv paper
win = 100;
figure
crossCorr2SNcStr(data6, goodX == 31, deconvolve, win,[1,2,1]);
subplot(1,2,2)
caxis([-0.1,0.7])

%% find raw examples
% Fig. 7C and Fig. S10A
bestMice = {'VGlut-2292-20220401-0002' 'Calb-1792-20210415-0002' 'Raldh-2071-20211201-0003'...
            'Raldh-1849-20210623-0004' 'Anxa-J075-20220308-0004' 'Dat-2111-20211026-0001'};
xlims = [149 134 571 868 300 504];
colors = [C.red{1}; C.orange{1}; C.green{1}; C.green{1}; C.aqua{1}; C.grey{1}];
colorsD = [C.darkRed{1}; C.darkOrange{1}; C.darkGreen{1}; C.darkGreen{1}; C.darkAqua{1}; [0 0 0]];

plotWin = 60; % in sec
clear x
bestIdx = nan(5,1);
% find miceIdx
for i = 1:6
    id = findMouse(data6,bestMice{i});
    if size(id,1) > 1
        cont = 1;
        j = 1;
        while cont == 1
            if (strcmp(data6.chG{id(j,1)},'snc') && snc == 1) || (~strcmp(data6.chG{id(j,1)},'snc') && snc == 0)
                cont = 0;
                bestIdx(i) = id(j,1);
            else
                j = j+1;
            end
        end
    else
        bestIdx(i) = id(1);
    end
end
% plot
vars = {'chGreen' 'chRed'};
figure
set(gcf, 'Position',  [406         291        1040         678])
for i = 1:5
    e = bestIdx(i);
    colors2 = [colorsD(i,:);colors(i,:)];
    for v = 1:length(vars)
        subplotMN(5,5,i,1:4);
        hold on
        temp = data6.data{e}.(vars{v}){1};
        temp = temp * data6.norm{e}(v); %de normalize
        temp(isnan(temp)) = temp(circshift(isnan(temp),1));
        temp2 = fastsmooth(temp,10);
        temp_405 = data6.data{e}.([vars{v} '405']){1};
        temp_405(isnan(temp_405)) = temp_405(circshift(isnan(temp_405),1));
        temp2_405 = fastsmooth(temp_405,10);
        if i == 5
            temp3 = (temp2-min(temp2))./(max(temp2)-min(temp2)); % anxa section huge trans
            temp3_405 = (temp2_405-min(temp2))./(max(temp2)-min(temp2));
        elseif i == 2
            temp3 = (temp2-min(temp2))./(max(temp2)-min(temp2))*1.6-0.1; % calb small trans
            temp3_405 = (temp2_405-min(temp2))./(max(temp2)-min(temp2))*1.6-0.1; 
        else
            temp3 = (temp2-min(temp2))./(max(temp2)-min(temp2))*1.2-0.1;
            temp3_405 = (temp2_405-min(temp2))./(max(temp2)-min(temp2))*1.2-0.1;
        end
        plot(x(temp3,100),temp3_405+v-1,'Color',C.lightBlue{1});
        plot(x(temp3,100),temp3+v-1,'Color',colors2(v,:));
    end
    title(bestMice{i})
    xlim([xlims(i) xlims(i)+plotWin])
    
    % plot cross corr
    subplotMN(5,5,i,5)
    plot((-win:win)/100,ccData405(e,:),'Color',C.lightBlue{1},'LineWidth',2);
    hold on
    plot((-win:win)/100,ccData(e,:),'Color',colors2(2,:),'LineWidth',2);
    xlabel('Lag (s)')
    ylim([-0.1 0.8])
    verticalLine(0);
end

% plot Dat example with velocity and acc
i = 6;
e = bestIdx(i);
vars = {'chGreen' 'chGreen405' 'chRed' 'chRed405' 'Acceleration' 'chMov'};
vals = cell(6,1);
for v = 1:length(vars)
    temp = data6.data{e}.(vars{v}){1};
    if ~ismember(vars{v},{'Acceleration' 'chMov'})
        temp = temp * data6.norm{e}(ceil(v/2)); %de normalize
    end
    temp(isnan(temp)) = temp(circshift(isnan(temp),1));
    temp2 = fastsmooth(temp,10);
    vals{v} = temp2;
end

figure
set(gcf, 'Position',  [54         531        1006         420])
subplot(1,5,1:4)
hold on
plot(x(temp2,100),vals{2},'Color',C.blue{1});
plot(x(temp2,100),vals{1},'Color',[0 0 0]);
plot(x(temp2,100),vals{4}+1,'Color',C.blue{1});
plot(x(temp2,100),vals{3}+1,'Color',C.grey{1});

plot(x(temp2,100),vals{5}/max(vals{5})-0.5,'Color',C.grey{1});
plot(x(temp2,100),vals{6}/max(vals{6})*0.9-1.8,'Color',[0 0 0]);

title(bestMice{i})
xlim([xlims(i) xlims(i)+plotWin])

% plot cross corr
subplot(1,5,5)
plot((-win:win)/100,ccData405(e,:),'Color',C.lightBlue{1},'LineWidth',2);
hold on
plot((-win:win)/100,ccData(e,:),'Color',C.grey{1},'LineWidth',2);
xlabel('Lag (s)')
ylim([-0.1 0.8])
verticalLine(0);


end










