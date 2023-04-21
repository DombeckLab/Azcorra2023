
function [ave,dataOrdNorm,temp405,dataOrd,dataNotOrd,sigs,sigsAV] = plotAverages_paper(data, depth, trig, RG, good, onoff, VA, neg, pos, idx, sp, smoothYN)
% necesary inputs: data, depth, trig, RG, good
hz = 100;

%% check RG variable is correct
if ~(strcmp(RG, 'R') || strcmp(RG, 'G') || strcmp(RG, 'L'))
    error("RG input must be 'R', 'G', or 'L'") 
end

%% check that good matches data
if length(good) ~= size(data,1)
    error("'good' input varaible must match size of data")
end


%% set defaults
if isempty(onoff)
    onoff = 1;
elseif ~(onoff == 1 || onoff == -1)
    error('onoff input must be 1, -1, or ~')
end
if isempty(VA)
    VA = 'A';
elseif ~(strcmp(VA, 'V') || strcmp(VA, 'A') || strcmp(VA, 'L'))
    error("VA input must be 'V', 'A', 'L', or ~")
end  
if isempty(neg)
    neg = 50;
end
if isempty(pos)
    pos = 100;
end

%% plot each mice averages (or just get values)
nThresh = 2; % at least this many triggers to accept average %%%%%%%%%%

good2 = find(good == 1);
bad = false(size(good2,1),1); % not enough n
if ~exist('sp','var')
    figure(99);
    clf
    if length(good2)/2 < 4
        sM = 2;
        sN = 4;
    else
        sM = 3;
        sN = ceil(length(good2)/3);
    end
end
means = nan(3,length(good2),neg+pos);
sigs = cell(length(good2),1);
sigsAV = cell(length(good2),1);
for e = 1:length(good2)
    if ~exist('sp','var') % never used
        subplot(sM,sN,e)
        [m,~,~,n] = plotSignals_paper(data.data{good2(e)},depth,trig,RG,1,VA,neg,pos);
        title(['n = ' num2str(n)])
    else
        [m,~,s,n] = alignSignals_paper(data.data{good2(e)},depth,trig,1,VA,neg,pos);
    end
    if n > nThresh
        if strcmp(RG,'R') && exist('sp','var')
            means(:,e,:) = m([1 4 5],:);
            sigs{e} = cell2mat(s(:,4));
        elseif strcmp(RG,'L') && exist('sp','var')
            means(:,e,:) = [m(1,:);m(6,:);m(6,:)];
            sigs{e} = cell2mat(s(:,6));
        elseif strcmp(RG,'G') && exist('sp','var')
            means(:,e,:) = m(1:3,:);
            sigs{e} = cell2mat(s(:,2));
        end
        sigsAV{e} = cell2mat(s(:,1));
    else
        if ~exist('sp','var')
            cla
            yyaxis left
            cla
        end
        bad(e) = true;
    end
end
ave = nan(3,neg+pos);
sem = nan(3,neg+pos);

%% smooth?
if exist('smoothYN','var') 
    if smoothYN == 1
        smoothWin = 5; 
        for i = 1:size(means,1)
            for j = 1:size(means,2)
                temp = means(i,j,:);
                temp = fastsmooth(temp,smoothWin);
                temp(temp == 0) = means(i,j,temp == 0); 
                means(i,j,:) = temp;
            end
        end
    end
end


%% plot average for all mice
if ~exist('sp','var')
    figure
    subplot(1,2,1);
else
    subplot(sp(1),sp(2),sp(3));
end

for i = 1:3
    a = mean(means(i,:,:),2,'omitnan');
    %reshape(a,1,neg+pos);
    ave(i,:) = a;
    s = std(means(i,:,:),[],2,'omitnan')/sqrt(length(good2));
    %reshape(s,1,neg+pos);
    sem(i,:) = s;
end
exp = data.Exp{find(good,1)};
plotSignals2_paper(ave,sem,'G',VA,neg,pos,hz,exp);
yyaxis left
ylim([-0.6 0.6])
if strcmp(VA,'A')
    horizontalLine(0);
end
yyaxis right
%ylim([-0.05 0.5])


if ~exist('sp','var')
    subplot(1,2,2)
else
    subplot(sp(1),sp(2),sp(3)+1)
    cla
end
    
% plot heatmap - normalize
peakHeatmap = {'peaks'};
if startsWith(trig,peakHeatmap) 
    dataNotOrd = reshape(means(1,:,:),length(good2),neg+pos);
    temp405 = [];
else
    dataNotOrd = reshape(means(2,:,:),length(good2),neg+pos);
    temp405 = reshape(means(3,:,:),length(good2),neg+pos);
end

%mark empty recordings to exclude from mice later
badMice = any(isnan(dataNotOrd),2);

% reorder if idx present
if exist('idx','var') && ~isempty(idx)
    dataOrd = dataNotOrd(idx,:);
else
    dataOrd = dataNotOrd;
end
% remove empty recordings
bad = any(isnan(dataOrd),2);
dataOrd(bad,:) = [];
% normalize heatmap - from 0 to 1
if startsWith(trig,{'Rew','Light','Air'}) % do not normalize!
    dataOrdNorm = dataOrd;
else
    %dataOrdNorm = dataOrd - min(dataOrd,[],2);
    %dataOrdNorm = dataOrdNorm ./ max(dataOrdNorm,[],2);
    
    dataOrdNorm = dataOrd - mean(dataOrd,2);
    dataOrdNorm = dataOrdNorm ./ max(abs(dataOrdNorm),[],2);
end
%plot
imagesc(dataOrdNorm);
xticks([])

% title and n
mice = data.Mouse(good,:);
mice = strip(mice,'right','L');
mice = mice(~badMice);
mice = unique(mice);
exp = unique(data.Exp(good,:));
% type = unique(data.Type(good,:));
% if length(type) > 1
%     type2 = type{1};
%     for i = 2:length(type)
%         type2 = [type2 ', ' type{i}];
%     end
% else 
%     type2 = type{1};
% end

title([exp{1} ' - ' num2str(length(mice)) 'm, ' num2str(sum(~bad)) ' r']);





    
    
    
