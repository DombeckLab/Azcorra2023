function RewardPuffAnalysis(data6,snc)
% This script plots many different figures for the analysis of REWARDS and
% AIR PUFF SIGNALING of different DA subtypes as presented in Azcorra et 
% al. 2023.
% Each section is labelled to show which figure and panel it generates
% within this manuscript. 

% The 'snc' variable indicates whether the analysis will be run on only
% axonal recordings (snc = 0) or only somatic recordings (snc = 1).

% This scripts was written for Windows at MATLAB R2021a

%% Set default variables - if snc isn't present, default is no
if ~exist('snc','var') || snc ~= 1
    sncYN = 'no';
    snc = 0;
else
    sncYN = 'yes';
end
disp(['    SNc? ' sncYN]);

% set type of heatmap normalization for dff-acceleration cross-corr
normHeat = 2; %2 = scale maintaining 0
hz = 100;

%% set DATA PROCESSING FOLDER (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
i = strfind(filePath, '\');
filePath = filePath(1:i(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder');

%% Switch channels for recordings where str/snc is in chR, and duplicate 
% and switch for those with double str 
% From here on only the chG will be included in analysis, and so any
% recordings of interest found in ChR need to be switched

data6(data6.dup,:) = []; %remove previous dups if present

% duplicate any recordings where chR is not empty
dupIdx = ~strcmp(data6.chR,'o'); 
cc = data6(dupIdx,:);
cc.dup = not(cc.dup); %label copies as dups
data6 = [data6; cc]; %append copies at end of data6
dupNum = sum(dupIdx); %number of copied rows for next steps

% flip chR and chG in duplicated recordings (so all non-empy recording are
% now in chG in one row)
flipIdx = false(size(data6,1),1);
flipIdx(end-dupNum+1:end) = true;
flipIdx2 = find(flipIdx); %which recordings to flip
data6 = flipRG(data6,flipIdx2); %separate script
data6.flip(flipIdx2) = true; % label as flipped

% sort rows
data6 = sortrows(data6,1);

%% get exclusion criteria for reward/air puff analysis
% get either all snc recordings or all str recordings
if snc == 0
    strSNc = data6.exStr(:,1) == 1; % only str
else
    strSNc = data6.exStr(:,1) == 0; % only snc
end
% include recordings in the correct location, with good sig2noise, and 
% during which rewards were delivered. These were calculated in the script 
% exclusionCriteria_allFigs.m
good = double(strSNc & data6.exSig2Noise(:,1) & strcmp(data6.RunRew,'rew'));

% export info - for summary table showing which mice were included for what
% figure/analysis
export = data6(logical(good),[1 2 3 4 5 6 7 9 10 11 12 16]);
export = [cell(size(export,1),1), export];
export.Properties.VariableNames{1} = 'Idx';
for i = 1:size(export,1)
    export.Idx{i} = export.data{i}.Properties.RowNames{1};
end
export.data = [];
if snc == 0
    exportRewStr = export;
    save([dataProcessingFolder '\exportRewStr.mat'], 'exportRewStr');
else
    exportRewSNc = export;
    save([dataProcessingFolder '\exportRewStr.mat'], 'exportRewSNc');
end

%% for each recording, exclude the bad 405 data points from data.Bad405G
for rec = 1:size(data6,1)
    if good(rec) == 1
        skip = find(strcmp(data6.data{rec}.Properties.VariableNames,'Bad405G'));
        for var = 2:size(data6.data{rec},2)
            if var ~= skip
                data6.data{rec}.(var){1}(data6.data{rec}.Bad405G{1}) = [];
            end
        end
        try
        data6.data{rec}.(skip){1}(data6.data{rec}.Bad405G{1}) = [];
        catch
            disp(rec)
        end
    end
end

%% also get goodRun - which recordings pass locomotion analysis criteria
% necessary for simultaneous analysis of locomotion and rew/puff signaling
goodRun = double(strSNc & data6.exSig2Noise(:,1) & data6.ex405Acc(:,1) & data6.exRun & data6.Bad405All(:,1) == 0 & strcmp(data6.RunRew,'rew'));

%% asign numbers for different subtypes (for both good and goodRun)
% this groups them for plotting.
% also set identifying colors for each
subtypes = {'VGlut' 'Calb' 'Raldh' 'Anxa' 'Dat'};
load colors.mat C
colors = [C.red{1}; C.orange{1}; C.green{1}; C.aqua{1}; C.grey{1}];
good(good == 1) = 6; % any subtypes different from above will be 6
goodRun(goodRun > 0) = 6;
for s = 1:length(subtypes)
    sIdx = strcmp(data6.Exp,subtypes{s});
    good(good > 0 & sIdx) = s;
    goodRun(goodRun > 0 & sIdx) = s;
end

%% Get reward trig ave to sort heatmaps 
% set parameters for trig averages - these will be used on all future trig ave
AV = 'A'; %show acceleration - change to 'V' to see velocity averages
neg0 = 100; % bins to show in trig ave BEFORE trigger point (100 bin = 1s)
pos0 = 100; % bins to show in trig ave AFTER trigger point (100 bin = 1s)
neg = neg0 + 5; % these are necessary because edge effect of fill function. Plot 5 extra bins but hide
pos = pos0 + 5;
% save trig ave to sort
rewData = nan(size(data6,1),neg0+pos0+1); %preallocate
% get rew trig ave
figure
for s = 1:length(subtypes)
	[~,~,~,~,trigAves] = plotAverages_paper(data6, '1', 'Reward', 'G', good == s, 1, AV, neg, pos, randperm(sum(good == i)),[1,2,1],1);
	rewData(good == s,:) = trigAves(:,neg-neg0:neg+pos0); %save only neg0:pos0
end
close() %this figure is unnecessary, we just need the data sorting. Will be replotted later

% Calculate reward response as baseline corrected integral under the curve
% after reward
% get baseline - mean of window before stimulus
base = mean(rewData(:,1:neg0),2);
% subtract baseline from traces
baseCorr = rewData - base;
% get response after stimulus as integral under curve
tt = (1:size(baseCorr,2))/hz; % time vector for integral calculation below
rewResp = trapz(tt(neg0+1:neg0+pos0), baseCorr(:,neg0+1:neg0+pos0),2);

% sort recordings by reward response size within subtype (both good and
% goodRun)
idxRew = nan(size(data6,1),1);
idxRewRun = nan(size(data6,1),1);
for s = 1:length(subtypes)
    rewRespSub = rewResp(good == s);
    [~,rewRespSubSort] = sort(rewRespSub);
    idxRew(good == s) = rewRespSubSort;
    
    rewRespSubRun = rewResp(goodRun == s);
    [~,rewRespSubRunSort] = sort(rewRespSubRun);
    idxRewRun(goodRun == s) = rewRespSubRunSort;
end


%% Crosscorr trans vs acc/vel 
% to check relationship between locomotion and reward/air puff signaling

% make variable to store cross-corr and trig averages for all recordings
% for later use
allData = cell(1,9);
% set parameters for cross-corr
win = 100; 
allData{1} = nan(size(data6,1),win*2+1); % preallocate
% plot
figure
set(gcf, 'Position',  [1 42 1590 954])
for s = 1:length(subtypes)
    ccDffAcc = crossCorr2(data6, goodRun == s, 'G', 'A', win,[length(subtypes),8,s*8-7],idxRewRun(goodRun == s),normHeat);
    ccDffAcc = cell2mat(ccDffAcc(:,1));
    allData{1}(goodRun == s,:) = ccDffAcc(goodRun == s,:);
end
for s = 2:8:length(subtypes)*8
    subplot(length(subtypes),8,s)
    if normHeat == 0
        caxis([-0.15 0.1])
    else
        caxis([-0.8,0.8])
    end
    subplot(length(subtypes),8,s-1)
    title('Cross corr with Acc')
end
disp('cross-corr done')

% PCA on cross corr - use PCs from striatal locomotion data (see p2_f1_Run)
load([dataProcessingFolder '\PCA_striatum.mat'], 'coeff','mu'); % mu will be zero if centered is false
ccData = allData{1};
scores = (ccData-mu)*coeff;
% Get integral of cross-corr at positive lags - to use later to plot vs
% reward response for Fig. S1I
ccInt = nan(size(ccData,1),1);
for rec = 1:size(ccData)
    if good(rec) ~= 0
        cc = ccData(rec,:);
        tt = (1:size(cc,1))/hz; % time vector for integral calculation below
        ccInt(rec) = trapz(tt,cc(ceil(length(cc)/2):end));
    end
end

%% Triggered averages on reward, air puff and light
% Fig. 4D,F and Fig. S1G, S8A,C

% set parameters for triggered averages - other parameters set above
trigs = {'Reward', 'AirPuff', 'Light'}; 
% save cross-corr and triggered averages for later use
% trig aves for each recording in allData, and individual events (470 and
% acceleration) in allDataEvents470 and allDataEventsAcc
allDataEvents470 = cell(size(data6,1),7);
allDataEventsAcc = cell(size(data6,1),7);
for t = 2:1+length(trigs)
    allData{t} = nan(size(data6,1),neg0+pos0+1); %preallocate
end
% plot - on same figure as cross-corr above
nn = cell(5,3); %save number of events
for s = 1:length(subtypes)
    for t = 1:length(trigs)
        [~,~,~,~,trigAves,events470,eventsAV] = plotAverages_paper(data6, '1', trigs{t}, 'G', good == s, 1, AV, neg, pos, idxRew(good == s),[length(subtypes),8,s*8-(7-t*2)],1);
        allData{t+1}(good == s,:) = trigAves(:,neg-neg0:neg+pos0);
        allDataEvents470(good == s,t) = events470;
        allDataEventsAcc(good == s,t) = eventsAV;
        % get number of events averaged per recording
        nn{s,t} = nan(size(events470,1),1);
        for rec = 1:length(nn{s,t})
            nn{s,t}(rec) = size(events470{rec},1);
        end
        nn{s,t}(isnan(nn{s,t}) | nn{s,t} == 0) = [];
    end
end
% set axis lims etc
ylimsR = [-0.02 0.65; -0.02 0.65; -0.02 0.65];
for s = 1:length(subtypes)
    for t = 1:length(trigs)
        subplotMN(length(subtypes),8,s,t*2+1)
        title(trigs{t})
        yyaxis left
        ylim([-2 2])
        yyaxis right
        ylim(ylimsR(t,:))
        xlim([-neg0 pos0]/hz)
        verticalLine(0);
        
        subplotMN(length(subtypes),8,s,t*2+2)
        xlim([5.5 pos0+neg0+5.5])
        caxis([0 0.6])
    end
end
disp('rew/puff/light done')

% get average number of events + STD (not SEM)
for t = 1:length(trigs)
    eventNum = cell2mat(nn(:,t));
    mm = mean(eventNum);
    ss = std(eventNum);
    disp(['average events ' trigs{t} ' = ' num2str(mm) ' +- ' num2str(ss) ' std']);
end


%% Trig ave on small rewards, large rewards, and rewards at rest
% Fig. S8E

% set parameters for triggered averages - other parameters set above
trigs2 = {'RewardShort', 'RewardLong', 'RewardRest'};
% plot (new figure)
figure
set(gcf, 'Position',  [1 42 1590 954])
for t = length(trigs)+2:length(trigs2)+length(trigs)+1
    allData{t} = nan(size(data6,1),neg0+pos0+1); %preallocate
end
for s = 1:length(subtypes)
    for t = 1:length(trigs2)
        [~,~,~,~,trigAve,events470,eventsAV] = plotAverages_paper(data6, '1', trigs2{t}, 'G', good == s, 1, AV, neg, pos, idxRew(good == s),[length(subtypes),8,s*8-(7-t*2)],1);
        allData{t+length(trigs)+1}(good == s,:) = trigAve(:,neg-neg0:neg+pos0);
        allDataEvents470(good == s,t+length(trigs)) = events470;
        allDataEventsAcc(good == s,t+length(trigs)) = eventsAV;
    end
end
% set axis lims etc
ylimsR = [-0.02 0.65; -0.02 0.65; -0.02 0.65];
for s = 1:length(subtypes)
    for t = 1:length(trigs2)
        subplotMN(length(subtypes),8,s,t*2+1)
        title(trigs2{t})
        yyaxis left
        ylim([-2 2])
        yyaxis right
        ylim(ylimsR(t,:))
        xlim([-neg0 pos0]/hz)
        verticalLine(0);
        
        subplotMN(length(subtypes),8,s,t*2+2)
        xlim([5.5 pos0+neg0+5.5])
        caxis([0 0.6])
    end
end
disp('rew sizes/rest done')


%% Plot average licking for all reward trig ave (in left column of same figure above)
% Fig. 4E and Fig. S8B
allData{8} = nan(size(data6,1),neg0+pos0+1); %preallocate
for s = 1:length(subtypes)
    [~,~,~,~,trigAve] = plotAverages_paper(data6, '1', 'Reward', 'L', good == s, 1, AV, neg, pos, idxRew(good == s),[length(subtypes),8,s*8-7],1);
    allData{8}(good == s,:) = trigAve(:,neg-neg0:neg+pos0);
end
% set axis lims etc
for s = 1:length(subtypes)
    subplotMN(length(subtypes),8,s,1)
    title('Licking')
    xlim([-neg0 pos0]/hz)
    yyaxis left
    ylim([-2 2])
    yyaxis right
    ylim([0 1])
    verticalLine(0);
    subplotMN(length(subtypes),8,s,2)
    xlim([5.5 pos0+neg0+5.5])
    colormap('gray') %set colormap for licking to gray. Will affect other panels, so change for figure
    caxis([-0.3 1])
end

%% remove rewards at rest from rewards 
% this will be used for comparison between rewards at rest and rewards not
% at rest
for rec = 1:size(data6,1)
    try
        rewAll = data6.data{rec}.Reward{1};
        rewRest = data6.data{rec}.RewardRest{1};
        rewNoRest = rewAll &~ rewRest;
        data6.data{rec}.RewardNoRest{1} = rewNoRest;
    catch 
    end
end
% trig ave on rewards not at rest
figure
allData{9} = nan(size(data6,1),neg0+pos0+1); %preallocate
for s = 1:length(subtypes)
    [~,~,~,~,trigAve,events470,eventsAV] = plotAverages_paper(data6, '1', 'RewardNoRest', 'G', good == s, 1, AV, neg, pos, idxRew(good == s),[length(subtypes),8,s*8-7],1);
    allData{9}(good == s,:) = trigAve(:,neg-neg0:neg+pos0);
    allDataEvents470(good == s,7) = events470;
    allDataEventsAcc(good == s,7) = eventsAV;
end
close() %don't need the figure


%% mouse ID heatmaps
% make a bar showing which recording in the above heatmaps is from what
% mouse by color-coding them for each subtype. 
figure
trigsAll = [trigs trigs2];
for s = 1:length(subtypes)
    for t = 2:7
        subplotMN(length(subtypes),6,s,t-1)
        mouseHeatmapIdx = idxRew(good == s); %get order of recordings in heatmap for subtype
        mice = data6.Mouse(good == s); % get mice
        mice = strip(mice,'right','L'); % one mouse was recorded from left hemisphere one day, but same mouse
        [miceUnique,~,recIdx] = unique(mice); % get unique mice and which recordings come from each
        %only the ones with data for this particular trigger event
        keep = ~isnan(allData{t}(good == s));         
        miceIdOrder = recIdx(mouseHeatmapIdx); % get which mice match the heatmap for that subtype
        keep = keep(mouseHeatmapIdx);
        miceIdOrder = miceIdOrder(keep);
        % plot
        colorsU = hsv(length(miceUnique))/4*3; %colormap - scale the hsv colormpa for the number of unique mice in subtype
        imagesc(miceIdOrder)
        colormap(colorsU)
        title([subtypes{s} ' - ' trigsAll{t-1}])
    end
end


%% Split rew/air puff events into two quartiles based on accelerations
% this is used to determine whether the amplitude of the movement that
% accompanies the stimulus affects the response
% Fig. 5C,D,E,F

% For each recording, split events into two sets based on the amplitude of
% the deceleration (calculated as the integral within a window after the 
% trigger point), using the median to split. Then average each set of 
% events to get two triggered averages per recording (average acc and DF/F
% traces). Then average together all recordings' average for small dec, and
% all averages for large dec, and plot.
% Also calculate the amplitude of the response using the integral of the
% DF/F trace within a window and plot the response for small dec vs large
% dec for each recording and subtype.

accWin = 75; %integral of acceleration within 0.75s of trigger point
accQuart = cell(size(data6,1),4,2); % save acc and 470 - 1 = DF/F for largest dec, 2 = DF/F for smallest dec, 3 = Acc for largest dec, 4 = Acc for smallest dec
for t = 1:2 %reward and air puff
    for rec = 1:size(allDataEvents470,1)
        if ~isempty(allDataEvents470{rec,t})
            event470 = allDataEvents470{rec,t}; 
            eventAcc = allDataEventsAcc{rec,t};
            % calculate integral of acc
            tt = (1:size(eventAcc,2))/hz; % time vector for integral calculation below
            tempA_d = trapz(tt(neg:neg+accWin),eventAcc(:,neg:neg+accWin),2); 
            % split acc in half using median
            mm = median(tempA_d);
            idxTopAcc = tempA_d <= mm; % greatest half of decelerations (smallest values)
            % save average of all events within each half of events based
            % on decelerations. Save both 470 trace and acc
            accQuart{rec,1,t} = mean(event470(idxTopAcc,:),1); %top DF/F
            accQuart{rec,2,t} = mean(event470(~idxTopAcc,:),1);
            accQuart{rec,3,t} = mean(eventAcc(idxTopAcc,:),1); %top Acc
            accQuart{rec,4,t} = mean(eventAcc(~idxTopAcc,:),1);
        end
    end
end
% plot triggered averages for small vs large dec, as well as amplitude of
% the responses for each
dffWin = 50; %integral of DF/F within 0.5s of trigger point
pureSubs = [1 2 4]; %only plot pure subtypes, Vglut, Calb and Anxa
colorsDark = [C.darkRed{1}; C.darkOrange{1}; C.darkGreen{1}; C.darkAqua{1}; 1 1 1];
tAxis = (-neg+1:pos)/hz; %time axis for plotting
figure
for ss = 1:3
    for t = 1:2 %rew and puff
        % plot acc averages
        subplotMN(3,6,ss,t*3-2)
        plot(tAxis,mean(cell2mat(accQuart(good == pureSubs(ss),3,t)),1,'omitnan'),'Color',colorsDark(pureSubs(ss),:),'LineWidth',1.5)  %big acc (Acc)
        hold on
        plot(tAxis,mean(cell2mat(accQuart(good == pureSubs(ss),4,t)),1,'omitnan'),'Color',colors(pureSubs(ss),:),'LineWidth',1.5)   %small acc (Acc)
        if t == 1
            ylim([-3 1.5])
        else
            ylim([-4 2])
        end
        xlim([-1 1])
        verticalLine(0);
        % plot DF/F averages
        subplotMN(3,6,ss,t*3-1)
        plot(tAxis,mean(cell2mat(accQuart(good == pureSubs(ss),1,t)),1,'omitnan'),'Color',colorsDark(pureSubs(ss),:),'LineWidth',1.5)  %big acc (DF)
        hold on
        plot(tAxis,mean(cell2mat(accQuart(good == pureSubs(ss),2,t)),1,'omitnan'),'Color',colors(pureSubs(ss),:),'LineWidth',1.5)   %small acc (DF)
        ylim([0 0.6])
        xlim([-1 1])
        verticalLine(0);
        
        % Amplitude of the response and plot for small vs large dec. 
        % Calculate the integral of the DF/F after the trigger and subtract
        % the integral of the DF/F before the trigger for each event and
        % average per recording to get the response amplitude.
        subplotMN(3,6,ss,t*3)
        aveLarge = cell2mat(accQuart(good == pureSubs(ss),1,t)); % big acc (DF) all events
        aveSmall = cell2mat(accQuart(good == pureSubs(ss),2,t)); % small acc (DF) all events
        tt = (1:size(aveLarge,2))/hz; % time vector for integral calculation below        
        baseLarge = trapz(tt(neg-dffWin:neg),aveLarge(:,neg-dffWin:neg),2); %get baseline (before trig)
        respLarge = trapz(tt(neg+5:neg+5+dffWin),aveLarge(:,neg+5:neg+5+dffWin),2); % get response (after trig)        
        changeLarge = respLarge-baseLarge; %big acc (DF)
        baseLarge = trapz(tt(neg-dffWin:neg),aveSmall(:,neg-dffWin:neg),2); %get baseline
        respSmall = trapz(tt(neg+5:neg+5+dffWin),aveSmall(:,neg+5:neg+5+dffWin),2); % get response        
        changeSmall = respSmall-baseLarge; %small acc (DF)   
        % plot responses for small vs large dec for each recording plus
        % mean and sem per subtype
        plot([1 2],[changeSmall changeLarge],'o-','Color',colors(pureSubs(ss),:))
        hold on
        errorbar([1 2],[mean(changeSmall,'omitnan') mean(changeLarge,'omitnan')],[std(changeSmall,'omitnan')/sqrt(sum(~isnan(changeSmall))) std(changeLarge,'omitnan')/sqrt(sum(~isnan(changeLarge)))],'k','LineWidth',1.5,'CapSize',20)
        ylim([-0.15 0.3])
        xlim([0 3])
        horizontalLine(0);
        % pvalue to compare amplitude of response for small vs large
        p = signrank(changeSmall,changeLarge)*6;
        title(['m1 = ' num2str(mean(changeSmall,'omitnan')) ', m2 = ' num2str(mean(changeLarge,'omitnan')) ', p = ' num2str(p) ', n = ' num2str(length(changeSmall))])    
    end
end


%% calculate responses to each stimulus as the integral of the DF/F
% Rather than calculate the response amplitude from the triggered average, 
% calculate the response event by event and then get the average for the
% session.
% The response is calculated as the integral of the DF/F in a window after 
% the event minus the integral of the DF/F in a window before the event,
% then averaged by session)
winInt = 50; %0.5 s
gap = 5; % ignore the first 5 bins after the event, as there is a bit of a delay
response = nan(size(allData{1},1),7);
for v = 1:7 %variables rew, puff, light, rew small, rew large, reward rest, reward no rest
    for rec = 1:size(allData{1},1)
        if size(allDataEvents470{rec,v},1) > 0
            eventDFF = allDataEvents470{rec,v};
            % replace nans with nearest bin
            if any(any(isnan(eventDFF))) 
                eventDFFshift = circshift(eventDFF,1,2);
                eventDFF(isnan(eventDFF)) = eventDFFshift(isnan(eventDFF));
            end
            % get response: integral after minus integral before event
            tt = (1:size(eventDFF,2))/hz;
            base = trapz(tt(neg0-winInt:neg0),eventDFF(:,neg0-winInt:neg0),2); %get baseline
            resp = trapz(tt(neg0+gap:neg0+winInt+gap),eventDFF(:,neg0+gap:neg0+winInt+gap),2); % get response
            response(rec,v) = mean(resp-base,'omitnan');
        end
    end
end


%% Plot bar graph of reward and airpuff responses by subtype
% Fig. 4G and Fig. S9C (for snc)
meanRewPuff = nan(5,2);
semRewPuff = nan(5,2);
n = nan(5,2);
for s = 1:length(subtypes)
    % get average response by subtype for bar graph
    meanRewPuff(s,1) = mean(response(good == s,1),'omitnan'); %rew
    meanRewPuff(s,2) = mean(response(good == s,2),'omitnan'); %puff
    % get s.e.m.
    n(s,1) = sum(~isnan(response(good == s,1)));
    n(s,2) = sum(~isnan(response(good == s,2)));
    semRewPuff(s,1) = std(response(good == s,1),'omitnan')/sqrt(n(s,1));
    semRewPuff(s,2) = std(response(good == s,2),'omitnan')/sqrt(n(s,2));
end
% plot bar graph for average responses, errorbars and dots for each
% recording
good2 = good;
good2(good2 == 0 | good2 >5) = nan; %make nan all idx not for subtypes
figure
for s = 1:length(subtypes)
    for i = 1:2
        subplot(1,2,i)
        hold on
        % plot bars for averages
        bar(s,meanRewPuff(s,i),'FaceColor',colors(s,:),'FaceAlpha',0.7)
        % plot errorbars
        errorbar(s,meanRewPuff(s,i),semRewPuff(s,i),'LineStyle','none','Color','k','LineWidth',1,'CapSize',10)
        % plot dots for each recording
        plotSpread(response(:,i),'distributionIdx',good2,'distributionColors',colors)
        % set lims and ticks etc
        ylim([-0.15 0.25])
        xticks(1:5)
        xticklabels(subtypes)
        if i == 1
            ylabel('Reward response (Int Norm \DeltaF/F)')
        else
            ylabel('Air Puff response (Int Norm \DeltaF/F)')
        end
    end
end

% get p-values for significance of average responses
pvalsRewPuff = nan(length(subtypes),2);
for s = 1:length(subtypes)
    for i = 1:2
        pvalsRewPuff(s,i) = signrank(response(good == s,i));
    end
end
% bonferroni correction (only for V, C A, D, not Raldh)
pvalsRewPuff = pvalsRewPuff*4;
disp('p-val REW (V,C,A,D) (BONF-4): ')
disp(['    ' num2str(pvalsRewPuff(1,1)) ', ' num2str(pvalsRewPuff(2,1)) ', ' num2str(pvalsRewPuff(4,1)) ', ' num2str(pvalsRewPuff(5,1))])
disp('p-val PUFF (V,C,A,D) (BONF-4): ') 
disp(['    ' num2str(pvalsRewPuff(1,2)) ', ' num2str(pvalsRewPuff(2,2)) ', ' num2str(pvalsRewPuff(4,2)) ', ' num2str(pvalsRewPuff(5,2))])


%% plot reward vs airpuff, as well as reward vs reward at rest
% to compare the amplitude the response to reward vs air puff for each
% recording, and to determine whether the amplitude of the response to
% reward changes when the mouse moves or not.
% Fig. 4H, Fig. S9D and Fig. S8F 

vars = ['acc crossCorr' trigs trigs2 'rewardNoRest'];
id = [1 2 4 5 3]; %which subtype gets plotted where
varIdxAll = [1 2; 7 6]; % rewards vs airpuff and rewards at rest vs not at rest
for v = 1:2
    varIdx = varIdxAll(v,:);
    figure
    set(gcf, 'Position',  [630, 550, 680, 420])
    pvals = nan(5,1);
    for s = 1:length(subtypes)
        subplot(2,3,id(s))
        % scatter points for each recording
        scatter(response(good == s,varIdx(1)),response(good == s,varIdx(2)),250,colors(s,:),'.');
        % lims, lines, labels
        xlim([-0.25 0.25])
        ylim([-0.25 0.25])
        verticalLine(0);
        horizontalLine(0);
        xlabel(vars{varIdx(1)+1})
        ylabel(vars{varIdx(2)+1})
        % identity line (diagonal)
        tAxis = refline(1,0);
        tAxis.Color = 'k';
        % plot center of mass
        com = [mean(response(good == s,varIdx(1)),'omitnan'),mean(response(good == s,varIdx(2)),'omitnan')];
        hold on
        scatter(com(1),com(2),200,'k','x')
        %p-val
        pvals(s) = signrank(response(good == s,varIdx(1)),response(good == s,varIdx(2)));
        pvals(s) = pvals(s) * 4; %bonf correction, for V,C,A,D
        pvals(s) = round(pvals(s),2,'significant');
    end
    subplot(2,3,2)
    if varIdx(2) == 2
        title('Reward vs Air puff response')
        disp('p-values (BONF) reward vs air puff:')
    else
        title('Reward vs Reward at rest')
        disp('p-values (BONF) reward not at rest vs reward at rest:')
    end
    disp(['    V,C,A,D = ' num2str(pvals(1)) ', ' num2str(pvals(2)) ', ' num2str(pvals(4)) ', ' num2str(pvals(5))])
end


%% plot paired small vs large rewards
% Fig. 4I and Fig. S9E
figure
hold on
% plot pair of connected dots per recording, plus average and sem per each
% subtype
for s = 1:length(subtypes)
    plot([s*2-1 s*2],[response(good == s,4) response(good == s,5)],'o-','Color',C.grey{1},'MarkerEdgeColor',colors(s,:))
    errorbar([s*2-1 s*2],[mean(response(good == s,4),'omitnan') mean(response(good == s,5),'omitnan')],[std(response(good == s,4),'omitnan')/sqrt(sum(~isnan(response(good == s,4)))) std(response(good == s,5),'omitnan')/sqrt(sum(~isnan(response(good == s,5))))],'LineStyle','none','Color','k','LineWidth',1.5,'CapSize',15)
    errorbar([s*2-1 s*2],[mean(response(good == s,4),'omitnan') mean(response(good == s,5),'omitnan')],[0 0],'LineStyle','none','Color','k','LineWidth',1.5,'CapSize',30)
end
xlim([0.5 10.5])
ylim([-0.1 0.25])
horizontalLine(0);
xticks(1:10)
xticklabels(repmat({'Small' 'Large'},1,5))
ylabel('Reward response (Int Norm \DeltaF/F)')

% pvalues (paired) - plus bonferroni correction for V,C,A,D
disp('p-vals small vs large rew (BONF-4): (V,C,A,D)')
pval = nan(1,5);
n = nan(1,5);
m = nan(1,5); %mice
for s = 1:5
    respSmall = response(good == s,4);
    respLarge = response(good == s,5);
    p = signrank(respSmall,respLarge);
    pval(s) = p *4; %bonferroni correction
    n(s) = min(sum(~isnan(respLarge)),sum(~isnan(respSmall)));
    goodTemp = find(good == s);
    goodTemp(isnan(respLarge) |  isnan(respSmall)) = [];
    m(s) = length(unique(data6.Mouse(goodTemp)));
end
disp(['    ' num2str(pval(1)) ', ' num2str(pval(2)) ', ' num2str(pval(4)) ', ' num2str(pval(5))])
disp(['    n = ' num2str(n(1)) ', ' num2str(n(2)) ', ' num2str(n(4)) ', ' num2str(n(5))])
disp(['    m = ' num2str(m(1)) ', ' num2str(m(2)) ', ' num2str(m(4)) ', ' num2str(m(5))])


%% plot reward response vs locomotion cross-corr integral
% to show that in Raldh (mixed subtype), locomotion response is linked to
% reward response amplitude
% Fig. S1I
id = [1 2 4 5 3]; %order of subtypes to plot
figure
set(gcf, 'Position',  [630, 550, 680, 420])
for s = 1:5
    subplot(2,3,id(s))
    scatter(response(good == s,1),ccInt(good == s),250,colors(s,:),'.');
    xlim([-0.25 0.25])
    ylim([-0.10 0.10])
    verticalLine(0);
    horizontalLine(0);
    xlabel('Reward response')
    ylabel('Acc cross-corr Int')
    
    % plot center of mass
    com = [mean(response(good == s,1),'omitnan'),mean(ccInt(good == s),'omitnan')];
    hold on
    scatter(com(1),com(2),200,'k','x')
    view([90 -90])
end
subplot(2,3,2)
title('Reward vs Locomotion Cross-corr')


%% cross-correlation between licking and acceleration and other variables
% only consider variables during running
%Fig. S7H
colorsLight = [C.lightRed{1};C.lightYellow{1};C.lightGreen{1}; C.lightAqua{1}; C.lightGrey{1}];
vars = {'chMov' 'Acceleration' 'chGreen' 'Licking'};
plotVarId = [4 1; 4 2; 4 3; 1 3]; %which comparisons to do in each plot
yl = [-0.5 0.5; -0.12 0.12; -0.2 0.2];
ylMovLick = [-0.5 0.5; -0.3 0.3; -0.2 0.2]; %ylims for each plot
lickAccCC = nan(size(data6,1),win*2+1,size(plotVarId,1)); %save cross-corr

figure
for s = [1 2 4] %only pure subtypes
    % get cross-correlations between variable pairs
    idx = find(good == s);
    for recIdx = 1:length(idx)
        rec = idx(recIdx);
        if data6.exRun(rec) == 1
            % get data
            mov = logical(data6.data{rec}.MovOnOff{1});
            varData = nan(length(vars),length(mov));
            for t = 1:length(vars)
                varData(t,:) = data6.data{rec}.(vars{t}){1};
            end
            varData(:,~mov) = []; % only use data during movement
            varData(:,any(isnan(varData),2)) = []; %remove nans
            % calculate cross-corr between pairs in variables as indicated
            % by plotVarId
            for t = 1:size(plotVarId,1)
                lickAccCC(rec,:,t) = crosscorr(varData(plotVarId(t,1),:),varData(plotVarId(t,2),:),win);
            end
        end        
    end 
    % plot
    for t = 1:size(plotVarId,1)
        cc = lickAccCC(good == s,:,t);
        cc(any(isnan(cc),2),:) = [];
        mem = mean(cc,1,'omitnan'); %nean
        sem = std(cc,[],1,'omitnan')/sqrt(size(cc,1));
        if s == 4
            ss = 3; %subplot number for Raldh is 3rd
        else
            ss = s;
        end
        subplotMN(3,size(plotVarId,1),ss,t)
        ran = (-win:win)/hz;
        fill([ran fliplr(ran)], [mem+sem fliplr(mem-sem)],colorsLight(s,:),'EdgeColor',colors(s,:),'LineStyle','-')
        % add lims etc
        if t ~= 4
            ylim(yl(t,:))
        else
            ylim(ylMovLick(ss,:))
        end
        if s == 1
            title([vars{plotVarId(t,1)} ' vs ' vars{plotVarId(t,2)}])
        end
        verticalLine(0);
        xlabel('Time (s)')
        ylabel('cross-corr')
    end
end 


%% compare Anxa and Calb in same region
% select an area of striatum where Anxa1 and Calb axons overlap and compare
% triggered averages for only recordings made within this area.
% This requires recordings to have XYZ coordinates, added using the 
% FiberLocations2 script.
% Select two areas in two most anterior brain slices (bregma +0.86 and +0.5)
if snc == 0  
    % load brain image
    load('brain_Empty.mat', 'brain') 
    scale = 35.6667; % see FiberLocations2, convert XYZ from pixesl to mm
    secEdge = 204.4; % in pixels, this indicates the start of each slice in the image
    % set area to select recordings
    radius = 0.5; % mm
    radiusPix = radius*scale;
    center1 = [91,101]; %center of circle in first brain slice, bregma +0.86
    center2 = [306,101]; %center of circle in first brain slice, bregma +0.5
    % If combine = 1, combine both slices into one. If = 0, keep
    % coordinates separate
    combine = 1;
    if combine == 1
        brain = [brain(:,round(secEdge+1):round(secEdge*2)) brain(:,round((secEdge+1)*5):round(secEdge*6))];
    end
    % get all recordings within these circles
    loc = cell2mat(data6.RecLocG);
    locYes1 = (loc(:,1)-center1(1)).^2 + (loc(:,2)-center1(2)).^2 <= radiusPix^2;
    locYes2 = (loc(:,1)-center2(1)).^2 + (loc(:,2)-center2(2)).^2 <= radiusPix^2;
    locYes = locYes1 | locYes2; %recordings within circle in either slice one or two
    % plot triggered averages only for Calb and Anxa recordings within circle
    figure
    set(gcf, 'Position',  [50 551 1340 420])
    pureSubs = [2 4]; % C and A
    for s = 1:2
        idx = good == pureSubs(s) & locYes; % subtype and lcoation
        % plot reward and air puff triggered average and heatmaps
        plotAverages_paper(data6, '1', 'Reward', 'G', idx, 1, AV, neg, pos, 1:sum(idx),[2 6 6*s-3],1);
	    subplotMN(2,6,s,3)
        ylim(ylimsR(1,:))
        yyaxis left
        ylim([-2 2])
        plotAverages_paper(data6, '1', 'AirPuff', 'G', idx, 1, AV, neg, pos, 1:sum(idx),[2 6 6*s-1],1);
        subplotMN(2,6,s,5)
        ylim(ylimsR(2,:))
        yyaxis left
        ylim([-2 2])
        % plot brain slices and selected circle
        subplotMN(2,6,s,[1,2])
        imshow(brain)
        xlim([1,secEdge*2])
        hold on
        if combine == 1 %combine both slices on one
            [xx,yy] = drawCircle(center1(1),center1(2),radiusPix);
            xx = xx + 11;
            plot(xx,yy,'r:','LineWidth',1.5)
            [xx,yy] = drawCircle(center2(1),center2(2),radiusPix);
            xx = xx - 204;
            plot(xx,yy,'r:','LineWidth',1.5)
        else
            [xx,yy] = drawCircle(center1(1),center1(2),radiusPix);
            plot(xx,yy,'r:','LineWidth',1.5)
            [xx,yy] = drawCircle(center2(1),center2(2),radiusPix);
            plot(xx,yy,'r:','LineWidth',1.5)
        end
        % plot recording locations
        idx = find(idx);
        for j = 1:length(idx)
            if combine == 1
                RecLoc = data6.RecLocGshift{idx(j)};
            else
                RecLoc = data6.RecLocG{idx(j)};
            end
            scatter(RecLoc(1),RecLoc(2),[],[0.8 0.8 0.8],'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceAlpha',0.5)
        end
    end
end          


%% 3D plot of locomotion, reward and airpuff signaling for all subtypes
% we cannot plot locomotion PC1, PO1, reward and air puff, as it would 
% require a 4D plot. Instead, we simplify the locomotion signaling to just 
% the angle between PC1 and PC2.
% Fig. 6A,B and Fig. S9F,G (snc)

% calculate locomotion PC1/PC2 angle
angleSplit = 96; %this angle was chosen because it best separated the subtypes - it was the sparsest angle
score2Norm = scores./std(scores,'omitnan');
angles = atan2d(score2Norm(:,2),score2Norm(:,1));
angles = mod(angles, 360);
angles(angles<angleSplit) = angles(angles<angleSplit)+360;
% plot
varLims = [angleSplit angleSplit+360; -0.10 0.25; -0.10 0.25];
figure
subplot(2,2,1)
for s = [1 2 4] %only V,C,A
    scatter3(angles(good == s), response(good == s,1), response(good == s,2), [], colors(s,:),'filled')
    hold on
    % ploy center of mass
    scatter3(mean(angles(good == s),'omitnan'), mean(response(good == s,1),'omitnan'), mean(response(good == s,2),'omitnan'), 200, [0 0 0],'x')
end
% add labels and lims etc
xlabel('locomotion PC1/PC2 angle'); xlim([angleSplit angleSplit+360])
ylabel('reward response')
zlabel('air puff response')
axis square
view([1 -1 1])
xlim(varLims(1,:))
ylim(varLims(2,:))
zlim(varLims(3,:))
% plot 3 additional plots showing each pair of the 3 variables 
varTit = {'locomotion PC1/PC2 angle', 'reward response', 'air puff response'};
varI = [1 2; 1 3; 2 3];
for s = [1 2 4]
    vars = {angles(good == s), response(good == s,1), response(good == s,2)};
    for i = 1:3
        subplot(2,2,i+1)
        scatter(vars{varI(i,1)}, vars{varI(i,2)}, [], colors(s,:),'filled')
        hold on
        xlim(varLims(varI(i,1),:))
        ylim(varLims(varI(i,2),:))
        xlabel(varTit{varI(i,1)})
        ylabel(varTit{varI(i,2)})
        verticalLine(0);
        horizontalLine(0);
    end
end


%% PCA on all variables to reduce dimensionality
% To visualize all 4 variables, locomotion PC1, PO1, reward and air puff, 
% we conduct PCA on these 4 variables and plot the best two PCs. This will
% not be used in the paper, but it is useful to understand the outcome of
% K-means classification below.

% get locomotion scores and reward/airpuff responses
allPeaksP = [scores(:,1:2) response(:,1:2)];
goodP = good(~any(isnan(allPeaksP),2));
allPeaksP = allPeaksP(~any(isnan(allPeaksP),2),:);
% only pure subtypes (V,C,A)
allPeaksP(goodP == 3 | goodP == 5,:) = [];
goodP(goodP == 3 | goodP == 5,:) = [];
% normalize data - this is important as otherwise variables with greater
% variance will be overrepresented
allPeaksPN = normalize(allPeaksP,1);
% run PCA
[~,scoreP,~,~,explainedP,] = pca(allPeaksPN);
% plot 
figure;
for s = [1 2 4]
    scatter(scoreP(goodP == s,2), scoreP(goodP == s,1), [], colors(s,:),'filled')
    hold on
end    
axis square
set(gca, 'YDir', 'reverse'); % to match best 3D representation in previous section
title(['PC1 = ' num2str(round(explainedP(1),1)) '%, PC2 = ' num2str(round(explainedP(2),1)) '% - total = ' num2str(round(sum(explainedP(1:2)),1)) '% of variance'])
xlabel('PC2')
ylabel('PC1')


%% K-means classification
% run K-means on locomotion PC1, PO1, reward and air puff for 3 pure
% subtypes (V,C,A) with 3 clusters, matching subtype number. 
% Then match each resulting cluster with one existing subtype and determine
% how similar the unsupervised clustering is to the real subtypes.
% Fig. 6C
clusterNum = 3;
if snc == 0
    % run PCA on striatal data
    [Kidx,Kcent] = kmeans(allPeaksPN,clusterNum);
    save([dataProcessingFolder '\k-meansCent.mat'], 'Kcent'); % save centroids for use with SNc data
else
    % for SNc data, load striatal k-mean centroids and determine which
    % cluster each recording corresponds to by calculating distance.
    load([dataProcessingFolder '\k-meansCent.mat'], 'Kcent');
    distKmeans = nan(size(allPeaksPN,1),size(Kcent,1));
    for p = 1:size(Kcent,1) %number of k-means centroids
        % calculate ecludian distance for each recording to each centroid
        dist = zeros(size(allPeaksPN,1),1);
        for d = 1:size(Kcent,2) %dimensions
            dist = dist + (allPeaksPN(:,d)-Kcent(p,d)).^2;
        end
        dist = sqrt(dist);
        distKmeans(:,p) = dist;
    end
    % get cluster identity by determining the closest centroid for each
    % recording (min distance)
    [~,Kidx] = min(distKmeans,[],2);
end

% match clusters to subtypes
match = [1 nan; 2 nan; 4 nan];  
Kidx2 = Kidx;
for i = 1:3
    match(i,2) = median(Kidx(goodP == match(i,1)));
end
for i = 1:3
    Kidx2(Kidx == match(i,2)) = match(i,1);
end
% get % correct (how many of the recordings of each subtype were classified
% in a single cluster
correct = zeros(5,1);
for i = [1 2 4]
    correct(i) = sum(Kidx2(goodP == i) == i)/sum(goodP == i);
end
correct(3) = [];
correct(4) = sum(Kidx2 == goodP)/length(goodP);
% bar and pie graphs of % correct classified
colorsK = colors([1 2 4 5],:);
f = figure;
f.Position = [1200 550 330 420];
tiledlayout(4,3,'TileSpacing','compact');
nexttile([4,2])
hold on
for i = 1:4
    bar(i,correct(i)*100,'FaceColor',colorsK(i,:),'FaceAlpha',0.7)
end
xticks(1:4)
xticklabels([subtypes([1:2 4]) 'Total'])
ylabel('% correctly classified')
if snc == 0
    title('k-means clustering')
else
    title('k-means clustering - str clusters')
end
horizontalLine(100/3); 
for i = 1:4
    nexttile
    p = pie([correct(i) 1-correct(i)]);
    p(1).FaceColor = colorsK(i,:);
    p(1).FaceAlpha = 0.7;
    p(3).FaceColor = [0.9 0.9 0.9];
    p(3).FaceAlpha = 0.7;
end
      

%% Triggered averages per mouse
% Fig. S8H,I

% identify recordings from same mouse and average
mice = data6.Mouse(good > 0); 
mice = strip(mice,'right','L'); % one mouse had recording from left hemisphere, but still same mouse
type = good(good > 0);
[mice,ii] = unique(mice);   % list of all good mice names
type = type(ii);            % mouse exp type
mouseAve = cell(5,2);
figure
for s = 1:5
    miceX = mice(type == s); % good, unique mice of subtype s
    miceI = data6.Mouse;     % all mice, in order like data6
    miceI = strip(miceI,'right','L'); 
    [~,miceIdx] = ismember(miceI,miceX); % idx of mouse for data6
    miceIdx(good == 0) = 0;              % keep only good

    %for each mouse, average triggered averages from all its recordings
    for i = 1:2 %rew and puff
        mouseAve{s,i} = nan(max(miceIdx),size(allData{1},2));
        for m = 1:max(miceIdx)
            mouseAve{s,i}(m,:) = mean(allData{i+1}(miceIdx == m,:),1,'omitnan');
        end
        mouseAve{s,i}(any(isnan(mouseAve{s,i}),2),:) = [];
    end
    % plot mouse averages and heatmaps
    colorsSet = [colorsLight(s,:); colors(s,:)];
    ran = (-win:win)/hz;
    for j = 1:2
        mem = mean(mouseAve{s,j},1);
        sem = std(mouseAve{s,j},1)/sqrt(size(mouseAve{s,j},1));
        subplotMN(5,4,s,2*j-1) %trig ave
        fill([ran fliplr(ran)], [mem+sem fliplr(mem-sem)],colorsSet(1,:),'EdgeColor',colorsSet(2,:))
        ylim([-0.02 0.65])
        verticalLine(0);
        subplotMN(5,4,s,2*j) %heatmap
        imagesc(mouseAve{s,j})
        caxis([0 0.6])
        title(['n = ' num2str(size(mouseAve{s,j},1))])
    end
end
subplot(5,4,1)
title('AVE REW PER MOUSE')
subplot(5,4,3)
title('AVE PUFF PER MOUSE')


%% get raw examples for rew and puff
% Fig. 4C
if snc == 0
    colorsRaw = [C.red{1}; C.orange{1}; C.aqua{1}]; %V,C,A
    bestMice = {'VGlut-B612-20200309-0003' 'Calb-2290-20220321-0004' 'Anxa-J075-20220310-0001'};
    xlims = [403 394 248;339 381 228]; %what part of each recording to show
    plotWin = 15; % in sec
    %set figure dimensions
    figure
    set(gcf, 'Position',  [406 291 660 678])
    % find idx for bestMice 
    bestIdx = nan(3,1);
    for s = 1:3
        id = findMouse(data6,bestMice{s});
        bestIdx(s) =  id(1);
    end
    % Plot - see function at bottom
    for s = 1:3
        rec = bestIdx(s);
        for i = 1:2
            plotRawData(data6,rec,plotWin,xlims(i,s),{3,2,s,i},colorsRaw(s,:))
        end
    end
end

%% get raw examples (with rew and run) for supplemental figure
% Fig. S1D
if snc == 0
    colorsRaw = [C.red{1}; C.orange{1}; C.green{1}; C.green{1}];
    bestMice = {'VGlut-B612-20200309-0003' 'Calb-2290-20220321-0004' 'Raldh-1573-20201015-0002' 'Raldh-1846-20210623-0002'};
    xlims = [403 393 111 152];
    plotWin = 30; % in sec
    %set figure dimensions
    figure
    set(gcf, 'Position',  [406 291 660 678])
    % find idx for bestMice 
    bestIdx = nan(4,1);
    for s = 1:4
        id = findMouse(data6,bestMice{s});
        bestIdx(s) =  id(1);
    end
    % Plot raw traces - see function at bottom
    for s = 1:4
        rec = bestIdx(s);
        plotRawData(data6,rec,plotWin,xlims(s),{4,5,s,1:4},colorsRaw(s,:))
    end
    
    % plot cross-corr DF/F and acc for each example
    for s = 1:4
        rec = bestIdx(s);
        % get cross-corr
        g = data6.data{rec}.chGreen{1};
        g405 = data6.data{rec}.chGreen405{1};
        acc = data6.data{rec}.Acceleration{1};
        run = data6.data{rec}.Run{1};
        bad = isnan(g) | isnan(g405) | isnan(acc);
        g(bad | ~run') = [];
        g405(bad | ~run') = [];
        acc(bad | ~run') = [];
        cc = crosscorr(acc,g,win); %470
        cc405 = crosscorr(acc,g405,win); %405
        % plot
        subplotMN(4,5,s,5);
        plot(-win/hz:1/hz:win/hz,cc405,'Color',C.lightBlue{1})
        hold on
        plot(-win/hz:1/hz:win/hz,cc,'Color',colorsRaw(s,:))
        ylim([-0.3 0.3])
        xlabel('Lag (s)')
        ylabel('cross-corr')
    end
end


%% save good and rew/puff responses for future access
temp = data6(:,[1:3 5:6 10:11]);
RewardResp = response(:,1);
AirPuffResp = response(:,2);
recID = cell(size(data6,1),1);
for rec = 1:size(data6,1)
    recID{rec} = data6.data{rec}.Properties.RowNames{1};
end
goodRew = [table(recID) temp table(RewardResp) table(AirPuffResp)];

if snc == 0
    exportRewStr = export;
    save([dataProcessingFolder '\good_Rew.mat'], 'goodRew');
    save([dataProcessingFolder '\exportRewStr.mat'], 'exportRewStr');
else
    exportRewSNc = export;
    save([dataProcessingFolder '\goodSNc_Rew.mat'], 'goodRew');
    save([dataProcessingFolder '\exportRewSNc.mat'], 'exportRewSNc');
end

end



%% FUNCTION TO PLOT RAW DATA
function plotRawData(data6,rec,plotWin,xl,sp,colorRaw)
load colors.mat C
% On the plot, we will need scales for each of the variables, and to
% plot each variable normalized to a similar amplitude. To do so,
% normalize each variable by its max and min, and normalize a standard
% scale bar (scales0) in the same way.
vars = {'chMov','Acceleration', 'chGreen405','chGreen','Reward','AirPuff','Licking'};
scales0 = [0.3 1 0.2 0.2]; % vel, acc, DF/F twice
% get variables
vals = cell(length(vars),2);
for v = 1:length(vars)
    %de normalize DF/F
    temp = data6.data{rec}.(vars{v}){1};
    if any(strcmp(vars{v},{'chGreen405','chGreen'}))
        temp = temp * data6.norm{rec}(1);
    end
    % replace nans with closest values
    temp(isnan(temp)) = temp(circshift(isnan(temp),1));
    % smooth vel, acc, and DF/F
    if v <= 4
        temp = fastsmooth(temp,10);
    end
    vals{v} = temp;
end
% normalize variables and scales for plotting
for j = 1:4 % add scales to vals table
    vals{j,2} = [0 scales0(j)];
end
for j = 2:-1:1 %normalize each scale then each variable by max/min
    vals{1,j} = (vals{1,j}-min(vals{1,1}))./(max(vals{1,1})-min(vals{1,1}))*1.2-0.1;
    vals{2,j} = (vals{2,j})./(max(abs(vals{2,1})))*1.2+0.5; %acc normalize by abs max value to keep zero
    vals{3,j} = (vals{3,j}-min(vals{4,1}))./(max(vals{4,1})-min(vals{4,1}))*1.2-0.1;
    vals{4,j} = (vals{4,j}-min(vals{4,1}))./(max(vals{4,1})-min(vals{4,1}))*1.2-0.1; %use chGreen to norm chGreen 405 too
end
% get yAxis ticks from scales
% we will plot each variable on top of each other, adding +1 for
% each variable - so add these on too
scaleTicks = cell2mat(vals(1:3,2));
scaleTicks = scaleTicks + [0 0; 1 1; 2 2]; % add shift for each variable
scaleTicks2 = reshape(scaleTicks',1,6); %get scale positions
scaleNums = fliplr(upsample(fliplr(scales0(1:3)),2)); %get scale values
% plot
vals2 = vals;
subplotMN(sp{1},sp{2},sp{3},sp{4});
hold on
% plot reward, air puff and licking to span whole ploy
plot(x(vals2{5},100),vals2{5}*10-5,'Color',C.blue{1}); %reward
plot(x(vals2{6},100),vals2{6}*10-5,'Color',C.red{1}); %air puff
plot(x(vals2{7},100),vals2{7}*10-5,'Color',C.lightBlue{1}); %licking
% plot vel at bottom, accel on top of that (+1) and DF/F (405
% and 470) on top of all (+2)
plot(x(vals2{1},100),vals2{1},'Color',[0 0 0]); % vek
plot(x(vals2{2},100),vals2{2}+1,'Color',C.grey{1}); %acc
horizontalLine(1.5);
plot(x(vals2{3},100),vals2{3}+2,'Color',C.lightBlue{1}); %405
plot(x(vals2{4},100),vals2{4}+2,'Color',colorRaw); %470
% set title, scale, and lims
title(data6.data{rec}.Properties.RowNames{1})
yticks(scaleTicks2)
yticklabels(scaleNums)
ylim([min(vals2{1}) max(vals2{4}+2)])
xlim([xl xl+plotWin])
xticks(xl:10:xl+plotWin)
xticklabels(0:10:plotWin)
xlabel('Time (s)')
end
    



