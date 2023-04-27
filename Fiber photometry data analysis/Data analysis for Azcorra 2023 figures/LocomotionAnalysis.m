function LocomotionAnalysis(data6,snc)
% This script plots many different figures for the analysis of LOCOMOTION
% SIGNALING of different DA subtypes as presented in Azcorra et al. 2023.
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
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
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


%% get exclusion criteria for locomotion analysis
% get either all snc recordings or all str recordings
if snc == 0
    strSNc = data6.exStr(:,1) == 1; % only str
else
    strSNc = data6.exStr(:,1) == 0; % only snc
end
% include recordings in the correct location, with good sig2noise, with 
% enough running, without high control cross-correlation 405-acceleration,
% and without too many movement artifacts in 405. These were calculated in
% the script exclusionCriteria_allFigs.m
good = double(strSNc & data6.exSig2Noise(:,1) & data6.ex405Acc(:,1) & data6.exRun & data6.Bad405All(:,1) == 0);

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
    exportMoveStr = export;
    save([dataProcessingFolder '\exportMoveStr.mat'], 'exportMoveStr');
else
    exportMoveSNc = export;
    save([dataProcessingFolder '\exportMoveSNc.mat'], 'exportMoveSNc');
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
        data6.data{rec}.(skip){1}(data6.data{rec}.Bad405G{1}) = [];
    end
end

%% asign numbers for different subtypes
% this groups them for plotting.
% also set identifying colors for each
subtypes = {'VGlut' 'Calb' 'Raldh' 'Anxa' 'Dat'};
load colors.mat C
colors = [C.red{1}; C.orange{1}; C.green{1}; C.aqua{1}; C.grey{1}];
s = length(subtypes);
good(good == 1) = 6; % any subtypes different from above will be 6
for i = 1:length(subtypes)
    sIdx = strcmp(data6.Exp,subtypes{i});
    good(good > 0 & sIdx) = i;
end

%% calculate crosscorr DF/F vs acc for PCA analysis
% it will be plotted again further down, sorted by PCA
win = 100;
ccData470 = nan(size(data6,1),win*2+1); % save cross-corr (cc) here
ccData405 = ccData470; % save control cross-corr here
figure
for i = 1:s
    cc = crossCorr2(data6, good == i, 'G', 'A', win,[1,2,1],1:sum(good == i),normHeat); %separate script
    cc470 = cell2mat(cc(:,1)); %
    cc405 = cell2mat(cc(:,2));
    ccData470(good == i,:) = cc470(good == i,:);
    ccData405(good == i,:) = cc405(good == i,:);
end
close() %this figure is unnecessary, we just need the data for PCA analysis

%% PCA analysis
% Run pca on cross-correlations between DF/F and acceleration to reduce the
% dimensionality of the cross-correlation shapes, which will allow
% comparison between subtypes. 
% Run only on data from VGlut, Calb, and Anxa!!! Not on mixed subtypes
% (however, running pca on all subtypes still separates subtypes, it just
% skewes the center of the PC coordinates)
% Do not center data ('Centered' = false) - this is to maintain zero in the
% cross-correlation, which has meaning
ccDataPCA = ccData470;
ccDataPCA(good == 3,:) = nan; %remove Raldh
ccDataPCA(good == 5,:) = nan; %remove Dat
ccPureSubs = ccDataPCA(any(~isnan(ccDataPCA),2),:); % remove Raldh and dat

% get PCs of all data
% for str save PCs, for SNc use str PCs
if snc == 0
    [coeff,~,~,~,explained,mu] = pca(ccPureSubs,'Centered',false,'Algorithm','eig'); %latent is how much variance each PC explains
    disp(['   First PC explains ' num2str(round(explained(1),1)) '% of variance'])
    disp(['   Second PC explains ' num2str(round(explained(2),1)) '% of variance'])
    save([dataProcessingFolder '\PCA_striatum.mat'], 'coeff','explained','mu'); % mu will be zero if centered is false
else
    load([dataProcessingFolder '\PCA_striatum.mat'], 'coeff', 'mu');
    % calculate the percent of the SNc variance explained by the striatal 
    % PCs but WITHOUT centering - method calculated by John and corroborated
    varianceMat = coeff'*(ccPureSubs'*ccPureSubs)*coeff;  
    latent = diag(varianceMat)/(size(ccPureSubs,1)-1); 
    explained = latent./sum(latent)*100;
    disp(['   First PC explains ' num2str(round(explained(1),1)) '% of variance'])
    disp(['   Second PC explains ' num2str(round(explained(2),1)) '% of variance'])
end

% Get PC scores for the orignal data (with lots of empty rows and DAT/Raldh)
% to keep the order of the recordings (and match with 'good')
scores = (ccData470-mu)*coeff;

%% Plot scores for each recording along PC1 and PC2 by subtype 
% Fig. 2K (Str) and Fig. S9K (SNc)

% get max scores for plotting (will use for x/y lims)
limsPC = [min(scores,[],1)-0.1;max(scores,[],1)+0.1]';
% plot scores for first 2 PCs
iPC = [1 2]; % which PCs to plot
subOrder = [1 2 4 5 3];
figure
set(gcf, 'Position',  [100, 550-500*snc, 680, 420])
for s = 1:length(subtypes)
    subplot(2,3,subOrder(s))
    scatter(scores(good == s,iPC(1)),scores(good == s,iPC(2)),150,colors(s,:),'.');
    title(subtypes{s})
    xlabel(['PC' num2str(iPC(1))])
    ylabel(['PC' num2str(iPC(2))])
    
    xlim([-max(abs(limsPC(iPC(1),:))) max(abs(limsPC(iPC(1),:)))]); %all plots use same lims = largest of all
    ylim([-max(abs(limsPC(iPC(2),:))) max(abs(limsPC(iPC(2),:)))]);
    horizontalLine(0);
    verticalLine(0);
    
    % plot center of mass (COM)
    COM = [mean(scores(good == s,iPC(1))),mean(scores(good == i,iPC(2)))];
    hold on
    scatter(COM(1),COM(2),200,'k','x')
end

%% plot loadings of first 2 PCs and their combinations per quadrants 
% only for str (snc will use striatal ones)
% Fig. 2J
if snc == 0
    xAxis = -1:0.01:1;
    plotNum = [1:4 6:9];
    pcMult = [-1 1; 0 1; 1 1; -1 0; 1 0; -1 -1; 0 -1; 1 -1];
    titles = {'Quad2' ['+PC' num2str(iPC(2))] 'Quad1' ['-PC' num2str(iPC(1))] ['+PC' num2str(iPC(1))] 'Quad3' ['-PC' num2str(iPC(2))] 'Quad4';};
    pc1 = coeff(:,iPC(1));
    pc2 = 0.7*coeff(:,iPC(2));
    
    figure
    set(gcf, 'Position',  [780, 550, 420, 420])
    for i = 1:length(plotNum)
        subplot(3,3,plotNum(i))
        plot(xAxis, pc1*pcMult(i,1) + pc2*pcMult(i,2),'k')
        ylim([-0.3 0.3])
        verticalLine(0);
        horizontalLine(0);
        title(titles{i})
        if i >= 7
            xlabel('Time (s)')
        end
    end
end

%% Calculate PC1/PC2 angle
% This will be used to sort recordings for heatmaps in triggered ave etc.,
% and to plot radial histograms below
% First normalize each PC scores - otherwise the angle will be skewed 
% towards the PC with greater scores

scoresNorm = scores./std(scores,'omitnan'); %normalize scores
angles = atan2d(scoresNorm(:,2),scoresNorm(:,1)); %calculate angles - 0 is right, positive angles go counter clockwise
angles = mod(angles, 360); % convert angles into all positives with 0 in right and increasing values counter clockwise

% For sorting recordings by a circular score, need to split the circle
% somewhere. For each subtype, split by the middle of the quadrant opposite
% to where it's center of mass is located
idxPC = nan(size(data6,1),1); %sorting order
angleSplit = [315 45 135 135 45]; % split angle for each subtype
for s = 1:length(subtypes)
    subAngles = angles(good == s);
    subAngles(subAngles < angleSplit(s)) = subAngles(subAngles < angleSplit(s))+360; %split
    [~,subIdx] = sort(subAngles);
    idxPC(good == s) = subIdx; %sorted recordings within subtypes
end


%% Made radial histograms
% use PC1/2 angles as calculated above
% Fig. 2L
figure;
set(gcf, 'Position',  [100, 50, 680, 420])
subOrder = [1 2 4 5 3];
edges = 0:pi/10:2*pi;
anglesBySub = cell(5,1);
samples = nan(5,1);
for s = 1:length(subtypes)
    subAngles = angles(good == s);
    anglesBySub{s} = subAngles;
    anglesRad = deg2rad(subAngles);
    samples(s) = length(anglesRad);
    subplot(2,3,subOrder(s))
    polarhistogram(anglesRad,edges,'FaceColor',colors(s,:)) %,'Normalization','probability')
end

% show mean angles
subplot(2,3,length(subtypes))
title(['Mean angles (V,C,A) = ' num2str(mean(anglesBySub{1})) ', ' num2str(mean(anglesBySub{2})) ', ' num2str(mean(anglesBySub{4}))]);

% statistical test for whether subtypes have different angles
% this is not ideal because it would need a circular statistic. But because
% almost no subtypes (or mixtures) have recordings at 45 degrees, we can
% split there and do an aproximation
angleSplit = 45; % make this 0
for s = 1:length(subtypes)
    anglesBySub{s}(anglesBySub{s} < angleSplit) = anglesBySub{s}(anglesBySub{s} < angleSplit) + 360;
    anglesBySub{s} = anglesBySub{s} - angleSplit;
end
% Mann-Whitney U test (non-parametric unpaired) with Bonf correction
% compare only pure subtypes
subsPure = [1 2 4];
pvals = nan(3,3);
for i = 1:length(subsPure)
    for j = 1:length(subsPure)
        if i < j
            pvals(i,j) = ranksum(anglesBySub{subsPure(i)},anglesBySub{subsPure(j)});
            pvals(i,j) = pvals(i,j)*3; %Bonferroni correction (3 comparisons)
            pvals(i,j) = round(pvals(i,j), 2, 'significant'); % round
        end
    end
end
subplot(2,3,2)
title(['BONF p-values loc angles (VC,VA,CA): ' num2str(pvals(1,2)) ', ' num2str(pvals(1,3)) ', ' num2str(pvals(2,3))])


%% for Aldh, plot PC1 and PC2 color coded by depth
% Fig. S6E
if snc == 0 % only for striatum
    iPC = [1 2];
    % get aldh data
    s = 3;
    pointsX = scores(good == s,iPC(1));
    pointsY = scores(good == s,iPC(2));
    depths = cell2mat(data6.depthG(good == s));
    
    %get color scale
    cL = [1 1 0.01]; % light = yellow
    cD = [0 0 0]; % dark = black
    cM = C.green{1}; %middle color = green
    nTotal = 100; % total number of steps in color scale
    nL = 10; % steps in scale that are all light color
    nLM = 20; % steps in gradient between light and medium color
    nDM = nTotal-nL-nLM; % steps in gradient beteween dark and medium color
    gLM = [linspace(cM(1),cL(1),nLM); linspace(cM(2),cL(2),nLM); linspace(cM(3),cL(3),nLM)]'; %light to medium gradient
    gDM = [linspace(cM(1),cD(1),nDM); linspace(cM(2),cD(2),nDM); linspace(cM(3),cD(3),nDM)]'; %dark to medium gradient
    g = [repmat(cL,nL,1);flipud(gLM);gDM]; % total color scale: light solid + light-to-medium + medium-to-dark
    
    % get depth scale
    gDepth = linspace(1.5,4,nTotal);
    % plot each recording in  PC1/PC2 space but color coded by depth
    figure
    set(gcf, 'Position',  [1000, 700 300, 220])
    hold on
    for p = 1:length(pointsX)
        depthIdx = interp1(gDepth,1:length(gDepth),depths(p),'nearest'); % get the depth index based on depth scale
        cc = g(depthIdx,:); %get color based on depth idx
        scatter(pointsX(p),pointsY(p),300,cc,'.'); %plot
    end
    title(subtypes{s})
    xlabel(['PC' num2str(iPC(1))])
    ylabel(['PC' num2str(iPC(2))])
    
    xlim(limsPC(iPC(1),:));
    ylim(limsPC(iPC(2),:));
    horizontalLine(0);
    verticalLine(0);
    
    colormap(g);
    h = colorbar;
    set( h, 'YDir', 'reverse' );
    ylabel(h,'Depth (mm)');
    caxis([1.5 4])
end

%% plot all triggered averages and cross-corr for all subtypes, with heatmaps sorted by PC1/PC2 angle
% Fig. 2F,G,H and Fig. S6B,C,D (str), as well as Fig. S9H,I,J (Snc)
subNum = length(subtypes); % number of subtypes
% set parameters for cross-corr
win = 100; 
normHeat = 2;
% set parameters for triggered averages
AV = 'A'; %show acceleration - change to 'V' to see velocity averages
trigs = {'AccOn', 'DecOn', 'peaksGRun'};
neg0 = 100; % bins to show in trig ave BEFORE trigger point (100 bin = 1s)
pos0 = 100; % bins to show in trig ave AFTER trigger point (100 bin = 1s)
neg = neg0 + 5; % these are necessary because edge effect of fill function. Plot 5 extra bins but hide
pos = pos0 + 5;
% save cross-corr and triggered averages for later use
allData = cell(1,4);

% Cross-correlation DF/F vs acceleration - same as used for PCA but now
% plot heatmap sorted by PC1/PC2 angle
figure
set(gcf, 'Position',  [1 42 1590 954])
allData{1} = nan(size(data6,1),win*2+1); %preallocate
for s = 1:subNum
    cc = crossCorr2(data6, good == s, 'G', 'A', win,[subNum,8,s*8-7],idxPC(good == s),normHeat);
    cc470 = cell2mat(cc(:,1)); %keep only cross-corr between 470 and acc
    allData{1}(good == s,:) = cc470(good == s,:);
end
for s = 2:8:subNum*8
    subplot(subNum,8,s)
    caxis([-1,1])
    subplot(subNum,8,s-1)
    title('Cross corr with Acc')
end

% Trig ave on acc/dec/trans
for t = 2:4
    allData{t} = nan(size(data6,1),neg0+pos0+1); %preallocate 
end
eachEvent = cell(5,3); %save trace for each event that goes into triggered averages per recording
nEvents = cell(5,3); % number of events that go to each triggered avera per recording

for s = 1:subNum
    for t = 1:length(trigs)
        subplotInput = [subNum,length(trigs)*2+2,s*(length(trigs)*2+2)-(length(trigs)*2+1-t*2)]; %[subNum,8,s*8-(7-t*2)]
        [~,~,~,~,trigAve,sigsDF,sigsAV] = plotAverages_paper(data6, '1', trigs{t}, 'G', good == s, 1, AV, neg, pos, idxPC(good == s),subplotInput,1); % plot triggered averages - separate script
        allData{t+1}(good == s,:) = trigAve(:,neg-neg0+1:neg+pos-(pos-pos0-1)); %save only between neg0 and pos0, not neg and pos (the edges are only needed for plotting)
        % save DF/F or acceleration for each event that goes into triggered averages per recording
        if t < 3
            eachEvent{s,t} = sigsDF; % DF/F for acc/dec trig ave
        else
            eachEvent{s,t} = sigsAV; % acceleration for transient trig ave
        end
        % save number of events that go to each triggered avera per recording
        nEvents{s,t} = nan(size(sigsDF,1),1); 
        for e = 1:length(nEvents{s,t})
            nEvents{s,t}(e) = size(sigsDF{e},1);
        end
        nEvents{s,t}(isnan(nEvents{s,t}) | nEvents{s,t} == 0) = [];
    end
end
% set axis limits, titles and legend for the figure
ylimsL = [-2 2; -2 2; -0.5 0.5];
ylimsR = [-0.02 0.2; -0.02 0.2; -0.02 0.35];
for s = 1:subNum
    for t = 1:length(trigs)
        subplotMN(subNum,8,s,t*2+1)
        title(trigs{t})
        yyaxis left
        ylim(ylimsL(t,:))
        yyaxis right
        ylim(ylimsR(t,:))
        xlim([-neg0 pos0]/100)
        verticalLine(0);
        
        subplotMN(subNum,8,s,t*2+2)
        xlim([5.5 pos0+neg0+5.5])
    end
end

% Get average number of events for each triggered average plus S.T.D. (NOT
% SEM as elsehwere in figures!!)
for t = 1:length(trigs)
    nEventsTrig = cell2mat(nEvents(:,t));
    Mean = mean(nEventsTrig);
    Std = std(nEventsTrig);
    disp(['average events ' trigs{t} ' = ' num2str(Mean) ' +- ' num2str(Std) ' std']);
end

%% make a bar showing which recording in the above heatmaps is from what
% mouse by color-coding them for each subtype. 
figure
for s = 1:subNum
    subplot(1,5,s)
    mouseHeatmapIdx = idxPC(good == s); %get order of recordings in heatmap for subtype
    mice = data6.Mouse(good == s); % get mice
    mice = strip(mice,'right','L'); % one mouse was recorded from left hemisphere one day, but same mouse
    [miceUnique,~,recIdx] = unique(mice); % get unique mice and which recordings come from each
    miceIdOrder = recIdx(mouseHeatmapIdx); % get which mice match the heatmap for that subtype
    colorsU = hsv(length(miceUnique))/4*3; %colormap - scale the hsv colormpa for the number of unique mice in subtype
    imagesc(miceIdOrder)
    colormap(colorsU)
    title(subtypes{s});
end


%% Calculate the probability that an acc/dec is followed by a DF/F increase
% in each subtype. Also calculate the probability that a DF/F transient is
% followed by an acceleration/deceleration
% Fig. S7F

% check pairing between Vglut/Calb transients and DECELERATION, and Anxa
% transients and ACCELERATION
testPairs = [1 2; 2 2; 4 1]; % vglut=dec, calb=dec, anxa=acc
winAcc = 75; %window of 0.75 s post acc/dec/trans
probTransAcc = cell(size(eachEvent));

% get the % of acc/dec with a positive DF/F integral after the acc/dec
for i = 1:size(testPairs,1)
    events = eachEvent{testPairs(i,1),testPairs(i,2)}; %get DF/F for each event for the given subtype and acc/dec pair
    probTransAcc{testPairs(i,1),testPairs(i,2)} = nan(size(events,1),1);
    for e = 1:size(events,1)
        dff = events{e}(:,neg:neg+winAcc) - events{e}(:,neg); %DF/F in window after acc/dec
        tt = (1:size(dff ,2))/hz; % time vector for integral calculation below
        dffInt = trapz(tt,dff ,2); %integral of the DF/F trace after acc/dec minus DF/F at t=0 for each recording
        perPos = sum(dffInt > 0)/length(dffInt)*100; % for each recording, % of positive
        probTransAcc{testPairs(i,1),testPairs(i,2)}(e) = perPos;
    end
end
% get the % of trans with a positive/negative accleration integral after the trans
for i = 1:size(testPairs,1)
    events = eachEvent{testPairs(i,1),3}; % get acc for each event for the given subtype transients
    probTransAcc{testPairs(i,1),3} = nan(size(events,1),1);
    for e = 1:size(events,1)
        acc = events{e}(:,50:neg) - events{e}(:,neg); %acc in window after trans
        tt = (1:size(acc,2))/hz; % time vector for integral calculation below
        accInt = trapz(tt,acc,2); %integral of the acc trace after trans minus DF/F at t=0 for each recording
        if testPairs(i,1) == 4 % anxa = % with pos acc
            perPosNeg = sum(accInt > 0)/length(accInt)*100;
        else %vglut and cal = % with neg acc
            perPosNeg = sum(accInt < 0)/length(accInt)*100;
        end
        probTransAcc{testPairs(i,1),3}(e) = perPosNeg;
    end
end
% plot histograms
figure
for i = 1:size(testPairs,1)
    subplotMN(2,size(testPairs,1),1,i) %trans 
    histogram(probTransAcc{testPairs(i,1),testPairs(i,2)},0:5:100,'FaceColor',colors(testPairs(i,1),:))
    m = mean(probTransAcc{testPairs(i,1),testPairs(i,2)}); %get subtype mean
    title(['mean = ' num2str(m) ' %'])
    
    subplotMN(2,size(testPairs,1),2,i)
    histogram(probTransAcc{testPairs(i,1),3},0:5:100,'FaceColor',colors(testPairs(i,1),:))
    m = mean(probTransAcc{testPairs(i,1),3});
    title(['mean = ' num2str(m) ' %'])
end
for i = 1:size(testPairs,1)*2
    subplot(2,size(testPairs,1),i)
    ylabel('# of recordings')
    yl = ylim;
    ylim([0 yl(2)+1])
    verticalLine(m);
    xlim([0 100])
end
subplot(2,size(testPairs,1),2)
xlabel('% acc/dec followed by DF/F increase')
subplot(2,size(testPairs,1),5)
xlabel('% trans followed by acc/dec')

%% Trig ave on MovOn and MovOff
% Fig. S7G
% use same AV, neg/pos and neg0/pos0 as above trig ave on acc/dec/trans
trigs = {'MovOn' 'MovOff'};
nEvents = cell(subNum,length(trigs));

figure
set(gcf, 'Position',  [412    42   718   954])
for s = 1:subNum
    for t = 1:2
        [~,~,~,~,~,sig] = plotAverages_paper(data6, '1', trigs{t}, 'G', good == s, 1, AV, neg, pos, idxPC(good == s),[subNum,4,s*4-(5-t*2)],1);
        numEventsTemp = nan(size(sig,1),1);
        for j = 1:size(sig,1)
            numEventsTemp(j) = size(sig{j},1);
        end
        nEvents{s,t} = numEventsTemp(numEventsTemp ~= 0);
    end
end
for s = 1:subNum
    for t = 1:2
        subplotMN(subNum,4,s,t*2-1)
        title(trigs{t})
        yyaxis left
        ylim([-1.3 1.3])
        verticalLine(0);
    end
end

% Get average number of events for each triggered average plus S.T.D. (NOT SEM)
for t = 1:length(trigs)
    nEventsTrig = cell2mat(nEvents(:,t));
    Mean = mean(nEventsTrig);
    Std = std(nEventsTrig);
    disp(['average events ' trigs{t} ' = ' num2str(Mean) ' +- ' num2str(Std) ' std']);
end


%% quantify "timing" of locomotion signaling by looking at trough or peak 
% of cross correlation and of triggered averages on acc/dec and trans
% Fig. 2I and Fig. S6F-G
vars = [1 2 3 4]; % cross-corr, accel, decel, trans - ref from allData
varNames = {'CrossCorr', 'Acc-trigAve', 'Dec-trigAve', 'Trans-trigAve'};
troughIdx = cell(length(subtypes),length(vars));
for s = subsPure
    for v = 1:length(vars)
        sig = allData{vars(v)}(good == s,:);
        if vars(v) == 4 % for trans-trig ave, look at trough value BEFORE 0
            sig = fliplr(sig); 
        elseif vars(v) == 3 % subtypes with peaks in cross-corr will have troughts in dec-trig ave, so mirror
            sig = -sig; 
        end
        if s == 4 %for anxa, flip so peaks are throughs (opposite to Calb and vglut)
            sig = -sig; 
        end
        % get troughs after trigger point
        % (which correspond to peaks in some cases or troughs before 
        % triggger point based on above transformations)
        troughIdx{s,v} = nan(size(sig,1),1);
        for i = 1:size(sig,1)
            [~,idx] = min(sig(i,neg0:end)); %find min
            if idx > neg0 %if min before trigger point
                idx = neg0;
            end
            troughIdx{s,v}(i) = idx/hz; %seconds
        end
    end
end

%pvals to see differences between Vglut and Calb
disp('-NO- BONF p-values timing difference VGlut vs Calb:')
for v = 1:length(vars)
    disp(['  ' varNames{v} ' = ' num2str(ranksum(troughIdx{1,v},troughIdx{2,v}))]);
end
% plot
figure;
set(gcf, 'Position',  [10 600 800 308])
for v = 1:length(vars)
    subplot(1,length(vars),v)
    plotSpread([troughIdx{1,v};troughIdx{2,v};troughIdx{4,v}],'distributionIdx',[ones(length(troughIdx{1,v}),1); ones(length(troughIdx{2,v}),1)*2; ones(length(troughIdx{4,v}),1)*3],'distributionColors',colors([1 2 4],:)) %Jonas (2023). plot spread points (beeswarm plot) (https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot), MATLAB Central File Exchange. Retrieved March 20, 2023.
    for i = 1:length(subsPure)
        s = subsPure(i);
        errorbar(i,mean(troughIdx{s,v}),std(troughIdx{s,v})/sqrt(size(troughIdx{s,v},1)),'LineStyle','none','Color','k','LineWidth',1,'CapSize',10)
        errorbar(i,mean(troughIdx{s,v}),0,'LineStyle','none','Color','k','LineWidth',1,'CapSize',20)
    end
    ylabel('Trough/peak lag (s)')
    xticks(1:3)
    xticklabels({'VGlut' 'Calb' 'Anxa'})
    ylim([0 1])
    xlim([0.5 3.5])
    title(varNames{v})
end

%% Get integral of cross-corr at positive lags and plot by depth
% Fig. S1H
ccPureSubs = allData{1};
ccInt = nan(size(ccPureSubs,1),1);
for e = 1:size(ccPureSubs)
    if good(e) ~= 0
        cc = ccPureSubs(e,:);
        tt = (1:length(cc))/hz;
        ccInt(e) = trapz(tt(ceil(length(cc)/2):end),cc(ceil(length(cc)/2):end));
    end
end

% plot by depth
bins = 1.15:0.5:4.3;
binCents = (bins(2:end)+bins(1:end-1))/2;
figure
set(gcf, 'Position',  [680   677   890   300])
for s = 1:length(subtypes)
    subplot(1,length(subtypes),s)
    depths = cell2mat(data6.depthG(good == s));
    depths2 = depths + rand(size(depths))*0.2-0.1;
    ccIntS = ccInt(good == s);
    scatter(ccIntS,depths2,[],colors(s,:),'filled');
    set(gca, 'YDir', 'reverse');
    ylim([1.4 4])
    xlim([-0.1 0.1])
    verticalLine(0); 
    
    %moving average
    movAve = nan(size(binCents));
    for j = 1:length(binCents)
        subIdx = depths >= bins(j) & depths < bins(j+1);
        if sum(subIdx) > 1
            movAve(j) = mean(ccIntS(subIdx));
        end
    end
    hold on
    plot(movAve,binCents,'k','LineWidth',2)
end
subplot(1,length(subtypes),1)
ylabel('Depth (mm)')
subplot(1,length(subtypes),ceil(5/2))
xlabel('Int cross-corr')



%% Cross-corr with Acc averaged PER MOUSE
% Fig. S6I
figure
set(gcf, 'Position',  [200 42 300 954])
colorsL = [C.lightRed{1};C.lightOrange{1};C.lightGreen{1}; C.lightAqua{1}; C.lightGrey{1}];
crossCorrMice = cell(5,1);
crossCorrMiceN = cell(5,1);
for s = 1:length(subtypes)
    recIdx = double(good == s); %get recordings for subtype
    % find which recordings are from same mice
    mice = cell(30,1); 
    goodIdx = find(recIdx);
    mouseCount = 0;
    for j = 1:length(goodIdx)
        e = goodIdx(j);
        mouseID = data6.Mouse{e};
        if ~any(strcmp(mice,mouseID)) % new mouse
            mouseCount = mouseCount+1; % assign 'mouseCount' as ID number
            mice{mouseCount} = mouseID; % save in 'mice' to check for other recordings from same micce
        end
        recIdx(e) = mouseCount; %if new mouse then new mouse Count. If not, use previous mouse count.
        % THIS REQUIRES THE RECORDINGS TO BE SORTED BY MOUSE!! Standard
    end
    
    % for each mouse, get average of averages
    crossCorrMice{s} = nan(mouseCount,size(allData{1},2));
    for j = 1:mouseCount
        eachRec = allData{1}(recIdx == j,:);
        crossCorrMice{s}(j,:) = mean(eachRec,1); %get mean
    end
    % get PCs to sort
    scoresMice = (crossCorrMice{s}-mu)*coeff;
    % sort by PC1/Pc2 angle in same way as previously
    scoresMiceNorm = scoresMice./std(scores,'omitnan');
    anglesMice = atan2d(scoresMiceNorm(:,2),scoresMiceNorm(:,1));
    anglesMice = mod(anglesMice, 360);
    angleSplit = [315 45 135 135 45];
    subAngles = anglesMice;
    subAngles(subAngles < angleSplit(s)) = subAngles(subAngles < angleSplit(s))+360;
    [~,idxPCMice] = sort(subAngles); %get indexes to sort
    crossCorrMice{s} = crossCorrMice{s}(idxPCMice,:); %sort
    
    % normalize (for plotting heatmap)
    minmax = max(abs(crossCorrMice{s}),[],2);
    crossCorrMiceN{s} = crossCorrMice{s} ./ minmax;
    % plot heatmap
    subplotMN(length(subtypes),2,s,2)
    imagesc(crossCorrMiceN{s})
    caxis([-1 1])
    title([subtypes{s} ' - ' num2str(mouseCount) 'm']);
    verticalLine(win+1);
    xticks([0 win win*2])
    xticklabels([-win/100 0 win/100])
    xlabel('Time (s)')
    % plot averages
    subplotMN(length(subtypes),2,s,1)
    Mean = mean(crossCorrMice{s},1,'omitnan');
    Sem = std(crossCorrMice{s},[],1,'omitnan')/sqrt(size(crossCorrMice{s},1));
    ran = (-win:win)/100;
    fill([ran fliplr(ran)], [Mean+Sem fliplr(Mean-Sem)],colorsL(s,:),'EdgeColor',colors(s,:))
    ylim([-0.11 0.11])
    verticalLine(0);
end


%% Plot raw examples plus cross-corr for those examples
% Fig. 2C
load colors.mat C
colorsRaw = [C.red{1}; C.orange{1}; C.aqua{1}];
if snc == 0
    bestMice = {'VGlut-2298-20220324-0003' 'Calb-1789-20210406-0003' 'Anxa-J075-20220308-0004'};
    xlims = [0 255 300]; % for each example, best Times
    
    plotWin = 35; % in sec
    clear x
    bestIdx = nan(length(bestMice),1);
    % find miceIdx
    for i = 1:length(bestMice)
        subOrder = findMouse(data6,bestMice{i});
        if size(subOrder,1) > 1 %when there are two fibers for that recording, get striatal fiber
            cont = 1;
            j = 1;
            while cont == 1
                if (strcmp(data6.chG{subOrder(j,1)},'snc') && snc == 1) || (~strcmp(data6.chG{subOrder(j,1)},'snc') && snc == 0)
                    cont = 0;
                    bestIdx(i) = subOrder(j,1);
                else
                    j = j+1;
                end
            end
        else
            bestIdx(i) = subOrder(1);
        end
    end
    % plot
    vars = {'chMov','Acceleration', 'chGreen405','chGreen'};
    figure
    set(gcf, 'Position',  [406 291 660 678])
    scales0 = [0.2 0.5 0.1 0.1]; % vel, acc, DF/F twice - this is for scalebars
    for i = 1:3
        e = bestIdx(i);
        vals = cell(4,2);
        for v = 1:length(vars)
            trace = data6.data{e}.(vars{v}){1};
            if any(strcmp(vars{v},{'chGreen405','chGreen'}))
                trace = trace * data6.norm{e}(1); %de normalize DF/F
            end
            trace(isnan(trace)) = trace(circshift(isnan(trace),1)); %remove nans
            traceS = fastsmooth(trace,10); % smooth for plotting
            vals{v} = traceS;
        end
        % scale for plotting - and get scalebars
        for j = 1:4
            vals{j,2} = [0 scales0(j)];
        end
        for j = 2:-1:1
            vals{1,j} = (vals{1,j}-min(vals{1,1}))./(max(vals{1,1})-min(vals{1,1}))*1.2-0.1;
            vals{2,j} = (vals{2,j})./(max(abs(vals{2,1})))*1.2+0.5;
            vals{3,j} = (vals{3,j}-min(vals{4,1}))./(max(vals{4,1})-min(vals{4,1}))*1.2-0.1; %use chGreen to scale
            vals{4,j} = (vals{4,j}-min(vals{4,1}))./(max(vals{4,1})-min(vals{4,1}))*1.2-0.1; %use chGreen to scale
        end
        
        % get scalebars
        scaleTicks = cell2mat(vals(1:3,2));
        scaleTicks = scaleTicks + [0 0; 1 1; 2 2];
        scaleTicks2 = reshape(scaleTicks',1,6);
        scaleNums = fliplr(upsample(fliplr(scales0(1:3)),2));
        
        % plot
        subplot(3,1,i);
        hold on
        plot(x(vals{1},100),vals{1},'Color',[0 0 0]);
        plot(x(vals{2},100),vals{2}+1,'Color',C.grey{1});
        horizontalLine(1.5);
        plot(x(vals{3},100),vals{3}+2,'Color',C.lightBlue{1});
        plot(x(vals{4},100),vals{4}+2,'Color',colorsRaw(i,:));
        
        title(bestMice{i})
        xlim([xlims(i) xlims(i)+plotWin])
        yticks(scaleTicks2)
        yticklabels(scaleNums)
    end
end


%% compare Anxa and Calb in same region
% see FiberLocations2 script for calculations of XYZ coordinates for recordings
% Fig. 3A
if snc == 0  
    % set area to include recordings within
    scale = 35.6667; %see FiberLocations2 script
    secEdge = 204.4; %in pixels
    radius = 0.5; %mm
    radiusPix = radius*scale;
    center1 = [91,101]; %slice bregma +0.86
    center2 = [306,101]; %slice bregma +0.50
    shift = 1; % if shift == 1, combine all slices onto bregma +0.86. Otherwise, plot separately
    %get recordings within these circles
    loc = cell2mat(data6.RecLocG);
    locYes1 = (loc(:,1)-center1(1)).^2 + (loc(:,2)-center1(2)).^2 <= radiusPix^2;
    locYes2 = (loc(:,1)-center2(1)).^2 + (loc(:,2)-center2(2)).^2 <= radiusPix^2;
    locYes = locYes1 | locYes2;
    % load brain images and keep only the ones of interest (see 'shift')
    load('brain_Empty.mat', 'brain')
    if shift == 1
        brain = [brain(:,round(secEdge+1):round(secEdge*2)) brain(:,round((secEdge+1)*5):round(secEdge*6))];
    end
    figure
    set(gcf, 'Position',  [50 551 720 420])
    subsPure = [2 4]; %Calb and Anxa
    for s = 1:2
        % plot cross-corr for only recordings within area - average and heatmap
        idx = good == subsPure(s) & locYes;
        crossCorr2(data6, idx, 'G', 'A', win,[2,4,s*4-1],1:sum(idx),normHeat); 
        % plot brain slices and area included
        subplotMN(2,4,s,[1,2])
        imshow(brain)
        xlim([1,secEdge*2])
        hold on
        if shift == 1
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
            if shift == 1
                RecLoc = data6.RecLocGshift{idx(j)};
            else
                RecLoc = data6.RecLocG{idx(j)};
            end
            scatter(RecLoc(1),RecLoc(2),[],[0.8 0.8 0.8],'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceAlpha',0.5)
        end
    end
end
             

%% plot cross-corr but heatmaps sorted by integral at pos lag
% Crosscorr DF/F vs acc
% Fig. S1F

%sort by cross-corr integral at pos lag
idxInt = nan(size(data6,1),1);
for i = 1:5
    cc = ccInt(good == i); %from a previous section when we calculated the integral of the cross-corr at positive lags
    [~,subIdx] = sort(cc); 
    idxInt(good == i) = subIdx;
end
% plot cross-corr sortex by int
figure
set(gcf, 'Position',  [1 42 300 954])
for i = 1:5
    crossCorr2(data6, good == i, 'G', 'A', win,[5,2,i*2-1],idxInt(good == i),normHeat);
end
for i = 2:2:5*2
    subplot(5,2,i)
    caxis([-1,1])
    subplot(5,2,i-1)
    title('Cross corr with Acc')
end

%% plot relationship between DF/F and vel
% bin velocity and calculate the average DF/F for each bin and subtype
% Fig. S7C
binVsize = 0.1;
binV = [-inf -0.15:binVsize:0.75 inf];
% only for Calb, Vglut, and Anxa
goodAny = find(good==1 | good==2 | good ==4);
% average all DF/F for velocities within iach bin
velDF = nan(size(data6,1),length(binV)-1);
for i = 1:length(goodAny)
    e = goodAny(i);
    v = data6.data{e}.chMov{1};
    g = data6.data{e}.chGreen{1};
    for j = 1:length(binV)-1
        vI = v >= binV(j) & v < binV(j+1);
        velDF(e,j) = mean(g(vI),'omitnan');
    end
end
% plot
ran = (binV(1:end-1)+binV(2:end))/2;
ran(1) = ran(2)-binVsize;
ran(end) = ran(end-1)+binVsize;
figure
for i = [1 2 4]
    m = mean(velDF(good==i,:),1,'omitnan');
    s = std(velDF(good==i,:),1,'omitnan')./sqrt(sum(~isnan(velDF(good==i,:))));
    if i == 1
        colors = [C.lightRed{1};C.red{1}];
    elseif i == 2
        colors = [C.lightYellow{1};C.orange{1}];
    elseif i == 3
        colors = [C.lightGreen{1};C.green{1}];
    elseif i == 4
        colors = [C.lightAqua{1};C.aqua{1}];
    elseif i == 5
        colors = [C.lightGrey{1};C.grey{1}];
    end
    hold on;plot(ran,m,'Color',colors(2,:))
    hold on;plot(ran,m-s,'Color',colors(1,:))
    hold on;plot(ran,m+s,'Color',colors(1,:))
end
xlabel('Velocity bin (m/s)')
ylabel('Mean Norm %\DeltaF/F')


%% save 'good' and 'angles' for future access
scoresNorm = scores./std(scores,'omitnan');
angles = atan2d(scoresNorm(:,2),scoresNorm(:,1));
LocomAngles = mod(angles, 360);
LocomPC1 = scoresNorm(:,1);
LocomPC2 = scoresNorm(:,2);
recID = cell(size(data6,1),1);
for e = 1:size(data6,1)
    recID{e} = data6.data{e}.Properties.RowNames{1};
end
LocomCrossCorr = num2cell(ccData470,2);

cc = data6(:,[1:3 5:6 10:11]);
goodRun = [table(recID) cc table(LocomPC1) table(LocomPC2) table(LocomAngles) table(LocomCrossCorr)];
if snc == 0
    save([dataProcessingFolder '\good_Run.mat'], 'goodRun');
else
    save([dataProcessingFolder '\goodSNc_Run.mat'], 'goodRun');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to find recording in data6
function idx = findMouse(data,mouse)
% script to find a recoding within a dataset
% mouse format: Exp-Mouse-Date-RecNum (ex: 'VGlut-2298-20220324-0003')

% get name of each recording in same format as input mouse
idxList = cell(size(data,1),2);
for rec = 1:size(data,1)
    idxList{rec,1} = [rec,1];
    idxList{rec,2} = data.data{rec}.Properties.RowNames{1};
end
% find mouse
idx1 = ismember(idxList(:,2),mouse);
if sum(idx1) == 0
    error('Mouse not found')
end
idx = nan(sum(idx1),2);
idx2 = find(idx1);
for i = 1:sum(idx1)
    idx(i,:) = idxList{idx2(i),1};
end
end



