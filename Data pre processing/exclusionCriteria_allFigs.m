function data6 = exclusionCriteria_allFigs(data5,newData)
% This scripts determines the exclusion criteria for each recording and
% analysis type. For this it:
    % Separates each session into individual recordings (in data5 each
        % session is one row, in data6 each recording is one row).
    % Crops the ends of recordings if there is a recording artifact in the
        % 405 channel (Picoscope error) - it automatically detects putative
        % problems and allows you to manually crop them off
    % Automatically finds artifacts in the 405 channel (see section for
        % details) and excludes them from the recording (it saves a new 
        % data variable for each fiber, 'Bad405G' and 'Bad405R'. If more 
        % than 5% of a recording has artifacts, it excludes the whole 
        % recording.
    % Calculates the time the mouse spent running in each recording 
        %('RunTime') 
    % Calculates the cross-correlation between the 405 channel for each 
        % fiber and get the max value ('Acc405'), used to determine
        % eligibility for locomotion analysis
    % Calculates the cross-correlation between the 405 in one fiber vs the
        % 470 of the other and vice versa, and averages them to get a max 
        % value ('RG405') uses to determine eligibility for dual
        % comparisons between simultaneous recordings.    
    % Adds empty 'dup' and 'flip' variables - these will be useful in the 
        % future, because of having two fibers: a recording with data for 
        % chG and chR will be duplicated and one of the copies with have 
        % chG/chR flipped

% It then creates variables showing whether each recording passes exclusion
% criteria:
    % 'exSig2Noise' - which recordings have high enough signal to noise 
        % ratio (1 = good, 0 = bad)
    % 'Bad405All' - which recordings have too many artifacts in the 405
        % channel (1 = bad, 0 = good, nan = not enough sig2noise)
    % 'exRun' - which recordings have enough running time to include for
        % movement analysis (1 = good, 0 = bad)
    % 'exStr' - which recordings are made in axons vs cell bodies (1 =
        % axons, 0 = cell bodies snc/vta, nan = empty)
    % 'ex405Acc' - which recordings must be excluded from locomotion 
        % analyhave due to high correlations between the 405 ch and
        % acceleration  (1 = good, 0 = bad, nan = not enough running)
    % 'exBothCh' - which recording pairs have both fibers (chG and chR)
        % with signal-to-noise ratios above thresh, for dual comparisons
        % (1 = good, 0 = bad)
    % 'exSNcStr' - which recording pairs have one fiber in SNc and one in
        % striatum, for simultaneous comparisons of somas and axons 
        % (1 = good, 0 = bad)
    % 'exRG405' - which recording must be EXCLUDED from simultaneous
        % fiber-to-fiber comparisons due to high correlation between the
        % 405 and 470 channels of opposite fibers (1 = good, 0 = bad)
        
% The input newData is not necessary - if missing, this script will analyze
% all recordings in data5. If present, it will only analyze the data listed
% there, load data6 from dataProcessingFolder and append the newly analyzed
% data to that existing data6. The input data5 isn't necessary either - if
% not provided, the script will load it from dataProcessingFolder.

% Output: data6 = master table of all analyzed recordings with exclusion
% criteria determined. It is automatically saved in dataProcessingFolder

% This scripts was written for Windows at MATLAB R2021a
        
%% set MAIN FOLDER and DATA PROCESSING FOLDER (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder');

%% set important variables
hz = 100;

%% get data to process
if ~exist('data5','var') || isempty(data5)
    load([dataProcessingFolder '\data5.mat'], 'data5');
end
if ~exist('newData','var') %if newData not provided, process all data in data5
    data6 = data5;
else
    % find data to process from newData/data5
    badIdx = zeros(size(newData,1),1);
    for n = 1:size(newData,1)
        mouse = newData{n,5};
        date = newData{n,6};
        badIdx(n) = find(strcmp(data5.Mouse,mouse) & strcmp(data5.Date,date));
    end
    data6 = data5(badIdx,:);
end

%% Separate each recording per row in with correct depths
% In data5, each session is one row. In data6, each recording is one row.
% Also set correct depths for chG and chR
newSetup = 20211203; %from this date on the two fibers could move 
    % independently, so there are two sets of depths. Before, we only had
    % the depth difference between the fiber tips.
for rec = size(data6,1):-1:1
    depths = size(data6.data{rec},1);
    data6 = [data6(1:rec-1,:); repmat(data6(rec,:),depths,1); data6(rec+1:end,:)];
    for d = 1:depths % separate not only data and depths but also other metadata
        data6.data{rec+d-1} = data6.data{rec+d-1}(d,:);
        data6.cropStart{rec+d-1} = data6.cropStart{rec+d-1}(d);
        data6.base{rec+d-1} = data6.base{rec+d-1}(d,:);
        data6.norm{rec+d-1} = data6.norm{rec+d-1}(d,:);
        data6.sig2noise{rec+d-1} = data6.sig2noise{rec+d-1}(d,:);
        data6.depthG{rec+d-1} = data6.data{rec+d-1}.Depth{1};
        if max(size(data6.DepthDiff{rec+d-1})) > 1 || (depths == 1 && str2double(data6.Date{rec+d-1}) > newSetup)
            data6.depthR{rec+d-1} = data6.DepthDiff{rec+d-1}(d);
        else
            data6.depthR{rec+d-1} = data6.data{rec+d-1}.Depth{1} - data6.DepthDiff{rec+d-1};
        end
    end
end
% reorder vars and remove depthDiff variable (no need for them anymore!)
if size(data6,2) ~= 17 || ~strcmp(data6.Properties.VariableNames{16},'depthG')
    error('too many variables??')
else
    data6 = [data6(:,1:4) data6(:,16:17) data6(:,5:10) data6(:,12:15)];
end


%% Crop bad ends in 405 channel
% for some reason, sometimes Picoscope stops recording one channel before
% the others and causes an edge artifact at the end. If this has ocurred,
% crop end manually
bad405thresh = 1.5;
figure
for rec = 1:size(data6,1)
    chGbad405 = any(strcmp(data6.data{rec}.Properties.VariableNames,'chGreen')) && any(data6.data{rec}.chGreen405{1} > bad405thresh);
    chRbad405 = any(strcmp(data6.data{rec}.Properties.VariableNames,'chRed')) && any(data6.data{rec}.chRed405{1} > bad405thresh);
    if  chGbad405 || chRbad405
        subplot(2,1,1)
        try
            plot(data6.data{rec}.chGreen405{1})
        catch
            cla
        end
        title([num2str(rec) ' - click left to not crop, click in plot to crop END'])
        try
            subplot(2,1,2)
            plot(data6.data{rec}.chGRed405{1})
        catch
            cla
        end
        x = ginput(1);
        x = round(x(1));
        if x < length(data6.data{rec}.chMov{1}) && x > 1
            data6.cropStart{rec} = [data6.cropStart{rec}(1), data6.cropStart{rec}(1)+x];
            for v = 2:size(data6.data{rec},2)
                if length(data6.data{rec}.(v){1}) > 2
                    data6.data{rec}.(v){1} = data6.data{rec}.(v){1}(1:x);
                end
            end
        end
    end
end
close()

%% Set exclusion criteria for transient sig2noise ('exSig2Noise')
% only use recordings with signal-to-noise ratios above 10
    % exSig2Noise: Good = 1; bad = 0
threshS2N = 10;
sig2noise = cell2mat(data6.sig2noise);
data6.sig2noise = sig2noise;
goodPeaksG = sig2noise(:,1) > threshS2N;
goodPeaksR = sig2noise(:,2) > threshS2N;
data6.exSig2Noise = double([goodPeaksG goodPeaksR]);


%% Exclude bad recordings with lots of 405 artifacts, or crop out bad 405 sections if small
% Fist identify the noise in the 405 channel for each fiber, and determine
% an upper bound (0.5 percentile of abs noise). Any sections of the 405
% channel more than 3 times above that noise max are automatically
% excluded, as well as any sections above the noise max for longer than 0.2
% seconds (less than half of a transient width). For any excluded times,
% exclude an additional 0.1 s on either side, just to be safe.
% If more than 5% of any recording needs to be excluded, exclude the entire
% recording.
    % Bad405All: Bad = 1, good = 0, nan = not enough sig2noise to process
    
noiseVeryBad = 3; % exclude any points above 3 times the max noise level in 405
noiseBadTimeThresh = 20; % exclude any sections of 405 longer than 0.2s 
% under half a transient in duration)with values above noise
noisePer = 0.005; % noise percentile to use to exclude values (top 0.5 percentile)
extraRemoveWin = 10; %remove 10 aditional bins on either side of bad sections just in case
maxBadBinPer = 0.05; %if more than 5 % of recording needs to be excluded due to artifacts, exclude whole recording

[recordings,fibers] = find(data6.sig2noise > threshS2N); %find all recordings (in both fibers) with sig2noise above thresh
data6.Bad405All = nan(size(data6,1),2);

for recID = 1:length(recordings)
    rec = recordings(recID);
    fib = fibers(recID);
    % get 405 channel for all recordings with 470 above sig2noise ratio
    if fib == 1
        dff = data6.data{rec}.chGreen405{1};
    elseif fib == 2
        dff = data6.data{rec}.chRed405{1};
    end
    dff_notNorm = dff * data6.norm{rec}(1,fib); %de-normalize
    % remove nan (without changing length!)
    bad = find(isnan(dff_notNorm));
    try
        dff_notNorm(bad) = dff_notNorm(bad+1);
    catch
        dff_notNorm(bad) = dff_notNorm([bad(1:end-1)+1;bad(end)-1]);
    end
    % get NOISE - use the top 8th percentile of noise to find anything that 
    % goes above it (artifacts)
    dff_smooth = fastsmooth(dff_notNorm,10);
    noiseAbs = abs(dff_notNorm-dff_smooth);
    noiseSort = sort(noiseAbs);
    maxNoise = noiseSort(round(length(noiseSort)*(1-noisePer)));
    % remove anything that is too far above or below maxNoise, even if short
    idxVeryBad = dff_notNorm > maxNoise*noiseVeryBad | dff_notNorm < -maxNoise*noiseVeryBad;
    % find points that go above maxNoise
    % these will only be removed if they are a certain width
    idxBad = dff_notNorm > maxNoise | dff_notNorm < -maxNoise;
    
    % keep only regions above maxNoise that are longer than 0.2 s (less than half a transient)
    % first find bad section onsets and offsets - and fix edges if necessary
    idxBadOn = find(diff(idxBad) == 1);
    idxBadOff = find(diff(idxBad) == -1);
    if ~isempty(idxBadOn) && ~isempty(idxBadOff)
        if idxBadOn(1) < idxBadOff(1) && length(idxBadOn) == length(idxBadOff)
            %do nothing: each onset has a matching offset
        elseif idxBadOn(1) < idxBadOff(1) && length(idxBadOn) > length(idxBadOff)
            % there is an onset at the end with no offset - remove last onset
            idxBadOn = idxBadOn(1:end-1);
        elseif idxBadOn(1) > idxBadOff(1) && length(idxBadOn) < length(idxBadOff)
            % there is an offset at the start with no onset - remove first offset
            idxBadOff = idxBadOff(2:end);
        elseif idxBadOn(1) > idxBadOff(1) && length(idxBadOn) == length(idxBadOff)
            % both the above conditions are true, so remove the first off
            % and last onset
            idxBadOff = idxBadOff(2:end);
            idxBadOn = idxBadOn(1:end-1);
        end
        % If some bad sections are too close to edge, remove even if short
        if idxBadOn(end) <= length(dff_notNorm) - noiseBadTimeThresh
                idxVeryBad(idxBadOn(end):end) = true;
        end
        if idxBadOff(1) >= noiseBadTimeThresh
            idxVeryBad(1:idxBadOff(1)) = true;
        end
        % Calculate duration of each bad section and exclude if too long
        % Also remove 10 extra bins (0.1s) on either side of bad sections,
        % just in case.
        dur = idxBadOff-idxBadOn;
        badIdx = find(dur >= noiseBadTimeThresh);
        for i = 1:length(badIdx)
            idxOn = max(idxBadOn(badIdx(i))-extraRemoveWin,1);
            idxOff = min(idxBadOff(badIdx(i))+extraRemoveWin,length(dff_notNorm));
            idxVeryBad(idxOn:idxOff) = true; 
        end
    end
    
    % if too much is bad (more than 5%) - exclude recording because of noise
    binsDeleted = sum(idxVeryBad)/length(idxVeryBad);
    if binsDeleted > maxBadBinPer
        data6.Bad405All(rec,fib) = 1; %bad = 1
    else
        data6.Bad405All(rec,fib) = 0; %indicate it was processed and is good
    end
    % add bad bins to table - these will not be removed yet because they
    % can be different for fiber 1 and fiber 2
    if fib == 1
        data6.data{rec}.Bad405G{1} = idxVeryBad;
    else
        data6.data{rec}.Bad405R{1} = idxVeryBad;
    end
end


%% Get time spent running ('timeRun', from 'Run' variable)
timeRun = zeros(size(data6,1),1);
for rec = 1:size(data6,1)
    timeRun(rec) = sum(data6.data{rec}.Run{1})/hz;
end
data6.timeRun = timeRun;

%% Set exclusion criteria for running time ('exRun')
% run at least 100s in session
    % exRun: good = 1, bad = 0
    
threshRun = 100;
timeRun = data6.timeRun;
goodRun = timeRun > threshRun;
data6.exRun = double(goodRun); % good = 1


%% Add exclusion criteria variable for str/nac vs snc/vta ('exStr')
% this indicates if recording was made in axons (exStr = 1) or cell bodies 
% (exStr = 0) or nothing (exStr = nan)

goodStrG = double(~strcmp(data6.chG,'snc') & ~strcmp(data6.chG,'vta'));
goodStrR = double(~strcmp(data6.chR,'snc') & ~strcmp(data6.chR,'vta'));
goodStrG(strcmp(data6.chG,'o')) = nan;
goodStrR(strcmp(data6.chR,'o')) = nan;
data6.exStr = [goodStrG goodStrR];


%% set exclusion criteria for Movement analysis ('Acc405' and 'ex405Acc')
% This is used to determine which recordings can be used for movment
% analysis and which ones have too great of a movement artifact
% Use the cross correlation between the 405 channel and acceleration,
% calculated in the same way as will be done in the paper
    % ex405Acc: 1 = good, 0 = bad, nan = non enough running
    
win = 50;
threshPeak405 = 0.1; % max cross-corr with 405 channel

maxAcc405G = nan(size(data6,1),1);
maxAcc405R = nan(size(data6,1),1);
for rec = 1:size(data6,1)
    if data6.exRun(rec) == 1
        % calculate cross-corr between each 405 (if present) and
        % acceleration, and save max value in 'Acc405'
        try
            cc = crossCorr(data6.data{rec},'GA',win); % see function below, same as crossCorr2 script
            maxAcc405G(rec,:) = max(abs(cc));
        end
        try
            cc = crossCorr(data6.data{rec},'RA',win);
            maxAcc405R(rec,:) = max(abs(cc));
        end
    end
end
data6.Acc405 = [maxAcc405G maxAcc405R];
% Check which recordings have high 405-acceleration max to exclude, and
% save in 'ex405Acc'
good405AccG = data6.Acc405(:,1) < threshPeak405;
good405AccR = data6.Acc405(:,2) < threshPeak405;
data6.ex405Acc = double([good405AccG good405AccR]); % 'ex405Acc' = movement artifacts below thresh


%% set exclusion criteria for both ch above sig2noise threshold ('exBothCh')
% necessary when comparing simultaneous recordings from both fibers
% threshS2N from above
    % exBothCh: good = 1, bad = 0
    
sig2noise = data6.sig2noise;
goodBothCh = ~(any(sig2noise < threshS2N,2)) & ~any(isnan(sig2noise),2);
data6.exBothCh = goodBothCh;


%% 'exSncStr' only pairs of recordings from SNc somas and axons in str
% check if one fiber is in snc and the other one isn't empty
% assumes that no pairs of recordings are both from SNc (never done)
    % exSncStr: good = 1; bad = 0
    
chGchR = [data6.chG data6.chR];
sncYN = double(contains(chGchR,'snc'));
empty = strcmp(chGchR,'o');
sncYN(empty) = nan;
bothSNcStr = sum(sncYN,2) == 1;
data6.exSNcStr = bothSNcStr;


%% set exclusion criteria for noise in chG vs chR cross corr (using 405)
% This is used to determine which pairs of chG/chR recorings can be used
% for simultaneous comparisons
% Check for cross-correlations between the noise in one ch and the DF/F in
% the other channel. Then get the mean of the two comparisons and get the
% max value.
    % exRG405: good = 1, bad = 0
    
win = 50;
threshPeak405 = 0.12;

maxRG405 = nan(size(data6,1),1);
for rec = 1:size(data6,1)
    if data6.exBothCh(rec) == 1 % only check recordings pairs where both have good sig2noise
        ccTemp1 = crossCorr(data6.data{rec},'GR405',win); % see function below, same as crossCorr2SNcStr script
        ccTemp2 = crossCorr(data6.data{rec},'G405R',win);
        cc = mean([ccTemp1,ccTemp2],2);
        
        maxRG405(rec) = max(abs(cc));
    end
end

data6.RG405 = maxRG405;
goodRG405 = data6.RG405 < threshPeak405;
data6.exRG405 = double(goodRG405);


%% add dup and flip variables
% these will be useful in the future, because of having two fibers: 
% a recording with data for chG and chR will be duplicated and one of the
% copies with have chG/chR flipped
    % flip = indicates the two fibers are flipped (chG is now chR and viceversa)
    % dup = indicates the recording is a duplicate
data6.flip = false(size(data6,1),1);
data6.dup = false(size(data6,1),1);


%% add data to data6 (if new) and save
%add empty columns for the recording locations - these will be added on a
%different script and need empty columns to combine old and new data tables
for rec = 1:size(data6,1)
    data6.RecLocR{rec} = [nan nan]; 
    data6.RecLocG{rec} = [nan nan];
end
% add to main data6
if exist('newData','var')
    try
        data6New = data6;
        load([dataProcessingFolder '\data6.mat'], 'data6');
        data6 = [data6; data6New];
    catch
    end
    % sort
    data6 = sortrows(data6,[1 2 7]); 
end
% save
save([dataProcessingFolder '\data6.mat'], 'data6', '-v7.3');
disp('data6 saved!')

end



%% FUNCTION to get cross correlations (same as CrossCorr2 function)
function cc = crossCorr(data,RGA,win)          
    % get data
    if strcmp(RGA,'RA')
        ch1 = data.chRed405{1};
        ch2 = data.Acceleration{1};
    elseif strcmp(RGA,'GA')
        ch1 = data.chGreen405{1};
        ch2 = data.Acceleration{1};
    elseif strcmp(RGA,'GR405')
        ch1 = data.chGreen{1};
        ch2 = data.chRed405{1};
    elseif strcmp(RGA,'G405R')
        ch1 = data.chGreen405{1};
        ch2 = data.chRed{1};
    else
        error('')
    end

    % remove nan
    bad = isnan(ch2) | isnan(ch1);
    if sum(bad) > 15
        warning('Please check, too many nans in data');
    end
    ch2(bad) = [];
    ch1(bad) = [];

    % keep only run times - but only for GA or RA
    if strcmp(RGA,'RA') || strcmp(RGA,'GA')
        run = logical(full(data.Run{1}));
        if any(run)
            ch2 = ch2(run);
            ch1 = ch1(run);
            try
                cc = crosscorr(ch2,ch1,win);
            catch
                cc = nan(win*2+1,1);
            end
        else
            cc = nan(win*2+1,1);
        end
    else
        try
            cc = crosscorr(ch2,ch1,win);
        catch
            cc = nan(win*2+1,1);
        end
    end
end











    
