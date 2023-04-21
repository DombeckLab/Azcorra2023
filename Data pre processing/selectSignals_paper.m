function data5 = selectSignals_paper(data4, newData)
% This script selects variables of interest for further analysis,
% particularly for triggered averages, by making new binnary vectors for
% each recording where ONESs mark the location of different variables:
    % 'Licking' - times where the mouse was licking.    
    % 'AirPuff' and 'Light' - times where the air puff or ligh stimuli were
        % delivered
    % 'RewardsAll' - all times where rewards were delivered
    % 'Rewards' - times where any rewards were delivered but only if the
        % mouse licks (consumes the reward) within a certain time. It also
        % excludes times were the lick sensor was not working properly.
    % 'RewardShort' and 'RewardLong' - rewards are also split based on the
        % duration of the reward delivery signal pulse.
    % 'RewardRest', 'RewardShortRest', and 'RewardLongRest' - only rewards
        % that were delivered while the mouse was at rest and for which the
        % mouse remained at rest after the reward was delivered (all
        % rewards in 'RewardRest' split into Short and Long)
    % 'AirPuffRest' and 'LightRest' - same as for rewards
    
    % 'MovOnOff' - times when the mouse was resting (no movement) vs all
        % other times when the mouse was moving 
    % 'Run' - times when the mouse was running, defined as moving bouts
        % longer than 0.5 s with a mean velocity greater than 0.2 m/s, and
        % excluding times 5 s after rewards/air puff/light stimuli were
        % delivered
    % 'RunRew' - times when the mouse was running defined as above, but
        % without excluding stimuli delivery times
    % 'AccOn' and 'DecOn' - onsets of accelerations and decelerations above
        % a threshold ocurring during 'Run' periods defined as above, and
        % not too close to other large accelerations or decelerations (to
        % prevent confounds)
    % 'MovOn' and 'MovOff' - clean movement onsets and offsets: for 'Run' 
        % periodslonger than 3s where the velocity reached at least 0.4 m/s
        % within the first 0.75 s, and with a first acceleration of at 
        % least 1 m/s2 and no negative velocity (backwards movements) at
        % the movement initiation - all reverse for offsets.
        
    % 'peaksG' and 'peaksR' - large fluorescence transient peaks that are 
        % well isolated from other transients, for both fibers if present
    % 'peaksGRun' and 'peaksRRun' - same as previous but only those
        % ocurring during running periods as defined above for 'Run' 
        % (stimuli times excluded)
        
% see each Section for more details on each of the above variables
    
% This script also removes the reward, air puff, and light variables from
% recordings where none were delivered, to reduce file size.

% The input newData is not necessary - if missing, this script will analyze
% all recordings in data4. If present, it will only analyze the data listed
% there, load data5 from dataProcessingFolder and append the newly analyzed
% data to that existing data5. The input data4 isn't necessary either - if
% not provided, the script will load it from dataProcessingFolder.

% Output: data5 = master table of all analyzed recordings with the above
% variables identified. It is automatically saved in dataProcessingFolder

% This scripts was written for Windows at MATLAB R2021a

%% set MAIN FOLDER and DATA PROCESSING FOLDER (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
load([filePath '\MatlabFolders.mat'], 'dataProcessingFolder');

%% get pre-processed data (see data_procesing.m script)
% if data4 was not provided, load
if ~exist('data4','var') || isempty(data4)
    load([dataProcessingFolder '\data4.mat'], 'data4');
end
% if newData is not provided, analyze all of data4. If it is, find in data4
% the recordings specified in newData and analyze only those.
if ~exist('newData','var')
    data5 = data4;
else
    % find newData recordings in data4 (using mouse and date ID)
    idx = zeros(size(newData,1),1);
    for rec = 1:size(newData,1)
        mouse = newData{rec,5};
        date = newData{rec,6};
        idx(rec) = find(strcmp(data4.Mouse,mouse) & strcmp(data4.Date,date));
    end
    data5 = data4(idx,:);
    data5 = sortrows(data5,1:3);
end


%% remove rew/light/puff variables in non-reward recordings
% they are useless (all zeros) and take up space
binaryThresh = 0.05;
for rec = 1:size(data5,1)
    if strcmp(data5.RunRew{rec},{'run'}) %these recordings should not have rewards
        remove = zeros(1,size(data5.data{rec},1));
        for depth = 1:size(data5.data{rec},1) %check just in case they were mislabelled
            if any(data5.data{rec}.Reward{depth} > binaryThresh) || any(data5.data{rec}.AirPuff{depth} > binaryThresh) || any(data5.data{rec}.Light{depth} > binaryThresh)
                disp(['Exp ' num2str(rec) '-' num2str(depth) ' has rewards but should not'])
            else
                remove(depth) = 1;
            end
        end
        % remove variables
        data5.data{rec}.Light = [];
        data5.data{rec}.Licking = [];
        data5.data{rec}.Reward = [];
        data5.data{rec}.AirPuff = [];
    end
end

%% Convert reward, air puff, light, and licking variables to binary signals
binaryThresh = 0.05; 
for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec},1)
        if strcmp(data5.RunRew{rec},{'rew'}) %these recordings only have rewards
            data5.data{rec}.Light{depth} = data5.data{rec}.Light{depth} > binaryThresh;
            data5.data{rec}.Licking{depth} = data5.data{rec}.Licking{depth} > binaryThresh;
            data5.data{rec}.Reward{depth} = data5.data{rec}.Reward{depth} > binaryThresh;
            data5.data{rec}.AirPuff{depth} = data5.data{rec}.AirPuff{depth} > binaryThresh;
        end
    end
end

%% Select reward signals
% First check that all recordings that are supposed to have rewards do
% Keep only rewards where mice licked after but not continuously (this
% indicates an error with the lickometer)
% Also exclude rewards too close to the edge of the recording.
% Then split rewards into large and small by the duration of the reward
% delivery signal.

rewThreshOld = 5.5; % duration threshold (bins) separating long/short rewards
rewThreshNew = 2.5; % The tubbing was changed and re-callibrated
dateThresh = 20220223; % last day of old reward was 20220222
lickWin = 100; % 1s - window within which mice must lick to consume reward
for rec = 1:size(data5,1)
    try
        if strcmp(data5.RunRew{rec},{'rew'}) % recordings with rewards only
            for depth = 1:size(data5.data{rec},1)
                % check that each recording actually has rewards
                if ~any(data5.data{rec}.Reward{depth})
                    disp(['Exp ' num2str(rec) '-' num2str(depth) ' does not have rewards etc. but should'])
                    data5.data{rec}.RewardLong{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                    data5.data{rec}.RewardShort{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                    data5.data{rec}.RewardAll{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                else
                    % keep only rewards where mice licked within a window 
                    % after but not continuously. A New variable 'RewardAll' 
                    % has all rewards, keeping only the good ones in
                    % 'Reward'
                    data5.data{rec}.RewardAll{depth} = data5.data{rec}.Reward{depth};
                    data5.data{rec}.Reward{depth} = zeros(size(data5.data{rec}.RewardAll{depth}));
                    rewOn = find(diff(data5.data{rec}.RewardAll{depth}) == 1);
                    rewOff = find(diff(data5.data{rec}.RewardAll{depth}) == -1);
                    lick = data5.data{rec}.Licking{depth};
                    for i = 1:length(rewOn)
                        try %if rew too close to edge
                            lickInWin = sum(lick(rewOn(i):rewOn(i)+lickWin));
                            if lickInWin < lickWin && lickInWin > 5 %if lickInWin = lickWin, it means mouse was 'licking' the entire time, which means the lickometer was not working 
                                data5.data{rec}.Reward{depth}(rewOn(i)+1) = 1;
                                data5.data{rec}.Reward{depth}(rewOff(i)+1) = -1;
                            end
                        catch
                        end
                    end
                    data5.data{rec}.Reward{depth} = cumsum(data5.data{rec}.Reward{depth});
                    data5.data{rec}.Reward{depth} = sparse(data5.data{rec}.Reward{depth}); %sparse to take less space
                    
                    % separate small and large rewards based on duration
                    rewOn = find(diff(data5.data{rec}.Reward{depth}) == 1);
                    rewOff = find(diff(data5.data{rec}.Reward{depth}) == -1);
                    dur = rewOff - rewOn;
                    % get thresh based on date (reward tubbing was changed 
                    % and callibrated and thus small/large durations changed)
                    try
                        date = data5.Date{rec}{1};
                    catch
                        date = data5.Date{rec};
                    end
                    date = str2double(date);
                    if date >= dateThresh
                        rewThresh = rewThreshNew;
                    else
                        rewThresh = rewThreshOld;
                    end
                    % add long/short rewards to table
                    data5.data{rec}.RewardLong{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                    data5.data{rec}.RewardShort{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                    data5.data{rec}.RewardLong{depth}(rewOn(dur > rewThresh)) = 1;
                    data5.data{rec}.RewardLong{depth}(rewOff(dur > rewThresh)) = -1;
                    data5.data{rec}.RewardShort{depth}(rewOn(dur < rewThresh)) = 1;
                    data5.data{rec}.RewardShort{depth}(rewOff(dur < rewThresh)) = -1;
                    data5.data{rec}.RewardLong{depth} = cumsum(data5.data{rec}.RewardLong{depth});
                    data5.data{rec}.RewardShort{depth} = cumsum(data5.data{rec}.RewardShort{depth});
                    data5.data{rec}.RewardLong{depth} = sparse(data5.data{rec}.RewardLong{depth});
                    data5.data{rec}.RewardShort{depth} = sparse(data5.data{rec}.RewardShort{depth});
                end
            end
        end
    catch
        disp(['problem in reward selection for ' num2str(rec)])
    end
end


%% select reward signals during REST (but also try air puffs and light)
% these are rewards delivered while the mouse was not moving beforehand, 
% but also for which the mouse did not move after the delivery
% also try for air puffs, but never happens (they always move upon delivery)
winNoMoveShort = 40; % very strict no movement in window 0.4 s after event
winNoMoveLong = 75; % no movement in 0.75 s before, during and after event
accThreshShort = 1.5; 
accThreshLong = 2.5; 
vars = {'Reward', 'RewardLong', 'RewardShort', 'AirPuff', 'Light'};

for rec = 1:size(data5,1)
    if strcmp(data5.RunRew{rec},{'rew'}) && any(data5.data{rec}.Reward{depth})
        for v = 1:length(vars)
            var = vars{v};
            % create new REST variable
            data5.data{rec}.([var 'Rest']) = data5.data{rec}.(var);
            
            for depth = 1:size(data5.data{rec},1)
                % create new RewardRest etc variables
                data5.data{rec}.RewardRest{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                data5.data{rec}.RewardLongRest{depth} = zeros(size(data5.data{rec}.Reward{depth}));
                data5.data{rec}.RewardShortRest{depth} = zeros(size(data5.data{rec}.Reward{depth}));              
                % get acceleration
                acc = data5.data{rec}.Acceleration{depth}; % get acceleration
                % get variable on/off
                idxOn = find(diff(data5.data{rec}.(var){depth}) == 1);
                idxOff = find(diff(data5.data{rec}.(var){depth}) == -1);
                for i = 1:length(idxOn)
                    on = idxOn(i);
                    off = idxOff(i);
                    if startsWith(var,'Rew') && off > on + 100
                        warning(['Problem with ' var ' at rest rec=' num2str(rec) ', depth=' num2str(depth)])
                    else
                        try
                            % check that acceleration within short/long window
                            % isn't above respective thresholds
                            if any(abs(acc(on-winNoMoveLong:off+winNoMoveLong)) > accThreshLong) || any(abs(acc(on:off+winNoMoveShort)) > accThreshShort)
                                data5.data{rec}.RewardRest{depth}(on:off) = 0;
                                data5.data{rec}.RewardLongRest{depth}(on:off) = 0;
                                data5.data{rec}.RewardShortRest{depth}(on:off) = 0;
                            end
                        catch % if rew/puff etc too close to edge of recording
                            data5.data{rec}.RewardRest{depth}(on:off) = 0;
                            data5.data{rec}.RewardLongRest{depth}(on:off) = 0;
                            data5.data{rec}.RewardShortRest{depth}(on:off) = 0;
                        end
                    end
                end
            end
        end
    end
end

%% Select conditional rewards and air puffs %don't need
for xx = [] 
for rec = 1:size(data5,1)
    try
        if strcmp(data5.RunRew{rec},{'rew'}) % recordings with rewards only
            for depth = 1:size(data5.data{rec},1)
                light = data5.data{rec}.Light{depth};
                rew = data5.data{rec}.Reward{depth};
                puff = data5.data{rec}.AirPuff{depth};

                omission = zeros(size(light));
                condRew = omission;
                condPuff = omission;
                
                if any(light) % check that each recording actually has lights
                    lightOn = find(diff(light) == 1);
                    lightOff = find(diff(light) == -1);
                    for i = 1:length(lightOn)
                        onIdx = lightOn(i);
                        offIdx = lightOff(find(lightOff>onIdx,1));
                        if rew(offIdx+1) == 1
                            condRew(onIdx) = 1;
                            condRew(offIdx) = -1;
                        elseif puff(offIdx+1) == 1
                            condPuff(onIdx) = 1;
                            condPuff(offIdx) = -1;
                        else
                            omission(onIdx) = 1;
                            omission(offIdx) = -1;
                        end
                    end
                    data5.data{rec}.Omission{depth} = [0;cumsum(omission)];
                    data5.data{rec}.CondRew{depth} = [0;cumsum(condRew)];
                    data5.data{rec}.CondPuff{depth} = [0;cumsum(condPuff)];
                end
            end
        end
    end
end
end

%% Select rest vs movement periods (MovOnOff)
% this selects strict rest periods by excluding any movements, even the 
% shorterst, but also short rest periods (<0.5s) within movement bouts, as
% these are not real rests
thresh1 = 0.0237; %thresholds are so strange because they were determined in volts then converted to m/s
thresh2 = 0.0101;
timeThresh = 50;  %0.5s

for rec = 1:size(data5,1)
     try
         for depth = 1:size(data5.data{rec},1)
             rawMov = data5.data{rec}.chMov{depth};
             % double thresh in both negative and positive directions to
             % find any movements even if small
             prestateRunRest = double_thresh_updown(rawMov,thresh1,thresh2);
             % remove 'movements' 1 frame long with no other near movements
             % (2 frames around) = artifacts
             stateRunRest = prestateRunRest;
             for j = 3:length(prestateRunRest)-2
                 if prestateRunRest(j) == 1 && sum(prestateRunRest(j-2:j+2)) == 1
                     stateRunRest(j) = 0;
                 end
             end
             % select true rest periods (longer than time threshold 0.5s)
             stateRunRestThresh = ones(size(stateRunRest)); %start assuming everthing is movement
             firstRun = find(stateRunRest == 1,1);
             if ~isempty(firstRun) %in case there are no movement periods at all
                 if firstRun ~= 1
                     firstRun = firstRun-1;
                     stateRunRestThresh(1:firstRun) = 0; % all time before fist run is rest
                 end
                 lastRun = find(stateRunRest == 1,1,'last')+1;
                 stateRunRestThresh(lastRun:end) = 0;
                 lastRun = min(lastRun, length(stateRunRest)-timeThresh); % to avoid issues with test windows                 
                 % keep rest periods only if longer than 0.5s
                 for timeBin = firstRun:lastRun
                     if stateRunRest(timeBin) == 0 %already rest in previous step
                         if mean(stateRunRest(timeBin:timeBin+timeThresh)) == 0 % no movement within window after
                             stateRunRestThresh(timeBin:timeBin+timeThresh) = zeros(1,timeThresh+1);
                         end
                     end
                 end
             else
                 stateRunRestThresh = zeros(size(stateRunRest));
             end
             data5.data{rec}.MovOnOff{depth} = stateRunRestThresh;
         end
     catch
         disp(['problem in run/rest selection for ' num2str(rec)])
     end
end


%% Select strict running (Run)
% Out of all non-rest periods as selected about, only consider 'Run' periods
% those bouts with meanVel > 0.2 m/s and duration of more than 0.5 s
velThresh = 0.2;
durThresh = 50; %0.5m/s
for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec})
        movOnOff = data5.data{rec}.MovOnOff{depth};
        data5.data{rec}.Run{depth} = sparse(zeros(size(movOnOff)));
        onIdx = find(diff(movOnOff) == 1);
        offIdx = [find(diff(movOnOff) == -1) length(movOnOff)];
        if ~isempty(onIdx)
            for o = 1:length(onIdx)
                % find onset and offset for each movement bout
                on = onIdx(o);
                off = offIdx(find(offIdx > on,1));
                dur = off-on;
                meanVel = mean(data5.data{rec}.chMov{depth}(on:off),'omitnan');
                if meanVel > velThresh && dur > durThresh % check that mean velocity and duration are abouve thresh
                    data5.data{rec}.Run{depth}(on+10) = 1;
                    data5.data{rec}.Run{depth}(off-10) = -1;
                end
            end
        end
        data5.data{rec}.Run{depth} = cumsum(data5.data{rec}.Run{depth}); %from on/off to binary
    end
end
disp('Finished strict run selection')

%% Exclude reward times (and other stimuli) from Run times
% exclude 5s window after reward, aipuff or light stimulus delivery to
% avoid confounds. 
% New variable 'RunRew' is all run periods without exclusion, 'Run' is now
% only non-rew/puff/light running. 
for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec})
        data5.data{rec}.RunRew{depth} = data5.data{rec}.Run{depth};
        data5.data{rec} = excludeRew(data5.data{rec},depth,'Run'); %separate function, see below
    end
end

%% get DF/F transient locations and calculate Signal-to-Noise ratio
% Find transient peaks, defined as transients with a rise of at least 0.03
% DF/F followed immediatelly by a decay of at least -0.005 DF/F
% Only calculate signal-to-noise ratio if there's a minimun number of
% transients per second, otherwise it's most probably a bad recording.
riseThresh = 0.03;
decayThresh = -0.005;
minTransFreq = 0.0222; %min transients per bin (200 trans in 1 min 30 s)
topTransPecent = 0.2; % top 20 percentile of transients for signal-to-noise

data5.sig2noise = cell(size(data5,1),1);
vars = {'chGreen','chRed'};
for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec})
        for v = 1:2
            var = vars{v};
            try
                % get data and de-normalize
                dff = data5.data{rec}.(var){depth};
                dff = dff * data5.norm{rec}(depth,v); %de-normalize
                
                % remove nan (without changing length!)
                bad = find(isnan(dff));
                try
                    dff(bad) = dff(bad+1);
                catch
                    dff(bad) = dff([bad(1:end-1)+1;bad(end)-1]);
                end
                
                % get NOISE - calculated by smoothing trace and subtracting
                % this from the unsmooth trace. What remain is the noise
                dff_Smooth = fastsmooth(dff,10);
                noise = std(dff - dff_Smooth,'omitnan');
                
                % get TRANSIENTS
                % first use high-pass filter to remove baseline fluctuations
                % (keep only sharp transients)
                dff_noFluct = removeFluctuations(dff); % separate function, see below
                % smooth trace
                dff_noFluct_smooth = fastsmooth(dff_noFluct,20);
                % get integral and smooth it too (to get chances in DF/F)
                int = diff(dff_noFluct_smooth);
                int_Smooth = fastsmooth(int,20);
                % get zero crossings of integral to separate rises and
                % decays in transients
                int_zeroX = ZeroX(int_Smooth); % separate fucntion - output is column
                
                % Get transient locations and characteristics:
                % for each point the integral crosses zero, if the next
                % period is a rise (abs max is positive), get the location
                % of the rise start, the peak of the transient (end of the
                % rise) value and location, the amplitude of the rise, and
                % the amplitude of the transient decay
                trans = nan(length(int_zeroX)-2,5); % 1=location, 2=peakLoc, 3=rise peak, 4=rise delta, 5=decay delta
                for i = 1:length(int_zeroX)-2
                    if max(int_Smooth(int_zeroX(i):int_zeroX(i+1))) > abs(min(int_Smooth(int_zeroX(i):int_zeroX(i+1)))) %rise
                        trans(i,1) = int_zeroX(i);
                        trans(i,2) = int_zeroX(i+1);
                        trans(i,3) = int_Smooth(int_zeroX(i+1));
                        trans(i,4) = dff_noFluct_smooth(int_zeroX(i+1)) - dff_noFluct_smooth(int_zeroX(i));
                        trans(i,5) = dff_noFluct_smooth(int_zeroX(i+2)) - dff_noFluct_smooth(int_zeroX(i+1));
                    end
                end
                trans(any(isnan(trans),2),:) = [];
                % keep only trans above thresh for rise and decay
                goodTransIdx = trans(:,4) > riseThresh & trans(:,5) < decayThresh;
                goodTrans = trans(goodTransIdx,:);
                peakLoc = goodTrans(:,2);
                % get transient amplitude size - from high-pass filtered signal! IMPORTANT!
                transAmplitude = zeros(size(dff));
                transAmplitude(peakLoc) = dff_noFluct(peakLoc);
                
                % get SIGNAL to calculate signal-to-noise ratio
                % make this 0 if number of trans below thresh. 
                transSort = sort(transAmplitude(peakLoc));
                if length(transSort)/length(dff) < minTransFreq
                    sig = 0;
                else
                    sig = transSort(round(length(transSort)*(1-topTransPecent)));
                end
                sig2noise = sig/noise; %calculate signal to noise ratio
                
                % transRun: keep only transients ocurring during Run periods
                transRun = transAmplitude;
                transRun(~data5.data{rec}.Run{depth}) = 0;
                
                % add to data (peaks,peaksRun and sig2noise)
                RG = var(3);
                data5.data{rec}.(['peaks' RG]){depth} = transAmplitude;
                data5.data{rec}.(['peaks' RG 'Run']){depth} = transRun;
                data5.sig2noise{rec}(depth,v) = sig2noise;
                
            catch %ch does not exist
                if v == 1
                    disp('whaaat?? no chGreen?')
                end
                data5.sig2noise{rec}(depth,v) = nan;
            end
        end
    end
end


%% Select accelerations and decelerations above thresh
accMinThresh = 2; %m/s2 (or negative for dec)
accMaxThresh = inf; %this is in case only acc/dec within a certain range want to be selected
for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec},1)
        acceleration = data5.data{rec}.Acceleration{depth};
        [acc,dec] = selectAccDec(acceleration,accMinThresh,accMaxThresh); %separate function below
        data5.data{rec}.AccOn{depth} = acc;
        data5.data{rec}.DecOn{depth} = dec;
        
        % exclude reward
        data5.data{rec} = excludeRew(data5.data{rec},depth,'AccOn');
        data5.data{rec} = excludeRew(data5.data{rec},depth,'DecOn');
    end
end


%% Select movement onsets and offsets
% criteria:
durThresh = 300; %run period at least 3s
velThresh = 0.4; %reach at least 0.4 m/s  win 0.75s window (below) from start or end of running
velWin = 75;
accMinThresh = 1; %first acc in onset at least 1 m/s2 (or last dec in offset negative this thresh)
edgeMinVel = -0.05;% no velocity below -0.05 m/s right when onset starts or when offset ends

durationThresh = 5; % min duration for acceleration/deceleration = 5 bins = 0.05s
peakThresh = 0.002; % min acc/dec peak (m/s2)

for rec = 1:size(data5,1)
    for depth = 1:size(data5.data{rec})
        % get vel
        vel = data5.data{rec}.chMov{depth};
        % get coarse movement onsets/offsets
        movOn = find(diff(data5.data{rec}.MovOnOff{depth}) == 1);
        movOff = find(diff(data5.data{rec}.MovOnOff{depth}) == -1);

        % get accelerations and decelerations
        acc = data5.data{rec}.Acceleration{depth};
        accSmooth = fastsmooth(acc,5);
        accZero = ZeroX(accSmooth); %get zero crossings - output is column
        accOnOff = nan(size(accZero,1)-1,4); %1=accOnIdx, 2=accOffIdx, 3=peakLoc, 4=peakValue
        decOnOff = nan(size(accZero,1)-1,4); %same for dec
        for z = 1:length(accZero)-1
            dur = accZero(z+1) - accZero(z);
            [~,peakLoc] = max(abs(accSmooth(accZero(z):accZero(z+1))));
            peak = accSmooth(peakLoc + accZero(z)-1);
            if dur > durationThresh && peak > peakThresh
                accOnOff(z,1) = accZero(z);
                accOnOff(z,2) = accZero(z+1);
                accOnOff(z,3) = peakLoc + accZero(z)-1;
                accOnOff(z,4) = peak;
            elseif dur > durationThresh && peak < -peakThresh
                decOnOff(z,1) = accZero(z);
                decOnOff(z,2) = accZero(z+1);
                decOnOff(z,3) = peakLoc + accZero(z)-1;
                decOnOff(z,4) = peak;
            end
        end
        accOnOff(any(isnan(accOnOff),2),:) = [];    
        decOnOff(any(isnan(decOnOff),2),:) = [];
        
        % check each automatic *ONSET* for criteria
        movOnFinal = zeros(size(vel));
        for onidx = 1:length(movOn)
            on = movOn(onidx);
            off = movOff(find(movOff > on,1)); % find offset to get duration below
            if isempty(off)
                dur2 = length(vel)-on;
            else
                dur2 = off - on;
            end
            if dur2 > durThresh % minimum bout length
                maxVelWin = max(vel(on:on+velWin)); % max velocity within window from onset
                if maxVelWin > velThresh % reach min velocity within window
                    firstAccOnIdx0 = find(accOnOff(:,2) > on ,1); %find first acc in movement bout
                    firstAccOnIdx = accOnOff(firstAccOnIdx0,1); % START of first acc
                    firstAccPeak = accOnOff(firstAccOnIdx0,4);
                    if firstAccPeak > accMinThresh % first acc min peak
                        velPre = vel(firstAccOnIdx); % get velocity at start of first acc
                        if velPre > edgeMinVel % first acc cannot ocurr at negative velocity
                            movOnFinal(firstAccOnIdx) = 1; % all criteria fulfilled!!
                        end
                    end
                end
            end
        end
        
        % check each automatic *OFFSET* for criteria
        movOffFinal = zeros(size(vel));
        for offidx = 1:length(movOff)
            off = movOff(offidx);
            on = movOff(find(movOff < off,1,'last')); % find prvious onset to get duration below
            if isempty(on)
                dur3 = off;
            else
                dur3 = off - on;
            end
            if dur3 > durThresh % minimum bout length
                maxVelWin = max(vel(off-velWin:off)); % max velocity within window before offset
                if maxVelWin > velThresh % reach min velocity within window
                    lastDecOffIdx0 = find(decOnOff(:,1) < off ,1,'last'); %find last dec in mov bout
                    lastDecOffIdx = decOnOff(lastDecOffIdx0,2); % END of last dec
                    lastDecPeak = decOnOff(lastDecOffIdx0,4);
                    if lastDecPeak < -accMinThresh % last dec min peak
                        velPost = vel(lastDecOffIdx); % get velocity at end of last dec
                        if velPost > edgeMinVel % last dec cannot end at negative velocity
                            movOffFinal(lastDecOffIdx) = 1; % all criteria fulfilled!!
                        end
                    end
                end
            end
        end
        
        % add MonOn and MovOff to data
        data5.data{rec}.MovOn{depth} = movOnFinal;
        data5.data{rec}.MovOff{depth} = movOffFinal;
    end
end


%% add data to data5 (if new) and save
if exist('newData','var')
    try
        temp = data5;
        load([dataProcessingFolder '\data5.mat'], 'data5'); 
        data5 = [data5;temp]; %add 
        data5 = table2cell(data5(:,1:2));
        for i = 1:size(data5,1)
            for j = 1:2
                try
                    data5{i,j} = data5{i,j}{1}; % sometimes saved as cell inside cell and it's annoying
                end
            end
        end
    catch
    end
end
data5 = sortrows(data5,1); % sort recordings
save([dataProcessingFolder '\data5.mat'], 'data5','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 
disp('saved data5')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTION to exclude reward (and other stimuli) times
% exclude all times at which a reward, air puff, or light stimulus was
% delivered, as well as a period of 5 s after the delivery of that stimulus
function data = excludeRew(data,depth,varName)
    win = 500; %5s
    if any(strcmp(data.Properties.VariableNames,'Reward'))
        rew = find(data.Reward{depth} | data.AirPuff{depth} | data.Light{depth});
        if ~isempty(rew)
            for i = 1:length(rew)
                try
                    data.(varName){depth}(rew(i):rew(i)+win) = 0;
                catch
                    data.(varName){depth}(rew(i):end) = 0; %if stimulus too close to end
                end
            end
        end
    end
end


%% FUNCTION to remove baseline fluctuations from DF/F to find transients
% subtract 8th percentile of small sliding window (high-pass filter)
function dff_noFluc = removeFluctuations(dff_raw) 
per = 0.08;
win = 250; %2.5s = over twice as large as CGaMP transients 

% for each point, get all values in win around, sort them and get base
len = length(dff_raw);
idx0 = 1:len;
idxWin = repmat(idx0,win+1,1) + repmat(-floor(win/2):floor(win/2),length(idx0),1)';
empty = idxWin < 1 | idxWin > len; %edges
idxWin(empty) = 1;
valsInWin = dff_raw(idxWin);
valsInWin(empty) = nan;

% get baseline as 8th percentile
sortedValsInWin = sort(valsInWin,1);
idxBaseMN = round(sum(~isnan(sortedValsInWin))'*per);
idxBase = sub2ind(size(sortedValsInWin),idxBaseMN,idx0(:));
base = sortedValsInWin(idxBase);

% subtract base
dff_noFluc = dff_raw(:) - base;
end


%% FUNCTION to select accelerations and decelerations within a range (accMin-accMax)
% For acc/dec greater than a min value, accMax = inf or omit variable
% Select acc/dec within a range that don't have any other large (> +-
% 2ms/2) accelerations or decelerations within a 0.25s window on either
% side.
function [acc2,dec2,peaks] = selectAccDec(acceleration,accMin,accMax)
    durationThresh = 5; % min duration for acceleration/deceleration = 5 bins = 0.05s
    threshRemove = 2; % acc threshold for nearby acc/dec - no large acc/dec nearby
    excludeWin = 50; %0.5s - 0.25 on each side

    if ~exist('accMax','var')
        accMax = inf;
    end
    accOn = ZeroX(acceleration); %get zero crossings of acceleration trace - output is column
    accOn(1) = []; %otherwise first is 0
    dur = [diff(accOn)' nan]; %acc durations

    peakValue = nan(1,length(accOn));
    peaks = nan(1,length(accOn));
    acc = false(size(acceleration));
    dec = false(size(acceleration));
    accDecRemove = false(size(acceleration));

    for i = 1:length(accOn)-1
        if dur(i) > durationThresh % min duration
            % get acc/dec peak
            [~,peakIdx] = max(abs(acceleration(accOn(i):accOn(i+1))));
            peakValue(i) = acceleration(accOn(i)+peakIdx-1);
            % check if acc/dec within range accMin-accMax
            if peakValue(i) > accMin && peakValue(i) < accMax
                acc(accOn(i)) = true;
                peaks(accOn(i)) = peakValue(i);
            elseif peakValue(i) < -accMin && peakValue(i) > -accMax
                dec(accOn(i)) = true;
                peaks(accOn(i)) = peakValue(i);
            end
            % save accelerations above threshRemove to exclude acc/dec too
            % near to other large acc/dec below
            if peakValue(i) > threshRemove || peakValue(i) < -threshRemove
                accDecRemove(accOn(i)) = true; 
            end
        end
    end

    % remove any acc/dec that are too close to each other - to exclude zigzag running
    allAccDec = acc | dec;
    allAboveRemoveThresh0 = [false(ceil(excludeWin/2),1);accDecRemove(:);false(floor(excludeWin/2),1)];
    allAboveRemoveThresh = nan(size(allAboveRemoveThresh0,1),excludeWin+1);
    for i = -ceil(excludeWin/2):floor(excludeWin/2)
        if i ~= 0
            % exclude acc/dec with either acc or dec next to it
            allAboveRemoveThresh(:,i+ceil(excludeWin/2)+1) = circshift(allAboveRemoveThresh0,i);
        end
    end
    allAboveRemoveThresh(:,ceil(excludeWin/2)+1) = [];
    allRemove = logical(sum(allAboveRemoveThresh,2) > 0);
    allRemove = allRemove(ceil(excludeWin/2)+1:end-floor(excludeWin/2));
    allAccDec(allRemove) = false;
    
    acc2 = acc(:) & allAccDec(:);
    dec2 = dec(:) & allAccDec(:);  
end





