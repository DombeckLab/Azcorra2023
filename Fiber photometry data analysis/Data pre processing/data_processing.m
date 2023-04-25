function data4 = data_processing(MiceList,newData)
% This script calculates DF/F from raw fluorescence traces by:
    % subtracting non-GCaMP fluorescence (autofluorescence from fibers and
        % brain tissue)
    % dividing fluorecence by baseline (F0) to get F/F0, the first step in
        % calculating DF/F
    % correcting for bleaching (by doing the two previous steps using a 
        % moving window, making sure this window is wide enough to not
        % remove fluctuations of interest (like transients)
    % manually cropping the start of each recording if necessary (due to an
        % edge effect of the wide moving window used for bleachig
        % correction the start of the recording is sometimes not flat and
        % needs to be excluded)
    % subtracting the baseline, step two to get DF/F - this baseline is
        % selected manually because different subtypes can have different
        % transient frequencies, which makes automatic baseline selection
        % systems biased
    % normalize transients from 0 to 1 - to allow better comparisons 
        % between recordings and subtypes    
    % convert velocity from Volts to m/s (using a callibrated conversion 
        % factor)
    % add the newly analyzed data to the data master, check there are no 
        % duplicate mice, and save
        
% For ADDITIONAL DETAILS see each section and coresponding functions below     

% The inputs MiceList and newData are not necessary - if they are not
% provided, the script with load the MiceList file from the
% dataProcessingFolder and process all the data in MiceList to generate
% data4. If newData is provided, it will only analyze the data listed
% there, load data4 from dataProcessingFolder and append the newly analyzed
% data to that existing data4.

% Output: data4 = master table of all analyzed recordings previously. It is
% automatically saved in dataProcessingFolder

% This scripts was written for Windows at MATLAB R2021a


%% set MAIN FOLDER and DATA PROCESSING FOLDER (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder','mainFolder');

%% set variables and load colors
load colors.mat C %color file for default plotting colors
hz = 100; %data sampling frequency


%% get data to process
if ~exist('MiceList','var') || isempty(MiceList) %if MiceList is not provided, load
    load([dataProcessingFolder '\MiceList.mat'], 'MiceList');
end
if isempty(newData)
    data = MiceList; %if newData not provided, process all data in MiceList
else
    data = MiceList(end-size(newData,1)+1:end,:); %get data to process from MiceList
    % Check that data matches newData matches - must be the last entries in MiceList
    recIDmiceList = table2cell(data(:,2:3));
    recIdnewData = newData(:,5:6);
    for i = 1:size(recIDmiceList,1)
        recIDmiceList_i = strcat(recIDmiceList{i,1},recIDmiceList{i,2});
        recIdnewData_i = strcat(recIdnewData{i,1},recIdnewData{i,2});
        if ~strcmp(recIDmiceList_i,recIdnewData_i)
            error('Problem with mice order')
        end
    end
end

%% find data (already concatenated and rebinned to 100hz - see concatPhotom405)
% This section:
    % loads the data for each recording (must be already concatenated and
        % rebinned to 100Hz - see concatPhotom405)
    % removes the end of each recording, which is 0 due to a Picoscope error
for rec = 1:size(data,1)
    try
        % find data
        file = dir([mainFolder '\BinnedData\**\' data.Exp{rec} '-' data.Mouse{rec}]);
        if isempty(file)
            file = dir([mainFolder '\BinnedData\**\' data.Exp{rec} '_' data.Mouse{rec}]);
        end 
        file = file(1).folder;
        day = data.Day(rec);
        fileList = dir([file '\*Binned*']);
        if size(fileList,1) >= day{1}
            % load data already concatenated and rebinned to 100hz - see concatPhotom405
            load([fileList(day{1}).folder '\' fileList(day{1}).name], 'T')
            
            % remove the end of the recording - it's 0 for some variables because
            % of a strange error in picoscope
            for depth = 1:size(T,1)
                remove = find(T.chGreen405{depth} == 0, 1);
                if ~isempty(remove)
                    for var = 2:size(T,2)
                        T.(var){depth} = T.(var){depth}(1:remove-2);
                    end
                end
            end            
        else % can't find the recording - data must not be concatendated and rebinned
            idx = strfind(file,'\');
            warning(['Binned data missing for mouse ' file(idx(end)+1:end) ' #' num2str(day)]);
            T = [];
        end
        data.data{rec} = T;
    catch
        warning(['error in file #' num2str(rec)])
    end
end
% reorder
data = movevars(data,'data','After',3);


%% remove cortex recordings 
% these were used to determine 85th percentile used in the cortexCorrect
% function, but are not useful otherwise
cortexThresh = 1.5; %mm
for rec = 1:size(data,1)
    depths = cell2mat(data.data{rec}.Depth);
    remove = depths < cortexThresh;
    data.data{rec}(remove,:) = []; 
end
   
%% Get DF/F STEP 1 - correct for autofluorescence, bleaning, and divide by baseline to get F/F0
    % subtracts non-GCaMP fluorescence (autofluorescence from fibers and
        % brain) - see cortexCorrect function at bottom for details
    % divides fluorecence by baseline (F0) to get F/F0, the first step in
        % calculating DF/F - see baselineCorrect function at bottom
    % corrects for bleaching (by doing the two previous steps using a 
        % moving window)
% for more details see the sub-functions cortexCorrect and baselineCorrect
% at the bottom

for rec = 1:size(data,1)
    try
        % find data
        file = dir([mainFolder '\**/' data.Exp{rec} '-' data.Mouse{rec}]);
        if isempty(file)
            file = dir([mainFolder '\**/' data.Exp{rec} '_' data.Mouse{rec}]);
        end
        file = file(1).folder;
        day = data.Day(rec);
        fileList = dir([file '\*Binned*']);
        if size(fileList,1) >= day{1}
            % load data already concatenated and rebinned to 100hz - see concatPhotom405
            load([fileList(day{1}).folder '\' fileList(day{1}).name], 'T')
            
            % remove the end of the recording - it's 0 for some variables because
            % of a strange error in picoscope
            for depth = 1:size(T,1)
                remove = find(T.chGreen405{depth} == 0, 1);
                if ~isempty(remove)
                    for var = 2:size(T,2)
                        T.(var){depth} = T.(var){depth}(1:remove-2);
                    end
                end
            end
            
            % subtract autofluorescence, correct for bleaching and get F/F0
            vars = {'chGreen','chGreen405','chRed','chRed405'};
            dataFF = T;
            for depth = 1:size(T,1)
                for v = 1:length(vars)
                    fluorRaw = T.(vars{v}){depth};
                    % correct for autofluorescence - see function at bottom for details
                    fluorNoAutofluor = cortexCorrect(fluorRaw);
                    % get DF/F step1: divide by baseline to get F/F0 - see function at bottom for details
                    fluorFF = baselineCorrect(fluorNoAutofluor);
                    dataFF.(vars{v}){depth} = fluorFF;
                end
            end
        else % can't find the recording - data must not be concatendated and rebinned
            idx = strfind(file,'\');
            warning(['Binned data missing for mouse ' file(idx(end)+1:end) ' #' num2str(day{1})]);
            dataFF = [];
        end
        data.data{rec} = dataFF;
    catch
        warning(['error in file #' num2str(rec)])
    end
    disp(['   rec ' num2str(rec) '/' num2str(size(data,1)) ' done'])
end

%save data in case next step fails
if ~isfolder([dataProcessingFolder '\Delete\'])
    mkdir([dataProcessingFolder '\Delete\'])
end
save([dataProcessingFolder '\Delete\dataTemp1.mat'], 'data','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later.
disp('Checkpoint 1 - DF/F step 1 completed and saved')


%% exclude bad recordings by looking at 405, and crop start based on both channels
% if 405 channel shows many movement artifacts (rare), exclude recording
% if the start of the recording is still decaying after the initial correction, crop it off
data2 = data;
data2.cropStart = cell(size(data2,1),1); %save the crop point
figure
for rec = 1:size(data2,1)
    bad = false(size(data2.data{rec},1),1); %recordings to exclude for 405 artifacts
    crop = data2.cropStart{rec};
    if isempty(crop)
        crop = ones(size(data2.data{rec},1),1);
    end
    
    for depth = 1:size(data2.data{rec},1)
        finish = false;
        while ~finish
            plotAll(data2,rec,depth,C,'405') %plot recording including 405 - see function at bottom
            subplot(2,1,1); ylim('auto')
            subplot(2,1,2); ylim('auto')
            title('Click right of plot = exclude recording; inside plot = crop START; left of plot = keep as is')
            
            subplot(2,1,1)
            xl = xlim();
            [inputX,~] = ginput(1);
            inputX = round(inputX*hz);
            if inputX > xl(2)*hz %click right of plot -> exclude recording
                disp(['Bad: e = ' num2str(rec) ', d = ' num2str(depth)])
                crop(depth) = 0;
                bad(depth) = true;
                finish = true;
            elseif inputX > 0 %click inside of plot -> crop start there and replot (in case more cropping is necessary)
                for var = 2:size(data2.data{rec},2)
                    data2.data{rec}.(var){depth} = data2.data{rec}.(var){depth}(inputX:end); %crop all variables
                end
                crop(depth) = crop(depth) + inputX - 1;
            else %click left of plot -> keep as is
                finish = true;
            end
        end
    end
    
    %exclude bad recordings
    data2.data{rec}(bad,:) = [];
    crop(bad) = [];
    data2.cropStart{rec} = crop;
end
close()

%save data in case next step fails
save([dataProcessingFolder '\Delete\dataTemp2.mat'], 'data2','-v7.3'); 
disp('Checkpoint 2 - cropping of start completed and saved')

%% Get DF/F - STEP 2
% Manually select baseline and subtract - go from F/F to DF/F
data3 = data2;
figure
for rec = 1:size(data3,1)
    try
        [DFF,base] = baseSelect(data3.data{rec},C); %manually select baseline and subtract
        %see baseSelect function at bottom
        data3.base{rec} = base; %save selected baseline
        data3.data{rec} = DFF; %save final DF/F 
    catch
        disp(['problem in DF/F step 2 for ' num2str(rec)])
    end
end
close

%save data in case next step fails
save([dataProcessingFolder '\Delete\dataTemp3.mat'], 'data3','-v7.3');
disp('Checkpoint 3 - DF/F step 2 completed and saved')


%% Normalize trace from 0 (baseline) to 1 (biggest transient)
% The goal of this step is to allow better comparisons between recordings
% from regions with higher or lower GCaMP expression, and especially
% between DA subtypes that might have different axonal densities or calcium
% dynamics. 
% This step can be skipped however and won't affect the differences between
% subtypes observed, it will just increase the variability.
% The 405 channel is normalized the same way as its corresponding 470
data4 = data3;
incThresh = 0.03; %max jump between DF/F bins - greater jumps are artifacts
for rec = 1:size(data4,1)
    data4.norm{rec} = nan(size(data4.data{rec},1),2); % save normalization in case it is needed for de-normalization
    for depth = 1:size(data4.data{rec},1)
        % Make sure the 'biggest transient' is not an artifact!: ignore if
        % max value is much greater than surounding points
        fibers = {'chGreen' 'chRed'};
        for f = 1:2
            if any(strcmp(data4.data{rec}.Properties.VariableNames,fibers{f}))
                ready = false;
                while ~ready
                    [maxPeak,idx] = max(data4.data{rec}.(fibers{f}){depth}); %get max value of trace
                    if data4.data{rec}.(fibers{f}){depth}(idx-1) > (maxPeak - incThresh) || data4.data{rec}.(fibers{f}){depth}(idx+1) > (maxPeak - incThresh) %check if surounding points are much lower: not a transient
                        ready = true; %all good!
                    else
                        data4.data{rec}.(fibers{f}){depth}(idx) = data4.data{rec}.(fibers{f}){depth}(idx-1); %clean artifact
                    end
                end
                data4.data{rec}.(fibers{f}){depth} = data4.data{rec}.(fibers{f}){depth} / maxPeak; %normalize by max peak
                data4.norm{rec}(depth,f) = maxPeak;
                % normalize 405 the same way
                data4.data{rec}.([fibers{f} '405']){depth} = data4.data{rec}.([fibers{f} '405']){depth} / maxPeak;
            end
        end
    end
end

%% convert velocity in m/s and get acceleration
% must callibrate for each treadmill and rotary encoder
convFactor = 0.6766; %ms-1/V - see Maite LabBook 2020-21 p63
for rec = 1:size(data4,1)
    try
        for depth = 1:size(data4.data{rec},1)
            % convert velocity to m/s
            data4.data{rec}.chMov{depth} = data4.data{rec}.chMov{depth}*convFactor;
            % get acceleration
            movS = smooth(data4.data{rec}.chMov{depth},6);
            data4.data{rec}.Acceleration{depth} = [0;diff(movS)*hz];    %multiply by hz to get acc in m/s2 and not m/(s*bin)
        end
    catch
        disp(['problem in vel/acc for ' num2str(rec)])
    end
end

%save data in case next step fails
save([dataProcessingFolder '\Delete\dataTemp4.mat'], 'data4','-v7.3');
disp('Checkpoint 4 - transient normalization and acceleration calculation completed and saved')


%% add newly analyzed data to existing data4 and save
data4_new = data4;
try % load existing data4
    load([dataProcessingFolder '\data4.mat'], 'data4'); 
catch
    clear data4
end
if ~exist('data4','var') %first time analysis -> create new data4
    data4 = data4_new; 
elseif exist('data4','var') && ~isempty(newData)
    data4 = [data4;data4_new]; % add new data to existing data4
elseif exist('data4','var') && isempty(newData)
    warning('You are about to OVERWRITE the existing saved data - are you sure?')
    overwrite = input("Input 'yes' to overwrite, any other input will save new as ERRORDATA.mat and quit: ");
    if strcmp(overwrite,'yes')
        data4 = data4_new; 
    else
        save([dataProcessingFolder '\ERRORDATA.mat'], 'data4','-v7.3');
    end
end

%sort rows by exp, mouse, and date for ease of finding data
recIdnewData = table2cell(data4(:,1:3));
for i = 1:size(recIdnewData,1)
    for j = 1:3
        try
            recIdnewData{i,j} = recIdnewData{i,j}{1}; %need to take out of cell for sorting
        catch
        end
    end
end
[~,idx] = sortrows(recIdnewData);
data4 = data4(idx,:);

%save final data4 
save([dataProcessingFolder '\data4.mat'], 'data4','-v7.3'); 
disp('Checkpoint 5 - SAVED FULL DATA 4!')
 

%% check for duplicate mice - just in case!
mice = data4.Mouse;
for rec = size(mice,1)-1:-1:1
    if strcmp(mice(rec),mice(rec+1)) == 1
        mice(rec) = [];
    end
end
[uni,idx] = unique(mice);
if length(uni) ~= length(mice)
    x = 1:length(mice);
    x(idx) = [];
    disp([num2str(length(x)) ' mice are repeated!'])
    for i = 1:length(x)
        disp(['     Mice number: ' mice{x(i)}])
    end
end

end

%%% FUNCTIONS %%%

%% Plot data to visualize - both fibers if available and 405 ch in black
function plotAll(data,e,d,C,x405)
dataToPlot = data.data{e};
fibers = {'chGreen','chRed'};
colors = [C.green{1}; C.red{1}];
clf
for f = 1:2
    subplot(2,1,f);
    hold on
    if any(strcmp(dataToPlot.Properties.VariableNames,fibers{f}))
        plot(x(dataToPlot.(fibers{f}){d},100),dataToPlot.(fibers{f}){d},'Color',colors(f,:));
        if exist('x405','var') && strcmp(x405,'405')
            plot(x(dataToPlot.([fibers{f},'405']){d},100),dataToPlot.([fibers{f},'405']){d},'k');
        end
        horizontalLine(1);
        ylim([0.9 1.3])
        title([num2str(e) '-' num2str(d)])
        xlim([1 length(dataToPlot.(fibers{f}){d})]/100)
    end
end
end


%% Correct for fiber autofluorescence and bleaching
% input: fluorIN = fluorescence trace to correct
% output: fluorOUT = corrected fluorescence trace

% In fiber photometry the fluorescence (F) measured is not all comming from 
% the fluorescent indicator (GCaMP); a large portion is autofluorescence of
% the fibers or other components of the setup, and autofluorescence of the
% brain tissue. To determine what percentage of the measured F is 'real',
% we made several recordings in cortex, where there is no GCaMP, and then
% moved to the striatum where there is abundant GCaMP. The F in cortex was
% on average 85% of the F in striatum, and thus we will subtract 85% of the
% base F in our recordings to get GCaMP dependent F only. 

% Furthermore, F bleaches throughout the recording (from both GCaMP and
% non-GCaMP). To correct for this, we subtact the 85% percentile explained
% above using a moving window. This window must be much wider that any
% fluctuations that we might be interested in - calcium transients are
% about 1s in width and running bouts are 5-10s, so a window of 20s will
% not filter these out. However, this wide window means that the start of
% each recording, whih bleaches more sharply, will often not be properly
% corrected for. We will cut it out in a subsequent step.

function fluorOUT = cortexCorrect(fluorIN) %subtract 85% of 8th percentile of large sliding window
nonGcampFpercent = 0.85; %percent of fluorescence estimated to come from non-GCaMP sources
basePercent = 0.08; %conservative percentile of the F assumed to be baseline
win = 2001; %20s moving window for bleaching correction

% for each time point, get all values in win around, sort them and get baseline value
% first make matrix where each colum contains each time point and points before and after within window
len = length(fluorIN);
idx0 = 1:len;
idx = repmat(idx0,win,1) + repmat(-floor(win/2):floor(win/2),length(idx0),1)';
empty = idx < 1 | idx > len; %edges
idx(empty) = 1;
vals = fluorIN(idx);
vals(empty) = nan;
% sort values in window around each time point = column
sortedFluor = sort(vals,1);
% get baseline value for each time point
idx2 = round(sum(~isnan(sortedFluor))'*basePercent);
idx3 = sub2ind(size(sortedFluor),idx2,idx0(:));
base = sortedFluor(idx3);

% subtract from fluorescence 85 % of base (non GCaMP fluorescence)
fluorOUT = fluorIN(:) - nonGcampFpercent*base(:);
end

%% DF/F step 1 - Divide Fluorescence by fluorescence basline to go from F to F/F0
% input: fluorIN = fluorescence trace to correct
% output: fluorOUT = corrected fluorescence trace

% to get DF/F from fluorescence, you first need to divide by the the
% baseline to get F/F0. To do this, we calculate the baseline as the 8th
% percentile of the fluorescence within a sliding window of 20s (much wider
% than any fluctuations of interest) and divide the fluorescence by this
% value.

function fluorOUT = baselineCorrect(fluorIN) %divide by 8th percentile of small sliding window (get DF/F step1)
basePercent = 0.08; %conservative percentile of the F assumed to be baseline
win = 2001; %20s moving window for bleaching correction

% for each point, get all values in win around, sort them and get base = per
len = length(fluorIN);
idx0 = 1:len;
idx = repmat(idx0,win,1) + repmat(-floor(win/2):floor(win/2),length(idx0),1)';
empty = idx < 1 | idx > len; %edges
idx(empty) = 1;
vals = fluorIN(idx);
vals(empty) = nan;
sortedVar = sort(vals,1);
idx2 = round(sum(~isnan(sortedVar))'*basePercent);  
idx3 = sub2ind(size(sortedVar),idx2,idx0(:));
base = sortedVar(idx3);

% divide fluorescence by base (get F/F0)
fluorOUT = fluorIN(:)./base(:);
end


%% get DF/F step 2 - subtract baseline to go from F/F to DF/F
% inputs:
    % dataIN = input data table (F/F)
    % C = color file
% outputs:
    % dataOUT = output data table (DF/F)
    % base = manually selected baseline values to save for reference
    
function [dataOUT,base] = baseSelect(dataIN,C)
base = nan(size(dataIN,1),2);
dataOUT = dataIN;
fibers = {'chGreen','chRed'};
colors = [C.green{1}; C.red{1}];

for d = 1:size(dataIN,1)
    % plot data
    clf
    for f = 1:2
        subplot(2,1,f)
        hold on
        if any(strcmp(dataIN.Properties.VariableNames,fibers{f}))
            plot(x(dataIN.(fibers{f}){d},100),dataIN.(fibers{f}){d},'Color',colors(f,:));
            xlim([1 length(dataIN.(fibers{f}){d})]/100)
        end
    end
    
    % select baseline manually for 470 channel
    % select the top ege of the baseline noise - for subtraction in the 
    % next step, the median of all data below this top edge will be used
    for f = 1:2
        if any(strcmp(dataIN.Properties.VariableNames,fibers{f}))
            subplot(2,1,1)
            title(['Select baseline ' fibers{f}])
            subplot(2,1,f)
            [~,manualBase] = ginput(1);
            base(d,f) = manualBase;
        end
    end    
    
    % subtract 470 baseline from 470 to get DF/F
    % also correct 405 channel by subtracting median of whole 405 recording (it is all noise=baseline)
    for f = 1:2
        if any(strcmp(dataIN.Properties.VariableNames,fibers{f}))
            baseMedian = median(dataIN.(fibers{f}){d}(dataIN.(fibers{f}){d} < base(d,f)),'omitnan');
            dataOUT.(fibers{f}){d} = dataIN.(fibers{f}){d} - baseMedian;
            dataOUT.([fibers{f} '405']){d} = dataIN.([fibers{f} '405']){d} - median(dataIN.([fibers{f} '405']){d},'omitnan');
        end
    end
end
end








