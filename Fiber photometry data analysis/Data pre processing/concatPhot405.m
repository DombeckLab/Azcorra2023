function T = concatPhot405(firstFolder,depths)
% This script takes raw fluorescence and behavioral data collected in 
% PicoScope at 4 kHz, concatenates the raw data for each recording and 
% separates the fluorescence from 470 nm vs 405 nm (isosbestic control), 
% illumination, saving it as a binned file in both the original folder and 
% a new backup folder

% inputs - optional, if not set the script will ask for them:
    % firstFolder = directory for first recording in session set
    % depths = vector of depths for fiber 1 recording locations (chG), one
        % depth per recording in set. This script was originally designed 
        % for two fibers located at the same depth, so only one depth is 
        % asigned to each recording at this stage. The second set of depths
        % will be linked in 'data_processing.m'
% outputs:
    % T = table. Automatically saved with below format in RawDataFolder 
        % and MainFolder (edit MatlabFolders.mat file in Desktop): 
        % 'Binned405_(experiment)-(mouseID)-(recording date yyyymmdd).mat'  
    
% This scripts was written for Windows at MATLAB R2021a

% The variables in Picoscope are as follows (can be edited):
    % A = velocity (chMov)
    % B = fiber 2 fluorescence (from here on called channel red - chR)
    % C = fiber 1 fluorescence (channel green - chG)
    % D = light stimulus trigger
    % E = waveforem generator output indicating illumination wavelength (1 = 470nm, 0 = 405nm)
    % F = reward delivery trigger
    % G = licking sensor ourput
    % H = air puff delivery trigger
    
% Picoscope data must be saved in MATLAB format and organized in the
% following way:
    % Each recording in a set from a single mouse and session should be
        % named as 'date-recording number', where recording number = 
        % number within the recording set and is 4 digits long: 
        % yyyymmdd-00##. this is Picoscope's default naming system
    % Each set of recording folders must be inside a folder named as
        % 'experiment type-mouse ID' - example 'Anxa-A005'. All recording
        % sets from the same mouse should be in this same folder
    % All mouse folders can be in any folder desired, but must ultimatedly
        % be inside a RAW DATA FOLDER for automatic saving.
% The RAW DATA FOLDER will be specified in a file called 'MatlabFolders.mat' 
% that should be saved in the Desktop. Additionally, this file should
% include a MAIN FOLDER where the processed data will be automatically
% saved in addition to in the RAW DATA FOLDER, as well as a DATA PROCESSING
% FOLDER where the final master dataset will be saved at.


%% set variables
binsToIgnore = 5; %how many timepoints after 405-470 boundaries to ignore (including boundary) - default = 5 at 4 kHz
rawDataFreq = 4000; % data recorded on picoscope at 4 kHz
finalDataFreq = 100; % binned data output at 100 hz
binsize = rawDataFreq/finalDataFreq; %Number of bins to be averaged together

if binsize ~= 40
    warning(['Careful!! binning set at ' num2str(binsize) '.']) %default
end

%% set RAW DATA FOLDER and MAIN FOLDER for backup
%(edit MatlabFolders.mat file in code folder). The output table T will be 
% saved in the original raw data folder and a new folder within MainFolder
filePath = matlab.desktop.editor.getActiveFilename;
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
load([filePath 'MatlabFolders.mat'], 'rawDataFolder', 'mainFolder');

currentFolder = pwd; % to return to at the end

%% Get list of folders for each recording in set
% data should be saved as .mat files from PicoScope in the following format
% MasterFolder > Experiment-mouseID > YYYYMMDD-000#
% inside there are many .mat folders for that recording

% select folder if not set
if isempty(firstFolder)
    firstFolder = uigetdir('','Select file'); % select first folder of the recording set (0001)
end

% get mouse and date
k = strfind(firstFolder,'\');
mouse = firstFolder(k(length(k)-1)+1:k(length(k))-1);
k = k(length(k));
folder = firstFolder(1:k);
cd(folder)
if ~strcmp(firstFolder(k+1:k+2),'20')
    error('File selected is not the right type')
end
date = firstFolder(k+1:k+8);

% get all recordings in set (list of files for the same date experiment)
fileList = dir([folder date '-*' '']);
dirFlags = [fileList.isdir];
fileList = fileList(dirFlags);
files = cell(2,size(fileList,1));
for i = 1:size(fileList,1)
    files{1,i} = fileList(i).name;
end

numFiles = size(fileList,1);
disp(['Mouse: ' mouse '  ExpDate: ' date '   NumFiles: ' num2str(numFiles)])


%% Get depths
if ~exist('depths','var') % if not set, input manually - must match number of recordings in folder with same date
    depths = [];
    while length(depths)~=numFiles
        if ~isempty(depths)
            warning(['Number of depths (' num2str(length(depths)) ')does not match number of files (' num2str(numFiles) ')']);
        end
        depths = input('Depths (as vector):  ');
    end
else
    depths = cell2mat(depths);
end
files(2,:) = num2cell(depths);

clear firstFolder k i depths date2 fileList


%% Concatenate files and put into table T

% link picoscope vars (A-H) with associated variables
vars = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
%           1:A     2:B      3:C      4:D     5:E     6:F       7:G       8:H
varNames = {'chMov' 'chRed' 'chGreen' 'Light' 'ch470' 'Reward' 'Licking' 'AirPuff'}; %original vars list

% make empty table
varNamesNew = {'chMov' 'chRed' 'chGreen' 'Light' 'Reward' 'Licking' 'AirPuff' 'chRed405' 'chGreen405'}; %final vars list
empty = cell(numFiles,1);
T = cell2table(repmat(empty,1,length(varNamesNew)),'VariableNames',varNamesNew); %T table with all variables

flist = files(1,:);
fdepths = files(2,:);
w = waitbar(0,'Concatenating data...');
for t = 1:numFiles
    missingVar = 0;
    clear A B C D E F G H
    waitbar(t/numFiles,w,'Concatenating data...');
    % open folder for recording t
    filename = [folder flist{t}];
    name = flist{t};
    cd(filename);
    % load first file in recording
    filelist = dir(fullfile(filename, '*.mat'));
    numSections = size(filelist,1);
    load(filelist(1).name);
    % check that all variables exist
    for v = 1:length(vars)
        if exist(vars{v},'var') == 0
            warning(['Variable does not exist: ' vars{v} ' (' varNames{v} ') - filenumber: ' num2str(t)])
            missingVar = 1;
        end
    end
    if missingVar == 1
        error(['Missing variables in file #' num2str(t) '- please check']); 
        %sometimes a variable isn't saved properly
    end
    
    fileSize = length(A); % Picoscope saves recordings in many 10s files 
    % containing 40007 or 40006 bins (should be 4000 = 10s*4kHz, but it 
    % isn't too precise).
    
    len = fileSize*numSections;
    unbinned = cell(len,length(vars));
    for tn = 1:numSections
        load(filelist(tn).name);
        H = [H;zeros(size(A,1)-size(H,1),1)];  %for some strange reason, 
        % the last file in each recording has a truncated H variable.
        
        % check that file lengths are correct
        try
            unbinned((tn-1)*fileSize+1:(tn)*fileSize,:) = num2cell([A B C D E F G H]);
        catch
            if tn ~= numSections
                error('something is wrong...')
            end
        end
    end
    clear A B C D E F G H
    
    
    %% remove infinities (from signals saturating the detector)
    unbinnedNoInf = cell2mat(unbinned);
    % just in case one channel is permanently saturated (error) - make it all 0
    if any(any(isfinite(unbinnedNoInf)))
        bad = find(~any(isfinite(unbinnedNoInf))); % bad variables
        for b = bad
            unbinnedNoInf(:,b) = zeros(size(unbinnedNoInf,1),1)+5;
        end
    end
    
    % for binnary variables (D,E,F,G,H) +Inf = 5 (max value) 
    tempDEFGH = unbinnedNoInf(:,4:8);
    tempDEFGH(tempDEFGH == Inf) = 5;
    unbinnedNoInf = [unbinnedNoInf(:,1:3) tempDEFGH];
    
    % for non binary variables (A,B,C), sustitute saturation for interpolated value
    temp = (1:size(unbinnedNoInf,1))';
    unbinnedNoInf = bsxfun(@(x,y) interp1(y(isfinite(x)),x(isfinite(x)),y),unbinnedNoInf,temp);
    
    
    
    %% separate fluorescence from 405nm and 470nm illumination
    % LEDs are alternated at 100Hz, as reported by variable E (output of
    % waveform generator). E (var5) ~= 0 for 405, E ~= 5 for 470.
    % use E to find boundaries, then exclude bins right after transition
    % (binsToIgnore var at top)
    
    ch405Idx = unbinnedNoInf(:,5) <= 0.5;
    ch470Idx = unbinnedNoInf(:,5) > 0.5;
    first405on = find(diff(ch405Idx) == 1,1)+1; %for next section
    
    % exclude bins after transitions
    boundaries = [abs(diff(ch405Idx));false(binsToIgnore+1,1)];
    boundariesSum = false(size(boundaries,1),binsToIgnore);
    for i = 0:binsToIgnore-1
        boundariesSum(:,i+1) = circshift(boundaries,i);
    end
    boundaries = logical(sum(boundariesSum,2));
    boundaries(end-binsToIgnore+1:end) = [];
    ch405Idx(boundaries) = false; % indexes for 405 excluding the boundary times
    ch470Idx(boundaries) = false;
    
    % add new chR and chG (470 in place of original and 405 at end)
    chR405 = unbinnedNoInf(:,2);
    chR405(~ch405Idx) = nan;
    chG405 = unbinnedNoInf(:,3);
    chG405(~ch405Idx) = nan;
    unbinnedNoInf = [unbinnedNoInf chR405 chG405];
    unbinnedNoInf(~ch470Idx,2) = nan;
    unbinnedNoInf(~ch470Idx,3) = nan;
    
    % remove ch470 (waveform generator) - no longer useful. Now variables 
    % correspond to varNamesNew 
    unbinnedNoInf(:,5) = [];
    
    %% Rebin all channels
    % remove data until first 405 start to make even bins
    unbinnedNoInf(1:first405on,:) = [];
    % rebin by averaging (omiting nan)
    binned = nan(floor(size(unbinnedNoInf,1)/binsize),size(unbinnedNoInf,2));
    for var = 1:size(unbinnedNoInf,2)
        temp = unbinnedNoInf(:,var);
        binned(:,var) = arrayfun(@(i) mean(temp(i:i+binsize-1),'omitnan'),1:binsize:length(temp)-binsize+1)';
    end
    
    
    %% Add to T table
    for var = 1:length(varNamesNew)
        T{t,var} = {binned(:,var)};
    end
end


%% Reorder variables and add depths and row names
orderedVars = {'chMov' 'chRed' 'chGreen' 'chRed405' 'chGreen405' 'Light' 'Reward' 'Licking' 'AirPuff'};
T = T(:,orderedVars);
T = [table(fdepths','VariableNames',{'Depth'}) T];
flist = strcat(mouse,'-',flist);
T.Properties.RowNames = flist;

%% save 
% save T as 'Binned405_(experiment)-(mouseID)-(recording date yyyymmdd).mat' 
% in both original rawDataFolder and mainFolder for automatic backup
% detect saving errors, and stop for manual save if neither save was successful 
rawPath = folder;
mainPath = [mainFolder '\' mouse];
err = 0;
name = [mouse '-' name(1:length(name)-5)];
if ~strcmp(folder,rawPath) &&~strcmp(folder,mainPath)
    save([folder,'Binned405_',name],'T');
    if strcmp(folder(1:length(rawDataFolder)), rawDataFolder)
        rawPath = [rawDataFolder '\' folder(length(rawDataFolder)+2:end)];
        mainPath = [mainFolder '\' folder(length(rawDataFolder)+2:end)];
    end
end
try
    save([mainPath,'Binned405_',name],'T');
catch
    try
        mkdir(mainPath)
        save([mainPath,'Binned405_',name],'T');
    catch
        err = 1;
    end
end
try
    save([rawPath,'Binned405_',name],'T');
catch
    if err == 1
        err = 3;
    elseif err == 0
        err = 2;
    end
    
end

cd(currentFolder);
if err == 1
    warning('COULD NOT SAVE DATA IN MAIN FOLDER')
elseif err == 2
    warning('COULD NOT SAVE DATA IN RAW DATA FOLDER')
elseif err == 3
    warning('COULD NOT SAVE DATA IN RAW DATA FOLDER NOR MAIN FOLDER - SAVE MANUALLY!') 
    dbstop in concatPhot405 at 306 % do not remove this breakpoint! allows manual saving to not loose data
end

close(w)
