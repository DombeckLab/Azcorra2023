function data6 = addData_paper(newData)
% This script takes raw fluorescence and behavioral data collected in 
% PicoScope at 4 kHz, concatenates the raw data for each recording and adds
% missing metadata, separates the fluorescence from 470 nm vs 405 nm 
% (isosbestic control) illumination, and adds these new recordings to a
% master dataset with all recordings

% The input newData is not necessary and is only useful when a step in the
% script fails to avoid re-running the entire script. This variable will
% include the manually inputed metadata for each recording, and will be
% saved and replaced in the DATA PROCESSING FOLDER.

% This script will also request important metadata for each recording that
% is not automatically recorded, particularly:
    % the depths at which each recording was made (for each fiber)
    % the location of each fiber during recording: 'ds' = body of striatum, 
        %'ts' = tail of the striatum, 'snc' 'nac' or 'vta'
    % whether the mouse was only running ('run') or received stimuli
        % (reward, air puff, light... - 'rew')
    % the mouse sex ('f' = female, 'm' = male)
    % recording info - any notes you want to add for the record, but that
        % will not be used for any further analysis.
% If you are bulk analyzing multiple recording sessions, 

% Output: data6 = master table of all recordings previously analyzed. 
% It is automatically saved in dataProcessingFolder

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

%% set RAW DATA FOLDER as current directory (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
k = strfind(filePath, '\');
filePath = filePath(1:k(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder', 'rawDataFolder');
cd(rawDataFolder)

%% select new raw data to analyze, add metadata and save info in 'newData'
% as many recording sets as desired can be bulk analyzed. Each set will be
% saved as a row in the 'newData' variable.
% newData colums = folder, depths chG, depths chR, experiment type, 
    % mouse ID, date, session number per mouse, info, run/reward, chG 
    % location, chR location, gender
if ~exist('newData','var')
    newData = cell(20,12); 
    rec = 0;
    firstFolder = 1; %initialize
    while firstFolder ~= 0 %as many recording sets as desidered
        rec = rec+1;
        % get first folder of recording set. 
        % select 'cancel' to end recording selection
        firstFolder = uigetdir('','Select file'); % select first folder of the recording set (0001)
        if firstFolder ~= 0
            k0 = strfind(firstFolder,'\');
            mouse = firstFolder(k0(length(k0)-1)+1:k0(length(k0))-1);
            kx = strfind(mouse,'-');
            mouse2 = mouse(kx+1:end);
            k = k0(length(k0));
            folder = firstFolder(1:k);
            cd(folder)
            if ~strcmp(firstFolder(k+1:k+2),'20') %will stop working in the year 2100
                error('File selected is not the right type')
            end
            date = firstFolder(k+1:k+8);

            % get all recordings in set (list of files for the same date and mouse)
            fileList = dir([folder date '-*' '']);
            dirFlags = [fileList.isdir];
            fileList = fileList(dirFlags);
            files = cell(2,size(fileList,1));
            for i = 1:size(fileList,1)
                files{1,i} = fileList(i).name;
            end
            numFiles = size(fileList,1);
            disp(['Mouse: ' mouse '  ExpDate: ' date '   NumFiles: ' num2str(numFiles)])

            % manually input recording depths as a vector for each fiber 
            % must match number of recordings in folder with same date
            % if fiber was not used, input 'o'
            depthsG = [];
            while length(depthsG) ~= numFiles && ~strcmp(depthsG,'o')
                if ~isempty(depthsG)
                    warning(['Number of depths (' num2str(length(depthsG)) ') does not match number of files (' num2str(numFiles) ')']);
                end
                depthsG = input('Depths chG (vector):  ');
            end
            depthsR = [];
            while length(depthsR) ~= numFiles && ~strcmp(depthsR,'o')
                if ~isempty(depthsR)
                    warning(['Number of depths (' num2str(length(depthsR)) ') does not match number of files (' num2str(numFiles) ')']);
                end
                depthsR = input("Depth chR:  ");
            end

            % add recording info and depths to newData file 
            newData{rec,1} = firstFolder; % address of first recording in set
            newData{rec,2} = num2cell(depthsG); %depths fiber G
            newData{rec,3} = num2cell(depthsR); %deptgs fiber R
            idx = strfind(firstFolder(k0(end-1):k0(end)),'-');
            exp = firstFolder(k0(end-1)+1:k0(end-1)+idx-2);
            newData{rec,4} = exp; %experiment type
            newData{rec,5} = mouse2; %mouse ID
            newData{rec,6} = date; %recording date

            fileList2 = dir(folder);
            dirFlags = [fileList2.isdir];
            fileList2 = fileList2(dirFlags);
            idxGoodFolders = startsWith({fileList2(1:end).name},'20');
            fileList2 = cell2mat({fileList2(idxGoodFolders).name}');
            fileList2 = cellstr(fileList2(:,1:8));
            expList = unique(fileList2);
            dayNum = find(strcmp(expList,date));
            newData{rec,7} = dayNum; % session number for each mouse

            % get rest of metadata
            copyPrevious = input("Recording info? - 1 to copy previous:  "); % notes to describe experiment. 
            % For reference, will not be used. Input 1 to copy all following 
            % metadata from previous recording selected.
            if copyPrevious == 1
                newData(rec,8:end) = newData(rec-1,8:end);
            else
                newData{rec,8} = copyPrevious;
                newData{rec,9} = input("Run or rew:  "); %did the mouse get rewards in this session or not? (input either 'run' or 'rew')
                newData{rec,10} = input("What is is chG?:  "); %location of chG fiber in brain (i.e. 'snc', 'ds', 'ts')
                newData{rec,11} = input("What is is chR?:  "); %as previous for chR
                newData{rec,12} = input("Mouse gender? ('f' or 'm'):  "); %mouse gender
            end
            disp('   ')    
        end
    end
    newData(rec:end,:) = [];

    save([dataProcessingFolder '\newData.mat'], 'newData');
end
    
    

%% add new datasets to miceList to check for duplicate recordings
% this is a list of all recordings analyzed to date without the data
if isempty(newData)
    error('No data selected')
end
try
    load([dataProcessingFolder '\MiceList.mat'], 'MiceList');
    vars = MiceList.Properties.VariableNames;
catch
    vars = {'Exp','Mouse','Date','Day','Type','RunRew','chG','chR','Gen','DepthDiff'};
    MiceList = array2table(cell(1,size(vars,2)));
    MiceList.Properties.VariableNames = vars;
end
newMice = cell2table([newData(:,4:end) mat2cell(newData(:,3),ones(size(newData,1),1),1)]); %#ok<MMTC>
newMice.Properties.VariableNames = vars;
newMice.Day = num2cell(newMice.Day);

% check that new recordings don't already exist!
mouseDate = table2cell([MiceList(:,2:3); newMice(:,2:3)]);
mouseDate(any(cellfun(@isempty,mouseDate),2),:) = []; % remove empty rows for first time MiceList
mouseDate2 = cell(size(mouseDate,1),1);
for i = 1:size(mouseDate,1)
    mouseDate2{i} = strcat(mouseDate{i,1},mouseDate{i,2});
end

MiceList = [MiceList;newMice];
MiceList(cellfun(@isempty,MiceList{:,1}),:) = []; % remove empty rows for first time MiceList
if size(unique(mouseDate2),1) ~= size(MiceList,1)
    warning('Duplicate mice! Fix (temp & newData, or MiceList) or abort');
    dbstop in addData_paper at 183 % DO NOT REMOVE THIS BREAKPOINT!!!
    % Allows you to fix issue by deleting repeated mouse
    if size(unique(mouseDate2),1) ~= size(MiceList,1)
        error('Something is still wrong...')
    end
end

save([dataProcessingFolder '\MiceList.mat'], 'MiceList');

%% concatenate new data and separate 405 and 470 channels
bad = [];
for rec = 1:size(newData,1) % for each recording set selected
    k0 = strfind(newData{rec,1},'\');
    % if recording has already been concatenated, don't re-do
    if ~isfile([newData{rec,1}(1:k0(end)) 'Binned405_' newData{rec,4} '-' newData{rec,5} '-' newData{rec,6},'.mat'])
        try
            concatPhot405(newData{rec,1},newData{rec,2}); % SEPARATE SCRIPT
        catch
            warning(['error with current file - rec #' num2str(rec)]);
            dbstop in addData_paper at 202 % DO NOT REMOVE THIS BREAKPOINT!!! Allows for fixing
            disp('fixing...')
            % normally this is due to one recording not being saved
            % properly - you will know because the .mat files within the
            % folder with be smaller than in good recordings. 
            % Rename that recording folder changing the start 'bad-2023...'
            % and re-run the script after this one finishes. The newData
            % for the bad recording sets will be automatically saved - but
            % remember to remove the depths for the bad recordings!
            
            % save bad recordings newData for later re-analyzing
            bad = [bad; newData(rec,:)]; %#ok<AGROW>
            newData(rec,:) = [];
            idx = strcmp(MiceList.Mouse,newData{rec,5}) & strcmp(MiceList.Date,newData{rec,6});
            MiceList(idx,:) = []; %remove bad recording from MiceList
        end
    end
end 
if ~isempty(bad)
    save([dataProcessingFolder '\badNewData.mat'], 'bad');
    save([dataProcessingFolder '\newData.mat'], 'newData');
    save([dataProcessingFolder '\MiceList.mat'], 'MiceList');
end


%% pre-process data to get DF/F and save as data4
data4 = data_processing(MiceList,newData); % SEPARATE SCRIPT

%% select signals (accelerations, rewards...) and save as data5
data5 = selectSignals_paper(data4,newData); % SEPARATE SCRIPT

%% separate individual recordings, get exclusion criteria, and save as FINAL data 6
data6 = exclusionCriteria_allFigs(data5,newData); % SEPARATE SCRIPT









