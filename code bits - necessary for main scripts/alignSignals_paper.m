function [means, sems, signals, n] = alignSignals_paper(data, depth, trig, onoff, VA, neg, pos)
% necesary inputs: data, depth, trig

%% check that variable is valid
if ~any(strcmp(data.Properties.VariableNames,trig))
    error('Trigger variable does not exist in data')
end

%% set defaults
if isempty(onoff)
    onoff = 1;
elseif ~(onoff == 1 || onoff == -1)
    error('onoff input must be 1, -1, or ~')
end
if isempty(VA)
    VA = 'A';
elseif ~(strcmp(VA, 'V') || strcmp(VA, 'A'))
    error("VA input must be 'V', 'A', or ~")
end  
if isempty(neg)
    neg = 75;
end
if isempty(pos)
    pos = 150;
end


%% decide depths
% options: 
% - can be a number of depth (first, second) > then it must be in string format ('1', '2')
% - can be depth in mm > then it must be a double
% - can be a range of two depths > vector [x y]
if isa(depth,'char') 
    if length(depth) == 1
        if str2double(depth) > size(data,1)
            error('data does not have that many depths')
        end
        depths = false(size(data,1),1);
        depths(str2double(depth)) = true; 
    else
        error('depth number must be single digit string')
    end 
elseif isa(depth,'double') && length(depth) == 1
    depths = cell2mat(data.Depth) == depth; 
    if ~any(depths)
        error('depth is not present in data')
    end
elseif isa(depth,'double') && length(depth) == 2
    depths = cell2mat(data.Depth) >= depth(1) & cell2mat(data.Depth) <= depth(2);
end

data(~depths,:) = [];
    

%% Get trigger points and aligned signals
signals = cell(1,6); % 1=A/V/L, 2=G, 3=G405, 4=R, 5=R405, 6=L

count = 0;
for row = 1:size(data,1)
    % get trigger points: 2 types of trigger signals:
        % - bistable (mostly zero, ones periods for signal on) > reward, light, puff, movOnOff...
        % - peaks (mostly zero, isolated peaks of any value
    % for all, do diff and select trigger on > 0 (on) or < 0 (off)
    sig = data.(trig){row};
    sig(sig < 0) = 0;
    sig = sig(:);
    if onoff == 1
    	idx = find([0;diff(sig)] > 0);
    elseif onoff == -1
        idx = find([0;diff(sig)] < 0);
    end
    
    % exclude points too close to start or end
    idx(idx < neg) = [];
    idx(idx > length(data.(trig){row})-pos) = [];
    
    % get aligned signals
    if ~isempty(idx)
        for i = 1:length(idx)
            try
                count = count+1;
                if strcmp(VA,'V')
                    signals{count,1} = data.chMov{row}(idx(i)-neg+1 : idx(i)+pos);
                elseif strcmp(VA,'A')
                    signals{count,1} = data.Acceleration{row}(idx(i)-neg+1 : idx(i)+pos);
                end

                signals{count,2} = data.chGreen{row}(idx(i)-neg+1 : idx(i)+pos);
                signals{count,3} = data.chGreen405{row}(idx(i)-neg+1 : idx(i)+pos);

                if any(strcmp(data.Properties.VariableNames,'chRed'))
                    signals{count,4} = data.chRed{row}(idx(i)-neg+1 : idx(i)+pos);
                    signals{count,5} = data.chRed405{row}(idx(i)-neg+1 : idx(i)+pos);
                end
                if any(strcmp(data.Properties.VariableNames,'Licking'))
                    signals{count,6} = data.Licking{row}(idx(i)-neg+1 : idx(i)+pos);
                end
            end
        end
    end
end


%% get mean and sem
means = nan(size(signals,2),neg+pos);
sems = nan(size(signals,2),neg+pos);
for i = 1:size(signals,2)   
    % make row vectors
    try
        for n = 1:size(signals,1)
            signals{n,i} = reshape(signals{n,i},1,neg+pos);
        end
        
        % get mean/sem
        sig = cell2mat(signals(:,i));
        means(i,:) = mean(sig,1,'omitnan');
        sems(i,:) = std(sig,1,'omitnan')/sqrt(size(sig,1));
    end
end

n = size(signals,1);

    
    
    
    
