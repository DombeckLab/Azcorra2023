function idx = findMouse(data,mouse)

%% check mouse imput is correct
k = strfind(mouse,'-');
% if ~((length(k) == 3 && length(mouse(k(1)+1:end)) == 18) || (length(k) == 2 && length(mouse(k(1)+1:end)) == 13))
%     error('Mouse input is not correct')
% end    
    
%% get number of recordings
count = 0;
multiDepths = false;
for e = 1:size(data,1)
    count = count + size(data.data{e},1);
    if size(data.data{e},1) > 1
        multiDepths = true;
        n = 0;
    end
end

%% get name of each recording
idxList = cell(count,2);
for e = 1:size(data,1)
    if ~multiDepths
        idxList{e,1} = [e,1];
        idxList{e,2} = data.data{e}.Properties.RowNames{1};
    else
        for d = 1:size(data.data{e},1)
            n = n+1;
            idxList{n,1} = [e d];
            idxList{n,2} = data.data{e}.Properties.RowNames{d};
        end
    end
end

%% normally don't need to do - check for correctness
for x = []
    % check that all recording names match info on data
    for i = 1:size(idxList,1)
        rec = idxList{i,2};
        e = idxList{i,1}(1);
        k = strfind(rec,'-');
        if ~strcmp(rec(1:k(1)-1),data.Exp{e}) || ~strcmp(rec(k(1)+1:k(2)-1),data.Mouse{e}) || ~strcmp(rec(k(2)+1:k(3)-1),data.Date{e})
            warning(['recording ' num2str(idxList{i,1}) ' does not match'])
        end
    end

    % check that all recording names are unique
    if size(idxList,1) ~= length(unique(idxList(:,2)))
        error('Some recordings are duplicated!')
    end
end

%% to find recording DAY in multidepth not caring about individual recordings
if length(k) == 2
    for i = 1:size(idxList,1)
        idxList{i,2} = idxList{i,2}(1:end-5);
    end
    [~,b] = unique(idxList(:,2));
    idxList = idxList(b,:);
end
    

%% find mouse
idx1 = ismember(idxList(:,2),mouse);
if sum(idx1) == 0 
    error('Mouse not found')
end
idx = nan(sum(idx1),2);
idx2 = find(idx1);
for i = 1:sum(idx1)
    idx(i,:) = idxList{idx2(i),1};
end

    

