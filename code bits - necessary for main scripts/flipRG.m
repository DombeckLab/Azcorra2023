function data = flipRG(data,recs)

for r = 1:length(recs)
    e = recs(r);
    %% switch vars within table
    vars = {'chGreen' 'chRed'; 'chGreen405' 'chRed405'; 'peaksG' 'peaksR'; 'peaksGRun' 'peaksRRun'};
    for v = 1:size(vars,1)
        % switch chG/chR
        try
            chG = data.data{e}.(vars{v,1}){1};
            data.data{e}.(vars{v,1}){1} = data.data{e}.(vars{v,2}){1};
            data.data{e}.(vars{v,2}){1} = chG;
        catch
            data.data{e} = renamevars(data.data{e},vars{v,2},vars{v,1});
        end
    end
    
    %% switch 'Bad405G' and 'Bad405G' if exist
    g = any(ismember(data.data{e}.Properties.VariableNames,'Bad405G'));
    r = any(ismember(data.data{e}.Properties.VariableNames,'Bad405R'));
    if g && r
        chG = data.data{e}.Bad405G{1};
        data.data{e}.Bad405G{1} = data.data{e}.Bad405R{1};
        data.data{e}.Bad405R{1} = chG;
    elseif g
        data.data{e}.Bad405R{1} = data.data{e}.Bad405G{1};
        data.data{e}.Bad405G = [];
    elseif r
        data.data{e}.Bad405G{1} = data.data{e}.Bad405R{1};
        data.data{e}.Bad405R = [];
    end
    

    %% switch other vars
    % switch rec loc and depths
    vars = {'chG' 'chR'; 'depthG' 'depthR'};
    for v = 1:size(vars,1)
        chG = data.(vars{v,1}){e};
        data.(vars{v,1}){e} = data.(vars{v,2}){e};
        data.(vars{v,2}){e} = chG;
    end

    % switch base, norm, sig2noise
    vars = {'base' 'norm' 'sig2noise' 'Acc405' 'RG405'};
    for v = 1:length(vars)
        try
            data.(vars{v}){e} = fliplr(data.(vars{v}){e});
        catch
            data.(vars{v})(e,:) = fliplr(data.(vars{v})(e,:));
        end
    end
    
    %% switch any exclusion criteria vars
    vars = find(startsWith(data.Properties.VariableNames,'ex'));
    for v = 1:length(vars)
        if length(data.(vars(v))(e,:)) == 2
            data.(vars(v))(e,:) = fliplr(data.(vars(v))(e,:));
        end
    end
    
    %% try to flip recording locations if they exist
    nam = 'RecLoc';
    recloc = startsWith(data.Properties.VariableNames,[nam 'G']);
    if any(recloc)
        id = find(recloc);
        for i = 1:length(id)
            varName = data.Properties.VariableNames{id(i)};
            extra = varName(length(nam)+2:end);
            temp = data.([nam 'G' extra]){e};
            data.([nam 'G' extra]){e} = data.([nam 'R' extra]){e};
            data.([nam 'R' extra]){e} = temp;
        end
    end
    
    
    %% save flip
    flip = data.flip;
    flip(recs) = not(flip(recs));
    data.flip = flip;
end

