%% FUNCTION to get cross correlations
function ccData = crossCorr2(data,good,RG,VA,win,sp,idx,norm)
    if ~exist('idx','var')
        idx = [];
    end
    good2 = find(good);
    
    ccData = {nan(1,win*2+1)};
    ccData =  repmat(ccData,length(good),3);
        
    for e2 = 1:length(good2)
        e = good2(e2);
        
        % get data
        if strcmp(RG,'R')
            chRG = data.data{e}.chRed{1};
            ch405 = data.data{e}.chRed405{1};
        elseif strcmp(RG,'G')
            chRG = data.data{e}.chGreen{1};
            ch405 = data.data{e}.chGreen405{1};
        elseif strcmp(RG,'Rd')
            chRG = data.data{e}.chRedDec{1};
            ch405 = data.data{e}.chRed405Dec{1};
        elseif strcmp(RG,'Gd')
            chRG = data.data{e}.chGreenDec{1};
            ch405 = data.data{e}.chGreen405Dec{1};
        elseif strcmp(RG,'A')
            chRG = data.data{e}.Acceleration{1};
            ch405 = data.data{e}.Acceleration{1};
        end       
        
        if strcmp(VA,'V')
            chMov = data.data{e}.chMov{1};
        elseif strcmp(VA,'A')
            chMov = data.data{e}.Acceleration{1};
        end
        run = logical(full(data.data{e}.Run{1}));
        
        % if exists 405 exclusion var, remove those
        if (strcmp(RG,'R') || strcmp(RG,'Rd')) && any(strcmp(data.data{e}.Properties.VariableNames,'Bad405R'))
            chRG(data.data{e}.Bad405R{1}) = [];
            ch405(data.data{e}.Bad405R{1}) = [];
            chMov(data.data{e}.Bad405R{1}) = [];
            run(data.data{e}.Bad405R{1}) = [];
        elseif (strcmp(RG,'G') || strcmp(RG,'Gd')) && any(strcmp(data.data{e}.Properties.VariableNames,'Bad405G'))
            chRG(data.data{e}.Bad405G{1}) = [];
            ch405(data.data{e}.Bad405G{1}) = [];
            chMov(data.data{e}.Bad405G{1}) = [];
            run(data.data{e}.Bad405G{1}) = [];
        end
            
            
        % remove nan
        bad = isnan(chRG) | isnan(ch405) | isnan(chMov);
        if sum(bad) > length(chRG)*0.01/100
            warning('Please check, too many nans in data');
        end
        chMov(bad) = [];
        chRG(bad) = [];
        ch405(bad) = [];

        % keep only run times
        run(bad) = [];
        chMov = chMov(run);
        chRG = chRG(run);
        ch405 = ch405(run);

        % get cross correlation
        try
            ccRG = crosscorr(chMov,chRG,win);
            cc405 = crosscorr(chMov,ch405,win);
            ccMov = crosscorr(chMov,chMov,win);

            % add to table
            ccData{e,1} = ccRG;
            ccData{e,2} = cc405;
            ccData{e,3} = ccMov;
        end
    end
    
    % smooth
    smoothWin = 5;
    for e2 = 1:length(good2)
        e = good2(e2);
        for i = 1:3
            temp = ccData{e,i};
            temp = fastsmooth(temp,smoothWin);
            temp(temp == 0) = ccData{e,i}(temp == 0);
            ccData{e,i} = temp(:)';
        end
    end
    
    % get means and sems
    means = {nan(1,win*2+1)};
    means =  repmat(means,2,3);
    for i = 1:3
        temp = cell2mat(ccData(:,i));
        means{1,i} = mean(temp,1,'omitnan');
        means{2,i} = std(temp,[],1,'omitnan')/sqrt(sum(any(~isnan(temp),2)));
    end

    
    %% plot
    load colors.mat C
    exp = data.Exp{find(good,1)};
    if strcmp(exp,'VGlut')
        colors = [C.lightRed{1};C.red{1}];
    elseif strcmp(exp,'Calb')
        colors = [C.lightYellow{1};C.orange{1}];
    elseif strcmp(exp,'Raldh')
        colors = [C.lightGreen{1};C.green{1}];
    elseif strcmp(exp,'Anxa')
        colors = [C.lightAqua{1};C.aqua{1}];
    elseif strcmp(exp,'Dat')
        colors = [C.lightGrey{1};C.grey{1}];
    else
        colors = [C.lightYellow{1};C.yellow{1}];
    end

    % plot means and sem
    subplot(sp(1),sp(2),sp(3))
    cla
    hold on
    hz = 100;
    ran = (-win:win)/hz;
    fill([ran fliplr(ran)], [means{1,2}+means{2,2} fliplr(means{1,2}-means{2,2})],C.lightBlue{1},'EdgeColor',C.blue{1})
    fill([ran fliplr(ran)], [means{1,1}+means{2,1} fliplr(means{1,1}-means{2,1})],colors(1,:),'EdgeColor',colors(2,:))
    
    if strcmp(VA,'V')
        ylim([-0.25 0.25])
    else
        ylim([-0.11 0.11])
    end
    verticalLine(0);
    ylabel('Cross correlation')
    xlabel('Lag (s)')
    
    % plot heatmaps
    subplot(sp(1),sp(2),sp(3)+1)
    cla
    
    
    heatmap = cell2mat(ccData(:,1));
    heatmap(any(isnan(heatmap),2),:) = [];
    % normalize heatmap
    if norm == 1
        % normalize heatmap - from 0 to 1
        heatmap = heatmap - min(heatmap,[],2);
        heatmap = heatmap ./ max(heatmap,[],2);
    elseif norm == 2
        % normalize heatmap - from -1 to 1 preserving 0
        minmax = max(abs(heatmap),[],2);
        heatmap = heatmap ./ minmax;
    end
    
    % order heatmap
    if ~isempty(idx)
        heatmap = heatmap(idx,:);
    end

    
    % plot
    imagesc(heatmap);
    verticalLine(win+1);
    xticks([])
    
    mice = length(unique(data.Mouse(good,:)));
    exp = unique(data.Exp(good,:));
        
    title([exp{1} ' - ' num2str(mice) 'm, ' num2str(sum(good)) 'r']);

end
