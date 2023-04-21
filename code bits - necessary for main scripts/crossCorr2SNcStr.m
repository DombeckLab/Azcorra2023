%% FUNCTION to get cross correlations
function ccData = crossCorr2SNcStr(data,good,deconvolve,win,sp,idx)
    if ~exist('idx','var')
        idx = [];
    end
    good2 = find(good);
    
    ccData = {nan(1,win*2+1)};
    ccData =  repmat(ccData,length(good),2);
        
    for e2 = 1:length(good2)
        e = good2(e2);
        
        % get data
        if deconvolve == 0
            chG = data.data{e}.chGreen{1};
            chG405 = data.data{e}.chGreen405{1};
            chR = data.data{e}.chRed{1};
            chR405 = data.data{e}.chRed405{1};
        elseif deconvolve == 1
            chG = data.data{e}.chGreenDec{1};
            chG405 = data.data{e}.chGreen405Dec{1};
            chR = data.data{e}.chRedDec{1};
            chR405 = data.data{e}.chRed405Dec{1};
        end
        
        % if exists 405 exclusion var, remove those
        if any(strcmp(data.data{e}.Properties.VariableNames,'Bad405R'))
            remove = data.data{e}.Bad405G{1} | data.data{e}.Bad405R{1};
            chG(remove) = [];
            chG405(remove) = [];
            chR(remove) = [];
            chR405(remove) = [];
        end
        
        
        % remove nan
        bad = isnan(chG) | isnan(chG405) | isnan(chR) | isnan(chR405);
        if sum(bad) > length(chG)*0.01/100
            warning('Please check, too many nans in data');
        end
        chG(bad) = [];
        chG405(bad) = [];
        chR(bad) = [];
        chR405(bad) = [];

       % get cross correlation
        try
            ccGR = crosscorr(chG,chR,win);
            ccGR405 = crosscorr(chG,chR405,win);
            ccG405R = crosscorr(chG405,chR,win);
            cc405 = mean([ccGR405 ccG405R],2);

            % add to table
            ccData{e,1} = ccGR;
            ccData{e,2} = cc405;
        end
    end
    
    % smooth
    smoothWin = 5;
    for e2 = 1:length(good2)
        e = good2(e2);
        for i = 1:2
            temp = ccData{e,i};
            temp = fastsmooth(temp,smoothWin);
            temp(temp == 0) = ccData{e,i}(temp == 0);
            ccData{e,i} = temp(:)';
        end
    end
    
    % get means and sems
    means = {nan(1,win*2+1)};
    means =  repmat(means,2,3);
    for i = 1:2
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
    end

    % plot means and sem
    subplot(sp(1),sp(2),sp(3))
    cla
    hold on
    hz = 100;
    ran = (-win:win)/hz;
    fill([ran fliplr(ran)], [means{1,2}+means{2,2} fliplr(means{1,2}-means{2,2})],C.lightBlue{1},'EdgeColor',C.blue{1})
    fill([ran fliplr(ran)], [means{1,1}+means{2,1} fliplr(means{1,1}-means{2,1})],colors(1,:),'EdgeColor',colors(2,:))
    
    ylim([-0.1 0.8])
    verticalLine(0);
    ylabel('Cross correlation')
    xlabel('Lag (s)')
    
    % plot heatmaps
    subplot(sp(1),sp(2),sp(3)+1)
    cla
    
    
    heatmap = cell2mat(ccData(:,1));
    heatmap(any(isnan(heatmap),2),:) = [];
        
    % order heatmap
    if isempty(idx)
        peaks = max(heatmap,[],2);
        [~,idx] = sort(peaks);
    end
    heatmap = heatmap(idx,:);

    
    % plot
    imagesc(heatmap);
    verticalLine(win+1);
    xticks([])
    
    mice = length(unique(data.Mouse(good,:)));
    exp = unique(data.Exp(good,:));
        
    title([exp{1} ' - ' num2str(mice) 'm, ' num2str(sum(good)) 'r']);

end
