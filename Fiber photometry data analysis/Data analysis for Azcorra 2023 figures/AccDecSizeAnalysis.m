function AccDecSizeAnalysis(data6)
% Fig. 5A,B


%% switch channels for recordings where str is in chG, duplicate and switch for those with double str 
data6(data6.dup,:) = []; %remove previous dups
dup = ~strcmp(data6.chR,'snc') & ~strcmp(data6.chR,'o');
temp = data6(dup,:);
temp.dup = not(temp.dup); %label as dups
data6 = [data6; temp];
dupNum = sum(dup);

flip = false(size(data6,1),1);
flip(end-dupNum+1:end) = true;
flip2 = find(flip);

data6 = flipRG(data6,flip2);
data6.flip(flip2) = true;

% sort
data6 = sortrows(data6,1);

%% get good
strSNc = data6.exStr(:,1) == 1; % only str
good = double(strSNc & data6.exSig2Noise(:,1) & data6.ex405Acc(:,1) & data6.exRun & data6.Bad405All(:,1) == 0);
%good = double(strSNc & data6.exSig2Noise(:,1) & data6.ex405Acc(:,1) & data6.exRun & data6.Bad405All(:,1) == 0 & strcmp(data6.RunRew,'rew'));
subtypes = {'VGlut' 'Calb' 'Anxa'};
load colors.mat C
good(good == 1) = 6;
for j = 1:length(subtypes)
    sIdx = strcmp(data6.Exp,subtypes{j});
    good(good > 0 & sIdx) = j;
end


%% Select acceleration and decelerations by percentiles
accThreshPer = 0:0.2:1; % percentiles of acc above thresh
accThreshMin = 2;

for e = 1:size(data6,1)
    if good(e) > 0 && good(e) < 6
        acceleration = data6.data{e}.Acceleration{1};
        [acc,dec,peaks] = selectAccDec(acceleration,accThreshMin);
        
        data6.data{e}.AccOn{1} = zeros(size(acc));
        data6.data{e}.DecOn{1} = zeros(size(acc));
        data6.data{e}.AccOn{1}(acc) = peaks(acc);
        data6.data{e}.DecOn{1}(dec) = peaks(dec);
        
        % exclude reward
        data6.data{e} = excludeRew(data6.data{e},'AccOn');
        data6.data{e} = excludeRew(data6.data{e},'DecOn');
    end
end

for e = 1:size(data6,1)
    if good(e) > 0 && good(e) < 6
        acc = data6.data{e}.AccOn{1};
        dec = data6.data{e}.DecOn{1};
        accSort = sort(acc(acc ~= 0));
        decSort = sort(dec(dec ~= 0),'descend');
        for i = 1:length(accThreshPer)-1
            accThreshMin = accSort(round(accThreshPer(i)*length(accSort))+1);
            accThreshMax = accSort(round(accThreshPer(i+1)*length(accSort)));
            decThreshMin = decSort(round(accThreshPer(i)*length(decSort))+1);
            decThreshMax = decSort(round(accThreshPer(i+1)*length(decSort)));
            
            accTemp = acc;
            accTemp(accTemp < accThreshMin | accTemp > accThreshMax) = 0;
            decTemp = dec;
            decTemp(decTemp > decThreshMin | decTemp < decThreshMax) = 0;
            
            data6.data{e}.(['AccOn' num2str(accThreshPer(i)*100) '-' num2str(accThreshPer(i+1)*100)]){1} = logical(accTemp);
            data6.data{e}.(['DecOn' num2str(accThreshPer(i)*100) '-' num2str(accThreshPer(i+1)*100)]){1} = logical(decTemp);
        end
    end
end


%% Trig ave on acc/dec/transAV = 'A';
trigs = {'AccOn', 'DecOn'};
neg = 100;
pos = 100;
AV = 'A';

ave = cell(3,length(accThreshPer)-1,2);
aveV = cell(3,length(accThreshPer)-1,2);
pop = cell(3,length(accThreshPer)-1,2);

for t = 1:2
    figure
    set(gcf, 'Position',  [1 42+470*(t-1) 1590 470])
    for i = 1:3
        for j = 1:length(accThreshPer)-1
            extra = [num2str(accThreshPer(j)*100) '-' num2str(accThreshPer(j+1)*100)];
            [~,~,~,~,temp,~,temp2] = plotAverages_paper(data6, '1', [trigs{t} extra], 'G', good == i, 1, AV, neg, pos, 1:sum(good == i),[3,10,i*10-11+2*j],1);
            ave{i,j,t} = mean(temp,'omitnan');
            for e = 1:size(temp2,1)
                if ~isempty(temp2{e})
                    temp2{e} = mean(temp2{e},'omitnan');
                else
                    temp2{e} = nan(1,neg+pos);
                end
            end
            aveV{i,j,t} = mean(cell2mat(temp2),'omitnan');
            
            vals = [max(temp(:,neg:end),[],2)-temp(:,neg), min(temp(:,neg:end),[],2)-temp(:,neg)];
            [~,idx] = max(abs(vals),[],2);
            vals2 =  nan(size(vals,1),1);
            for v = 1:length(vals2)
                vals2(v) =  vals(v,idx(v));
            end
            pop{i,j,t} = vals2;
        end
    end
    close()
end

diffVal = nan(3,5,2);
figure
colors = [C.red{1}; C.orange{1}; C.aqua{1}];
for t = 1:2
    for s = 1:3
        c = colors(s,:);
        colors2 = [c+(1-c)/3*2; c+(1-c)/3; c; c/3*2; c/3];
        subplotMN(4,2,s,t)
        hold on
        for j = 1:length(accThreshPer)-1
            temp1 = ave{s,j,t};
            plot((-neg+1)/100:1/100:neg/100,temp1,'Color',colors2(j,:),'LineWidth',2)
            
            vals = [max(temp1(neg:end))-temp1(neg), min(temp1(neg:end))-temp1(neg)];
            if abs(vals(2)) > abs(vals(1))
                diffVal(s,j,t) = vals(2);
            else
                diffVal(s,j,t) = vals(1);
            end
        end
        verticalLine(0);
        xlabel('Time (s)')
        ylabel('Norm \DeltaF/F')
    end
    
    subplotMN(4,2,4,t)
    hold on
    for s = 1:3
        vAve = nan(5,1);
        vSem = nan(5,1);
        for j = 1:length(accThreshPer)-1
            vals = pop{s,j,t};
            vAve(j) = mean(vals,'omitnan');
            vSem(j) = std(vals,'omitnan')/sqrt(sum(~isnan(vals)));
        end
        errorbar(1:5,vAve,vSem,'LineStyle','none','Color',colors(s,:),'LineWidth',1.5,'CapSize',20)
        plot(1:5,vAve,'Color',colors(s,:));
        %plot(1:5,diffVal(s,:,t),'o-','Color',colors(s,:));
        xlim([0.5 5.5])
    end
end
disp('p-vals smalles vs largest acc/dec - with BONF:')
p(1) = signrank(pop{1,1,2},pop{1,5,2})*3;
p(2) = signrank(pop{2,1,2},pop{2,5,2})*3;
p(3) = signrank(pop{3,1,1},pop{3,5,1})*3;
disp(p)
n = [length(pop{1,1,2}), length(pop{2,1,2}), length(pop{3,1,1})];
disp(['   n = ' num2str(n)])

disp('increase in amplitude from smalles to largest acc/dec (%):')
k(1) = mean(pop{1,5,2},'omitnan')/mean(pop{1,1,2},'omitnan')*100;
k(2) = mean(pop{2,5,2},'omitnan')/mean(pop{2,1,2},'omitnan')*100;
k(3) = mean(pop{3,5,1},'omitnan')/mean(pop{3,1,1},'omitnan')*100;
disp(k)

% for acc
figure
colors = [C.red{1}; C.orange{1}; C.aqua{1}];
for t = 1:2
    for s = 1:3
        c = colors(s,:);
        colors2 = [c+(1-c)/3*2; c+(1-c)/3; c; c/3*2; c/3];
        subplotMN(3,2,s,t)
        hold on
        for j = 1:length(accThreshPer)-1
            temp1 = aveV{s,j,t};
            plot((-neg+1)/100:1/100:neg/100,temp1,'Color',colors2(j,:),'LineWidth',2)
        end
        xlabel('Time (s)')
        ylabel('Norm \DeltaF/F')
        ylim([-2.6 2.6])
        verticalLine(0);
    end
end

end


%% FUNCTION select acc/dec
function [acc2,dec2,peaks] = selectAccDec(acceleration,accMin,accMax)

threshRemove = 2;

accOn = ZeroX(acceleration); %get zero crossings - output is column
accOn(1) = []; %otherwise first is 0
dur = [diff(accOn)' nan]; %acc durations

peakLoc = nan(1,length(accOn));
peaks = nan(size(acceleration));
acc = false(size(acceleration));
dec = false(size(acceleration));
accDecRemove = false(size(acceleration));

for i = 1:length(accOn)-1
    if dur(i) > 5 % min duration
        [~,peakIdx] = max(abs(acceleration(accOn(i):accOn(i+1))));
        peakLoc(i) = acceleration(accOn(i)+peakIdx-1);
        if ~exist('accMax','var')
            if peakLoc(i) > accMin
                acc(accOn(i)) = true;
                peaks(accOn(i)) = peakLoc(i);
            elseif peakLoc(i) < -accMin
                dec(accOn(i)) = true;
                peaks(accOn(i)) = peakLoc(i);
            end
        else
            if peakLoc(i) > accMin && peakLoc(i) < accMax
                acc(accOn(i)) = true;
            elseif peakLoc(i) < -accMin && peakLoc(i) > -accMax
                dec(accOn(i)) = true;
            end
        end
        if peakLoc(i) > threshRemove || peakLoc(i) < -threshRemove
            accDecRemove(accOn(i)) = true;
        end
    end
end

% remove any acc/dec that are too close to each other - exclude zigzag running
timeMinDist = 50; %0.5s - 0.25 on each side
all = acc | dec;
all0 = [false(ceil(timeMinDist/2),1);accDecRemove(:);false(floor(timeMinDist/2),1)];
allMat = nan(size(all0,1),timeMinDist+1);
for i = -ceil(timeMinDist/2):floor(timeMinDist/2)
    if i ~= 0
        % exclude acc/dec with either acc or dec next to it
        allMat(:,i+ceil(timeMinDist/2)+1) = circshift(all0,i);
    end
end
allMat(:,ceil(timeMinDist/2)+1) = [];
allRemove = logical(sum(allMat,2) > 0);
allRemove = allRemove(ceil(timeMinDist/2)+1:end-floor(timeMinDist/2));
all(allRemove) = false;

acc2 = acc(:) & all(:);
dec2 = dec(:) & all(:);
end


%% FUNCTION exclude reward times
function data = excludeRew(data, varName)
if any(strcmp(data.Properties.VariableNames,'Reward'))
    win = 500; %5s
    rew = find(data.RewardAll{1} | data.AirPuff{1} | data.Light{1});
    if ~isempty(rew)
        for j = 1:length(rew)
            try
                data.(varName){1}(rew(j):rew(j)+win) = 0;
            catch
                data.(varName){1}(rew(j):end) = 0;
            end
        end
    end
end
end


