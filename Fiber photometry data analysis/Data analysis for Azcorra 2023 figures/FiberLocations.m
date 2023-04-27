function FiberLocations2(data6)
%% set DATA PROCESSING FOLDER (edit MatlabFolders.mat file in code folder)
filePath = matlab.desktop.editor.getActiveFilename;
i = strfind(filePath, '\');
filePath = filePath(1:i(end));
load([filePath 'MatlabFolders.mat'], 'dataProcessingFolder');

%% set permanent variables
secEdge = 204.4; %(in pix)
scale = 35.6667; %because the size of the image is 6 mm and 214 pixels tall

%% add new fiber locations per mouse and line
% this section was used to add fiber locations for new mice. It requires a
% brain image with the fiber traced onto it. To use, change the first line
% from [] to 1.
for xx = []
sub = 'Anxa'; %change for each
clear fibers
try
    load([dataProcessingFolder '\fibers_' sub '.mat'], 'fibers')
    newMice = input('Input NEW mice numbers as cell: ');
    fibers = [fibers; newMice(:) cell(length(newMice),1)];
    n = length(newMice);
catch
    mice = input('Input mice numbers as cell: ');
    fibers = [mice(:) cell(length(mice),1)];
    n = length(mice);
end

load([dataProcessingFolder '\brain_' sub '.mat'], 'brain')
figure;
imshow(brain)
n2 = size(fibers,1);
for i = n-1:-1:0
    title(fibers{n2-i,1})
    roi = drawline;
    input([fibers{n2-i,1} ' - done?']);
    fibers{n2-i,2} = roi.Position;
end
save([dataProcessingFolder '\fibers_' sub '.mat'], 'fibers')
end


%% get recording locations onto data6
% same as previous section, used to add fiber locations for new mice
for xx = []
if ~exist('data6','var')
    load([dataProcessingFolder '\data6.mat'], 'data6')
end
subtypes = {'Dat' 'Raldh' 'VGlut' 'Calb' 'Anxa'};
noHist = {};
for s = 1:5
    sub = subtypes{s};
    load([dataProcessingFolder '\fibers_' sub '.mat'], 'fibers')
    good = find(strcmp(data6.Exp,sub));

    mice2 = fibers(:,1);
    for i = 1:length(mice2)
        k = strfind(mice2{i},'-');
        if ~isempty(k)
            mice2{i} = mice2{i}(1:k-1);
        end
    end
    scale = 35.6667; %because the size of the image is 6mm
    
    for ee = 1:length(good)
        e = good(ee);
        mouse = data6.Mouse{e};
        try
            locG = data6.chG{e};
            locR = data6.chR{e};
            idxG = nan;
            idxR = nan;
            idx = find(strcmp(mice2,mouse));
            if isempty(idx) % no histology for that mouse
                if ~ismember(mouse,noHist)
                    noHist = [noHist; mouse];
                end
            elseif length(idx) == 1 % only one fiber track for that mouse
                if ~strcmp(locG,'o') && ~strcmp(locG,'snc')
                    idxG = idx; 
                end
                if ~strcmp(locR,'o') && ~strcmp(locR,'snc')
                    idxR = idx;
                end
            else % more than one fiber track per mouse
                details = cell(size(fibers,1),1);
                for j = 1:size(details,1)
                    if ismember(j,idx)
                        k = strfind(fibers{j,1},'-');
                        details{j} = fibers{j,1}(k+1:end);
                    else
                        details{j} = '';
                    end
                end
                if ~any(contains(details,'1')) %for ds/ts, dms/dls
                    idxG = find(strcmp(details,locG));
                    idxR = find(strcmp(details,locR));
                else %more than one fibertrack per mouse & loc - go case by case
                    if strcmp(mouse,'D092') && strcmp(data6.Date{e},'20200717')
                        idxG = find(strcmp(details,[locG '1']));
                        idxR = find(strcmp(details,[locR '1']));
                    elseif strcmp(mouse,'D092') && ~strcmp(data6.Date{e},'20200717')
                        idxG = find(strcmp(details,[locG '2']));
                        idxR = find(strcmp(details,[locR '2']));
                    elseif strcmp(mouse,'1849') && strcmp(data6.Date{e},'20210623')
                        idxG = find(strcmp(details,'2'));
                    elseif strcmp(mouse,'1849') && ~strcmp(data6.Date{e},'20210623')
                        idxG = find(strcmp(details,'1'));
                    elseif strcmp(mouse,'1790') && strcmp(data6.Date{e},'20210412')
                        idxG = find(strcmp(details,'2'));
                    elseif strcmp(mouse,'1790') && ~strcmp(data6.Date{e},'20210412')
                        idxG = find(strcmp(details,'1'));
                    else
                        warning(['error 1 with mouse ' num2str(e) ' - ' mouse]) 
                    end
                end    
            end
            
            if ~isnan(idxG)
                depth = data6.depthG{e};
                fiberTrack = fibers{idxG,2};
                v = fiberTrack(2,:) - fiberTrack(1,:);
                u = v/(sqrt(v(1)^2+v(2)^2));

                RecLoc = fiberTrack(1,:) + depth*scale*u;
                data6.RecLocG{e} = RecLoc;
            else
                data6.RecLocG{e} = nan(1,2);
            end
            if ~isnan(idxR)
                depth = data6.depthR{e};
                fiberTrack = fibers{idxR,2};
                v = fiberTrack(2,:) - fiberTrack(1,:);
                u = v/(sqrt(v(1)^2+v(2)^2));

                RecLoc = fiberTrack(1,:) + depth*scale*u;
                data6.RecLocR{e} = RecLoc;
            else
                data6.RecLocR{e} = nan(1,2);
            end
        catch
            warning(['error 2 with mouse ' num2str(e) ' - ' mouse]) 
        end
    end
end
save([dataProcessingFolder '\data6.mat'], 'data6','-v7.3')
end


%% convert coordinates to mm in XYZ
% same as previous sections, used to add fiber locations for new mice
for xx = []
midLine = [37 36 37 46 32 45 45]; %distance from sectionEdge to midline (in pix)
rosCau = [0.86, 0.50, 0.14, -0.22, -0.58, -0.94, -1.34]; %slice mm rostro-ventral from bregma
dorVen = [0.67 0.68 0.36 0.47 0.32 0.37 0.22]; %top mm from bregma

vars = {'RecLocR','RecLocG'};
for e = 1:size(data6,1)
    for v = 1:2
        if ~any(isnan(data6.(vars{v}){e}))
            temp = data6.(vars{v}){e};
            sec = ceil(temp(1)/secEdge);
            
            X = (temp(1)-(sec-1)*secEdge-midLine(sec))/scale; % in mm
            Y = rosCau(sec);
            Z = temp(2)/scale + dorVen(sec); % in mm
            RecLocNew = [X,Y,Z]; % in mm
            
            data6.([vars{v} 'mm']){e} = RecLocNew;
        else
            data6.([vars{v} 'mm']){e} = [nan nan nan];
        end
    end
end

save([dataProcessingFolder '\data6.mat'], 'data6','-v7.3')
end

%% shift coordinates to match to TS or DS mid-slice
% same as previous section, used to add fiber locations for new mice
for xx = []
secEdge = 204.4; %(all in pixels)
shiftX = [11 -204 -409 0 -598 -821 -1028];
shiftY = [0 0 -5 0 3 0 -6];


vars = {'RecLocR','RecLocG'};
for e = 1:size(data6,1)
    for v = 1:2
        if ~any(isnan(data6.(vars{v}){e}))
            temp = data6.(vars{v}){e};
            sec = ceil(temp(1)/secEdge);
            
            temp(1) = temp(1) + shiftX(sec);
            temp(2) = temp(2) + shiftY(sec);
            
            data6.([vars{v} 'shift']){e} = [temp(1) temp(2)];
            scatter(temp(1),temp(2))
        else
            data6.([vars{v} 'shift']){e} = [nan nan];
        end
    end
end

save([dataProcessingFolder '\data6.mat'], 'data6','-v7.3')
end


%% map locations to rew/loc values
subtypes = {'VGlut' 'Calb' 'Raldh' 'Anxa' 'Dat' };
load colors.mat C

if ~exist('data6','var')
    load([dataProcessingFolder '\data6.mat'], 'data6')
end
miceG = [data6.Exp data6.chG];
miceR = [data6.Exp data6.chR];
for e = 1:size(data6,1)
    miceG{e,1} = data6.data{e}.Properties.RowNames{1};
    miceR{e,1} = data6.data{e}.Properties.RowNames{1};
end
miceG = strcat(miceG(:,1),'-',miceG(:,2));
miceR = strcat(miceR(:,1),'-',miceR(:,2));

load([dataProcessingFolder '\good_Rew.mat'], 'goodRew');
load([dataProcessingFolder '\good_Run.mat'], 'goodRun');
load('brain_Empty.mat', 'brain')
goodAll = [{goodRew};{goodRun}];

for i = 1:2
    for e = 1:size(goodAll{i},1)
        mouse = strcat(goodAll{i}.recID(e),'-',goodAll{i}.chG(e));
        idG = find(strcmp(mouse,miceG));
        idR = find(strcmp(mouse,miceR));
        if ~isempty(idG) && ~isempty(idR) && any(~isnan(data6.RecLocG{idG})) && any(~isnan(data6.RecLocR{idR}))
            warning(['repeated mouse ' num2str(e) ' on i = ' num2str(i) ': ' mouse{1}])
        elseif isempty(idG) && isempty(idR)
            warning(['no hist mouse ' num2str(e) ' on i = ' num2str(i) ': ' mouse{1}])
        elseif ~isempty(idG)
            goodAll{i}.Loc{e} = data6.RecLocG{idG};
            goodAll{i}.LocShift{e} = data6.RecLocGshift{idG};
        elseif ~isempty(idR)
            goodAll{i}.Loc{e} = data6.RecLocR{idR};
            goodAll{i}.LocShift{e} = data6.RecLocRshift{idR};
        end
    end
end

%% plot dots for recording locations colorcoded by locomotion, reward, and airpuff
% Fig. S7A and S8J,K
varIdx = [2 1 1];
varNam = {'LocomAngles', 'RewardResp', 'AirPuffResp'};
cX = cell(3,1);
cScaleX = cX;
%cScaleX{1} = 0:360;
%cX{1} = circshift(customcolormap([0 0.25 0.5 0.75 1],[0 0 0; C.aqua{1}; C.orange{1}; C.red{1};0 0 0],360),45);
cScaleX{2} = linspace(-0.05,0.2,250);
cScaleX{3} = cScaleX{2};
cX{2} = customcolormap([0 0.3 1],[1 0 0; 1 0.8 0; 0 0 1],length(cScaleX{3}));
cX{3} = cX{2};

% make 2d colormap for velocity
[X,Y] = meshgrid(-2:1:2);
V = nan(size(X,1),size(X,2),3);
V([1,2],[1,2],:) = repmat(reshape(C.orange{1},1,1,3),2,2,1);
V([4,5],[1,2],:) = repmat(reshape(C.red{1},1,1,3),2,2,1);
V([4,5],[4,5],:) = repmat(reshape([0 0 0],1,1,3),2,2,1);
V([1,2],[4,5],:) =  repmat(reshape(C.aqua{1},1,1,3),2,2,1);
V(3,3,:) = [0.8 0.8 0.8];

V([1,2],3,:) = repmat(reshape(mean([C.orange{1};C.aqua{1}]),1,1,3),2,1,1);
V(3,[1,2],:) = repmat(reshape(mean([C.red{1};C.orange{1}]),1,1,3),1,2,1);
V([4,5],3,:) = repmat(reshape(mean([C.red{1};0 0 0]),1,1,3),1,2,1);
V(3,[4,5],:) = repmat(reshape(mean([C.aqua{1};0 0 0]),1,1,3),2,1,1);

[Xq,Yq] =  meshgrid(-2:0.01:2);
Vq = nan(size(Xq,1),size(Xq,2),3);
for i = 1:3
    Vq(:,:,i) = interp2(X,Y,V(:,:,i),Xq,Yq,'cubic');
end

Vq(Vq<0) = 0;
cX{1} = Vq;
cScaleX{1} = linspace(-4,4,size(Vq,1));

for s = 1:5
    sub = subtypes{s};
    % plot dots for recording locations colorcoded by locomotion, reward, and airpuff
    figure
    set(gcf, 'Position',  [260 41 1500 960])
    for vv = 1:3     
        subplot(3,1,vv)
        imshow(brain)
        hold on
        
        if vv > 1
            values1 = goodAll{varIdx(vv)}.(varNam{vv});
        else
            values1 = goodAll{varIdx(vv)}.LocomPC1;
            values2 = goodAll{varIdx(vv)}.LocomPC2;
        end
        c = cX{vv};
        cScale = cScaleX{vv};
        
        for i = 1:size(goodAll{varIdx(vv)},1)
            if strcmp(goodAll{varIdx(vv)}.Exp{i},sub) && ~isnan(values1(i))
                try
                    RecLoc = goodAll{varIdx(vv)}.Loc{i};
                    RecLoc(1) = RecLoc(1);
                    if vv == 1
                        [~,cIdxX] = min(abs(values2(i)-cScale));
                        [~,cIdxY] = min(abs(values1(i)-cScale));
                        scatter(RecLoc(1),RecLoc(2),[],reshape(c(cIdxX,cIdxY,:),1,3,1),'filled','MarkerEdgeColor',reshape(c(cIdxX,cIdxY,:),1,3,1),'MarkerFaceAlpha',0.5)
                    else
                        [~,cIdx] = min(abs(values1(i)-cScale));
                        scatter(RecLoc(1),RecLoc(2),[],c(cIdx,:),'filled','MarkerEdgeColor',c(cIdx,:),'MarkerFaceAlpha',0.5)
                    end
                catch
                end
            end
        end
        title([varNam{vv} ' - ' sub])
    end
end
figure
imagesc(Vq)
len = floor(size(Vq,1)/2);
xlim([len-len/max(cScaleX{1})*2.5, len+len/max(cScaleX{1})*2.5])
ylim([len-len/max(cScaleX{1})*2.5, len+len/max(cScaleX{1})*2.5])
verticalLine(len);
horizontalLine(len);
set(gca, 'YDir', 'normal');
axis equal



%% plot values shifted to match to TS or DS mid-slice plus jitter
% add jitter (than can be removed and matches each experiment type)
% Fig. 3B and Fig. 4J,K
jitterVal = 0.4; %mm
jitterValPix = round(jitterVal*scale*2)/2;
try
    load([dataProcessingFolder '\jitter.mat'], 'jitter'); %saved for consistency, generated as below
catch
    jitter = rand(size(goodAll{1},1),1)*(jitterValPix*2)-jitterValPix;
    save([dataProcessingFolder '\jitter.mat'], 'jitter');
end
temp = cell2mat(goodAll{1}.LocShift)+[jitter(:) zeros(size(jitter(:)))];
goodAll{1}.LocShift = mat2cell(temp,ones(size(temp,1),1),2);
goodAll{2}.LocShift = mat2cell(temp,ones(size(temp,1),1),2);

brain2 = [brain(:,round(secEdge+1):round(secEdge*2)) brain(:,round((secEdge+1)*5):round(secEdge*6))];
figure
set(gcf, 'Position',  [-1198 36 2683 960])
for s = 1:5   
    sub = subtypes{s};
    for vv = 1:3 
        subplotMN(3,5,vv,s)
        imshow(brain2)
        hold on
        
        if vv > 1
            values1 = goodAll{varIdx(vv)}.(varNam{vv});
        else
            values1 = goodAll{varIdx(vv)}.LocomPC1;
            values2 = goodAll{varIdx(vv)}.LocomPC2;
        end
        c = cX{vv};
        cScale = cScaleX{vv};
        
        for i = 1:size(goodAll{varIdx(vv)},1)
            if strcmp(goodAll{varIdx(vv)}.Exp{i},sub) && ~isnan(values1(i))
                try
                    RecLoc = goodAll{varIdx(vv)}.LocShift{i};
                    RecLoc(1) = RecLoc(1);
                    
                    if vv == 1
                        [~,cIdxX] = min(abs(values2(i)-cScale));
                        [~,cIdxY] = min(abs(values1(i)-cScale));
                        scatter(RecLoc(1),RecLoc(2),[],reshape(c(cIdxX,cIdxY,:),1,3,1),'filled','MarkerEdgeColor',reshape(c(cIdxX,cIdxY,:),1,3,1),'MarkerFaceAlpha',0.5)
                    else
                        [~,cIdx] = min(abs(values1(i)-cScale));
                        scatter(RecLoc(1),RecLoc(2),[],c(cIdx,:),'filled','MarkerEdgeColor',c(cIdx,:),'MarkerFaceAlpha',0.5)
                    end
                catch
                end
            end
        end
        title([varNam{vv} ' - ' sub])
    end
end

% remove jitter (than can be removed and matches each experiment type)
temp = cell2mat(goodAll{1}.LocShift)-[jitter(:) zeros(size(jitter(:)))];
goodAll{1}.LocShift = mat2cell(temp,ones(size(temp,1),1),2);
goodAll{2}.LocShift = mat2cell(temp,ones(size(temp,1),1),2);


%% Plot recording locations (no values) by slice with SMALL jitter and fiber tracks
% Fig. 2E
% add smaller jitter - 0.2 mm
jitterVal2 = 0.2; %mm
jitterValPix2 = round(jitterVal2*scale*2)/2;

jitter2 = jitter/jitterValPix*jitterValPix2;
temp = cell2mat(goodAll{1}.Loc)+[jitter2(:) zeros(size(jitter2(:)))];
goodAll{1}.Loc = mat2cell(temp,ones(size(temp,1),1),2);
goodAll{2}.Loc = mat2cell(temp,ones(size(temp,1),1),2);

figure
set(gcf, 'Position',  [260 41 1500 960]) 
for s = 1:5
    sub = subtypes{s};
    if s == 3
        figure
        set(gcf, 'Position',  [260 41 1500 960])
    end
    if s >= 3
        subplot(3,1,s-2)
    else
        subplot(3,1,s)
    end
    imshow(brain)
    hold on
    vv = 1; %locomotion
    values1 = goodAll{varIdx(vv)}.LocomPC1;
    for i = 1:size(goodAll{varIdx(vv)},1)
        if strcmp(goodAll{varIdx(vv)}.Exp{i},sub) && ~isnan(values1(i))
            try
                RecLoc = goodAll{varIdx(vv)}.Loc{i};
                RecLoc(1) = RecLoc(1);
                scatter(RecLoc(1),RecLoc(2),[],[0.8 0.8 0.8],'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceAlpha',0.5)
            catch
            end
        end
    end
    title(sub)
end

temp = cell2mat(goodAll{1}.Loc)-[jitter2(:) zeros(size(jitter2(:)))];
goodAll{1}.Loc = mat2cell(temp,ones(size(temp,1),1),2);
goodAll{2}.Loc = mat2cell(temp,ones(size(temp,1),1),2);


%% Calculate distance in PC1/2 space based on distance beetwen recording locations
% Fig. S6J
load([dataProcessingFolder '\good_Run.mat'], 'goodRun');
goodRun2 = goodRun;

miceG = [data6.Exp data6.chG];
miceR = [data6.Exp data6.chR];
for e = 1:size(data6,1)
    miceG{e,1} = data6.data{e}.Properties.RowNames{1};
    miceR{e,1} = data6.data{e}.Properties.RowNames{1};
end
miceG = strcat(miceG(:,1),'-',miceG(:,2));
miceR = strcat(miceR(:,1),'-',miceR(:,2));


for e = 1:size(goodRun2,1)
    mouse = strcat(goodRun2.recID(e),'-',goodRun2.chG(e));
    idG = find(strcmp(mouse,miceG));
    idR = find(strcmp(mouse,miceR));
    if ~isempty(idG) && ~isempty(idR) && any(~isnan(data6.RecLocG{idG})) && any(~isnan(data6.RecLocR{idR}))
        warning(['repeated mouse ' num2str(e) ' on i = ' num2str(i) ': ' mouse{1}])
    elseif isempty(idG) && isempty(idR)
        warning(['no hist mouse ' num2str(e) ' on i = ' num2str(i) ': ' mouse{1}])
    elseif ~isempty(idG)
        goodRun2.LocMm{e} = data6.RecLocGmm{idG};
    elseif ~isempty(idR)
        goodRun2.LocMm{e} = data6.RecLocRmm{idR};
    end
end
goodRun2(isnan(goodRun2.LocomPC1),:) = [];
remove = false(size(goodRun2,1),1);
for e = 1:size(goodRun2,1)
    if isempty(goodRun2.LocMm{e})
        remove(e) = true;
    end
end
goodRun2(remove,:) = [];

good = zeros(size(goodRun2,1),1);
subtypes = {'VGlut' 'Calb' 'Raldh' 'Anxa' 'Dat'};
for i = 1:length(subtypes)
    sIdx = strcmp(goodRun2.Exp,subtypes{i});
    good(sIdx) = i;
end


% plot
load colors.mat C
angDist = 1; %1 = angles, 2 = dist;
pairsXYZ = cell(length(subtypes),1);
pairsPC = cell(length(subtypes),1);
for s = 1:length(subtypes)
    pairsXYZ{s} = nan(sum(good == s));
    pairsPC{s} = nan(sum(good == s));
    idx = find(good == s);
    for e1x = 1:length(idx)-1
        for e2x = e1x+1:length(idx)
            e1 = idx(e1x);
            e2 = idx(e2x);
            
            pairsXYZ{s}(e1x,e2x) = pdist([goodRun2.LocMm{e1};goodRun2.LocMm{e2}]);
            if angDist == 2
                pairsPC{s}(e1x,e2x) = pdist([goodRun2.LocomPC1(e1) goodRun2.LocomPC2(e1);goodRun2.LocomPC1(e2) goodRun2.LocomPC2(e2)]);
            elseif angDist == 1
                normDeg = mod(goodRun2.LocomAngles(e1)-goodRun2.LocomAngles(e2),360);
                pairsPC{s}(e1x,e2x) = min(360-normDeg, normDeg);
            end
        end
    end
    pairsXYZ{s} = reshape(pairsXYZ{s},size(pairsXYZ{s},1)^2,1);
    pairsPC{s} = reshape(pairsPC{s},size(pairsPC{s},1)^2,1);
    emp = isnan(pairsXYZ{s}) | isnan(pairsPC{s});
    pairsXYZ{s}(emp) = [];
    pairsPC{s}(emp) = [];
end

colors = [C.red{1}; C.orange{1}; C.green{1}; C.aqua{1}; C.grey{1}];
figure
hold on
bins = 0:0.3:3.6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binCents = bins(1:end-1) + (bins(2)-bins(1))/2;
ave = nan(5,length(bins)-1);
sem = nan(5,length(bins)-1);
minN = 5;
for s = 1:length(subtypes)
    XYZ = pairsXYZ{s};
    PC = pairsPC{s};
    for i = 1:length(bins)-1
        idx = XYZ >= bins(i) & XYZ < bins(i+1);
        if sum(idx) >= minN
            ave(s,i) = mean(PC(idx));
            sem(s,i) = std(PC(idx))/sqrt(sum(idx));
        end
    end
    errorbar(binCents,ave(s,:),sem(s,:),'Color',colors(s,:),'LineStyle','none')
    plot(binCents,ave(s,:),'Color',colors(s,:),'LineWidth',2)
end
if angDist == 2 
    ylabel('PC1/PC2 euclidian distance')
else
    ylabel('PC1/PC2 Angle diff')
end
xlabel('Distance between recording pairs (mm)')


% add comparison between pure subtypes
pairsXYZ_2 = cell(3,1);
pairsPC_2 = cell(3,1);
s1 = [1 1 2];
s2 = [2 3 3];
for c = 1:3
    pairsXYZ_2{c} = nan(sum(good == s1(c)),sum(good == s2(c)));
    pairsPC_2{c} = nan(sum(good == s1(c)),sum(good == s2(c)));
    idx1 = find(good == s1(c));
    idx2 = find(good == s2(c));
    for e1x = 1:length(idx1)
        for e2x = 1:length(idx2)
            e1 = idx1(e1x);
            e2 = idx2(e2x);
            
            pairsXYZ_2{c}(e1x,e2x) = pdist([goodRun2.LocMm{e1};goodRun2.LocMm{e2}]);
            if angDist == 2
                pairsPC_2{c}(e1x,e2x) = pdist([goodRun2.LocomPC1(e1) goodRun2.LocomPC2(e1);goodRun2.LocomPC1(e2) goodRun2.LocomPC2(e2)]);
            elseif angDist == 1
                normDeg = mod(goodRun2.LocomAngles(e1)-goodRun2.LocomAngles(e2),360);
                pairsPC_2{c}(e1x,e2x) = min(360-normDeg, normDeg);
            end 
        end
    end
    pairsXYZ_2{c} = reshape(pairsXYZ_2{c},size(pairsXYZ_2{c},1)*size(pairsXYZ_2{c},2),1);
    pairsPC_2{c} = reshape(pairsPC_2{c},size(pairsPC_2{c},1)*size(pairsPC_2{c},2),1);
    pairsXYZ_2{c}(isnan(pairsXYZ_2{c})) = [];
    pairsPC_2{c}(isnan(pairsPC_2{c})) = [];
end
pairsXYZ_2 = cell2mat(pairsXYZ_2);
pairsPC_2 = cell2mat(pairsPC_2);

for i = 1:length(bins)-1
    idx = pairsXYZ_2 >= bins(i) & pairsXYZ_2 < bins(i+1);
    if sum(idx) >= minN
        ave(s,i) = mean(pairsPC_2(idx));
        sem(s,i) = std(pairsPC_2(idx))/sqrt(sum(idx));
    end
end
errorbar(binCents,ave(s,:),sem(s,:),'Color','k','LineStyle','none')
plot(binCents,ave(s,:),'Color','k','LineWidth',2)


%% Statistics on previous section - compare subtypes to mismatch
pairsXYZ{6} = pairsXYZ_2;
pairsPC{6} = pairsPC_2;

pvalsVSmismatch = nan(4,length(bins)-1);
n = nan(5,length(bins)-1);
for b = 1:length(bins)-1
    idx2 = pairsXYZ{6} >= bins(b) & pairsXYZ{6} < bins(b+1); %mismatch subtypes
    for i = [1 2 4] % no Raldh, no DAT
        idx1 = pairsXYZ{i} >= bins(b) & pairsXYZ{i} < bins(b+1);
        if sum(idx1) > 5
            pvalsVSmismatch(i,b) = ranksum(pairsPC{6}(idx2),pairsPC{i}(idx1));
        end
        n(i,b) = sum(idx1);
    end
    n(5,b) = sum(idx2);
end
pvalsVSmismatch = pvalsVSmismatch * sum(sum(~isnan(pvalsVSmismatch))); % BONF correction
pvalsVSmismatch = round(pvalsVSmismatch,2,'significant');
pvalsVSmismatch(3,:) = [];
disp('pvals each subtype (v,c,a) vs mismatch')
disp(num2str(pvalsVSmismatch))

% compare all subtypes together to mismatch
pvalsVSmismatch2 = nan(1,length(bins)-1);
for b = 1:length(bins)-1
    idx2 = pairsXYZ{6} >= bins(b) & pairsXYZ{6} < bins(b+1); %mismatch subtypes
    temp = [];
    for i = [1 2 4] % no Raldh, no DAT
        idx1 = pairsXYZ{i} >= bins(b) & pairsXYZ{i} < bins(b+1);
        temp = [temp; pairsPC{i}(idx1)];
    end
    if ~isempty(temp)
        pvalsVSmismatch2(b) = ranksum(pairsPC{6}(idx2),temp);
    end
end
pvalsVSmismatch2 = pvalsVSmismatch2 * sum(sum(~isnan(pvalsVSmismatch2))); % BONF correction
pvalsVSmismatch2 = round(pvalsVSmismatch2,2,'significant');
disp('pvals all subtypes together vs mismatch')
disp(num2str(pvalsVSmismatch2))


end
    






