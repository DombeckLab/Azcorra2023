function plotSignals2_paper(means,sems,RG,VA,neg,pos,hz,exp)
% necesary inputs: data, depth, trig, RG

mainHz = 100;

if isempty(hz)
    hz = mainHz;
    bin = 1;
else
    bin = mainHz/hz;
end

%% plot
load colors.mat C
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

ran = (-neg+1:bin:pos)/mainHz;
yyaxis left
hold off
fill([ran fliplr(ran)], [means(1,:)+sems(1,:) fliplr(means(1,:)-sems(1,:))],C.lightGrey{1},'EdgeColor',C.grey{1})
if strcmp(VA,'V')
    ylabel('Velocity (m/s)')
    ylim([-0.1 0.6])
elseif strcmp(VA,'A')
    ylabel('Acceleration (m/s^2)')
    ylim([-1 1])
else
    ylabel('Licking')
end

yyaxis right
ax = gca;
ax.YColor = colors(2,:);
hold off
if strcmp(RG, 'G')
    fill([ran fliplr(ran)], [means(3,:)+sems(3,:) fliplr(means(3,:)-sems(3,:))],C.lightBlue{1},'EdgeColor',C.blue{1},'LineStyle','-')
    hold on
    fill([ran fliplr(ran)], [means(2,:)+sems(2,:) fliplr(means(2,:)-sems(2,:))],colors(1,:),'EdgeColor',colors(2,:),'LineStyle','-')
else    
    fill([ran fliplr(ran)], [means(5,:)+sems(5,:) fliplr(means(5,:)-sems(5,:))],C.lightBlue{1},'EdgeColor',C.blue{1},'LineStyle','-')
    hold on
    fill([ran fliplr(ran)], [means(4,:)+sems(4,:) fliplr(means(4,:)-sems(4,:))],colors(1,:),'EdgeColor',colors(2,:),'LineStyle','-')
end
ylabel('Norm \DeltaF/F')

xlim([-neg pos]/100)
xticks((floor(-neg/50)*50:50:pos)/mainHz)
xlabel('Time (s)')


    
    
    
    
