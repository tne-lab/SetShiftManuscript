addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))
%% Load Data
load('Data/allBehav.mat')
bn = load('Data/allBehav2.mat');
behavior = rmfield(behavior, 'drug');
behavior=[behavior, bn.behavior(2:5)];
behav_sub = behavior;
behav_sub([31])=[];
behav_sub([15 16 34])=[];
sub_rats=find(~contains([behav_sub.subject],'DN'));
ntrials=0;
for i=1:length(sub_rats)
    for j=1:length(behav_sub(sub_rats(i)).sessions)
        ntrials=ntrials+length(behav_sub(sub_rats(i)).sessions{j});
    end
end
full_tbl = table('Size',[ntrials,8],'VariableTypes',{'double','double','double','logical','double','double','logical','double'});
full_tbl.Properties.VariableNames = {'RT','Rat','Session','Stim','Site','Light','Correct','Choice'};
sites=["mid","dorsomedial","ventral","dorsolateral"];
ind=1;
ses=1;
bad_sessions = [231, 243]
for i=1:length(sub_rats)
    for j=1:length(behav_sub(sub_rats(i)).sessions)
        if (sub_rats(i)~=13||j~=4)
            len=length(behav_sub(sub_rats(i)).sessions{j});
            full_tbl.RT(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.RT];
            full_tbl.Rat(ind:ind+len-1)=i;
            full_tbl.Session(ind:ind+len-1)=ses;
            full_tbl.Stim(ind:ind+len-1)=behav_sub(sub_rats(i)).stimOn(j);
            full_tbl.Site(ind:ind+len-1)=find(strcmp(sites,behav_sub(sub_rats(i)).site{j}));
            full_tbl.Light(ind:ind+len-1)=-([behav_sub(sub_rats(i)).sessions{j}.light]-2);
            full_tbl.Correct(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.performance];
            full_tbl.Choice(ind:ind+len-1)=-([behav_sub(sub_rats(i)).sessions{j}.frontChoice]-2);
            full_tbl.Rule(ind:ind+len-1)=extractBefore([behav_sub(sub_rats(i)).sessions{j}.rule],2);
            full_tbl.Trial(ind:ind+len-1)=1:len;
            full_tbl.RespTime(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.respTime];
            ind=ind+len;
            ses=ses+1;
        end
    end
end
full_tbl.RT(isnan(full_tbl.Choice),:)=-1;
full_tbl.Choice(isnan(full_tbl.Choice),:)=-1;
full_tbl(abs(full_tbl.RT)<=0.25,:)=[];

%% Stats
full_tbl.SL = strcmp(full_tbl.Rule,"L");
full_tbl.sSite = string(full_tbl.Site);
acc_mdl = fitglme(full_tbl,'Correct~1+SL+(1|Rat)','distribution','binomial');
rt_mdl = fitglme(full_tbl(full_tbl.RT>0,:),'RT~1+SL+(1|Rat)','distribution','gamma','link','identity');
rt_mdl = fitglme(full_tbl(full_tbl.RT>0,:),'RT~1+sSite+Stim:sSite+SL+(1|Rat)+(1|Session)','distribution','gamma','link','identity');

%% setup
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% b & c
S = max(full_tbl.Session);
errorsSide = zeros(S,1);
errorsLight = zeros(S,1);

for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    errorsSide(i) = nanmean(1-subtbl.Correct(strcmp(subtbl.Rule,'F')|strcmp(subtbl.Rule,'R')))*100;
    errorsLight(i) = nanmean(1-subtbl.Correct(strcmp(subtbl.Rule,'L')))*100;
end
errorsSide = rmmissing(errorsSide); errorsLight = rmmissing(errorsLight);

subplot_tight(3,12,[5,6])
hold on
handles = al_goodplot2({errorsSide,errorsLight}', 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.526,0.165,0.859;1,0.3,0]);
yticks(0:20:80)
xticks([])
ax = gca;
set(ax,'fontsize',18)
set(ax,'color','none')
set(ax,'xcolor','none')
ylabel('Error Percentage')
ylim([0,80])
set(gca,'linewidth',2)

subplot_tight(3,12,[7,8])
hold on
handles = al_goodplot2({full_tbl.RT(full_tbl.RT>0&(strcmp(full_tbl.Rule,'F')|strcmp(full_tbl.Rule,'R'))),...
    full_tbl.RT(full_tbl.RT>0&strcmp(full_tbl.Rule,'L'))}', 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.526,0.165,0.859;1,0.3,0]);
yticks(0.3:0.3:1.2)
ylim([0.2,1.3])
xticks()
ax = gca;
set(ax,'color','none')
set(ax,'xcolor','none')
set(ax,'fontsize',18)
ylabel('RT (s)')
set(gca,'linewidth',2)

%% c
shifts = nan(S*8,25);
light_rule = false(S*8,1);
for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    shiftInds = [1;find(~strcmp(subtbl.Rule(1:end-1),subtbl.Rule(2:end)))+1;height(subtbl)];
    for j=1:length(shiftInds)-1
        light_rule(8*(i-1)+j) = strcmp(subtbl.Rule(shiftInds(j)),"L");
        shifts(8*(i-1)+j,1:(shiftInds(j+1)-shiftInds(j)))=100*(1-subtbl.Correct(shiftInds(j):shiftInds(j+1)-1));
    end
end

subplot_tight(3,12,9:12)
hold off
ciL = bootci(1000,@(x) nanmean(x),shifts(light_rule,1:25));
ciS = bootci(1000,@(x) nanmean(x),shifts(~light_rule,1:25));
patch('XData',[1:25,flip(1:25)],'YData',[ciL(1,1:25),flip(ciL(2,1:25))],'FaceColor',[1,0.3,0],'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[1:25,flip(1:25)],'YData',[ciS(1,1:25),flip(ciS(2,1:25))],'FaceColor',[0.526,0.165,0.859],'EdgeColor','none','FaceAlpha',0.3)
hold on
plot(nanmean(shifts(light_rule,1:25)),'Color',[1,0.3,0],'LineWidth',3)
plot(nanmean(shifts(~light_rule,1:25)),'Color',[0.526,0.165,0.859],'LineWidth',3)
yticks(20:20:80)
ylabel('Error Percentage')
xlabel('Trials Since Rule Shift')
ylim([20,80])
xlim([1,15])
xticks([1,5,10,15])
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% d
addpath('F:\AtlasPlotter')
coords=readtable('Data/SSProbeCoords.xlsx');
coords2=readtable('Data/FCProbeCoords.xlsx');
coords2.Group=strcat(coords2.Group,"FC");
coords=[coords;coords2];
subplot_tight(3,12,13:16)
slice = Slice(1.4, 'c');
slice.plot()
locs=["midSTR","dmSTR","vSTR","dlSTR","midSTRFC"];
col=[0.949,0.631,0.008;0.102,0.102,0.961;0.792,0,0;0.0118,0.6275,0.3843;0.949,0.8,0.008];

hold on
for i=1:length(locs)
    subtbl = coords(strcmp(coords.Group,locs(i)),:);
    pts = [subtbl.AP_left,subtbl.ML_left,subtbl.DV_left;subtbl.AP_right,subtbl.ML_right,subtbl.DV_right];
    scatter(pts(:,2),-pts(:,3),60,"filled",'MarkerFaceColor',col(i,:),'MarkerEdgeColor',[0.3,0.3,0.3])
end
axis equal
axis off
%% f
subplot_tight(3,12,17:20)
rtChange = zeros(size(sub_rats));
sites = zeros(size(sub_rats));
for i=1:length(rtChange)
    subtbl = full_tbl(full_tbl.Rat==i,:);
    sites(i) = subtbl.Site(i);
    rtChange(i) = 100*mean(subtbl.RT(subtbl.RT>0&subtbl.Stim==1,:))/mean(subtbl.RT(subtbl.RT>0&subtbl.Stim==0,:)) - 100;
end
hold on
rng(623)
handles = al_goodplot2({rtChange(sites==1)',rtChange(sites==2)',rtChange(sites==3)',rtChange(sites==4)'}', 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.631,0.008;0.102,0.102,0.961;0.792,0,0;0.0118,0.6275,0.3843]);
for i=1:4
    scatter(i+(rand(nnz(sites==i),1)-0.5)*0.4,rtChange(sites==i),'filled','k')
end
yline(0,'k--')
ylim([-25,25])
xticks(1:4)
xlim([0.5,4.5])
xticklabels(["Mid", "Dorsomedial", "Ventral", "Dorsolateral"])
yticks(-25:25:25)
ylabel('% Change From Stim OFF')
title('RT')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% g
subplot_tight(3,12,21:24)
errorChange = zeros(size(sub_rats));
sites = zeros(size(sub_rats));
for i=1:length(errorChange)
    subtbl = full_tbl(full_tbl.Rat==i,:);
    noff = length(unique(subtbl.Session(subtbl.Stim==0)));
    non = length(unique(subtbl.Session(subtbl.Stim==1)));
    sites(i) = subtbl.Site(i);
    errorChange(i) = 100*(sum(1-subtbl.Correct(subtbl.Stim==1,:))/non)/(sum(1-subtbl.Correct(subtbl.Stim==0,:))/noff) - 100;
end
hold on
rng(623)
handles = al_goodplot2({errorChange(sites==1)',errorChange(sites==2)',errorChange(sites==3)',errorChange(sites==4)'}', 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.631,0.008;0.102,0.102,0.961;0.792,0,0;0.0118,0.6275,0.3843]);
for i=1:4
    scatter(i+(rand(nnz(sites==i),1)-0.5)*0.4,errorChange(sites==i),'filled','k')
end
yline(0,'k--')
ylim([-50,50])
xticks(1:4)
xlim([0.5,4.5])
xticklabels(["Mid", "Dorsomedial", "Ventral", "Dorsolateral"])
yticks([-50,0,50])
title('Errors')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% i
fc_data = readtable('Data/FCData.csv');
sessions=unique(fc_data.File_Name);
fc_accuracy=zeros(size(sessions));
fc_premature=zeros(size(sessions));
fc_omissions=zeros(size(sessions));
fc_rt=fc_data.Reaction_Time;
fc_on=zeros(size(sessions));
for i=1:length(sessions)
    subtbl = fc_data(strcmp(fc_data.File_Name,sessions{i}),:);
    fc_on(i) = strcmp(subtbl.Stim(1),"on");
    fc_accuracy(i)=nnz(strcmp(subtbl.Accuracy,"Reward"))/(nnz(strcmp(subtbl.Accuracy,"Reward"))+nnz(strcmp(subtbl.Accuracy,"Error")))*100;
    fc_premature(i)=nnz(strcmp(subtbl.Premature_or_not,"Premature"))/height(subtbl)*100;
    fc_omissions(i)=nnz(strcmp(subtbl.Omission_or_not,"Omission"))/height(subtbl)*100;
end
rats = unique(fc_data.Animal);
rt = zeros(size(rats));
pm = zeros(size(rats));
acc = zeros(size(rats));
om = zeros(size(rats));
for i=1:length(rats)
    sub_tbl = fc_data(strcmp(fc_data.Animal,rats(i)),:);
    on = strcmp(sub_tbl.Stim,'on');
    off = strcmp(sub_tbl.Stim,'off');
    pm(i) = nanmean(sub_tbl.Premature_binary(on))/nanmean(sub_tbl.Premature_binary(off))*100 - 100;
    accon = nnz(strcmp(sub_tbl.Accuracy(on),'Reward'))/(nnz(strcmp(sub_tbl.Accuracy(on,:),'Reward'))+nnz(strcmp(sub_tbl.Accuracy(on,:),'Error')));
    accoff = nnz(strcmp(sub_tbl.Accuracy(off),'Reward'))/(nnz(strcmp(sub_tbl.Accuracy(off,:),'Reward'))+nnz(strcmp(sub_tbl.Accuracy(off,:),'Error')));
    acc(i) = accon/accoff*100 - 100;
    rt(i) = nanmean(sub_tbl.Reaction_Time(on))/nanmean(sub_tbl.Reaction_Time(off))*100 - 100;
    om(i) = nanmean(sub_tbl.Omission_binary(on))/nanmean(sub_tbl.Omission_binary(off))*100 - 100;
end

subplot_tight(3,12,31:36)
hold on
rng(626)
al_goodplot2({pm', acc', om', rt'}, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.8,0.008;0.949,0.8,0.008;0.949,0.8,0.008;0.949,0.8,0.008])
scatter(1+(rand(size(rats))-0.5)*0.4,pm,'filled','k')
scatter(2+(rand(size(rats))-0.5)*0.4,acc,'filled','k')
scatter(4+(rand(size(rats))-0.5)*0.4,rt,'filled','k')
scatter(3+(rand(size(rats))-0.5)*0.4,om,'filled','k')
xticks(1:4)
xticklabels(["Premature Responses", "Accuracy", "Omissions", "Reaction Time"])
set(gca,'fontsize',18)
yticks(-40:40:120)
ylim([-40,120])
xlim([0.5,4.5])
ylabel('% Change From Stim Off')
%% Supplement Main
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% S1a & b
S = max(full_tbl.Session);
errorsLeft = zeros(S,1);
errorsRight = zeros(S,1);

for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    errorsLeft(i) = nanmean(1-subtbl.Correct(strcmp(subtbl.Rule,'F')))*100;
    errorsRight(i) = nanmean(1-subtbl.Correct(strcmp(subtbl.Rule,'R')))*100;
end
errorsLeft = rmmissing(errorsLeft); errorsRight = rmmissing(errorsRight);

subplot_tight(3,12,[1,2])
hold on
handles = al_goodplot2({errorsLeft,errorsRight}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.31,0.73,0.90;0.83,0.07,0.82]);
yticks(0:20:80)
xticks([])
ax = gca;
set(ax,'fontsize',18)
set(ax,'color','none')
set(ax,'xcolor','none')
ylabel('Error Percentage')
ylim([0,80])
set(gca,'linewidth',2)

subplot_tight(3,12,[3,4])
hold on
handles = al_goodplot2({full_tbl.RT(full_tbl.RT>0&(strcmp(full_tbl.Rule,'F'))),...
    full_tbl.RT(full_tbl.RT>0&strcmp(full_tbl.Rule,'R'))}', 'pos', [1,1], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.31,0.73,0.90;0.83,0.07,0.82]);
yticks(0.3:0.3:1.2)
ylim([0.2,1.3])
xticks()
ax = gca;
set(ax,'color','none')
set(ax,'xcolor','none')
set(ax,'fontsize',18)
ylabel('RT (s)')
set(gca,'linewidth',2)

%% S1c
subplot_tight(3,12,5:8)
accChange = zeros(size(sub_rats));
sites = zeros(size(sub_rats));
for i=1:length(rtChange)
    subtbl = full_tbl(full_tbl.Rat==i,:);
    sites(i) = subtbl.Site(i);
    accChange(i) = 100*mean(subtbl.Correct(subtbl.Stim==1,:))/mean(subtbl.Correct(subtbl.Stim==0,:)) - 100;
end
hold on
rng(623)
handles = al_goodplot2({accChange(sites==1)',accChange(sites==2)',accChange(sites==3)',accChange(sites==4)'}', 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.631,0.008;0.102,0.102,0.961;0.792,0,0;0.0118,0.6275,0.3843]);
for i=1:4
    scatter(i+(rand(nnz(sites==i),1)-0.5)*0.4,accChange(sites==i),'filled','k')
end
yline(0,'k--')
ylim([-10,10])
xticks(1:4)
xlim([0.5,4.5])
xticklabels(["Mid", "Dorsomedial", "Ventral", "Dorsolateral"])
yticks([-10,0,10])
ylabel('% Change from Stim OFF')
title('Accuracy')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% S1d
subplot_tight(3,12,9:12)
ttfChange = zeros(size(sub_rats));
sites = zeros(size(sub_rats));
for i=1:length(rtChange)
    subtbl = full_tbl(full_tbl.Rat==i,:);
    sites(i) = subtbl.Site(i);
    on_ttf = 0;
    on_ses = unique(subtbl.Session(subtbl.Stim==1));
    for j=on_ses'
        on_ttf = on_ttf + max(subtbl.RespTime(subtbl.Session==j));
    end
    on_ttf = on_ttf / length(on_ses);
    off_ttf = 0;
    off_ses = unique(subtbl.Session(subtbl.Stim==0));
    for j=off_ses'
        off_ttf = off_ttf + max(subtbl.RespTime(subtbl.Session==j));
    end
    off_ttf = off_ttf / length(off_ses);
    ttfChange(i) = 100*on_ttf/off_ttf - 100;
end
hold on
rng(623)
handles = al_goodplot2({ttfChange(sites==1)',ttfChange(sites==2)',ttfChange(sites==3)',ttfChange(sites==4)'}', 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.631,0.008;0.102,0.102,0.961;0.792,0,0;0.0118,0.6275,0.3843]);
for i=1:4
    scatter(i+(rand(nnz(sites==i),1)-0.5)*0.4,ttfChange(sites==i),'filled','k')
end
yline(0,'k--')
ylim([-30,30])
xticks(1:4)
xlim([0.5,4.5])
xticklabels(["Mid", "Dorsomedial", "Ventral", "Dorsolateral"])
yticks([-30,0,30])
ylabel('% Change from Stim OFF')
title('Time to Finish')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% S1e
rats = unique(full_tbl.Rat(full_tbl.Site==1));
subplot_tight(3,12,15:18)
hold on
for rat=rats'
    on_rt = mean(full_tbl.RT(full_tbl.Rat==rat&full_tbl.Stim==1&full_tbl.RT>0));
    off_rt = mean(full_tbl.RT(full_tbl.Rat==rat&full_tbl.Stim==0&full_tbl.RT>0));
    plot([1,2],[off_rt,on_rt],'-ok','LineWidth',2,'MarkerFaceColor','k')
end
ylim([0.25,1])
yticks(0.25:0.25:1)
xticks([1,2])
xticklabels(["OFF", "ON"])
xlim([0.6,2.4])
ylabel('Mean RT (s)')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
%% S1f
S = max(full_tbl.Session);
rts=zeros(S,1);
stims=zeros(S,1);
sites=zeros(S,1);

for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    rts(i) = median(subtbl.RT(subtbl.RT>0))/median(full_tbl.RT(full_tbl.Rat==subtbl.Rat(1)&full_tbl.RT>0&~full_tbl.Stim));
    stims(i) = subtbl.Stim(1);
    sites(i) = subtbl.Site(1);
end

subplot_tight(3,12,19:22)
hold on
on = rts(stims==1&sites==1)*100;
off = rts(stims==0&sites==1)*100;
al_goodplot2({off,on}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.3,0.3,0.3;0.949,0.631,0.008]);
scatter(0.99-0.19*rand(size(off)),off,'filled','MarkerFaceColor',[0.3,0.3,0.3])
scatter(1.01+0.19*rand(size(on)),on,'filled','MarkerFaceColor',[0.949,0.631,0.008])
xticks([])
ylabel('% Change Relative to Sham Median')
legend(["Sham","Mid"])
ax = gca;
set(ax,'color','none')
set(ax,'xcolor','none')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% SS Raw Data Plots
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

all_sessions = unique(full_tbl.Session);
ses_tbl = table('Size', [length(all_sessions),6], 'VariableNames',["Rat","Site", "EC", "Acc", "TTF", "Stim"], 'VariableTypes',["double","double", "double", "double", "double" "logical"]);
for i=1:length(all_sessions)
    sub_tbl = full_tbl(full_tbl.Session==i,:);
    ses_tbl.Rat(i) = sub_tbl.Rat(1);
    ses_tbl.Site(i) = sub_tbl.Site(1);
    ses_tbl.Stim(i) = sub_tbl.Stim(1);
    ses_tbl.EC(i) = sum(1-sub_tbl.Correct);
    ses_tbl.Acc(i) = mean(sub_tbl.Correct)*100;
    ses_tbl.TTF(i) = max(sub_tbl.RespTime);
end

rtOn = zeros(size(sub_rats));
rtOff = zeros(size(sub_rats));
ecOn = zeros(size(sub_rats));
ecOff = zeros(size(sub_rats));
accOn = zeros(size(sub_rats));
accOff = zeros(size(sub_rats));
ttfOn = zeros(size(sub_rats));
ttfOff = zeros(size(sub_rats));
sites = zeros(size(sub_rats));
for i=1:length(sub_rats)
    subtbl = full_tbl(full_tbl.Rat==i,:);
    sites(i) = subtbl.Site(1);
    nsesOn = length(unique(subtbl.Session(subtbl.Stim==1)));
    nsesOff = length(unique(subtbl.Session(subtbl.Stim==0)));
    rtOn(i) = mean(subtbl.RT(subtbl.Stim==1&subtbl.RT>0));
    rtOff(i) = mean(subtbl.RT(subtbl.Stim==0&subtbl.RT>0));
    ecOn(i) = sum(1-subtbl.Correct(subtbl.Stim==1)) / nsesOn;
    ecOff(i) = sum(1-subtbl.Correct(subtbl.Stim==0)) / nsesOff;
    accOn(i) = mean(subtbl.Correct(subtbl.Stim==1)) * 100;
    accOff(i) = mean(subtbl.Correct(subtbl.Stim==0)) * 100;
    ttfOn(i) = mean(ses_tbl.TTF(ses_tbl.Rat==i&ses_tbl.Stim==1));
    ttfOff(i) = mean(ses_tbl.TTF(ses_tbl.Rat==i&ses_tbl.Stim==0));
end

%% a
subplot_tight(4,1,1)
hold on
handles = al_goodplot2({full_tbl.RT(full_tbl.Site==1&full_tbl.Stim==0&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==1&full_tbl.Stim==1&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==2&full_tbl.Stim==0&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==2&full_tbl.Stim==1&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==3&full_tbl.Stim==0&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==3&full_tbl.Stim==1&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==4&full_tbl.Stim==0&full_tbl.RT>0),...
                        full_tbl.RT(full_tbl.Site==4&full_tbl.Stim==1&full_tbl.RT>0)}',...
                        'pos', [1,1,2,2,3,3,4,4], 'type', {'left','right','left','right','left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.631,0.008;0.4902,0.4902,0.6588;0.102,0.102,0.961;0.5216,0.2980,0.2980;0.792,0,0;0.2235,0.4510,0.3608;0.0118,0.6275,0.3843]);
subs = [rtOff', rtOn', nan(size(rtOn))']';
for i=1:4
    x = (i + [-0.1*ones(nnz(sites==i), 1), 0.1*ones(nnz(sites==i), 1), zeros(nnz(sites==i), 1)])';
    y = subs(:, sites==i);
    if i==1
        plot(x(:), y(:), 'Color',[0.6,0.6,0.6])
    else
        plot(x(:), y(:), 'Color',[0.8,0.8,0.8])
    end
end
xticks(1:4)
xlim([0.5,4.5])
xticklabels([])
ylim([0.25,1.5])
yticks(0.5:0.5:1.5)
ylabel('RT (s)')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% b
subplot_tight(4,1,2)
handles = al_goodplot2({ses_tbl.EC(ses_tbl.Site==1&ses_tbl.Stim==0),...
                        ses_tbl.EC(ses_tbl.Site==1&ses_tbl.Stim==1),...
                        ses_tbl.EC(ses_tbl.Site==2&ses_tbl.Stim==0),...
                        ses_tbl.EC(ses_tbl.Site==2&ses_tbl.Stim==1),...
                        ses_tbl.EC(ses_tbl.Site==3&ses_tbl.Stim==0),...
                        ses_tbl.EC(ses_tbl.Site==3&ses_tbl.Stim==1),...
                        ses_tbl.EC(ses_tbl.Site==4&ses_tbl.Stim==0),...
                        ses_tbl.EC(ses_tbl.Site==4&ses_tbl.Stim==1)}',...
                        'pos', [1,1,2,2,3,3,4,4], 'type', {'left','right','left','right','left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.631,0.008;0.4902,0.4902,0.6588;0.102,0.102,0.961;0.5216,0.2980,0.2980;0.792,0,0;0.2235,0.4510,0.3608;0.0118,0.6275,0.3843]);
subs = [ecOff', ecOn', nan(size(ecOn))']';
for i=1:4
    x = (i + [-0.1*ones(nnz(sites==i), 1), 0.1*ones(nnz(sites==i), 1), zeros(nnz(sites==i), 1)])';
    y = subs(:, sites==i);
    if i==1
        plot(x(:), y(:), 'Color',[0.6,0.6,0.6])
    else
        plot(x(:), y(:), 'Color',[0.8,0.8,0.8])
    end
end
ylim([0 120])
xticks(1:4)
xlim([0.5,4.5])
xticklabels([])
% yticks([50,100,150])
ylabel('Errors')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% c
subplot_tight(4,1,3)
handles = al_goodplot2({ses_tbl.Acc(ses_tbl.Site==1&ses_tbl.Stim==0),...
                        ses_tbl.Acc(ses_tbl.Site==1&ses_tbl.Stim==1),...
                        ses_tbl.Acc(ses_tbl.Site==2&ses_tbl.Stim==0),...
                        ses_tbl.Acc(ses_tbl.Site==2&ses_tbl.Stim==1),...
                        ses_tbl.Acc(ses_tbl.Site==3&ses_tbl.Stim==0),...
                        ses_tbl.Acc(ses_tbl.Site==3&ses_tbl.Stim==1),...
                        ses_tbl.Acc(ses_tbl.Site==4&ses_tbl.Stim==0),...
                        ses_tbl.Acc(ses_tbl.Site==4&ses_tbl.Stim==1)}',...
                        'pos', [1,1,2,2,3,3,4,4], 'type', {'left','right','left','right','left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.631,0.008;0.4902,0.4902,0.6588;0.102,0.102,0.961;0.5216,0.2980,0.2980;0.792,0,0;0.2235,0.4510,0.3608;0.0118,0.6275,0.3843]);
subs = [accOff', accOn', nan(size(accOn))']';
for i=1:4
    x = (i + [-0.1*ones(nnz(sites==i), 1), 0.1*ones(nnz(sites==i), 1), zeros(nnz(sites==i), 1)])';
    y = subs(:, sites==i);
    if i==1
        plot(x(:), y(:), 'Color',[0.6,0.6,0.6])
    else
        plot(x(:), y(:), 'Color',[0.8,0.8,0.8])
    end
end
% ylim([50,80])
% yticks(50:10:80)
xticks(1:4)
xlim([0.5,4.5])
xticklabels([])
% yticks([50,100,150])
ylim([40 80])
ylabel('Accuracy')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% d
subplot_tight(4,1,4)
handles = al_goodplot2({ses_tbl.TTF(ses_tbl.Site==1&ses_tbl.Stim==0),...
                        ses_tbl.TTF(ses_tbl.Site==1&ses_tbl.Stim==1),...
                        ses_tbl.TTF(ses_tbl.Site==2&ses_tbl.Stim==0),...
                        ses_tbl.TTF(ses_tbl.Site==2&ses_tbl.Stim==1),...
                        ses_tbl.TTF(ses_tbl.Site==3&ses_tbl.Stim==0),...
                        ses_tbl.TTF(ses_tbl.Site==3&ses_tbl.Stim==1),...
                        ses_tbl.TTF(ses_tbl.Site==4&ses_tbl.Stim==0),...
                        ses_tbl.TTF(ses_tbl.Site==4&ses_tbl.Stim==1)}',...
                        'pos', [1,1,2,2,3,3,4,4], 'type', {'left','right','left','right','left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.631,0.008;0.4902,0.4902,0.6588;0.102,0.102,0.961;0.5216,0.2980,0.2980;0.792,0,0;0.2235,0.4510,0.3608;0.0118,0.6275,0.3843]);
subs = [ttfOff', ttfOn', nan(size(ttfOn))']';
for i=1:4
    x = (i + [-0.1*ones(nnz(sites==i), 1), 0.1*ones(nnz(sites==i), 1), zeros(nnz(sites==i), 1)])';
    y = subs(:, sites==i);
    if i==1
        plot(x(:), y(:), 'Color',[0.6,0.6,0.6])
    else
        plot(x(:), y(:), 'Color',[0.8,0.8,0.8])
    end
end
% ylim([50,80])
% yticks(50:10:80)
xticks(1:4)
xlim([0.5,4.5])
xticklabels(["Mid", "Dorsomedial", "Ventral", "Dorsolateral"])
% yticks([50,100,150])
ylim([0.66 3] * 1e3)
ylabel('Time to Finish')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% Five Choice Raw Data Plots
fc_data = readtable('Data/FCData.csv');
sessions=unique(fc_data.File_Name);
fc_accuracy=zeros(size(sessions));
fc_premature=zeros(size(sessions));
fc_omissions=zeros(size(sessions));
fc_rt=fc_data.Reaction_Time;
fc_on=false(size(sessions));
for i=1:length(sessions)
    subtbl = fc_data(strcmp(fc_data.File_Name,sessions{i}),:);
    fc_on(i) = strcmp(subtbl.Stim(1),"on");
    fc_accuracy(i)=(1-nnz(strcmp(subtbl.Accuracy,"Reward"))/(nnz(strcmp(subtbl.Accuracy,"Reward"))+nnz(strcmp(subtbl.Accuracy,"Error"))))*100;
    fc_premature(i)=nnz(strcmp(subtbl.Premature_or_not,"Premature"))/height(subtbl)*100;
    fc_omissions(i)=nnz(strcmp(subtbl.Omission_or_not,"Omission"))/height(subtbl)*100;
end
rats = unique(fc_data.Animal);
rt_off = zeros(size(rats));
rt_on = zeros(size(rats));
pm_off = zeros(size(rats));
pm_on = zeros(size(rats));
acc_off = zeros(size(rats));
acc_on = zeros(size(rats));
om_off = zeros(size(rats));
om_on = zeros(size(rats));
for i=1:length(rats)
    sub_tbl = fc_data(strcmp(fc_data.Animal,rats(i)),:);
    on = strcmp(sub_tbl.Stim,'on');
    off = strcmp(sub_tbl.Stim,'off');
    pm_off(i) = nanmean(sub_tbl.Premature_binary(off))*100;
    pm_on(i) = nanmean(sub_tbl.Premature_binary(on))*100;
    acc_on(i) = (1-nnz(strcmp(sub_tbl.Accuracy(on),'Reward'))/(nnz(strcmp(sub_tbl.Accuracy(on,:),'Reward'))+nnz(strcmp(sub_tbl.Accuracy(on,:),'Error'))))*100;
    acc_off(i) = (1-nnz(strcmp(sub_tbl.Accuracy(off),'Reward'))/(nnz(strcmp(sub_tbl.Accuracy(off,:),'Reward'))+nnz(strcmp(sub_tbl.Accuracy(off,:),'Error'))))*100;
    rt_off(i) = nanmean(sub_tbl.Reaction_Time(off&sub_tbl.Reaction_Time>0&sub_tbl.Reaction_Time<5));
    rt_on(i) = nanmean(sub_tbl.Reaction_Time(on&sub_tbl.Reaction_Time>0&sub_tbl.Reaction_Time<5));
    om_off(i) = nanmean(sub_tbl.Omission_binary(off))*100;
    om_on(i) = nanmean(sub_tbl.Omission_binary(on))*100;
end

figure
subplot_tight(1,4,1)
hold on
al_goodplot2({fc_premature(~fc_on)', fc_premature(fc_on)'}, 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008])
xticks([])
set(gca,'fontsize',18)
ylabel('Premature Percentage')
x = 1+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [pm_off, pm_on, nan(size(pm_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])
subplot_tight(1,4,2)
hold on
al_goodplot2({fc_accuracy(~fc_on)', fc_accuracy(fc_on)'}, 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008])
xticks([])
set(gca,'fontsize',18)
ylabel('Error Percentage')
x = 1+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [acc_off, acc_on, nan(size(acc_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])
subplot_tight(1,4,3)
hold on
al_goodplot2({fc_omissions(~fc_on)', fc_omissions(fc_on)'}, 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008])
xticks([])
set(gca,'fontsize',18)
ylabel('Omission Percentage')
x = 1+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [om_off, om_on, nan(size(om_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])
subplot_tight(1,4,4)
hold on
al_goodplot2({fc_data.Reaction_Time(fc_data.Reaction_Time>0&~isnan(fc_data.Reaction_Time)&strcmp(fc_data.Stim, 'off')&fc_data.Reaction_Time<5)',...
    fc_data.Reaction_Time(fc_data.Reaction_Time>0&~isnan(fc_data.Reaction_Time)&strcmp(fc_data.Stim, 'on')&fc_data.Reaction_Time<5)'}, 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008])
xticks([])
set(gca,'fontsize',18)
ylabel('RT (s)')
ylim([0, 1.5])
x = 1+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [rt_off, rt_on, nan(size(rt_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])

%% S3 (20 vs 130)

freq_tbl = readtable('Data/SS20_130.csv');
subjects = unique(freq_tbl.Subject);
rt_effect = zeros(length(subjects),2);
errors_effect = zeros(length(subjects),2);
accuracy_effect = zeros(length(subjects),2);
ttf_effect = zeros(length(subjects),2);
groups = zeros(size(ttf_effect)) + [0,1];

for i=1:length(subjects)
    subtbl = freq_tbl(strcmp(freq_tbl.Subject,subjects(i)),:);
    rt_off = median(subtbl.RT(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==130&subtbl.RT<3));
    rt_effect(i,1) = median(subtbl.RT(subtbl.Stim_frequency==130&subtbl.RT<3))/rt_off;
    rt_off = median(subtbl.RT(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==20&subtbl.RT<3));
    rt_effect(i,2) = median(subtbl.RT(subtbl.Stim_frequency==20&subtbl.RT<3))/rt_off;
    errors_off = sum(1-subtbl.Correctness(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==130)) / length(unique(subtbl.Session(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==130)));
    errors_effect(i,1) = sum(1-subtbl.Correctness(subtbl.Stim_frequency==130)) / length(unique(subtbl.Session(subtbl.Stim_frequency==130))) / errors_off;
    errors_off = sum(1-subtbl.Correctness(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==20)) / length(unique(subtbl.Session(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==20)));
    errors_effect(i,2) = sum(1-subtbl.Correctness(subtbl.Stim_frequency==20)) / length(unique(subtbl.Session(subtbl.Stim_frequency==20)))  / errors_off;
    accuracy_off = mean(subtbl.Correctness(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==130));
    accuracy_effect(i,1) = mean(subtbl.Correctness(subtbl.Stim_frequency==130)) / accuracy_off;
    accuracy_off = mean(subtbl.Correctness(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==20));
    accuracy_effect(i,2) = mean(subtbl.Correctness(subtbl.Stim_frequency==20)) / accuracy_off;
    sessions = unique(subtbl.Session(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==130));
    ttf_off = 0;
    for j=1:length(sessions)
        ttf_off = ttf_off + max(subtbl.Time_in_sec(subtbl.Session==sessions(j)));
    end
    ttf_off = ttf_off / length(sessions);
    sessions = unique(subtbl.Session(subtbl.Stim_frequency==130));
    for j=1:length(sessions)
        ttf_effect(i,1) = ttf_effect(i,1) + max(subtbl.Time_in_sec(subtbl.Session==sessions(j)));
    end
    ttf_effect(i,1) = ttf_effect(i,1) / length(sessions) / ttf_off;
    sessions = unique(subtbl.Session(subtbl.Stim_frequency==0&subtbl.Stim_frequency_group==20));
    ttf_off = 0;
    for j=1:length(sessions)
        ttf_off = ttf_off + max(subtbl.Time_in_sec(subtbl.Session==sessions(j)));
    end
    ttf_off = ttf_off / length(sessions);
    sessions = unique(subtbl.Session(subtbl.Stim_frequency==20));
    for j=1:length(sessions)
        ttf_effect(i,2) = ttf_effect(i,2) + max(subtbl.Time_in_sec(subtbl.Session==sessions(j)));
    end
    ttf_effect(i,2) = ttf_effect(i,2) / length(sessions) / ttf_off;
end
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

subplot_tight(3,12,1:12)
hold on
rng(626)
boxchart(groups(:),rt_effect(:)*100 - 100, 'JitterOutliers', 'on', 'BoxFaceColor', [0.949,0.631,0.008])
scatter((rand(size(rt_effect(:))) - 0.5)*0.4 + groups(:), rt_effect(:)*100 - 100, 'filled','MarkerFaceColor','k')
boxchart(groups(:)+2.5,errors_effect(:)*100 - 100, 'JitterOutliers', 'on', 'BoxFaceColor', [0.949,0.631,0.008])
scatter((rand(size(errors_effect(:))) - 0.5)*0.4 + 2.5 + groups(:), errors_effect(:)*100 - 100, 'filled','MarkerFaceColor','k')
boxchart(groups(:)+5,accuracy_effect(:)*100 - 100, 'JitterOutliers', 'on', 'BoxFaceColor', [0.949,0.631,0.008])
scatter((rand(size(accuracy_effect(:))) - 0.5)*0.4 + 5 + groups(:), accuracy_effect(:)*100 - 100, 'filled','MarkerFaceColor','k')
boxchart(groups(:)+7.5,ttf_effect(:)*100 - 100, 'JitterOutliers', 'on', 'BoxFaceColor', [0.949,0.631,0.008])
scatter((rand(size(ttf_effect(:))) - 0.5)*0.4 + 7.5 + groups(:), ttf_effect(:)*100 - 100, 'filled','MarkerFaceColor','k')
yline(0,'--')
xline(1.75,':')
xline(4.25,':')
xline(6.75,':')
ylabel('% of Stim OFF')
xticks([0,1,2.5,3.5,5,6,7.5,8.5])
xticklabels(["130 Hz", "20 Hz", "130 Hz", "20 Hz", "130 Hz", "20 Hz", "130 Hz", "20 Hz"])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
ylim([-40,100])
yticks(-40:40:100)

%% Sex Differences Plots
MF_coords=readtable('Data/MF_coords.xlsx');
load('Data/MF_Data.mat')
rats = unique(MF_tbl.ratto);
rt_off = zeros(size(rats));
rt_on = zeros(size(rats));
errors_off = zeros(size(rats));
errors_on = zeros(size(rats));
sex = zeros(size(rats));
for i=1:length(rats)
    sub_tbl = MF_tbl(strcmp(MF_tbl.ratto,rats(i)),:);
    on = strcmp(sub_tbl.Stim,'on');
    off = strcmp(sub_tbl.Stim,'off');
    rt_off(i) = median(sub_tbl.RT(off));
    rt_on(i) = median(sub_tbl.RT(on));
    errors_off(i) = sum(1-sub_tbl.performance(off))/length(unique(sub_tbl.subses(off)));
    errors_on(i) = sum(1-sub_tbl.performance(on))/length(unique(sub_tbl.subses(on)));
    sex(i) = strcmp(sub_tbl.sex(1),"F");
end

sessions=unique(MF_tbl.subses);
sx_errors=zeros(size(sessions));
sx_on=false(size(sessions));
sx_sex=zeros(size(sessions));
sx_rat=strings(size(sessions));
for i=1:length(sessions)
    subtbl = MF_tbl(strcmp(MF_tbl.subses,sessions{i}),:);
    sx_on(i) = strcmp(subtbl.Stim(1),"on");
    sx_errors(i)=sum(1-subtbl.performance);
    sx_sex(i)=strcmp(subtbl.sex(1),"F");
    sx_rat(i)=subtbl.ratto(1);
end

err_tbl=table('Size',[length(sessions),4], 'VariableNames',["Stim","Rat","sex","errors"],'VariableTypes',{'logical','string','logical','double'});
err_tbl.Stim=sx_on;
err_tbl.Rat=sx_rat;
err_tbl.sex=sx_sex;
err_tbl.errors=sx_errors;
figure('Renderer', 'painters', 'Position', [100 1 1080 300])
subplot_tight(1,3,1)
addpath('F:\AtlasPlotter\src')
slice = Slice(1.4, 'c');
slice.plot()
hold on
pts = [MF_coords.AP_left,MF_coords.ML_left,MF_coords.DV_left;MF_coords.AP_right,MF_coords.ML_right,MF_coords.DV_right];
scatter(pts(:,2),pts(:,3),60,"filled",'MarkerFaceColor',[0.949,0.631,0.008],'MarkerEdgeColor',[0.3,0.3,0.3])
axis equal
axis off

subplot_tight(1,3,2)
hold on
al_goodplot2({MF_tbl.RT(strcmp(MF_tbl.Stim,'off')&strcmp(MF_tbl.sex,"M"))',...
              MF_tbl.RT(strcmp(MF_tbl.Stim,'on')&strcmp(MF_tbl.sex,"M"))',...
              MF_tbl.RT(strcmp(MF_tbl.Stim,'off')&strcmp(MF_tbl.sex,"F"))',...
              MF_tbl.RT(strcmp(MF_tbl.Stim,'on')&strcmp(MF_tbl.sex,"F"))'}, 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008;0.9294,0.8431,0.6706;0.949,0.8,0.008],'useMedian',true)
xticks([])
x = 1+repmat(sex',3,1)+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [rt_off, rt_on, nan(size(rt_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])
xticks([1,2])
xticklabels(["Males", "Females"])
ylabel("RT (s)")
ylim([0,1.5])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
xlim([0.5,2.5])
subplot_tight(1,3,3)
hold on
al_goodplot2({sx_errors(~sx_sex&~sx_on)',...
              sx_errors(~sx_sex&sx_on)',...
              sx_errors(sx_sex&~sx_on)',...
              sx_errors(sx_sex&sx_on)'}, 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.9294,0.8431,0.6706;0.949,0.8,0.008;0.9294,0.8431,0.6706;0.949,0.8,0.008],'useMedian',true)
xticks([])
x = 1+repmat(sex',3,1)+([-0.1*ones(length(rats), 1), 0.1*ones(length(rats), 1), zeros(length(rats), 1)])';
subs = [errors_off, errors_on, nan(size(errors_on))]';
plot(x(:), subs(:), 'Color',[0.6,0.6,0.6])
xticks([1,2])
xticklabels(["Males", "Females"])
ylabel("Errors")
set(gca,'fontsize',18)
set(gca,'linewidth',2)
ylim([10,90])
yticks(10:20:90)
xlim([0.5,2.5])