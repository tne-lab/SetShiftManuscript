%% Load data
motor_data = readtable('OF_data.xlsx');
load('allBehav.mat')
bn = load('allBehav2.mat');
behavior = rmfield(behavior, 'drug');
behavior=[behavior, bn.behavior(2:5)];
behav_sub = behavior;
behav_sub([31])=[];
behav_sub([15 16 34])=[];
valid_rats=["DN1" "S16" "S17" "S18" "S19" "S20" "S23" "S24" "S62" "DN2" "DN3" "DN4" "DN5"];
sub_rats=find(contains([behav_sub.subject],valid_rats));
ntrials=0;
for i=1:length(sub_rats)
    for j=1:length(behav_sub(sub_rats(i)).sessions)
        ntrials=ntrials+length(behav_sub(sub_rats(i)).sessions{j});
    end
end
full_tbl = table('Size',[ntrials,8],'VariableTypes',{'double','double','double','double','double','double','double','logical'});
full_tbl.Properties.VariableNames = {'RT','Rat','DT','Choice','immobility','distance','speed','stim'};
sites=["mid","dorsomedial","ventral","dorsolateral"];
ind=1;
ses=1;
for i=1:length(sub_rats)
    for j=1:length(behav_sub(sub_rats(i)).sessions)
        if (sub_rats(i)~=13||j~=4)
            sub_motor = motor_data(strcmp(motor_data.Subject,valid_rats(i)),:);
            len=length(behav_sub(sub_rats(i)).sessions{j});
            full_tbl.RT(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.RT];
            full_tbl.Rat(ind:ind+len-1)=i;
            full_tbl.Choice(ind:ind+len-1)=-([behav_sub(sub_rats(i)).sessions{j}.frontChoice]-2);
            full_tbl.DT(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.initTime]-[behav_sub(sub_rats(i)).sessions{j}.cueTime];
            full_tbl.immobility(ind:ind+len-1)=sub_motor.middle_body_time_stationary(j);
            full_tbl.distance(ind:ind+len-1)=sub_motor.traveled_distance_cm(j);
            full_tbl.speed(ind:ind+len-1)=sub_motor.speed_moving_cm_s(j);
            full_tbl.Stim(ind:ind+len-1)=behav_sub(sub_rats(i)).stimOn(j);
            ind=ind+len;
            ses=ses+1;
        end
    end
end
full_tbl.RT(isnan(full_tbl.Choice),:)=-1;
full_tbl.Choice(isnan(full_tbl.Choice),:)=-1;
full_tbl(full_tbl.RT<=0.25,:)=[];
%% Setup figure
figure('Renderer', 'painters', 'Position', [100 100 800 900])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
%% a
subplot_tight(3,2,1)
sub_tbl = full_tbl(~isnan(full_tbl.speed),:);
immmdl = fitglme(sub_tbl, 'RT ~ 1 + immobility + (1|Rat)', 'Distribution', 'gamma', 'link', 'identity');
im = 0:5:600;
rt = -0.5:0.01:0.5;
temp3 = normalized_ks(immmdl,sub_tbl.immobility,sub_tbl.RT,im,rt);
pcolor(im,rt,temp3)
hold on
plot(im, immmdl.Coefficients.Estimate(2) * im, 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
shading interp
colormap(viridis)
xlabel('Immobility (s)')
ylabel('Change in RT from intercept (s)')
caxis([0,0.05])
xticks([0,300,600])
set(gca,'fontsize',18)
yticks([-0.5,0,0.5])

%% b
subplot_tight(3,2,2)
sub_tbl = full_tbl(~isnan(full_tbl.distance),:);
dsmdl = fitglme(sub_tbl, 'RT ~ 1 + distance + (1|Rat)', 'Distribution', 'gamma', 'link', 'identity');
ds = 1000:5:2000;
rt = -0.5:0.01:0.5;
temp3 = normalized_ks(dsmdl,sub_tbl.distance,sub_tbl.RT,ds,rt);
pcolor(ds,rt,temp3)
hold on
plot(ds, dsmdl.Coefficients.Estimate(2) * ds, 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
shading interp
colormap(viridis)
xlabel('Distance Travelled (cm)')
yticks([])
caxis([0,0.035])
xticks([1000,1500,2000])
set(gca,'fontsize',18)
%% c
subplot_tight(3,2,3)
sub_tbl = full_tbl(~isnan(full_tbl.speed),:);
spmdl = fitglme(sub_tbl, 'RT ~ 1 + speed + (1|Rat)', 'Distribution', 'gamma', 'link', 'identity');
sp = 0.3:0.01:0.9;
rt = -0.5:0.01:0.5;
temp3 = normalized_ks(spmdl,sub_tbl.speed,sub_tbl.RT,sp,rt);
pcolor(sp,rt,temp3)
hold on
plot(sp, spmdl.Coefficients.Estimate(2) * sp, 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
shading interp
colormap(viridis)
xlabel('Speed (cm/s)')
ylabel('Change in RT from intercept (s)')
yticks([-0.5,0,0.5])
caxis([0,0.035])
xticks([0.3,0.6,0.9])
set(gca,'fontsize',18)
%% d
sub_tbl = full_tbl(~isnan(full_tbl.DT),:);
itimdl = fitglme(sub_tbl, 'RT ~ 1 + DT + (1|Rat)', 'Distribution', 'gamma', 'link', 'identity');
dt = 0:0.1:16;
rt = -0.5:0.01:0.5;
temp3 = normalized_ks(itimdl,sub_tbl.DT,sub_tbl.RT,dt,rt);
subplot_tight(3,2,4)
pcolor(dt,rt,temp3)
hold on
plot(dt, itimdl.Coefficients.Estimate(2) * dt, 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
shading interp
colormap(viridis)
xlabel('Delay Time (s)')
set(gca,'fontsize',18)
caxis([0,0.035])
yticks([])
xticks([0,8,16])
%% e
subjects = unique(motor_data.Subject);
changes = zeros(length(subjects),4);
for i=1:length(subjects)
    sub_tbl = motor_data(strcmp(subjects(i), motor_data.Subject),:);
    sub_tbl2 = full_tbl(find(strcmp(subjects(i),valid_rats)) == full_tbl.Rat,:);
    changes(i,1) = 100 * (mean(sub_tbl.middle_body_time_stationary(strcmp(sub_tbl.Stim,'ON'))) / ...
        mean(sub_tbl.middle_body_time_stationary(strcmp(sub_tbl.Stim,'OFF'))));
    changes(i,2) = 100 * (mean(sub_tbl.traveled_distance_cm(strcmp(sub_tbl.Stim,'ON'))) / ...
        mean(sub_tbl.traveled_distance_cm(strcmp(sub_tbl.Stim,'OFF'))));
    changes(i,3) = 100 * (mean(sub_tbl.speed_moving_cm_s(strcmp(sub_tbl.Stim,'ON'))) / ...
        mean(sub_tbl.speed_moving_cm_s(strcmp(sub_tbl.Stim,'OFF'))));
    changes(i,4) = 100 * (mean(sub_tbl2.DT(sub_tbl2.Stim)) / ...
        mean(sub_tbl2.DT(~sub_tbl2.Stim)));
end

changes2 =cell(4,1);
for i=1:4
    changes2{i} = rmmissing(changes(:,i));
end

subplot_tight(3,2,[5,6])
hold on
rng(623)
handles = al_goodplot2(changes2, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4,'col',[0.949,0.631,0.008;0.949,0.631,0.008;0.949,0.631,0.008;0.949,0.631,0.008]);
for i=1:4
    scatter(i+(rand(size(changes2{i},1),1)-0.5)*0.4,changes2{i},'filled','k')
end
yline(100,'k--')
xticks(1:4)
xticklabels(["Immobility", "Distance Travelled", "Speed", "Initiation Delay"])
ylabel("% change from stim OFF")
ylim([50,200])
yticks(50:50:200)
set(gca,'fontsize',18)

function hm = normalized_ks(mdl, x, y, a, b)
    D = designMatrix(mdl, 'random');
    B = randomEffects(mdl);
    B2 = fixedEffects(mdl);
    intercepts = D * B + B2(1);
    y = y - intercepts;
    
    [A,B] = meshgrid(a,b);
    temp = ksdensity([x, y], [A(:),B(:)]);
    temp2 = reshape(temp, size(A));
    hm = temp2 ./ repmat(sum(temp2), size(temp2,1),1);
end