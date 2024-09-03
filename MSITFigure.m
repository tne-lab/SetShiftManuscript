%% Load data
traces = readtable('Data/msit_traces.csv');

%% setup
figure('Renderer', 'painters', 'Position', [100 1 1080 320])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% c
subplot_tight(1,3,2)
effect_size=0.1;
hold on
dist=2*traces.dbs_coeff./(traces.dbs_coeff_std+traces.v_std);
disp(median(dist))
disp(nnz(abs(dist)<0.1)/length(dist))
disp(nnz(dist>0)/length(dist))
plot_95kd(dist,'c',[247,119,110]/255)
patch('XData',effect_size*[-1,1,1,-1],'YData',[0,0,0.8,0.8],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
xlim([-2,4])
set(gca,'ycolor','none')
set(gca,'fontsize',18)
xlabel('Group Effect of DBS on Drift Rate')
set(gca,'linewidth',2)

%% d
subplot_tight(2,3,3)
subject_effects=2*table2array(traces(:,48:52))./(traces.dbs_coeff_std+traces.v_std);
pts=-2:0.01:5;
hold on
C=colororder;
for i=1:size(subject_effects,2)
    kd=ksdensity(subject_effects(:,i),pts);
    patch('XData',[pts(1),pts,pts(end)],'YData',[0,kd,0],'FaceColor',C(i,:),'FaceAlpha',0.3,'EdgeColor',C(i,:),'linewidth',2)
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[0,0,0.8,0.8],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
xlim([-2,5])
xticks(-2:5)
legend(["P1","P2","P3","P4","P5"])
xlabel('Individual Effect of DBS on Drift Rate')
set(gca,'ycolor','none')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% e
data_tbl = readtable('Data/DBS_data.csv');
subjects = unique(data_tbl.subj_idx);
rt_diffs = zeros(size(subjects));

for i=1:length(rt_diffs)
    sub_tbl = data_tbl(strcmp(data_tbl.subj_idx, subjects(i)),:);
    rt_diffs(i) = mean(sub_tbl.rt(sub_tbl.dbs==1)) - mean(sub_tbl.rt(sub_tbl.dbs==0));
end

pts=linspace(-1.5,5.5,100)';
reg_mat=zeros(4000,100);
r2=zeros(4000,1);
for i=1:4000
    B=[ones(length(rt_diffs),1),subject_effects(i,:)']\rt_diffs;
    r2(i) = 1 - sum(([ones(length(rt_diffs),1),subject_effects(i,:)']*B-rt_diffs).^2)/sum((rt_diffs-mean(rt_diffs)).^2);
    reg_mat(i,:)=[ones(100,1),pts]*B;
end

subplot_tight(2,3,6)
sorted_mat = sort(reg_mat);
scatter(median(subject_effects),rt_diffs,'filled','MarkerFaceColor',[0.6,0.6,0.6],'MarkerFaceAlpha',0.3)
patch('xdata',[pts;flip(pts)],'ydata',[sorted_mat(100,:),flip(sorted_mat(3900,:))],'EdgeColor','none','facecolor',[0.949,0.631,0.008],'facealpha',0.3)
hold on
plot(pts,median(reg_mat),'linewidth',2,'Color',[0.949,0.631,0.008])
text(0.25*(max(pts)-min(pts))+min(pts), -0.2, strcat("R^2=", num2str(median(r2),"%.3f")),'FontSize',12)
legend('off')
xlim([min(pts),max(pts)])
ylim([-0.25,0.05])
xticks([-1.5, 0, 5.5])
yticks([-0.25,0,0.05])
ylabel('dRT (s)')
xlabel('dV')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% Supplement
sup_traces = readtable('Data/msit_traces2.csv');
figure('Renderer', 'painters', 'Position', [100 1 1080 320])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% S10b
subplot(2,2,2)
effect_size=0.1;
hold on
dist=2*sup_traces.dbs_coeff./(sup_traces.dbs_coeff_std+sup_traces.v_std);
disp(median(dist))
disp(nnz(abs(dist)<0.1)/length(dist))
disp(nnz(dist>0)/length(dist))
plot_95kd(dist,'c',[247,119,110]/255)
patch('XData',effect_size*[-1,1,1,-1],'YData',[0,0,2,2],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.5)
xlim([-1,2])
set(gca,'ycolor','none')
set(gca,'fontsize',18)
xlabel('Group Effect of DBS on Drift Rate')
set(gca,'linewidth',2)

%% S10c
subplot(2,2,3)
subject_effects=2*table2array(sup_traces(:,102:115))./(sup_traces.dbs_coeff_std+sup_traces.v_std);
subject_effects = [subject_effects(:,1), subject_effects(:,7:end), subject_effects(:,2:6)];
pts=-2:0.01:5;
hold on
C=colororder(turbo(size(subject_effects,2)));
for i=1:size(subject_effects,2)
    kd=ksdensity(subject_effects(:,i),pts);
    patch('XData',[pts(1),pts,pts(end)],'YData',[0,kd,0],'FaceColor',C(i,:),'FaceAlpha',0.3,'EdgeColor',C(i,:),'linewidth',2)
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[0,0,2,2],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.5)
xlim([-1,4])
xticks(-1:4)
xlabel('Individual Effect of DBS on Drift Rate')
set(gca,'ycolor','none')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% S10d
data_tbl2 = readtable('Data/Widge_Zorowitz2019_Data.csv');
subjects = unique(data_tbl2.subj_idx);
rt_diffs = zeros(size(subjects));

for i=1:length(rt_diffs)
    sub_tbl = data_tbl2(data_tbl2.subj_idx==subjects(i),:);
    rt_diffs(i) = mean(sub_tbl.rt(sub_tbl.dbs==1)) - mean(sub_tbl.rt(sub_tbl.dbs==0));
end

pts=linspace(-1,3,100)';
reg_mat=zeros(4000,100);
r2=zeros(4000,1);
for i=1:4000
    B=[ones(length(rt_diffs),1),subject_effects(i,:)']\rt_diffs;
    r2(i) = 1 - sum(([ones(length(rt_diffs),1),subject_effects(i,:)']*B-rt_diffs).^2)/sum((rt_diffs-mean(rt_diffs)).^2);
    reg_mat(i,:)=[ones(100,1),pts]*B;
end

subplot(2,2,4)
sorted_mat = sort(reg_mat);
scatter(median(subject_effects),rt_diffs,'filled','MarkerFaceColor',[0.6,0.6,0.6],'MarkerFaceAlpha',0.3)
patch('xdata',[pts;flip(pts)],'ydata',[sorted_mat(100,:),flip(sorted_mat(3900,:))],'EdgeColor','none','facecolor',[0.949,0.631,0.008],'facealpha',0.3)
hold on
plot(pts,median(reg_mat),'linewidth',2,'Color',[0.949,0.631,0.008])
text(0.25*(max(pts)-min(pts))+min(pts), -0.2, strcat("R^2=", num2str(median(r2),"%.3f")),'FontSize',12)
legend('off')
xlim([min(pts),max(pts)])
ylim([-0.3,0.05])
xticks([-1, 0, 3])
yticks([-0.25,0,0.05])
ylabel('dRT (s)')
xlabel('dV')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

function plot_95kd(data,c1,c2,warp)
    [ad, ax] = ksdensity(data,'NumPoints',1000);
    sorted=sort(data);
    l=sorted(0.025*length(sorted));
    h=sorted(0.975*length(sorted));
    if nargin > 3
        patch('XData',[ax(1),ax,ax(end)],'YData',warp([0,ad,0], max(ad)),'EdgeColor','none','FaceColor',c1)
        patch('XData',[l,ax(ax>l&ax<h),h],'YData',warp([0,ad(ax>l&ax<h),0], max(ad)),'EdgeColor','none','FaceColor',c2)
        [~,I] = min(abs(ax-median(data)));
        plot(median(data)*[1,1],warp(ad(I)*[0,1],max(ad)),'Color',0.8*c2,'linewidth',2)
    else
        patch('XData',[ax(1),ax,ax(end)],'YData',[0,ad,0],'EdgeColor','none','FaceColor',c1)
        patch('XData',[l,ax(ax>l&ax<h),h],'YData',[0,ad(ax>l&ax<h),0],'EdgeColor','none','FaceColor',c2)
    end
end