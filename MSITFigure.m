%% Load data
traces = readtable('Data/msit_traces.csv');

%% setup
figure('Renderer', 'painters', 'Position', [100 1 1080 320])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% b
subplot_tight(1,3,2)
effect_size=0.1;
hold on
plot_95kd(2*traces.dbs_coeff./(traces.dbs_coeff_std+traces.v_std),'c',[247,119,110]/255)
patch('XData',effect_size*[-1,1,1,-1],'YData',[0,0,0.8,0.8],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
xlim([-2,4])
set(gca,'ycolor','none')
set(gca,'fontsize',18)
xlabel('Group Effect of DBS on Drift Rate')
set(gca,'linewidth',2)

%% c
subplot_tight(1,3,3)
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