addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))
addpath(genpath('F:\AtlasPlotter'))
cfos_tbl = readtable('Data/cfos.csv');
cfos_tbl.Site(strcmp(cfos_tbl.Site,""))={'No Implant'};
cfos_mat = table2array(cfos_tbl(:,3:end));
g = repmat(cfos_tbl.Properties.VariableNames(3:end),height(cfos_tbl),1);
g2 = repmat(cfos_tbl.Site,1,6);
cfos_mat=cfos_mat(:);
g = g(:);
g = g(~isnan(cfos_mat));
g2 = g2(:);
g2 = g2(~isnan(cfos_mat));
cfos_mat = cfos_mat(~isnan(cfos_mat));

labels = repelem(unique(g),5);
labels(mod(1:length(labels),5) < 3) = {''};
labels(mod(1:length(labels),5) > 3) = {''};
%% Setup figure
figure('Renderer', 'painters', 'Position', [100 100 900 900])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(8,8,[1,2,3,9,10,11])
slice = Slice(3, 'c');
slice.plot()
axis off
axis equal

%% B
subplot_tight(8,6,13:24)
boxplot(cfos_mat',{g, categorical(g2,{'No Implant', 'mid','dorsomedial','ventral','dorsolateral'}, 'Ordinal', true)},'ColorGroup',g2','PlotStyle','compact','FactorGap',[5,0.2],'LabelOrientation','horizontal','Labels',labels,'Colors',[0.0118,0.6275,0.3843;0.102,0.102,0.961;0.949,0.631,0.008;0.792,0,0;0.3,0.3,0.3],'Widths',1.2)
box off
ticks=get(gca,'XTick');
set(gca,'XTickLabel',extractAfter(unique(g),'_'))
set(gca,'XTick',ticks(3:5:end))
set(gca,'linewidth',2)
set(gca,'fontsize',18)
ylim([0,300])
ylabel('Cell Count')

%% C
addpath('Utilities')
full_tbl = HDDM_setup(true);
rats = unique(full_tbl.subj_idx);
deltas = zeros(size(rats));
for i=1:length(rats)
    on_rt = mean(abs(full_tbl.rt(full_tbl.subj_idx==rats(i)&full_tbl.stim==1)));
    off_rt = mean(abs(full_tbl.rt(full_tbl.subj_idx==rats(i)&full_tbl.stim==0)));
    deltas(i) = on_rt - off_rt;
end
full_tbl = full_tbl(full_tbl.subj_idx~=8,:);
delta_bt = bootstrp(4000,@(x) bt(x,full_tbl), 1:height(full_tbl));

model_name = "full_var_stim_fRLDDM_NL";
chains = dir(strcat("Data/TraceData/", model_name));
chains = chains(3:end);
traces = readtable(fullfile(chains(1).folder, chains(1).name));
for i=2:length(chains)
    traces = [traces; readtable(fullfile(chains(i).folder, chains(i).name))];
end
cfos=cfos_tbl(strcmp(cfos_tbl.Site,'mid'),:);
cfos=cfos([9,1:6,8,10:13],:);
valid=[1:7,9:13];
vars = cell(3,1);
vars{1} = table2array(traces(:,82+valid));
vars{2} = table2array(traces(:,52+valid));
vars{3} = table2array(traces(:,112+valid));
vars{4} = delta_bt;
ylabels = ["Bound Sep.", "Drift Rate", "Bias", "RT"];
ylims = [-0.3,-0.4,-0.4,-0.25;0.1,0.8,0.6,0.1];

for v = 1:length(vars)
    for j=1:6
        fos_var=cfos{:,2+j};
        fos_valid = ~isoutlier(fos_var,'quartiles');
        fos_var=fos_var(fos_valid);
        pts=linspace(min(fos_var)-10,max(fos_var)+10,100)';
        reg_mat=zeros(4000,100);
        r2=zeros(4000,1);
        subplot_tight(8,6,24+j+(v-1)*6)
        for i=1:4000
            B=[ones(length(fos_var),1),fos_var]\vars{v}(i,fos_valid)';
            r2(i) = 1 - sum(([ones(length(fos_var),1),fos_var]*B-vars{v}(i,fos_valid)').^2)/sum((vars{v}(i,fos_valid)'-mean(vars{v}(i,fos_valid))).^2);
            reg_mat(i,:)=[ones(100,1),pts]*B;
        end
        sorted_mat = sort(reg_mat);
        scatter(fos_var,median(vars{v}(:,fos_valid)),'filled','MarkerFaceColor',[0.6,0.6,0.6],'MarkerFaceAlpha',0.3)
        patch('xdata',[pts;flip(pts)],'ydata',[sorted_mat(100,:),flip(sorted_mat(3900,:))],'EdgeColor','none','facecolor',[0.949,0.631,0.008],'facealpha',0.3)
        hold on
        plot(pts,median(reg_mat),'linewidth',2,'Color',[0.949,0.631,0.008])
        text(0.25*(max(pts)-min(pts))+min(pts), -0.2, strcat("R^2=", num2str(median(r2),"%.3f")),'FontSize',12)
        legend('off')
        xlim(round([min(pts),max(pts)]))
        ylim(ylims(:,v))
        xticks(round([min(pts),max(pts)]))
        xticklabels([])
        if j > 1
            yticklabels([])
        else
            ylabel(ylabels(v))
        end
        if v==1
                yticks([-0.3,0,0.1])
        elseif v==2
            yticks([-0.4,0,0.8])
        elseif v==3
            yticks([-0.4,0,0.6])
        else
            yticks([-0.2,0,0.1])
        end
        set(gca,'linewidth',2)
        set(gca,'fontsize',18)
    end
end

function deltas = bt(inds, full_tbl)
    sub_tbl = full_tbl(inds,:);
    rats = unique(sub_tbl.subj_idx);
    deltas = zeros(size(rats));
    for i=1:length(rats)
        on_rt = mean(abs(sub_tbl.rt(sub_tbl.subj_idx==rats(i)&sub_tbl.stim==1)));
        off_rt = mean(abs(sub_tbl.rt(sub_tbl.subj_idx==rats(i)&sub_tbl.stim==0)));
        deltas(i) = on_rt - off_rt;
    end
end