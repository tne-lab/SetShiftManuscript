addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))
addpath('Utilities')
addpath('Data')
%% Load Data
[full_tbl, light_seq, rule_seq] = HDDM_setup(false);
bad_sessions=[109,121];
full_tbl = full_tbl(~ismember(full_tbl.split_by,bad_sessions),:);
ndr = nnz(full_tbl.response==-1)/height(full_tbl);
lapse_rate = 0.01;
full_tbl.rt = abs(full_tbl.rt);
full_tbl.rt(full_tbl.response==-1)=-1;
full_tbl.response = -(full_tbl.response - 2);
full_tbl.light = -(full_tbl.light - 2);
model_name = "full_var_stim_fRLDDM_NL";
chains = dir(strcat("Data/TraceData/", model_name));
chains = chains(3:end);
traces = readtable(fullfile(chains(1).folder, chains(1).name));
for i=2:length(chains)
    traces = [traces; readtable(fullfile(chains(i).folder, chains(i).name))];
end
names = traces.Properties.VariableNames;
a_trace = table2array(traces(:, startsWith(names, 'a_subj')));
v_trace = table2array(traces(:, startsWith(names, 'v_subj')));
n_trace = table2array(traces(:, startsWith(names, 't_subj')));
l_trace = table2array(traces(:, startsWith(names, 'alpha_subj')));
b_trace = table2array(traces(:, startsWith(names, 'zt_subj')));
f_trace = table2array(traces(:, startsWith(names, 'forg_subj')));
s_trace = table2array(traces(:, startsWith(names, 'surp_subj')));
st_trace = table2array(traces(:, startsWith(names, 'st')));
a_stim = table2array(traces(:, startsWith(names, 'a_stim_subj')));
v_stim = table2array(traces(:, startsWith(names, 'v_stim_subj')));
n_stim = table2array(traces(:, startsWith(names, 't_stim_subj')));
l_stim = table2array(traces(:, startsWith(names, 'alpha_stim_subj')));
b_stim = table2array(traces(:, startsWith(names, 'zt_stim_subj')));
f_stim = table2array(traces(:, startsWith(names, 'forg_stim_subj')));
s_stim = table2array(traces(:, startsWith(names, 'surp_stim_subj')));
z_trace = table2array(traces(:, strcmp(names, 'z')));
qS_trace = table2array(traces(:, startsWith(names, 'qS')));
qL_trace = table2array(traces(:, startsWith(names, 'qL')));

%% Simulate
T = height(traces);
N = height(full_tbl);
S = length(unique(full_tbl.split_by));

slices = 32;
slice = reshape(1:T, [], slices)';
RT = cell(slices,1);
v_full = cell(slices,1);
a_full = cell(slices,1);
choices = cell(slices,1);
correct_full = cell(slices,1);
rule = cell(slices,1);
rule_num = cell(slices,1);
Vs = cell(slices,1);
Vs_hat = cell(slices,1);
b_full = cell(slices,1);

includes = [1,1,1];
RTMats = cell(size(includes,1),1);
correctMats = cell(size(includes,1),1);
rng(626)
% parpool(16)
for I=1:size(includes,1)
    include = includes(I,:);
    parfor s=1:slices
        inds = slice(s,:);
        a_temp = a_trace(inds,:);
        l_temp = l_trace(inds,:);
        n_temp = n_trace(inds,:);
        v_temp = v_trace(inds,:);
        b_temp = b_trace(inds,:);
        f_temp = f_trace(inds,:);
        s_temp = s_trace(inds,:);
        st_temp = st_trace(inds);
        z_temp = z_trace(inds);
        qs_temp = qS_trace(inds);
        ql_temp = qL_trace(inds);
        a_stim_temp = a_stim(inds,:);
        l_stim_temp = l_stim(inds,:);
        n_stim_temp = n_stim(inds,:);
        v_stim_temp = v_stim(inds,:);
        b_stim_temp = b_stim(inds,:);
        f_stim_temp = f_stim(inds,:);
        s_stim_temp = s_stim(inds,:);
        Vs_temp = zeros(N,T/slices,2);
        Vs_hat_temp = zeros(N,T/slices,2);
        rn_temp = ones(T/slices,S);
        V = [0,0];
        V_hat = [0,0];
        v_t = zeros(N,T/slices);
        a_t = zeros(N,T/slices);
        b_t = zeros(N,T/slices);
        v_that = zeros(N,T/slices);
        b_that = zeros(N,T/slices);
        n_t = zeros(N,T/slices);
        choices_temp = zeros(T/slices,N);
        rule_temp = false(T/slices,N);
        rt_temp = zeros(T/slices, N);
        correct_temp = false(T/slices, N);
        Q = [qs_temp(1), qs_temp(1), ql_temp(1)];
        Q_hat = [qs_temp(1), qs_temp(1), ql_temp(1)];
        blockNum = 1;
        for j=1:T/slices
            ses = 0;
            for i=1:N
                if i == 1 || full_tbl.split_by(i) ~= full_tbl.split_by(i-1)
                    if i~= 1
                        rn_temp(j,ses) = rn_temp(j,ses) + (blockNum - 1) / 5;
                    end
                    Q = [qs_temp(j), qs_temp(j), ql_temp(j)];
                    Q_hat = [qs_temp(j), qs_temp(j), ql_temp(j)];
                    blockNum = 1;
                    ses = ses + 1;
                end
                if full_tbl.stim(i)
                    alpha_t = a_temp(j, full_tbl.subj_idx(i)) + a_stim_temp(j, full_tbl.subj_idx(i)) * include(1);
                    lr_t = logit(l_temp(j, full_tbl.subj_idx(i)) + l_stim_temp(j, full_tbl.subj_idx(i)));
                    delta_t = v_temp(j, full_tbl.subj_idx(i)) + v_stim_temp(j, full_tbl.subj_idx(i)) * include(2);
                    ndt_t = n_temp(j, full_tbl.subj_idx(i)) + n_stim_temp(j, full_tbl.subj_idx(i));
                    bias_t = b_temp(j, full_tbl.subj_idx(i)) + b_stim_temp(j, full_tbl.subj_idx(i)) * include(3);
                    forg_t = logit(f_temp(j, full_tbl.subj_idx(i)) + f_stim_temp(j, full_tbl.subj_idx(i)));
                    surp_t = exp(s_temp(j, full_tbl.subj_idx(i)) + s_stim_temp(j, full_tbl.subj_idx(i)));
                else
                    alpha_t = a_temp(j, full_tbl.subj_idx(i));
                    lr_t = logit(l_temp(j, full_tbl.subj_idx(i)));
                    delta_t = v_temp(j, full_tbl.subj_idx(i));
                    ndt_t = n_temp(j, full_tbl.subj_idx(i));
                    bias_t = b_temp(j, full_tbl.subj_idx(i));
                    forg_t = logit(f_temp(j, full_tbl.subj_idx(i)));
                    surp_t = exp(s_temp(j, full_tbl.subj_idx(i)));
                end
                n_t(i,j) = ndt_t;
                a_t(i,j) = alpha_t;
                if full_tbl.rt(i) ~= -1
                    V(1) = Q(1) + (full_tbl.light(i) == 1) * (Q(3));
                    V(2) = Q(2) + (full_tbl.light(i) == 2) * (Q(3));
                    Vs_temp(i,j,:) = V;
                    v_t(i,j) = (V(1) - V(2)) * delta_t;
                    b_t(i,j) = logit((Q(1) - Q(2)) * bias_t + z_temp(j));
                    pe = full_tbl.feedback(i) - Q(full_tbl.response(i));
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q(full_tbl.response(i)) = Q(full_tbl.response(i)) + adj_lr * pe;
                    Q(abs(full_tbl.response(i)-3)) = Q(abs(full_tbl.response(i)-3)) * (1-forg_t);
                    if full_tbl.response(i) == full_tbl.light(i)
                        pe = full_tbl.feedback(i) - Q(3);
                        adj_lr = lr_t * abs(pe) ^ surp_t;
                        Q(3) = Q(3) + adj_lr * pe;
                    else
                        Q(3) = Q(3) * (1-forg_t);
                    end
                end
                if rand < ndr
                    choices_temp(j,i) = -1;
                    rt_temp(j,i) = -1;
                    correct_temp(j,i) = 0;
                else
                    V_hat(1) = Q_hat(1) + (light_seq(ses,mod(rn_temp(j,ses)-1,8)+1,blockNum) == 1) * (Q_hat(3));
                    V_hat(2) = Q_hat(2) + (light_seq(ses,mod(rn_temp(j,ses)-1,8)+1,blockNum) == 2) * (Q_hat(3));
                    Vs_hat_temp(i,j,:) = V_hat;
                    v_that(i,j) = (V_hat(1) - V_hat(2)) * delta_t;
                    b_that(i,j) = logit((Q_hat(1) - Q_hat(2)) * bias_t + z_temp(j));
                    if rand < lapse_rate
                        choices_temp(j,i) = -(round(rand)-2);
                        rt_temp(j,i) = unifrnd(0,3.005);
                    else
                        rt_temp(j,i) = wienerrng(a_t(i,j), n_t(i) + st_temp(j) * (rand - 0.5),b_that(i,j)*a_t(i,j),v_that(i,j));
                        choices_temp(j, i) = -((rt_temp(j,i)>0)-2);
                        rt_temp(j,i) = abs(rt_temp(j,i));
                    end
                    if choices_temp(j,i) == rule_seq(ses, mod(rn_temp(j,ses)-1,8)+1) || (choices_temp(j,i) == light_seq(ses, mod(rn_temp(j,ses)-1,8)+1, blockNum) && 3 == rule_seq(ses, mod(rn_temp(j,ses)-1,8)+1))
                        correct_temp(j,i) = 1;
                    else
                        correct_temp(j,i) = 0;
                    end
                    pe = correct_temp(j,i) - Q_hat(choices_temp(j,i));
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q_hat(choices_temp(j,i)) = Q_hat(choices_temp(j,i)) + adj_lr * pe;
                    Q_hat(abs(choices_temp(j,i)-3)) = Q_hat(abs(choices_temp(j,i)-3)) * (1-forg_t);
                    if choices_temp(j,i) == light_seq(ses, mod(rn_temp(j,ses)-1, 8) + 1, blockNum)
                        pe = correct_temp(j,i) - Q_hat(3);
                        adj_lr = lr_t * abs(pe) ^ surp_t;
                        Q_hat(3) = Q_hat(3) + adj_lr * pe;
                    else
                        Q_hat(3) = Q_hat(3) * (1-forg_t);
                    end
                end
                rule_temp(j,i) = rule_seq(ses, mod(rn_temp(j,ses)-1,8)+1) > 2;
                if correct_temp(j,i)
                    blockNum = blockNum + 1;
                    if blockNum == 6
                        rn_temp(j,ses) = rn_temp(j,ses) + 1;
                    end
                else
                    blockNum = 1;
                end
                blockNum = mod(blockNum - 1, 5) + 1;
            end
            rn_temp(j,end) = rn_temp(j,end) + (blockNum - 1) / 5;
            disp(j)
        end
        RT{s} = rt_temp;
        rule_num{s} = rn_temp;
        choices{s} = choices_temp;
        correct_full{s} = correct_temp;
        v_full{s} = v_t;
        a_full{s} = a_t;
        b_full{s} = b_t;
        Vs{s} = Vs_temp;
        Vs_hat{s} = Vs_hat_temp;
        rule{s} = rule_temp;
    end
    RTMats{I} = cell2mat(cellfun(@(x) transpose(x),RT,'UniformOutput',false)');
    correctMats{I} = cell2mat(cellfun(@(x) transpose(x),correct_full,'UniformOutput',false)');
end
delete(gcp('nocreate'))
b_full_mat = cell2mat(b_full');
RT_mat=cell2mat(cellfun(@(x) transpose(x),RT,'UniformOutput',false)');
dV_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs,'UniformOutput',false));
dV_hat_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs_hat,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs_hat,'UniformOutput',false));
choices_mat = cell2mat(cellfun(@(x) transpose(x),choices,'UniformOutput',false)');
rule_num_mat = cell2mat(rule_num);
correct_full_mat = cell2mat(cellfun(@(x) transpose(x),correct_full,'UniformOutput',false)');
rule_mat = cell2mat(cellfun(@(x) transpose(x),rule,'UniformOutput',false)');

%% setup
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% a
subplot_tight(4,4,[2,3])
rng(627)
trials = cell(5,1);
hold on
vhigh=4;
vlow=1;
for i=1:length(trials)
    [t, y] = wienerrng2(1.4, 0.18, 0.6, vhigh);
    tsteps=0.18:0.001:abs(t);
    y(y>1.4)=1.4;
    y(y<0)=0;
    plot(tsteps,y(1:length(tsteps)),'Color',[1,0,0,0.3])
    [t, y] = wienerrng2(1.4, 0.18, 0.6, vlow);
    tsteps=0.18:0.001:abs(t);
    y(y>1.4)=1.4;
    y(y<-1.4)=-1.4;
    plot(tsteps,y(1:length(tsteps)),'Color',[0,0,1,0.3])
end
dists = zeros(100000,2);
for i=1:length(dists)
    dists(i,1) = wienerrng(1.4, 0.18, 0.6, vhigh);
    dists(i,2) = wienerrng(1.4, 0.18, 0.6, vlow);
end
steps=0:0.001:1;
kd1=ksdensity(dists(dists(:,1)>0,1),steps)*mean(dists(:,1)>0);
kd2=ksdensity(abs(dists(dists(:,1)<0,1)),steps)*mean(dists(:,1)<0);
kd3=ksdensity(dists(dists(:,2)>0,2),steps)*mean(dists(:,2)>0);
kd4=ksdensity(abs(dists(dists(:,2)<0,2)),steps)*mean(dists(:,2)<0);
plot(steps,0.1*kd1+1.4,'r')
plot(steps,0.1*kd3+1.4,'b')
plot(steps,-0.1*kd2,'r')
plot(steps,-0.1*kd4,'b')
xline(0.18)
yline(0)
yline(1.4)
text(0.18/2, 0.5*1.4,"\tau",'horizontalalignment','center','fontsize',18)
text(0.35, 0.9*1.4,"\nu_t",'horizontalalignment','center','fontsize',18)
plot([0,0.18],[0.6,0.6]*1.4,'k--')
quiver(0.18, 0.6*1.4, 1.4*(1-0.6)/vhigh, 0.4*1.4, 'r')
quiver(0.18, 0.6*1.4, 1.4*(1-0.6)/vlow, 0.4*1.4, 'b')
xlim([0,1])
yticks([0,0.6*1.4,1.4])
yticklabels(["0","\beta_t\alpha","\alpha"])
xticks([])
set(gca,'xcolor','none')
set(gca,'fontsize',18)

%% b
subplot_tight(4,4,6)
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_choice=full_tbl.response(full_tbl.rt>0);
    [temp,I]=sort(dV_full_mat(i,full_tbl.rt>0)+rand(1,nnz(full_tbl.rt>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_choice(I),0.1,'samplepoints',temp),points);
    hat_rt=choices_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = 100*nanmedian(-choice_curves+2,2);
sor_curve = 100*(sort(-choice_curves+2,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=100*nanmedian(-choice_curves_hat+2,2);
sor_curve2 = 100*(sort(-choice_curves_hat+2,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('% Left Choice')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0:50:100)
ylim([0,100])
set(gca,'fontsize',18)

%% c
subplot_tight(4,4,7)
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_rt=full_tbl.rt(full_tbl.rt>0);
    [temp,I]=sort(dV_full_mat(i,full_tbl.rt>0)+rand(1,nnz(full_tbl.rt>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_rt(I),0.1,'samplepoints',temp),points);
    hat_rt=RT_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = nanmedian(choice_curves,2);
sor_curve = (sort(choice_curves,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=nanmedian(choice_curves_hat,2);
sor_curve2 = (sort(choice_curves_hat,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('RT (s)')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0.55:0.1:0.75)
ylim([0.55,0.75])
set(gca,'fontsize',18)

%% d
subplot_tight(4,4,8)
points = 0.35:0.001:0.65;
choice_curves = zeros(length(points),T);
rng(626)
valid=full_tbl.rt>0&~isnan(full_tbl.frontMostLast)&(full_tbl.frontMostLast~=full_tbl.rearMostLast);
var=full_tbl.frontMostLast(valid);
for i=1:T
    [temp,I]=sort(b_full_mat(valid,i)+rand(nnz(valid),1)*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(var(I),0.1,'samplepoints',temp),points);
    disp(i)
end
choice_curves(isnan(choice_curves))=0;
mcurve = median(choice_curves,2);
sor_curve = (sort(choice_curves,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'Color',[180,0,255]/255,'LineWidth',2)
xlabel('Bias')
ylabel('Proportion of Left ITI Pokes')
axis tight
ylim([0,1])
yticks(0:0.5:1)
xticks([0.35,0.5,0.65])
xlim([0.35,0.65])
set(gca,'fontsize',18)

%%
full_tbl.choiceTotal = full_tbl.frontTotal;
full_tbl.choiceTotal(full_tbl.response==2) = full_tbl.rearTotal(full_tbl.response==2);
full_tbl.nonchoiceTotal = full_tbl.frontTotal;
full_tbl.nonchoiceTotal(full_tbl.response==1) = full_tbl.rearTotal(full_tbl.response==1);
% subplot_tight(4,4,8)
points = 0:0.001:0.15;
choice_curves = zeros(length(points),T);
rng(626)
valid=full_tbl.rt>0&~isnan(full_tbl.frontMostLast)&full_tbl.stim==0;
var=full_tbl.choiceTotal(valid)>0;
for i=1:T
    [temp,I]=sort(abs(b_full_mat(valid,i)-0.5)+rand(nnz(valid),1)*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(var(I),0.1,'samplepoints',temp),points);
    disp(i)
end
choice_curves(isnan(choice_curves))=0;
mcurve = median(choice_curves,2);
sor_curve = (sort(choice_curves,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'Color',[180,0,255]/255,'LineWidth',2)
xlabel('Bias')
ylabel('Proportion of Left ITI Pokes')
axis tight
% ylim([0,1])
% yticks(0:0.5:1)
% xticks([0.35,0.5,0.65])
xlim([0,0.15])
set(gca,'fontsize',18)
%% E
mid_params=[traces.a_stim,traces.v_stim,traces.zt_stim, traces.t_stim,traces.alpha_stim,traces.forg_stim,traces.surp_stim];
mid_stds=[traces.a_std,traces.v_std,traces.zt_std,traces.t_std,traces.alpha_std,traces.forg_std,traces.surp_std];
mid_stds=sqrt((mid_stds.^2+[traces.a_stim_std,traces.v_stim_std,traces.zt_stim_std,traces.t_stim_std,traces.alpha_stim_std,traces.forg_stim_std,traces.surp_stim_std])/2);
effect_size=0.1;
subplot_tight(4,4,[9,13])
hold on
for i=1:size(mid_params,2)
    normed=mid_params(:,i)./mid_stds(:,i);
    sorted=sort(normed);
    l=sorted(0.025*length(sorted));
    h=sorted(0.975*length(sorted));
    plot_95kd(normed,'c',[247,119,110]/255,@(x, M)0.9*x/M-(i-1))
    disp(nnz(abs(normed)<effect_size)/length(normed)*100)
    disp(nnz(normed>0)/length(normed)*100)
    disp(median(normed))
    disp(' ')
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[1,1,1-size(mid_params,2),1-size(mid_params,2)],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
yticks(1-size(mid_params,2):0)
yticklabels(flip(["Boundary Separation","Drift Rate","Bias","Non-Decision Time","Learning Rate", "Forgetfulness", "Surprise"]))
ylim([0.5-size(mid_params,2), 1])
xlim([-1,1])
xticks(-1:1)
xlabel('Effect Size')
set(gca,'fontsize',18)
% view([-90,90])

%% e
includes = [0,0,0;1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1;1,1,1];
s=full_tbl.stim;
means = cell2mat(cellfun(@(x) sum(x(s,:).*(x(s,:)~=-1),1)./sum((x(s,:)~=-1)) -...
                              sum(x(~s,:).*(x(~s,:)~=-1),1)./sum((x(~s,:)~=-1)), RTMats, 'uniformoutput',false));
rat_effect = mean(full_tbl.rt(s&full_tbl.rt>-1))-mean(full_tbl.rt(~s&full_tbl.rt>-1));
C = zeros(8, size(means, 2));
C(8, :) = means(2, :) + means(3, :) + means(4, :) - means(7, :) - means(5, :) - means(6, :) + means(8, :) - means(1, :);
C(5, :) = means(7, :) - means(4, :) + means(6, :) - means(8, :) - means(1, :);
C(6, :) = means(5, :) - means(2, :) + means(7, :) - means(8, :) - means(1, :);
C(7, :) = means(5, :) - means(3, :) + means(6, :) - means(8, :) - means(1, :);
C(2, :) = means(2, :) - means(1, :);
C(3, :) = means(3, :) - means(1, :);
C(4, :) = means(4, :) - means(1, :);
C(1, :) = means(1, :);
colors=brewermap(size(C,1)+2,'YlOrRd');
upsetplot(100 * C / rat_effect, logical(includes)',["\alpha","\nu","\beta"], false,flip(colors(2:end-1,:)))
nexttile(3,[3,5])
ylabel("% Stim RT Effect Explained")
set(gca,'fontsize',18)
yticks([-20,0,50,100])
ylim([-20,120])
nexttile(22,[2,2])
xlabel("% Stim RT Effect Explained")
set(gca,'fontsize',18)
xticks([0,50,100])
xlim([-10,100])

%% Supplemental Figure
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
diq = [14248.134496890065, 14251.325314016845, 14252.059369716324, 14251.481475015062;
       14586.649038762325, 14580.88372054762, 14581.65040956457, 14590.874423258996;
       15226.503353013615, 15228.863468216718, 15229.57077875421, 15227.075240536517;
       15044.694554565267, 15046.391303344506, 15048.499382123087, 15047.234368966876;
       15522.173204067321, 15520.447252815724, 15526.486009003931, 15524.739081415728;
       15662.549464514876, 15665.276412214913, 15666.003955551743, 15663.13747151322;
       15603.20163034262, 15603.645197478056, 15601.546380514454, 15607.846798123324;
       16016.784385997948, 16015.49341265982, 16019.568097128431, 16017.608076759381];

subplot_tight(4,3,1)
barh(mean(diq,2), 'EdgeColor', 'none')
xlim([1.4e4, 1.65e4])
yticklabels(["Base Model", "No B/G", "No G/F", "No B/F", "No B", "No G", "No F", "Full Model"])
box off
set(gca,'fontsize',18)
xlabel('DIC')
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% B
subplot_tight(4,3,2)
sessions = unique(full_tbl.split_by);
acc_side = zeros(length(sessions), 4000);
acc_light = zeros(length(sessions), 4000);
rat_side = zeros(size(sessions));
rat_light = zeros(size(sessions));
for i=1:length(sessions)
    rule_seg = rule_mat(full_tbl.split_by==i,:);
    correct_seg = correct_full_mat(full_tbl.split_by==i,:);
    acc_side(i,:) = sum(correct_seg.*~rule_seg)./sum(~rule_seg);
    acc_light(i,:) = sum(correct_seg.*rule_seg)./sum(rule_seg);
    tbl_seg = full_tbl(full_tbl.split_by==i,:);
    rat_side(i) = mean(tbl_seg.feedback(~tbl_seg.rule));
    rat_light(i) = mean(tbl_seg.feedback(tbl_seg.rule));
end
hold on
pts = 0:0.005:1;
for i=1:50:4000
    plot(pts,ksdensity(acc_side(:,i),pts), 'Color',[0.31,0.73,0.90,0.05])
end
for i=1:50:4000
    plot(pts,ksdensity(acc_light(:,i),pts), 'Color', [0.83,0.07,0.82,0.05])
end
plot(pts,ksdensity(rat_side, pts),'Color',[0.31,0.73,0.90],'LineWidth',2)
plot(pts,ksdensity(rat_light, pts),'Color',[0.83,0.07,0.82],'LineWidth',2)
set(gca,'fontsize',18)
xlabel('Accuracy')
yticks([])
ylim([0,8])
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% C
subplot_tight(4,3,3)
hold on
pts = 0:0.001:5;
stims = full_tbl.stim==1;
all_vals = zeros(size(pts));
bw = 0.05;
for s = 1:slices
    for i=1:50:T/slices
        temp = RT{s}(i,:);
        temp1 = squeeze(temp(stims));
        temp1=rmmissing(temp1(temp1>0));
        temp2 = squeeze(temp(~stims));
        temp2=rmmissing(temp2(temp2>0));
        if ~isempty(temp1)
            vals = ksdensity(temp2,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.6,0.6,0.6,0.05])
            vals = ksdensity(temp1,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.949,0.631,0.008,0.05])
            all_vals = all_vals + vals;
        end
    end
end
hold on
vals=ksdensity(squeeze(full_tbl.rt(~stims&full_tbl.rt>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.6,0.6,0.6],'LineWidth',1.5)
vals=ksdensity(squeeze(full_tbl.rt(stims&full_tbl.rt>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.949,0.631,0.008],'LineWidth',1.5)
xticks(0:3)
xlim([0,3])
xlabel('RT (s)')
yticks([])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% D
subplot_tight(4,6,[8,9])
hold on
plot_95kd(reshape(sum(RT_mat(full_tbl.stim==0,:).*(choices_mat(full_tbl.stim==0,:)~=-1),1),1,[])./reshape(sum((choices_mat(full_tbl.stim==0,:)~=-1),1),1,[]), 'c',[0.6,0.6,0.6])
plot_95kd(reshape(sum(RT_mat(full_tbl.stim==1,:).*(choices_mat(full_tbl.stim==1,:)~=-1),1),1,[])./reshape(sum((choices_mat(full_tbl.stim==1,:)~=-1),1),1,[]), 'c',[0.949,0.631,0.008])
xline(mean(full_tbl.rt(full_tbl.stim==0&full_tbl.rt>-1)),'--','LineWidth',2,'Color','k')
xline(mean(full_tbl.rt(full_tbl.stim==1&full_tbl.rt>-1)),'--','LineWidth',2,'Color','k')
set(gca,'fontsize',18)
yticks([])
xlabel('RT (s)')
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% E
subplot_tight(4,6,[10,11])
plot_95kd(mean(rule_num_mat,2)-1, 'c',[247,119,110]/255)
yticks([])
xline(8,'k--','LineWidth',2)
xlabel('Rules Completed')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%%
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