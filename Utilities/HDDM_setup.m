function [full_tbl, light_seq, rule_seq] = HDDM_setup(remove_omissions)
    load('Data/allBehav.mat')
    bn = load('Data/allBehav2.mat');
    behavior = rmfield(behavior, 'drug');
    behavior=[behavior, bn.behavior(2:5)];
    behav_sub = behavior;
    behav_sub([31])=[];
    behav_sub([15 16 34])=[];

    files_off = dir('Data/ITITablesOff');
    files_on = dir('Data/ITITablesOn');
    
    site = "mid";
    
    sub_rats=1:length(behav_sub);
    site_rats = false(size(sub_rats));
    for i=1:length(site_rats)
        site_rats(i) = strcmp(behav_sub(i).site{1},site);
    end
    sub_rats=sub_rats(site_rats);
    
    ntrials=0;
    for i=1:length(sub_rats)
        for j=1:length(behav_sub(sub_rats(i)).sessions)
            ntrials=ntrials+length(behav_sub(sub_rats(i)).sessions{j});
        end
    end
    full_tbl = table('Size',[ntrials,7],'VariableTypes',{'double','double','double','logical','double','logical','double'});
    full_tbl.Properties.VariableNames = {'rt','subj_idx','split_by','stim','light','feedback','response'};
    ind=1;
    ses=1;
    for i=1:length(sub_rats)
        foff=files_off(contains(string({files_off.name}),behav_sub(sub_rats(i)).subject));
        fon=files_on(contains(string({files_on.name}),behav_sub(sub_rats(i)).subject));
        rat_files = [foff;fon];
        names = [extractAfter(string({foff.name}),'GLMSTIMOFF'), extractAfter(string({fon.name}),'GLMSTIMON')];
        [~,I]=sort(names);
        rat_files=rat_files(I);
        for j=1:length(behav_sub(sub_rats(i)).sessions)
            if (sub_rats(i)~=13||j~=4)
                len=length(behav_sub(sub_rats(i)).sessions{j});
                if j <= length(rat_files)
                    load(fullfile(rat_files(j).folder, rat_files(j).name))
                else
                    full_tbl.frontMostLast(ind:ind+len-1)=NaN;
                    full_tbl.rearMostLast(ind:ind+len-1)=NaN;
                end
                full_tbl.response(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.frontChoice];
                full_tbl.rt(ind:ind+len-1)=2*(full_tbl.response(ind:ind+len-1)-0.5).*[behav_sub(sub_rats(i)).sessions{j}.RT]';
                full_tbl.subj_idx(ind:ind+len-1)=i;
                full_tbl.split_by(ind:ind+len-1)=ses;
                full_tbl.stim(ind:ind+len-1)=behav_sub(sub_rats(i)).stimOn(j);
                full_tbl.light(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.light];
                full_tbl.feedback(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.performance];
                full_tbl.rule(ind:ind+len-1)=contains([behav_sub(sub_rats(i)).sessions{j}.rule],"L");
                full_tbl.dt(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.initTime]'-[behav_sub(sub_rats(i)).sessions{j}.cueTime]';
                if height(setshift) == len
                    full_tbl.frontMostLast(ind:ind+len-1)=double(setshift.lasttwofront>setshift.lasttworear);
                    full_tbl.rearMostLast(ind:ind+len-1)=double(setshift.lasttwofront<setshift.lasttworear);
                    full_tbl.frontTotal(ind:ind+len-1) = setshift.lasttwofront;
                    full_tbl.rearTotal(ind:ind+len-1) = setshift.lasttworear;
                else
                    full_tbl.frontMostLast(ind:ind+len-1)=NaN;
                    full_tbl.rearMostLast(ind:ind+len-1)=NaN;
                    full_tbl.frontTotal(ind:ind+len-1)=NaN;
                    full_tbl.rearTotal(ind:ind+len-1)=NaN;
                end
                full_tbl.ruleFull(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.rule];
                ind=ind+len;
                ses=ses+1;
            end
        end
    end
    light_seq = zeros(ses-1, 8, 5);
    rule_seq = zeros(ses-1, 8);
    num = 1;
    for r=1:length(sub_rats)
        for i=1:length(behav_sub(sub_rats(r)).sessions)
            if (sub_rats(r)~=13||i~=4)
                ss = SetShift(behav_sub(sub_rats(r)).sessions{i}, false);
                light_seq(num,:,:) = -ss.lightSeq+2;
                rule_seq(num,:) = ss.ruleSeq;
                num = num + 1;
            end
        end
    end
    if remove_omissions
        full_tbl(isnan(full_tbl.response),:)=[];
    end
    full_tbl.rt(isnan(full_tbl.response),:)=-1;
    full_tbl.response(isnan(full_tbl.response),:)=-1;
    full_tbl(full_tbl.subj_idx==0,:)=[];
    full_tbl.q_init = ones(height(full_tbl),1)*0.33;
    full_tbl.q1 = ones(height(full_tbl),1)*0.33;
    full_tbl.q2 = ones(height(full_tbl),1)*0.67;
end