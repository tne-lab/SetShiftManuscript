addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))
addpath('Utilities')
%% Load Data
load('Data/allBehav.mat')
bn = load('Data/allBehav2.mat');
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
files_off = dir('Data/ITITablesOff');
files_on = dir('Data/ITITablesOn');
full_tbl = table('Size',[ntrials,8],'VariableTypes',{'double','double','double','logical','double','double','logical','double'});
full_tbl.Properties.VariableNames = {'RT','Rat','Session','Stim','Site','Light','CorrectTrial','Choice'};
sites=["mid","dorsomedial","ventral","dorsolateral"];
ind=1;
ses=1;
for i=1:length(sub_rats)
    foff=files_off(contains(string({files_off.name}),valid_rats(i)));
    fon=files_on(contains(string({files_on.name}),valid_rats(i)));
    rat_files = [foff;fon];
    names = [extractAfter(string({foff.name}),'GLMSTIMOFF'), extractAfter(string({fon.name}),'GLMSTIMON')];
    [~,I]=sort(names);
    rat_files=rat_files(I);
    for j=1:length(behav_sub(sub_rats(i)).sessions)
        if (sub_rats(i)~=13||j~=4)
            if j <= length(rat_files)
                if (sub_rats(i)==13 && j > 4) || (i==8 && j > 8)
                    load(fullfile(rat_files(j-1).folder, rat_files(j-1).name))
                else
                    load(fullfile(rat_files(j).folder, rat_files(j).name))
                end
            else
                full_tbl.frontMostLast(ind:ind+len-1)=NaN;
                full_tbl.rearMostLast(ind:ind+len-1)=NaN;
                full_tbl.midPokeLast(ind:ind+len-1)=NaN;
                full_tbl.noPokeLast(ind:ind+len-1)=NaN;
            end
            len=length(behav_sub(sub_rats(i)).sessions{j});
            full_tbl.RT(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.RT];
            full_tbl.Rat(ind:ind+len-1)=i;
            full_tbl.Session(ind:ind+len-1)=ses;
            full_tbl.Stim(ind:ind+len-1)=behav_sub(sub_rats(i)).stimOn(j);
            full_tbl.Site(ind:ind+len-1)=find(strcmp(sites,behav_sub(sub_rats(i)).site{j}));
            full_tbl.Light(ind:ind+len-1)=-([behav_sub(sub_rats(i)).sessions{j}.light]-2);
            full_tbl.CorrectTrial(ind:ind+len-1)=[behav_sub(sub_rats(i)).sessions{j}.performance];
            full_tbl.Choice(ind:ind+len-1)=-([behav_sub(sub_rats(i)).sessions{j}.frontChoice]-2);
            full_tbl.Rule(ind:ind+len-1)=extractBefore([behav_sub(sub_rats(i)).sessions{j}.rule],2);
            full_tbl.Trial(ind:ind+len-1)=1:len;
            if height(setshift) == len
                lasttwochange = setshift.lasttwofront .* [0;1-setshift.frontChoice(1:end-1)] + setshift.lasttworear .* [0;setshift.frontChoice(1:end-1)];
                firsttwochange = setshift.firsttwofront .* [0;1-setshift.frontChoice(1:end-1)] + setshift.firsttworear .* [0;setshift.frontChoice(1:end-1)];
                full_tbl.prechoiceMostFirst(ind:ind+len-1)=double(setshift.firsttwocompulsive>firsttwochange);
                full_tbl.PreChoiceLateITI(ind:ind+len-1)=double(setshift.lasttwocompulsive>0);
                full_tbl.NonPreChoiceLateITI(ind:ind+len-1)=double(lasttwochange>0);
                full_tbl.nonprechoiceMostFirst(ind:ind+len-1)=double(setshift.firsttwocompulsive<firsttwochange);
                full_tbl.MorePreChoiceLateITI(ind:ind+len-1)=double(setshift.lasttwocompulsive>lasttwochange);
                full_tbl.MoreNonPreChoiceLateITI(ind:ind+len-1)=double(setshift.lasttwocompulsive<lasttwochange);
                full_tbl.notequalEarly(ind:ind+len-1)=double(setshift.firsttwofront~=setshift.firsttworear);
                full_tbl.notequalLast(ind:ind+len-1)=double(setshift.lasttwofront~=setshift.lasttworear);
                full_tbl.MoreLeftEarlyITI(ind:ind+len-1)=double(setshift.firsttwofront>setshift.firsttworear);
                full_tbl.MoreLeftLateITI(ind:ind+len-1)=double(setshift.lasttwofront>setshift.lasttworear);
                full_tbl.rearMostLast(ind:ind+len-1)=double(setshift.lasttwofront<setshift.lasttworear);
                full_tbl.MidPokeLateITI(ind:ind+len-1)=double(setshift.lasttwomiddle>0);
                full_tbl.noPokeLast(ind:ind+len-1)=double(setshift.lasttwomiddle==0&setshift.lasttwofront==0&setshift.lasttworear==0);
                full_tbl.frontPoke(ind:ind+len-1)=double(setshift.preititotalfront>0);
                full_tbl.rearPoke(ind:ind+len-1)=double(setshift.preititotalrear>0);
                full_tbl.frontPokeFirst(ind:ind+len-1)=double(setshift.firsttwofront>0);
                full_tbl.rearPokeFirst(ind:ind+len-1)=double(setshift.firsttworear>0);
                full_tbl.frontPokeMid(ind:ind+len-1)=double(setshift.midthreefront>0);
                full_tbl.rearPokeMid(ind:ind+len-1)=double(setshift.midthreerear>0);
                full_tbl.frontPokeLast(ind:ind+len-1)=double(setshift.lasttwofront>0);
                full_tbl.rearPokeLast(ind:ind+len-1)=double(setshift.lasttworear>0);
                full_tbl.sidePoke(ind:ind+len-1)=full_tbl.frontPoke(ind:ind+len-1)|full_tbl.rearPoke(ind:ind+len-1);
                full_tbl.NumCorrectPokeEarlyITI(ind:ind+len-1)=setshift.firsttwocorrective.*(1-full_tbl.CorrectTrial(ind:ind+len-1))+...
                    setshift.firsttwocompulsive.*full_tbl.CorrectTrial(ind:ind+len-1);
                full_tbl.NumCorrectPokeLateITI(ind:ind+len-1)=setshift.lasttwocorrective.*(1-full_tbl.CorrectTrial(ind:ind+len-1))+...
                    setshift.lasttwocompulsive.*full_tbl.CorrectTrial(ind:ind+len-1);
            else
                full_tbl.prechoiceMostFirst(ind:ind+len-1)=NaN;
                full_tbl.nonprechoiceMostFirst(ind:ind+len-1)=NaN;
                full_tbl.MorePreChoiceLateITI(ind:ind+len-1)=NaN;
                full_tbl.MoreNonPreChoiceLateITI(ind:ind+len-1)=NaN;
                full_tbl.MoreLeftLateITI(ind:ind+len-1)=NaN;
                full_tbl.rearMostLast(ind:ind+len-1)=NaN;
                full_tbl.midPokeLast(ind:ind+len-1)=NaN;
                full_tbl.noPokeLast(ind:ind+len-1)=NaN;
            end
            ind=ind+len;
            ses=ses+1;
        end
    end
end
full_tbl=full_tbl(full_tbl.Session>0,:);
full_tbl.RT(isnan(full_tbl.Choice),:)=-1;
full_tbl.Choice(isnan(full_tbl.Choice),:)=-1;
%%
load('..\behaviorData\DN_entries.mat')
load('..\behaviorData\allBehavDN.mat')
kstd = 0.2;
times = 0:0.001:7.1;
timesFull = -7.1:0.001:14.2;
ratsCorrectMatch = zeros(length(sub_rats),length(times));
ratsCorrectNonMatch = zeros(length(sub_rats),length(times));
ratsCorrectMid = zeros(length(sub_rats),length(times));
ratsIncorrectMatch = zeros(length(sub_rats),length(times));
ratsIncorrectNonMatch = zeros(length(sub_rats),length(times));
ratsIncorrectMid = zeros(length(sub_rats),length(times));
for i=1:length(sub_rats)
    ctrials = 0;
    itrials = 0;
    if contains(valid_rats(i), "DN")
        r = find(strcmp([all_entries.subject], valid_rats(i)));
        for j=1:length(all_entries(r).sessions)
            setshift = behavior(r).sessions{j};
            setshiftall = all_entries(r).sessions{j};
            for k=1:length(setshift)-1
                endoftrial = find([setshiftall.taskTrial]==k,1,'last');
                startofnext = find([setshiftall.taskTrial]==k+1,1);
                itipokes = setshiftall(endoftrial+1:startofnext-1);
                poketimes = [itipokes.entryTime] - setshift(k).respTime;
                poketypes = [itipokes.location];
                if ~isnan(setshift(k).frontChoice)
                    if setshift(k).performance
                        ctrials = ctrials + 1;
                        if nnz(poketypes==3)
                            ratsCorrectMid(i,:) = ratsCorrectMid(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==3)')/kstd).^2),1),7101,14201);
                        end
                        if nnz(poketypes==1)
                            if setshift(k).frontChoice == 1
                                ratsCorrectNonMatch(i,:) = ratsCorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1)')/kstd).^2),1),7101,14201);
                            else
                                ratsCorrectMatch(i,:) = ratsCorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1)')/kstd).^2),1),7101,14201);
                            end
                        end
                        if nnz(poketypes==2)
                            if setshift(k).frontChoice == 0
                                ratsCorrectNonMatch(i,:) = ratsCorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2)')/kstd).^2),1),7101,14201);
                            else
                                ratsCorrectMatch(i,:) = ratsCorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2)')/kstd).^2),1),7101,14201);
                            end
                        end
                    else
                        itrials = itrials + 1;
                        if nnz(poketypes==3)
                            ratsIncorrectMid(i,:) = ratsIncorrectMid(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==3)')/kstd).^2),1),7101,14201);
                        end
                        if nnz(poketypes==1)
                            if setshift(k).frontChoice == 1
                                ratsIncorrectNonMatch(i,:) = ratsIncorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1)')/kstd).^2),1),7101,14201);
                            else
                                ratsIncorrectMatch(i,:) = ratsIncorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1)')/kstd).^2),1),7101,14201);
                            end
                        end
                        if nnz(poketypes==2)
                            if setshift(k).frontChoice == 0
                                ratsIncorrectNonMatch(i,:) = ratsIncorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2)')/kstd).^2),1),7101,14201);
                            else
                                ratsIncorrectMatch(i,:) = ratsIncorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2)')/kstd).^2),1),7101,14201);
                            end
                        end
                    end
                end
            end
        end
    else
        offfiles = dir(fullfile('Data/RawBehaviorOff',valid_rats(i)));
        onfiles = dir(fullfile('Data/RawBehaviorMidOn',valid_rats(i)));
        files = [offfiles(3:end);...
                 onfiles(3:end)];
        for j=1:2:length(files)
            filename=fullfile(files(j).folder,files(j).name);
            setshiftall = readtable(filename);
            filenametwo = fullfile(files(j+1).folder,files(j+1).name);
            setshift =  readtable(filenametwo);
            for k=1:height(setshift)-1
                endoftrial = find(setshiftall.taskTrial==k,1,'last');
                startofnext = find(setshiftall.taskTrial==k+1,1);
                itipokes = setshiftall(endoftrial+1:startofnext-1,:);
                poketimes = itipokes.entryTime - setshift.respTime(k);
                poketypes = itipokes.location;
                if ~isnan(setshift.frontChoice(k))
                    if setshift.performance(k)
                        ctrials = ctrials + 1;
                        if nnz(poketypes==3)
                            ratsCorrectMid(i,:) = ratsCorrectMid(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==3))/kstd).^2),1),7101,14201);
                        end
                        if nnz(poketypes==1)
                            if setshift.frontChoice(k) == 1
                                ratsCorrectNonMatch(i,:) = ratsCorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1))/kstd).^2),1),7101,14201);
                            else
                                ratsCorrectMatch(i,:) = ratsCorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1))/kstd).^2),1),7101,14201);
                            end
                        end
                        if nnz(poketypes==2)
                            if setshift.frontChoice(k) == 0
                                ratsCorrectNonMatch(i,:) = ratsCorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2))/kstd).^2),1),7101,14201);
                            else
                                ratsCorrectMatch(i,:) = ratsCorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2))/kstd).^2),1),7101,14201);
                            end
                        end
                    else
                        itrials = itrials + 1;
                        if nnz(poketypes==3)
                            ratsIncorrectMid(i,:) = ratsIncorrectMid(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==3))/kstd).^2),1),7101,14201);
                        end
                        if nnz(poketypes==1)
                            if setshift.frontChoice(k) == 1
                                ratsIncorrectNonMatch(i,:) = ratsIncorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1))/kstd).^2),1),7101,14201);
                            else
                                ratsIncorrectMatch(i,:) = ratsIncorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==1))/kstd).^2),1),7101,14201);
                            end
                        end
                        if nnz(poketypes==2)
                            if setshift.frontChoice(k) == 0
                                ratsIncorrectNonMatch(i,:) = ratsIncorrectNonMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2))/kstd).^2),1),7101,14201);
                            else
                                ratsIncorrectMatch(i,:) = ratsIncorrectMatch(i,:) + boundary_correction(sum(1/(kstd*sqrt(2*pi))*exp(-1/2*((timesFull-poketimes(poketypes==2))/kstd).^2),1),7101,14201);
                            end
                        end
                    end
                end
            end
        end
    end
    ratsCorrectMid(i,:) = ratsCorrectMid(i,:) / ctrials;
    ratsCorrectMatch(i,:) = ratsCorrectMatch(i,:) / ctrials;
    ratsCorrectNonMatch(i,:) = ratsCorrectNonMatch(i,:) / ctrials;
    ratsIncorrectMid(i,:) = ratsIncorrectMid(i,:) / itrials;
    ratsIncorrectMatch(i,:) = ratsIncorrectMatch(i,:) / ctrials;
    ratsIncorrectNonMatch(i,:) = ratsIncorrectNonMatch(i,:) / ctrials;
    disp(i)
end

%% Tables
iti_tbl=full_tbl(full_tbl.Choice~=-1&~isnan(full_tbl.MoreLeftLateITI),:);
addpath(genpath('C:\Users\Evan\Documents\GitHub\py-behav-box-v2\source\reg2latex'))
iti_tbl.CorrectTrial=double(iti_tbl.CorrectTrial);
iti_tbl.LeftChosen=double(iti_tbl.Choice==1);
iti_tbl.LightChosen=double(iti_tbl.Choice==iti_tbl.Light);
iti_tbl.LeftCorrect=double((iti_tbl.Light==1&strcmp(iti_tbl.Rule,"L"))|strcmp(iti_tbl.Rule,"F"));
sessions=unique(iti_tbl.Session);
shifted_iti_tbl=iti_tbl(1:height(iti_tbl)-length(sessions),:);
ind=1;
for i=1:length(sessions)
    subtbl=iti_tbl(iti_tbl.Session==sessions(i),:);
    shifted_iti_tbl(ind:ind+height(subtbl)-2,:)=subtbl(2:end,:);
    shifted_iti_tbl.CorrectTrial(ind:ind+height(subtbl)-2,:)=subtbl.CorrectTrial(1:end-1,:);
    shifted_iti_tbl.LeftCorrect(ind:ind+height(subtbl)-2,:)=subtbl.LeftCorrect(1:end-1,:);
    ind = ind + height(subtbl)-1;
end
shifted_iti_tbl=shifted_iti_tbl(~isnan(shifted_iti_tbl.NumCorrectPokeEarlyITI),:);

mdl1=mdltostruct(fitglme(shifted_iti_tbl,'NumCorrectPokeEarlyITI~1+CorrectTrial+(1|Rat)','Distribution','Poisson','link','identity'));
mdl2=mdltostruct(fitglme(shifted_iti_tbl,'NumCorrectPokeLateITI~1+CorrectTrial+(1|Rat)','Distribution','Poisson','link','identity'));
mdl3=mdltostruct(fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalEarly==1,:),'MoreLeftEarlyITI~1+CorrectTrial*LeftCorrect+(1|Rat)','Distribution','Binomial'));
mdl4=mdltostruct(fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalLast==1,:),'MoreLeftLateITI~1+CorrectTrial*LeftCorrect+(1|Rat)','Distribution','Binomial'));
mdl5=mdltostruct(fitglme(iti_tbl(iti_tbl.notequalLast==1,:),'LeftChosen~1+MoreLeftLateITI*LightChosen+(1|Rat)','Distribution','Binomial'));
mdl6=mdltostruct(fitglme(iti_tbl,'RT~1+MidPokeLateITI+MorePreChoiceLateITI+MoreNonPreChoiceLateITI+MidPokeLateITI:MorePreChoiceLateITI+(1|Rat)','Distribution','Gamma','link','identity'));

fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalEarly==1&shifted_iti_tbl.CorrectTrial,:),'MoreLeftEarlyITI~1+LeftCorrect+(1|Rat)','Distribution','Binomial')
fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalLast==1&shifted_iti_tbl.CorrectTrial,:),'MoreLeftLateITI~1+LeftCorrect+(1|Rat)','Distribution','Binomial')
fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalEarly==1&~shifted_iti_tbl.CorrectTrial,:),'MoreLeftEarlyITI~1+LeftCorrect+(1|Rat)','Distribution','Binomial')
fitglme(shifted_iti_tbl(shifted_iti_tbl.notequalLast==1&~shifted_iti_tbl.CorrectTrial,:),'MoreLeftLateITI~1+LeftCorrect+(1|Rat)','Distribution','Binomial')
fitglme(shifted_iti_tbl,'NumCorrectPokeEarlyITI~1+CorrectTrial+(1|Rat)','Distribution','Poisson','link','identity')
fitglme(shifted_iti_tbl,'NumCorrectPokeLateITI~1+CorrectTrial+(1|Rat)','Distribution','Poisson','link','identity')
fitglme(iti_tbl(iti_tbl.notequalLast==1&~iti_tbl.LightChosen,:),'LeftChosen~1+MoreLeftLateITI+(1|Rat)','Distribution','Binomial')
fitglme(iti_tbl(iti_tbl.notequalLast==1&iti_tbl.LightChosen,:),'LeftChosen~1+MoreLeftLateITI+(1|Rat)','Distribution','Binomial')
fitglme(iti_tbl,'RT~1+MidPokeLateITI+MorePreChoiceLateITI+MoreNonPreChoiceLateITI+MidPokeLateITI:MorePreChoiceLateITI+(1|Rat)','Distribution','Gamma','link','identity')
fitglme(iti_tbl(~isnan(iti_tbl.MidPokeLateITI),:),'MidPokeLateITI~1+Stim+(1|Rat)','Distribution','Binomial')
fitglme(iti_tbl(~isnan(iti_tbl.MorePreChoiceLateITI),:),'MorePreChoiceLateITI~1+Stim+(1|Rat)','Distribution','Binomial')
fitglme(iti_tbl(~isnan(iti_tbl.MoreNonPreChoiceLateITI),:),'MoreNonPreChoiceLateITI~1+Stim+(1|Rat)','Distribution','Binomial')
%% setup
figure('Renderer', 'painters', 'Position', [100 100 1600 900])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
%% a
icorrChoice=zeros(S,1);
icorrNChoice=zeros(S,1);
icorrChoiceLate=zeros(S,1);
icorrNChoiceLate=zeros(S,1);
for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    shiftbl = subtbl(2:end,:);
    icorrChoice(i)=nanmean(shiftbl.prechoiceMostFirst(~subtbl.CorrectTrial(1:end-1)&(subtbl.prechoiceMostFirst(2:end)==1|subtbl.nonprechoiceMostFirst(2:end)==1)));
    icorrNChoice(i)=nanmean(shiftbl.nonprechoiceMostFirst(~subtbl.CorrectTrial(1:end-1)&(subtbl.prechoiceMostFirst(2:end)==1|subtbl.nonprechoiceMostFirst(2:end)==1)));
    icorrChoiceLate(i)=nanmean(shiftbl.MorePreChoiceLateITI(~subtbl.CorrectTrial(1:end-1)&(subtbl.MorePreChoiceLateITI(2:end)==1|subtbl.MoreNonPreChoiceLateITI(2:end)==1)));
    icorrNChoiceLate(i)=nanmean(shiftbl.MoreNonPreChoiceLateITI(~subtbl.CorrectTrial(1:end-1)&(subtbl.MorePreChoiceLateITI(2:end)==1|subtbl.MoreNonPreChoiceLateITI(2:end)==1)));
end
icorrChoice = rmmissing(icorrChoice); icorrNChoice = rmmissing(icorrNChoice);
icorrChoiceLate = rmmissing(icorrChoiceLate); icorrNChoiceLate = rmmissing(icorrNChoiceLate);

subplot_tight(2,6,[1,2,3])
hold on
plot(times,ratsIncorrectMid','Color',[1,0.604,0,0.2],'LineWidth',1.5)
plot(times,ratsIncorrectMatch','Color',[0.133,0.545,0.133,0.2],'LineWidth',1.5)
plot(times,ratsIncorrectNonMatch','Color',[0.796,0.377,0.082,0.2],'LineWidth',1.5)
plot(times,mean(ratsIncorrectMid),'Color',[1,0.604,0],'LineWidth',2)
plot(times,mean(ratsIncorrectMatch),'Color',[0.133,0.545,0.133],'LineWidth',2)
plot(times,mean(ratsIncorrectNonMatch),'Color',[0.796,0.377,0.082],'LineWidth',2)
xlabel('Time Since Start of ITI (s)')
title('Post-Incorrect Trial')
set(gca,'fontsize',18)
yticks(0:0.3:0.9)
xlim([0,7])
ylim([0,0.9])
yticklabels([])
axes('Position',[.65 .75 .1 .1])
box on
al_goodplot2({icorrChoice,icorrNChoice}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.133,0.545,0.133;0.796,0.377,0.082]);
xticks([])
xlabel('Early Window')
ylabel(["Trial Proportion With","More of Poke Type"])
axes('Position',[.825 .75 .1 .1])
box on
al_goodplot2({icorrChoiceLate,icorrNChoiceLate}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.133,0.545,0.133;0.796,0.377,0.082]);
xticks([])
xlabel('Late Window')
%% b
S = max(full_tbl.Session);
corrChoice=zeros(S,1);
corrNChoice=zeros(S,1);
corrChoiceLate=zeros(S,1);
corrNChoiceLate=zeros(S,1);
for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    shiftbl = subtbl(2:end,:);
    corrChoice(i)=nanmean(shiftbl.prechoiceMostFirst(subtbl.CorrectTrial(1:end-1)&(subtbl.prechoiceMostFirst(2:end)==1|subtbl.nonprechoiceMostFirst(2:end)==1)));
    corrNChoice(i)=nanmean(shiftbl.nonprechoiceMostFirst(subtbl.CorrectTrial(1:end-1)&(subtbl.prechoiceMostFirst(2:end)==1|subtbl.nonprechoiceMostFirst(2:end)==1)));
    corrChoiceLate(i)=nanmean(shiftbl.MorePreChoiceLateITI(subtbl.CorrectTrial(1:end-1)&(subtbl.MorePreChoiceLateITI(2:end)==1|subtbl.MoreNonPreChoiceLateITI(2:end)==1)));
    corrNChoiceLate(i)=nanmean(shiftbl.MoreNonPreChoiceLateITI(subtbl.CorrectTrial(1:end-1)&(subtbl.MorePreChoiceLateITI(2:end)==1|subtbl.MoreNonPreChoiceLateITI(2:end)==1)));
end
corrChoice = rmmissing(corrChoice); corrNChoice = rmmissing(corrNChoice);
corrChoiceLate = rmmissing(corrChoiceLate); corrNChoiceLate = rmmissing(corrNChoiceLate);
subplot_tight(2,6,[4,5,6])
hold on
plot(times,ratsCorrectMid','Color',[1,0.604,0,0.2],'LineWidth',1.5)
plot(times,ratsCorrectMatch','Color',[0.133,0.545,0.133,0.2],'LineWidth',1.5)
plot(times,ratsCorrectNonMatch','Color',[0.796,0.377,0.082,0.2],'LineWidth',1.5)
plot(times,mean(ratsCorrectMid),'Color',[1,0.604,0],'LineWidth',2)
plot(times,mean(ratsCorrectMatch),'Color',[0.133,0.545,0.133],'LineWidth',2)
plot(times,mean(ratsCorrectNonMatch),'Color',[0.796,0.377,0.082],'LineWidth',2)
xlabel('Time Since Start of ITI (s)')
title('Post-Correct Trial')
ylabel('Poke Density')
set(gca,'fontsize',18)
yticks(0:0.3:0.9)
xlim([0,7])
ylim([0,0.9])
axes('Position',[.11 .75 .1 .1])
box on
al_goodplot2({corrChoice,corrNChoice}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.133,0.545,0.133;0.796,0.377,0.082]);
xticks([])
xlabel('Early Window')
ylabel(["Trial Proportion With","More of Poke Type"])
axes('Position',[.275 .75 .1 .1])
box on
al_goodplot2({corrChoiceLate,corrNChoiceLate}', 'pos', [1,1], 'type', {'left','right'},'boxw',0.4,'col',[0.133,0.545,0.133;0.796,0.377,0.082]);
xticks([])
xlabel('Late Window')
%% c
S = max(full_tbl.Session);
lcFront = zeros(S,1);
lcRear = zeros(S,1);
nlcFront = zeros(S,1);
nlcRear = zeros(S,1);
for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    lcFront(i) = nanmean(subtbl.Choice((subtbl.Light==subtbl.Choice)&(subtbl.frontMostLast==1)));
    lcRear(i) = nanmean(subtbl.Choice((subtbl.Light==subtbl.Choice)&(subtbl.rearMostLast==1)));
    nlcFront(i) = nanmean(subtbl.Choice((subtbl.Light~=subtbl.Choice)&(subtbl.Choice~=-1)&(subtbl.frontMostLast==1)));
    nlcRear(i) = nanmean(subtbl.Choice((subtbl.Light~=subtbl.Choice)&(subtbl.Choice~=-1)&(subtbl.rearMostLast==1)));
end
lcFront=-(rmmissing(lcFront)-2);
lcRear=-(rmmissing(lcRear)-2);
nlcFront=-(rmmissing(nlcFront)-2);
nlcRear=-(rmmissing(nlcRear)-2);
subplot_tight(2,6,[7,8])
handles = al_goodplot2({lcFront,lcRear,nlcFront,nlcRear}', 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4,'col',[0.31,0.73,0.90;0.83,0.07,0.82;0.31,0.73,0.90;0.83,0.07,0.82]);
xticks([1,2])
xticklabels(["Light Chosen", "Light Not Chosen"])
ylabel('Proportion Front Choice')
yticks([0,0.5,1])
xlim([0.5,2.5])
legend([handles{1},handles{2}],["Front Most Late", "Rear Most Late"],'Location','southwest')
set(gca,'fontsize',18)
box off
title('Pre-Trial')
%% d
sub_tbl = full_tbl(full_tbl.RT>0,:);
sub_tbl.choiceMostLast = (sub_tbl.frontMostLast==1&sub_tbl.Choice==1)|(sub_tbl.rearMostLast==1&sub_tbl.Choice==2);
sub_tbl.nonchoiceMostLast = (sub_tbl.frontMostLast==1&sub_tbl.Choice==2)|(sub_tbl.rearMostLast==1&sub_tbl.Choice==1);
itimdl = fitglme(sub_tbl, 'RT ~ midPokeLast + midPokeLast:choiceMostLast + midPokeLast:nonchoiceMostLast + (1|Rat)', 'Distribution', 'gamma', 'link', 'identity');
D = designMatrix(itimdl, 'random');
B = randomEffects(itimdl);
B2 = fixedEffects(itimdl);
intercepts = D * B + B2(1);
sub_tbl.RT = sub_tbl.RT - intercepts;
rtmid = sub_tbl.RT(sub_tbl.midPokeLast==1&sub_tbl.frontMostLast==0&sub_tbl.rearMostLast==0);
rtchoice = sub_tbl.RT(sub_tbl.midPokeLast==1&sub_tbl.choiceMostLast);
rtnonchoice = sub_tbl.RT(sub_tbl.midPokeLast==1&sub_tbl.nonchoiceMostLast);
subplot_tight(2,6,[9,10])
handles = al_goodplot2({rtmid,rtchoice,rtnonchoice}','boxw',0.4,'col',[1,0.604,0;0.133,0.545,0.133;0.796,0.377,0.082]);
ylim([-0.4,0.4])
xticks(1:3)
yticks(-0.4:0.2:0.4)
XTickString = cell(3,1);
labels = {'Mid Poke', 'Choice +', 'Non-choice +';'','Mid Poke','Mid Poke'};
labels = strjust(pad(labels),'center');
tickLabels = strtrim(sprintf('%s\\newline%s\n', labels{:}));
set(gca,'TickLabelInterpreter','tex','XTickLabel',tickLabels)
xlim([0.5,3.5])
set(gca,'fontsize',18)
ylabel('\Delta RT (s)')
title('Pre-Trial')
%% e
S = max(full_tbl.Session);
Mid = nan(S,1);
Choice = nan(S,1);
NonChoice = nan(S,1);
stim = false(S,1);
for i=1:S
    subtbl = full_tbl(full_tbl.Session==i,:);
    stim(i) = subtbl.Stim(1) == 1;
    Mid(i) = nanmean(subtbl.MidPokeLateITI);
    Choice(i) = nanmean(subtbl.MorePreChoiceLateITI);
    NonChoice(i) = nanmean(subtbl.MoreNonePreChoiceLateITI);
end
valid = ~isnan(Mid);
subplot_tight(2,6,[11,12])
handles=al_goodplot2({Mid(valid&~stim),Mid(valid&stim),Choice(valid&~stim),Choice(valid&stim),NonChoice(valid&~stim),NonChoice(valid&stim)}', 'pos', [1,1,2,2,3,3], 'type', {'left','right','left','right','left','right'},'boxw',0.4,'col',[1,0.604,0;1,0.604,0;0.133,0.545,0.133;0.133,0.545,0.133;0.796,0.377,0.082;0.796,0.377,0.082]);
H = hatchfill(handles{2},'cross',45,5,[1,0.604,0],1.5);
H.Color = [0,0,0,0.2];
H = hatchfill(handles{4},'cross',45,5,[0.133,0.545,0.133],1.5);
H.Color = [0,0,0,0.2];
H = hatchfill(handles{6},'cross',45,5,[0.796,0.377,0.082],1.5);
H.Color = [0,0,0,0.2];
ylabel('Proportion of Trials with Poke')
ylim([0,1])
yticks([0,0.5,1])
xticks(1:3)
xticklabels(["Mid", "Choice", "Non-Choice"])
xlim([0.5,3.5])
set(gca,'fontsize',18)
ax = gca;
ax.Children = circshift(ax.Children,-3);

function corrected = boundary_correction(uncorrected,segstart,segend)
    corrected = uncorrected(segstart:segend);
    corrected(1:end) = corrected(1:end)+flip(uncorrected(1:segstart));
    corrected(1:end) = corrected(1:end)+flip(uncorrected(segend:end));
end