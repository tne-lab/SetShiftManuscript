%% SetShift
% Object defining an instance of the SetShift task
% Describes a sequence of rules (F/R/L) and lights
% Tracks the position of a player in the task

classdef SetShift < handle
    
    properties
        complete % Boolean indicating if the task is finished
        ruleNum % Index of current rule (1-8)
        blockNum % Index of current block (1-5)
        lightSeq % Boolean matrix (8,5) indicating if the light was on in front
        ruleSeq % Vector indicating sequence of rules(F=1,R=2,L=3)
        trialsCompleted % Counts the number of trials that have been done
    end
    
    methods
        
        % Create a new setshift task
        % Input:
        %   behavior = processed SetShift behavioral data (CSV)
        function setshift = SetShift(behavior, model)
            if nargin <2
                model=false;
            end
            % Initialize task tracking data
            setshift.complete=false;
            setshift.ruleNum=1;
            setshift.blockNum=1;
            setshift.trialsCompleted=0;
            % If sample data are provided
            if nargin < 1
                % Create a random task object
                % Initialize the random sequence of light activations
                setshift.lightSeq=rand(8,5);
                setshift.lightSeq=setshift.lightSeq<=median(setshift.lightSeq(:));
                % Initialize the rule sequence with alternating Light and Front/Rear
                rules=[1,1,2,2];
                lFirst=rand<0.5;
                setshift.ruleSeq=zeros(8,1);
                if lFirst
                    setshift.ruleSeq(1:2:end)=3;
                    setshift.ruleSeq(2:2:end)=rules(randsample(4,4));
                else
                    setshift.ruleSeq(2:2:end)=3;
                    setshift.ruleSeq(1:2:end)=rules(randsample(4,4));
                end
            elseif ~model
                % Create a task object corresponding to the data
                % Identify the rule sequence
                ruleNames=string({behavior.rule});
                rules=3*ones(size(ruleNames));
                rules(startsWith(ruleNames,'F'))=1;
                rules(startsWith(ruleNames,'R'))=2;
                changes=find(diff(rules)~=0);
                setshift.ruleSeq=[rules(1),rules(changes+1)];
                setshift.lightSeq=false(8,5);
                % Identify the light sequence
                ind=1;
                for i=1:length(behavior)
                    ruleInd=ind+find(setshift.ruleSeq(ind:end)==rules(i),1)-1;
                    if ruleInd > ind
                        ind=ind+1;
                    end
                    if ~isnan(behavior(i).frontChoice)
                        if behavior(i).light
                            setshift.lightSeq(ruleInd,behavior(i).block)=true;
                        else
                            setshift.lightSeq(ruleInd,behavior(i).block)=false;
                        end
                    end
                end
            else
                rules=[behavior.rule];
                changes=find(diff(rules)~=0);
                setshift.ruleSeq=[rules(1),rules(changes+1)];
                setshift.lightSeq=false(8,5);
                % Identify the light sequence
                ind=1;
                for i=1:length(behavior)
                    ruleInd=ind+find(setshift.ruleSeq(ind:end)==rules(i),1)-1;
                    if ruleInd > ind
                        ind=ind+1;
                    end
                    if ~isnan(behavior(i).frontChoice)
                        if behavior(i).light
                            setshift.lightSeq(ruleInd,behavior(i).blockNum)=true;
                        else
                            setshift.lightSeq(ruleInd,behavior(i).blockNum)=false;
                        end
                    end
                end
            end
        end
        
        % Choose a side
        % Input:
        %   front = boolean indicating whethere the front option was chosen
        %   (NaN indicates non-decision)
        function correct = choose(setshift,front)
            setshift.trialsCompleted=setshift.trialsCompleted+1;
            if setshift.trialsCompleted==1000
                setshift.complete=true;
            end
            % If a decision was made
            if ~isnan(front)
                % If front rule and choose front -> correct
                if setshift.ruleSeq(setshift.ruleNum)==1&&front
                    correct=true;
                % If rear rule and choose rear -> correct
                elseif setshift.ruleSeq(setshift.ruleNum)==2&&~front
                    correct=true;
                % If light rule and choose light -> correct
                elseif setshift.ruleSeq(setshift.ruleNum)==3&&setshift.getLightDir()==front
                    correct=true;
                else
                % Otherwise -> incorrect
                    correct=false;
                end
                
                % If decision was correct either complete the task, move to
                % next rule, or increment block number
                if correct
                    if setshift.blockNum==5&&setshift.ruleNum==8
                        setshift.complete=true;
                    elseif setshift.blockNum<5
                        setshift.blockNum=setshift.blockNum+1;
                    else
                        setshift.ruleNum=setshift.ruleNum+1;
                        setshift.blockNum=1;
                    end
                else % If incorrect, reset block
                    setshift.blockNum=1;
                end
            else
                % If no decision, return incorrect, reset block
                correct=false;
                setshift.blockNum=1;
            end
        end
        
        % Return side of light for the current trial
        function front = getLightDir(setshift)
            front = setshift.lightSeq(setshift.ruleNum,setshift.blockNum);
        end
        
        % Return rule for the current trial
        function rule = getRule(setshift)
            rule = setshift.ruleSeq(setshift.ruleNum);
        end
        
        % Restart the task
        function reset(setshift)
            setshift.complete=false;
            setshift.ruleNum=1;
            setshift.blockNum=1;
        end
    end
end