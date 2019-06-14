function RT_regression_model(alltbt,metadata,trialTypes)

% RT(n) = linear regression
%
% Regressors
%
% RT(n-i)
% touched pellet
% consumed pellet

doInteractions=1; % will add interaction terms to regression if 1
throwOutLongRTs=1; % will throw out all long reaction times if 1

% get reaction times
[rts,cueInd]=getRT(alltbt,'all_reachBatch','cueZone_onVoff',metadata,0,1);
reactionTimes.rts=rts';
f{1}='rts';
reactionTimes=getNback(reactionTimes,f);

% settings for regression
settings=regressionSettings2(alltbt,metadata,trialTypes,reactionTimes,cueInd);

% plot correlations
% plotCorrelationStructure(settings,alltbt,metadata,trialTypes,reactionTimes.rts);

if throwOutLongRTs==1
    reactionTimes.rts(reactionTimes.rts>5)=nan;
end

% calculate regression
calcRegress(settings,reactionTimes.rts,doInteractions);

end

function calcRegress(settings,dependent,doInteractions)

% Y_i = X_i * beta + epsilon_i;

rg=settings.rg;

for i=1:length(rg)
    if size(rg{i},2)>1
        % make it a column vector
        X(:,i)=rg{i}';
    else
        X(:,i)=rg{i};
    end
end

if doInteractions==1
    mdl=fitlm(X,dependent,'interactions');
    rg_names{1}='Intercept';
    j=2;
    for i=1:length(settings.rg_names)
        rg_names{j}=settings.rg_names{i};
        j=j+1;
    end
    % add interaction terms
    for i=1:length(settings.rg_names)-1
        for k=i+1:length(settings.rg_names)
            rg_names{j}=[settings.rg_names{i} ' : ' settings.rg_names{k}];
            j=j+1;
        end
    end           
else
    mdl=fitlm(X,dependent);
    rg_names{1}='Intercept';
    j=1;
    for i=1:length(settings.rg_names)
        rg_names{j}=settings.rg_names{i};
        j=j+1;
    end
end

mdl

% Sort by beta
betas=mdl.Coefficients.Estimate;
se=mdl.Coefficients.SE;
tStat=mdl.Coefficients.tStat;
pValue=mdl.Coefficients.pValue;

[~,si]=sort(abs(betas),'descend');
T=table(rg_names(si)',betas(si),se(si),tStat(si),pValue(si),'VariableNames',{'Coeff_names','betas','SE','tStat','pValue'});
disp(T);

end

function out=getNback(out,f)

for i=1:length(f)
    if ~(isnumeric(out.(f{i})) || islogical(out.(f{i})))
        continue
    end
    temp=out.(f{i});
    % 1 back
    newfieldname=[f{i} '_1back'];
    out.(newfieldname)=[nan; temp(1:end-1)];
    % 2 back
    newfieldname=[f{i} '_2back'];
    out.(newfieldname)=[nan; nan; temp(1:end-2)];
    % 3 back
    newfieldname=[f{i} '_3back'];
    out.(newfieldname)=[nan; nan; nan; temp(1:end-3)];
    % 4 back
    newfieldname=[f{i} '_4back'];
    out.(newfieldname)=[nan; nan; nan; nan; temp(1:end-4)];
    % 5 back
    newfieldname=[f{i} '_5back'];
    out.(newfieldname)=[nan; nan; nan; nan; nan; temp(1:end-5)];
    
    % 1 forward
    newfieldname=[f{i} '_1forward'];
    out.(newfieldname)=[temp(2:end); nan];
    % 2 forward
    newfieldname=[f{i} '_2forward'];
    out.(newfieldname)=[temp(3:end); nan; nan];
    % 3 forward
    newfieldname=[f{i} '_3forward'];
    out.(newfieldname)=[temp(4:end); nan; nan; nan];
    % 4 forward
    newfieldname=[f{i} '_4forward'];
    out.(newfieldname)=[temp(5:end); nan; nan; nan; nan];
    % 5 forward
    newfieldname=[f{i} '_5forward'];
    out.(newfieldname)=[temp(6:end); nan; nan; nan; nan; nan];
end

end

function [reactionTimes,maind]=getRT(tbt,whichReach,useAsCue,metadata,zscore_RTs,longRT_ifNoReach)

% cue ind
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

% which reaches to get
temp=tbt.(whichReach);
newtemp=nan(size(tbt.(whichReach)));

% get only first reaches
firstreaches=nan(1,size(temp,1));
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    if longRT_ifNoReach==1 && isempty(fi)
        fi=size(temp,2);
        firstreaches(i)=fi;
    end
    if ~isempty(fi)
        newtemp(i,:)=zeros(size(newtemp(i,:)));
        newtemp(i,fi)=1;
    end
end
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));
tbt.firstReachesAfterCueOnset=newtemp;

if zscore_RTs==1
    reactionTimes=Zscore_by_session(reactionTimes,metadata.sessid);
end

end

function settings=regressionSettings1(alltbt,metadata,trialTypes,reactionTimes)

% choose regressors
rg{1}=trialTypes.touched_pellet_1back;
rg{2}=trialTypes.consumed_pellet_1back;
rg{3}=reactionTimes.rts_1back;
rg{4}=reactionTimes.rts_2back;
rg{5}=reactionTimes.rts_3back;
rg{6}=reactionTimes.rts_4back;
rg{7}=reactionTimes.rts_5back;

rg_names{1}='touched pellet 1 trial back';
rg_names{2}='consumed pellet 1 trial back';
rg_names{3}='reaction time 1 trial back';
rg_names{4}='reaction time 2 trials back';
rg_names{5}='reaction time 3 trials back';
rg_names{6}='reaction time 4 trials back';
rg_names{7}='reaction time 5 trials back';

settings.rg=rg;
settings.rg_names=rg_names;

end

function settings=regressionSettings2(alltbt,metadata,trialTypes,reactionTimes,cueInd)

% choose regressors
rg{1}=trialTypes.touched_pellet_1back;
rg{2}=trialTypes.consumed_pellet_1back;
rg{3}=reactionTimes.rts_1back;
% rg{4}=reactionTimes.rts_2back;
% rg{5}=reactionTimes.rts_3back;
% rg{6}=reactionTimes.rts_4back;
% rg{7}=reactionTimes.rts_5back;
rg{4}=trialTypes.chewing_at_trial_start;
% get rid of reaches during chewing
% temp=alltbt.all_reachBatch;
% temp(alltbt.isChewing>0.5)=0;
% rg{5}=sum(temp(:,1:cueInd-1),2); % baseline reach rate
rg{5}=sum(alltbt.all_reachBatch(:,1:cueInd-1),2); % baseline reach rate
rg{6}=trialTypes.led_1back;

rg_names{1}='touched pellet 1 trial back';
rg_names{2}='consumed pellet 1 trial back';
rg_names{3}='reaction time 1 trial back';
% rg_names{4}='reaction time 2 trials back';
% rg_names{5}='reaction time 3 trials back';
% rg_names{6}='reaction time 4 trials back';
% rg_names{7}='reaction time 5 trials back';
rg_names{4}='chewing at this trial start';
rg_names{5}='number of reaches before cue';
rg_names{6}='led on 1 trial back';

settings.rg=rg;
settings.rg_names=rg_names;

end

function settings=regressionSettings3(alltbt,metadata,trialTypes,reactionTimes,cueInd)

% choose regressors
rg{1}=trialTypes.touched_pellet_1back;
rg{2}=reactionTimes.rts_1back;
rg{3}=trialTypes.chewing_at_trial_start;
rg{4}=sum(alltbt.all_reachBatch(:,1:cueInd-1),2); % baseline reach rate
rg{5}=trialTypes.led_1back;

rg_names{1}='touched pellet 1 trial back';
rg_names{2}='reaction time 1 trial back';
rg_names{3}='chewing at this trial start';
rg_names{4}='number of reaches before cue';
rg_names{5}='led on 1 trial back';

settings.rg=rg;
settings.rg_names=rg_names;

end

function accumX_perSess=getAccumulatedXinSession(metadata,x)



end

function plotCorrelationStructure(settings,alltbt,metadata,trialTypes,dependent)

% get regressors
rg=settings.rg;

% plot correlation of each with dependent variable
for i=1:length(rg)
    figure();
    scatter(rg{i},dependent);
    xlabel(['regressor ' num2str(i)]);
    ylabel(['dependent']);
end

end