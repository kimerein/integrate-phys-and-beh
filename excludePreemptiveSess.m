function [alltbt,metadata,trialTypes]=excludePreemptiveSess(alltbt,metadata,trialTypes,preCueWindow,postPreCueWindow,threshCutoff)

usess=unique(metadata.sessid);
fraction_of_precue=nan(1,length(usess));
exclude=nan(1,length(usess));
for i=1:length(usess)
    % filt to usess(i)
    alltbt.mouseid=metadata.mouseid;
    alltbt.sessid=metadata.sessid;
    trialTypes.sessid=metadata.sessid;
    tbt_filter.sortField='sessid';
    tbt_filter.range_values=[usess(i)-0.5 usess(i)+0.5];
    tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
    temp=tbt_filter.name;
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    tbt_filter.name=temp;
    tbt_filter.clock_progress=true;
    [filt_alltbt,filt_trialTypes,filt_metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

    [exclude(i),fraction_of_precue(i)]=testWhetherSessPreempt(filt_alltbt,filt_metadata,filt_trialTypes,preCueWindow,postPreCueWindow,threshCutoff);
end
excludeSess=nan(1,length(metadata.sessid));
for i=1:length(usess)
    excludeSess(metadata.sessid==usess(i))=exclude(i);
end

figure();
histogram(fraction_of_precue,100);
hold on;
line([threshCutoff threshCutoff],[0 10],'Color','r');

% filter for no preempt
alltbt.excludeSess=excludeSess;
trialTypes.excludeSess=excludeSess;
metadata.excludeSess=excludeSess;
tbt_filter.sortField='excludeSess';
tbt_filter.range_values=[-0.5 0.5];
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

end

function [exclude,fraction_of_precue]=testWhetherSessPreempt(alltbt,metadata,trialTypes,preCueWindow,postPreCueWindow,threshCutoff)

test.nInSequence=[2]; 
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
tbt_filter.name='any';
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/anything']; 
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; skipCorrected=true;
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 

returnThis=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching',true);
temp=returnThis.data2_mean{1};
precue=mean(temp(returnThis.time_for_x>=preCueWindow(1) & returnThis.time_for_x<=preCueWindow(2)));
postprecue=mean(temp(returnThis.time_for_x>=postPreCueWindow(1) & returnThis.time_for_x<=postPreCueWindow(2)));
if postprecue>threshCutoff*precue
    exclude=1;
else
    exclude=0;
end
fraction_of_precue=postprecue/precue;

end