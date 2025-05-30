function testingvigorchange(lowspeed_tbt,highspeed_tbt,alltbt,trialTypes,metadata)

cueatind=80;

%% Get kinematic tracking
addtotrial2=' & (trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1)';
% After success
whichEventName='cued success'; timeWindowToCountAsEventReach=[0 1]; whichReachName='reachBatch_success_reachStarts';
[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n successes:'); 
disp(nansum(logictogether));
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([~logictogether; true],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[X,Y,trial1Z_success]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([true; ~logictogether],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[X,Y,trial2Z_success]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
close all;
% After failure
whichEventName='all cued failures'; timeWindowToCountAsEventReach=[0 1]; whichReachName='anyFail';
[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n failures:'); 
disp(nansum(logictogether));
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([~logictogether; true],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[X,Y,trial1Z_failure]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([true; ~logictogether],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[X,Y,trial2Z_failure]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
close all;
figure(); plot(nanmean(trial1Z_success,1),'Color','k'); hold on; plot(nanmean(trial2Z_success,1)+15,'Color','r'); title('trial 1 then 2 after success');
figure(); plot(nanmean(trial1Z_failure,1),'Color','k'); hold on; plot(nanmean(trial2Z_failure,1)+15,'Color','r'); title('trial 1 then 2 after failure');


%% Get time to dislodge pellet
addtotrial2=' & (trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1)';
% After success
whichEventName='cued success'; timeWindowToCountAsEventReach=[0 1]; whichReachName='reachBatch_success_reachStarts';
[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n successes:'); 
disp(nansum(logictogether));
f=find(logictogether);
[reachtopelletdislodge_aftersucc,issuccess_aftersucc,isdrop_aftersucc,reachtopelletdislodge_secondtrial_aftersucc,issuccess_secondtrial_aftersucc,isdrop_secondtrial_aftersucc]=gettimetopelletdislodge(f,alltbt,trialTypes,cueatind)
% Just reaction times
[rt_aftersucc,issuccess_aftersucc,isdrop_aftersucc,rt_secondtrial_aftersucc,issuccess_secondtrial_aftersucc,isdrop_secondtrial_aftersucc]=getreactiontime(f,alltbt,trialTypes,cueatind)
% After failure
whichEventName='all cued failures'; timeWindowToCountAsEventReach=[0 1]; whichReachName='anyFail';
[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n failures:'); 
disp(nansum(logictogether));
f=find(logictogether);
[reachtopelletdislodge_afterfail,issuccess_afterfail,isdrop_afterfail,reachtopelletdislodge_secondtrial_afterfail,issuccess_secondtrial_afterfail,isdrop_secondtrial_afterfail]=gettimetopelletdislodge(f,alltbt,trialTypes,cueatind)
% Just reaction times
[rt_afterfail,issuccess_afterfail,isdrop_afterfail,rt_secondtrial_afterfail,issuccess_secondtrial_afterfail,isdrop_secondtrial_afterfail]=getreactiontime(f,alltbt,trialTypes,cueatind)

disp(rt_aftersucc)
disp(rt_afterfail)

end

function [rt,issuccess,isdrop,rt_secondtrial,issuccess_secondtrial,isdrop_secondtrial]=getreactiontime(f,alltbt,trialTypes,cueatind)

rt=nan(1,length(f));
issuccess=nan(1,length(f));
isdrop=nan(1,length(f));
rt_secondtrial=nan(1,length(f));
issuccess_secondtrial=nan(1,length(f));
isdrop_secondtrial=nan(1,length(f));
for i=1:length(f)
    % first trial in pair
    reachtemp=alltbt.all_reachBatch(f(i),:);
    f1=find(reachtemp(cueatind:end)>0.5,1,'first');
    if isempty(f1)
        continue
    end
    rt(i)=f1*0.035;
    issuccess(i)=trialTypes.after_cue_success(f(i));
    isdrop(i)=trialTypes.after_cue_drop(f(i));
    % second trial in pair
    reachtemp=alltbt.all_reachBatch(f(i)+1,:);
    f1=find(reachtemp(cueatind:end)>0.5,1,'first');
    if isempty(f1)
        continue
    end
    rt_secondtrial(i)=f1*0.035;
    issuccess_secondtrial(i)=trialTypes.after_cue_success(f(i)+1);
    isdrop_secondtrial(i)=trialTypes.after_cue_drop(f(i)+1);
end

end

function [reachtopelletdislodge,issuccess,isdrop,reachtopelletdislodge_secondtrial,issuccess_secondtrial,isdrop_secondtrial]=gettimetopelletdislodge(f,alltbt,trialTypes,cueatind)

reachtopelletdislodge=nan(1,length(f));
issuccess=nan(1,length(f));
isdrop=nan(1,length(f));
reachtopelletdislodge_secondtrial=nan(1,length(f));
issuccess_secondtrial=nan(1,length(f));
isdrop_secondtrial=nan(1,length(f));
for i=1:length(f)
    % first trial in pair
    reachtemp=alltbt.all_reachBatch(f(i),:);
    pellettemp=alltbt.pelletPresent(f(i),:);
    f1=find(pellettemp(cueatind:end)<0.5,1,'first');
    if isempty(f1)
        continue
    end
    %disp(f1)
    f2=find(reachtemp(cueatind:end)>0.5,1,'first');
    if isempty(f2)
        continue
    end
    %disp(f2)
    %size(reachtopelletdislodge(i))
    reachtopelletdislodge(i)=f2-f1;
    issuccess(i)=trialTypes.after_cue_success(f(i));
    isdrop(i)=trialTypes.after_cue_drop(f(i));
    % second trial in pair
    reachtemp=alltbt.all_reachBatch(f(i)+1,:);
    pellettemp=alltbt.pelletPresent(f(i)+1,:);
    f1=find(pellettemp(cueatind:end)<0.5,1,'first');
    if isempty(f1)
        continue
    end
    %disp(f1)
    f2=find(reachtemp(cueatind:end)>0.5,1,'first');
    if isempty(f2)
        continue
    end
    %disp(f2)
    %size(reachtopelletdislodge(i))
    reachtopelletdislodge_secondtrial(i)=f2-f1;
    issuccess_secondtrial(i)=trialTypes.after_cue_success(f(i)+1);
    isdrop_secondtrial(i)=trialTypes.after_cue_drop(f(i)+1);
end

end