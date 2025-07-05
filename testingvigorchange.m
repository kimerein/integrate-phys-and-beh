function testingvigorchange(lowspeed_tbt,highspeed_tbt,alltbt,trialTypes,metadata,saveTo)

% see Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\high and low speed aligned\first video\processed lowspeed_tbt to trialTypes
% for processing lowspeed_tbt to trialTypes

cueatind=80;

save(fullfile(saveTo,'alltbt.mat'),"alltbt");
save(fullfile(saveTo,'trialTypes.mat'),"trialTypes");
save(fullfile(saveTo,'metadata.mat'),"metadata");

%% Get kinematic tracking
addtotrial2=' & (trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1)';
% After success
whichEventName='cued success'; timeWindowToCountAsEventReach=[0 1]; whichReachName='reachBatch_success_reachStarts';
% [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
% trial2=[trial2 addtotrial2];
trial1=['(trialTypes.optoGroup~=1 & trialTypes.after_cue_success==1 & trialTypes.led==0 & trialTypes.led_1forward==0)'];        
trial2=['(trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.optoGroup~=1'];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n successes:'); 
disp(nansum(logictogether));
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([~logictogether; true],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[trial1X_success,trial1Y_success,trial1Z_success]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([true; ~logictogether],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[trial2X_success,trial2Y_success,trial2Z_success]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
save(fullfile(saveTo,'trial1X_success.mat'),"trial1X_success");
save(fullfile(saveTo,'trial1Y_success.mat'),"trial1Y_success");
save(fullfile(saveTo,'trial1Z_success.mat'),"trial1Z_success");
save(fullfile(saveTo,'trial2X_success.mat'),"trial2X_success");
save(fullfile(saveTo,'trial2Y_success.mat'),"trial2Y_success");
save(fullfile(saveTo,'trial2Z_success.mat'),"trial2Z_success");
close all;
% After failure
whichEventName='all cued failures'; timeWindowToCountAsEventReach=[0 1]; whichReachName='anyFail';
% [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
% trial2=[trial2 addtotrial2];
trial1=['trialTypes.optoGroup~=1 & (trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.led==0 & trialTypes.led_1forward==0']; 
trial2=['(trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.optoGroup~=1'];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n failures:'); 
disp(nansum(logictogether));
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([~logictogether; true],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[trial1X_failure,trial1Y_failure,trial1Z_failure]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
lowspeed_tbt.currentReachTypeToAnalyze=lowspeed_tbt.all_reachBatch;
lowspeed_tbt.currentReachTypeToAnalyze([true; ~logictogether],:)=0;
[allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,'currentReachTypeToAnalyze',1,'Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190815\March_A\DLC vids','vid_2019-08-15-145458-0000DLC_resnet50_Testing2DJan4shuffle1_500003D.mat',255,'all',0);
[trial2X_failure,trial2Y_failure,trial2Z_failure]=plotReachTrajectories(allX,allY,allZ,allX_from_under,reachTrajTimes,12,255);
save(fullfile(saveTo,'trial1X_failure.mat'),"trial1X_failure");
save(fullfile(saveTo,'trial1Y_failure.mat'),"trial1Y_failure");
save(fullfile(saveTo,'trial1Z_failure.mat'),"trial1Z_failure");
save(fullfile(saveTo,'trial2X_failure.mat'),"trial2X_failure");
save(fullfile(saveTo,'trial2Y_failure.mat'),"trial2Y_failure");
save(fullfile(saveTo,'trial2Z_failure.mat'),"trial2Z_failure");
close all;
figure(); plot(nanmean(trial1Z_success,1),'Color','k'); hold on; plot(nanmean(trial2Z_success,1)+15,'Color','r'); title('trial 1 then 2 after success');
figure(); plot(nanmean(trial1Z_failure,1),'Color','k'); hold on; plot(nanmean(trial2Z_failure,1)+15,'Color','r'); title('trial 1 then 2 after failure');


%% Get time to dislodge pellet
% addtotrial2=' & (trialTypes.after_cue_success==1 | trialTypes.after_cue_drop==1)';
addtotrial2=' & (trialTypes.after_cue_success==1)';
% After success
whichEventName='cued success'; timeWindowToCountAsEventReach=[0 1]; whichReachName='reachBatch_success_reachStarts';
% [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial1=['(trialTypes.optoGroup~=1 & trialTypes.after_cue_success==1 & trialTypes.led==0 & trialTypes.led_1forward==0)'];        
trial2=['trialTypes.after_cue_success==1 & trialTypes.optoGroup~=1'];
% trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n successes:'); 
disp(nansum(logictogether));
f=find(logictogether);
[reachtopelletdislodge_aftersucc,issuccess_aftersucc,isdrop_aftersucc,reachtopelletdislodge_secondtrial_aftersucc,issuccess_secondtrial_aftersucc,isdrop_secondtrial_aftersucc]=gettimetopelletdislodge(f,alltbt,trialTypes,cueatind)
save(fullfile(saveTo,'reachtopelletdislodge_aftersucc.mat'),"reachtopelletdislodge_aftersucc");
save(fullfile(saveTo,'issuccess_aftersucc.mat'),"issuccess_aftersucc");
save(fullfile(saveTo,'isdrop_aftersucc.mat'),"isdrop_aftersucc");
save(fullfile(saveTo,'reachtopelletdislodge_secondtrial_aftersucc.mat'),"reachtopelletdislodge_secondtrial_aftersucc");
save(fullfile(saveTo,'issuccess_secondtrial_aftersucc.mat'),"issuccess_secondtrial_aftersucc");
save(fullfile(saveTo,'isdrop_secondtrial_aftersucc.mat'),"isdrop_secondtrial_aftersucc");
% Just reaction times
[rt_aftersucc,issuccess_aftersucc,isdrop_aftersucc,rt_secondtrial_aftersucc,issuccess_secondtrial_aftersucc,isdrop_secondtrial_aftersucc]=getreactiontime(f,alltbt,trialTypes,cueatind)
save(fullfile(saveTo,'rt_aftersucc.mat'),"rt_aftersucc");
save(fullfile(saveTo,'issuccess_aftersucc.mat'),"issuccess_aftersucc");
save(fullfile(saveTo,'isdrop_aftersucc.mat'),"isdrop_aftersucc");
save(fullfile(saveTo,'rt_secondtrial_aftersucc.mat'),"rt_secondtrial_aftersucc");
save(fullfile(saveTo,'issuccess_secondtrial_aftersucc.mat'),"issuccess_secondtrial_aftersucc");
save(fullfile(saveTo,'isdrop_secondtrial_aftersucc.mat'),"isdrop_secondtrial_aftersucc");
% After failure
whichEventName='all cued failures'; timeWindowToCountAsEventReach=[0 1]; whichReachName='anyFail';
% [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,[],whichEventName,timeWindowToCountAsEventReach,whichReachName);
trial1=['trialTypes.optoGroup~=1 & (trialTypes.after_cue_drop==1 | trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.led==0 & trialTypes.led_1forward==0']; 
trial2=['trialTypes.after_cue_success==1 & trialTypes.optoGroup~=1'];
% trial2=[trial2 addtotrial2];
logic1=eval(trial1); logic2=eval(trial2); % these are only no LED trial pairs
logictogether=logic1(1:end-1) & logic2(2:end);
disp('n failures:'); 
disp(nansum(logictogether));
f=find(logictogether);
[reachtopelletdislodge_afterfail,issuccess_afterfail,isdrop_afterfail,reachtopelletdislodge_secondtrial_afterfail,issuccess_secondtrial_afterfail,isdrop_secondtrial_afterfail]=gettimetopelletdislodge(f,alltbt,trialTypes,cueatind)
save(fullfile(saveTo,'reachtopelletdislodge_afterfail.mat'),"reachtopelletdislodge_afterfail");
save(fullfile(saveTo,'issuccess_afterfail.mat'),"issuccess_afterfail");
save(fullfile(saveTo,'isdrop_afterfail.mat'),"isdrop_afterfail");
save(fullfile(saveTo,'reachtopelletdislodge_secondtrial_afterfail.mat'),"reachtopelletdislodge_secondtrial_afterfail");
save(fullfile(saveTo,'issuccess_secondtrial_afterfail.mat'),"issuccess_secondtrial_afterfail");
save(fullfile(saveTo,'isdrop_secondtrial_afterfail.mat'),"isdrop_secondtrial_afterfail");
% Just reaction times
[rt_afterfail,issuccess_afterfail,isdrop_afterfail,rt_secondtrial_afterfail,issuccess_secondtrial_afterfail,isdrop_secondtrial_afterfail]=getreactiontime(f,alltbt,trialTypes,cueatind)
save(fullfile(saveTo,'rt_afterfail.mat'),"rt_afterfail");
save(fullfile(saveTo,'issuccess_afterfail.mat'),"issuccess_afterfail");
save(fullfile(saveTo,'isdrop_afterfail.mat'),"isdrop_afterfail");
save(fullfile(saveTo,'rt_secondtrial_afterfail.mat'),"rt_secondtrial_afterfail");
save(fullfile(saveTo,'issuccess_secondtrial_afterfail.mat'),"issuccess_secondtrial_afterfail");
save(fullfile(saveTo,'isdrop_secondtrial_afterfail.mat'),"isdrop_secondtrial_afterfail");
disp(rt_secondtrial_aftersucc)
disp(rt_secondtrial_afterfail)

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