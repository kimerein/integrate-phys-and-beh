function pickReachTrajectories(lowspeed_tbt,highspeed_tbt,whichReachField,fromWhichArg,DLCoutput_location,fps)

% fps is frames per second of high speed movie

timeBeforeReach=0.2;
timeAfterReach=1;
nReachesFromEachTrial='all';
noReachesBeforeTime=-0.3; % wrt cue, to be sure have enough frames around reach
noReachesAfterTime=7.5; % wrt cue, to be sure have enough frames around reach

if fromWhichArg==1
    reaches=lowspeed_tbt.(whichReachField);
    reachtimes=lowspeed_tbt.times_wrt_trial_start;
    temp=nanmean(reachtimes,1);
    [~,ma]=nanmax(nanmean(lowspeed_tbt.cueZone_onVoff,1));
    avcuetime=temp(ma);
elseif fromWhichArg==2
    reaches=highspeed_tbt.(whichReachField);
    reachtimes=highspeed_tbt.times;
    temp=nanmean(reachtimes,1);
    tempcue=nanmean(highspeed_tbt.cue,1);
    ma=find(tempcue>0.5,1,'first');
    avcuetime=temp(ma);
else
    error('fromWhichArg must be 1 or 2');
end

% how many frames before reach, after reach
timestep_hs=1/fps;
framesBefore=ceil(timeBeforeReach/timestep_hs);
framesAfter=ceil(timeAfterReach/timestep_hs);

% for each reach, see if DLC made paw tracking output
% if yes, read in paw position during this reach
for i=1:size(reaches,1)
    f=find(reaches(i,:)>0.5);
    if strcmp(nReachesFromEachTrial,'all')
        n=length(f);
    else
        n=nReachesFromEachTrial;
    end
    for j=1:n
        if reachtimes(i,f(j))-avcuetime<noReachesBeforeTime
            continue
        end
        if reachtimes(i,f(j))-avcuetime>noReachesAfterTime
            continue
        end
        % for each reach, each trial
        [vid,frame]=getVidAndFrame_fromTimeAfterCue(reachtimes(i,f(j))-avcuetime,highspeed_tbt,i);
    end
end

end

function getPawTrajectory(DLCoutput_location,vid,frame,framesBefore,framesAfter)



end

function [vid,frame]=getVidAndFrame_fromTimeAfterCue(timeAfterCue,highspeed_tbt,whichTrial)

cuethistrial=highspeed_tbt.cue(whichTrial,:);
vidthistrial=highspeed_tbt.whichVid(whichTrial,:);
framethistrial=highspeed_tbt.whichFrame(whichTrial,:);
timethistrial=highspeed_tbt.times(whichTrial,:);

cuef=find(cuethistrial>0.5,1,'first');
cuetime=timethistrial(cuef);
times=timethistrial-cuetime;

[~,mi]=nanmin(abs(timeAfterCue-times));
vid=vidthistrial(mi);
frame=framethistrial(mi);

end