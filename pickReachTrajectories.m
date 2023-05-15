function [allX,allY,allZ,allX_from_under,reachTrajTimes]=pickReachTrajectories(lowspeed_tbt,highspeed_tbt,whichReachField,fromWhichArg,DLCoutput_location,vidName,fps,onlyLED)

% fps is frames per second of high speed movie

timeBeforeReach=2;
timeAfterReach=2.5;
nReachesFromEachTrial='all';
spaceOutReaches=true;
spaceOutByTime=1; % in seconds
noReachesBeforeTime=-0.3; % wrt cue, to be sure have enough frames around reach
noReachesAfterTime=9; % wrt cue, to be sure have enough frames around reach

switch onlyLED
    case 'noLED'
        % zero out the reach field for any trial where LED was on
        if fromWhichArg==1
            reaches=lowspeed_tbt.(whichReachField);
            reaches(any(lowspeed_tbt.optoOn(:,10:200)>0.5,2),:)=0;
            lowspeed_tbt.(whichReachField)=reaches;
        elseif fromWhichArg==2
            error('Not yet implemented for pick reach from high speed');
        end
        noReachesBeforeTime=0; % wrt cue, to be sure have enough frames around reach
        noReachesAfterTime=1; % wrt cue, to be sure have enough frames around reach
    case 'LED'
        % zero out the reach field for any trial where LED was NOT on
        if fromWhichArg==1
            reaches=lowspeed_tbt.(whichReachField);
            reaches(~any(lowspeed_tbt.optoOn(:,10:200)>0.5,2),:)=0;
            lowspeed_tbt.(whichReachField)=reaches;
        elseif fromWhichArg==2
            error('Not yet implemented for pick reach from high speed');
        end
        noReachesBeforeTime=0; % wrt cue, to be sure have enough frames around reach
        noReachesAfterTime=1; % wrt cue, to be sure have enough frames around reach
    case'all'
        % do nothing
    otherwise 
        error('Do not recognize onlyLED parameter value');
end

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
maxNReaches=10*size(reaches,1);
allX=nan(maxNReaches,framesBefore+framesAfter+1);
allY=nan(maxNReaches,framesBefore+framesAfter+1);
allZ=nan(maxNReaches,framesBefore+framesAfter+1);
allX_from_under=nan(maxNReaches,framesBefore+framesAfter+1);
reachTrajTimes=0:timestep_hs:(size(allX,2)-1)*timestep_hs;
reachcounter=1;
withinNInds_forSpaceOut=floor(spaceOutByTime/mode(diff(nanmean(reachtimes,1))));
for i=1:size(reaches,1)
    f=find(reaches(i,:)>0.5);
    % Discard reaches outside of range 
    toDiscard=reachtimes(i,f)-avcuetime<noReachesBeforeTime;
    toDiscard2=reachtimes(i,f)-avcuetime>noReachesAfterTime;
    f=f(~toDiscard & ~toDiscard2);
    if spaceOutReaches
        donottake=zeros(size(f));
        for j=1:length(f)-1
            donttake=f(j+1:end)-f(j)<withinNInds_forSpaceOut;
            donottake(j+1:end)=donottake(j+1:end)+donttake;
        end
        f=f(donottake<0.5);
    end
    if strcmp(nReachesFromEachTrial,'all')
        n=length(f);
    else
        n=nReachesFromEachTrial;
    end
    if isempty(f)
        continue
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
        [X,Y,Z,X_from_under]=getPawTrajectory(DLCoutput_location,vidName,vid,frame,framesBefore,framesAfter);
        if ~isempty(X)
            allX(reachcounter,:)=X;
            allY(reachcounter,:)=Y;
            allZ(reachcounter,:)=Z;
            allX_from_under(reachcounter,:)=X_from_under;
            reachcounter=reachcounter+1;
        end
    end
end

allX=allX(1:reachcounter-1,:);
allY=allY(1:reachcounter-1,:);
allZ=allZ(1:reachcounter-1,:);
allX_from_under=allX_from_under(1:reachcounter-1,:);

end

function [X,Y,Z,X_from_under]=getPawTrajectory(DLCoutput_location,vidName,vid,frame,framesBefore,framesAfter)

X=[]; Y=[]; Z=[]; X_from_under=[];

r=regexp(vidName,'0000');
firstHalf=vidName(1:r-1);
secondHalf=vidName(r+4:end);
vid_s=num2str(vid);
if length(vid_s)==1
    vid_s=['000' vid_s];
elseif length(vid_s)==2
    vid_s=['00' vid_s];
elseif length(vid_s)==3
    vid_s=['0' vid_s];
elseif length(vid_s)==4
    vid_s=[vid_s];
elseif length(vid_s)>4
    error('Need to fix code to accomodate more than 1000 highspeed videos');
end
currvidname=[firstHalf vid_s secondHalf];
if exist([DLCoutput_location '\' currvidname],'file')
    a=load([DLCoutput_location '\' currvidname]);
    disp(['Adding reach from ' DLCoutput_location '\' currvidname]);
    % expect X, Y, Z, X_from_under
    indstotake=frame-framesBefore:frame+framesAfter;
    padAtFront=0;
    padAtBack=0;
    if indstotake(1)<1
        padAtFront=nansum(indstotake<1);
        indstotake=1:indstotake(end);
    end
    if indstotake(end)>length(a.X)
        padAtBack=nansum(indstotake>length(a.X));
        indstotake=indstotake(1):length(a.X);
    end
    X=[nan(1,padAtFront) a.X(indstotake)' nan(1,padAtBack)];
    Y=[nan(1,padAtFront) a.Y(indstotake)' nan(1,padAtBack)];
    Z=[nan(1,padAtFront) a.Z(indstotake)' nan(1,padAtBack)];
    X_from_under=[nan(1,padAtFront) a.X_from_under(indstotake)' nan(1,padAtBack)];
else
    return
end

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