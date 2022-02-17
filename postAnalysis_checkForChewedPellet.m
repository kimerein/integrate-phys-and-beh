function [tbt,finaldata]=postAnalysis_checkForChewedPellet(tbt,finaldata,savehandles,zoneVals,eat)

overweightFP=6.5; 
% overweightFP=1; 
% will overweight false positives if not equal to 1, FP are real drops called a success
% false positives will count overweightFP times more than true positives
% useThreshFromNoPawOnWheel=true; % if true, apply threshold from no paw on wheel reaches to paw on wheel reaches

settings=autoReachAnalysisSettings();
disp(['priorToReach_chewWindow is ' num2str(settings.chew.priorToReach_chewWindow)]);
disp(['minTimeToChew_afterReach is ' num2str(settings.chew.minTimeToChew_afterReach)]);
disp(['withinXSeconds is ' num2str(settings.chew.withinXSeconds)]);
disp(['minTimeToChewPellet is ' num2str(settings.chew.minTimeToChewPellet)]);

minTimePelletChew=settings.chew.minTimeToChewPellet;
withinXSeconds=settings.chew.withinXSeconds;
priorSeconds=settings.chew.priorToReach_chewWindow;
dropIfChewingBefore=settings.chew.dropIfChewingBefore;
minTimeMoreStringent=settings.chew.minTimeToChew_afterReach;
fps=settings.movie_fps;
% Convert to inds
minIndToPelletChew=floor(minTimePelletChew/(1/fps));
withinXInds=floor(withinXSeconds/(1/fps));
priorXInds=floor(priorSeconds/(1/fps));
minIndMoreStringent=floor(minTimeMoreStringent/(1/fps));

% Set up fields for storing some ambiguity about drop vs. success
tbt.maybeDrop_reachStarts_pawOnWheel=zeros(size(tbt.cue));
tbt.maybeDrop_reachStarts=zeros(size(tbt.cue));

% Fix any reach values below 1
tbt=setAllAbove0to1(tbt,'reach');

% Check that for each reach classified as a success, there is at least this
% much chewing time after the cue, 

% Take backups
finaldata.flippedThese=finaldata.success_reachStarts~=finaldata.success_reachStarts_backup;
finaldata.flippedThese_pawOnWheel=finaldata.success_reachStarts_pawOnWheel~=finaldata.success_reachStarts_pawOnWheel_backup;
finaldata.success_reachStarts=finaldata.success_reachStarts_backup;
finaldata.drop_reachStarts=finaldata.drop_reachStarts_backup;
finaldata.success_reachStarts_pawOnWheel=finaldata.success_reachStarts_pawOnWheel_backup;
finaldata.drop_reachStarts_pawOnWheel=finaldata.drop_reachStarts_pawOnWheel_backup;

% First for reaches where paw does not start on wheel
finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
[finaldata.success_reachStarts,newDrops,tbt]=checkForSufficientChewing(finaldata.success_reachStarts,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese,false);
finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
finaldata.drop_reachStarts(newDrops==1)=1;

% For reaches where paw does start on wheel
finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
[finaldata.success_reachStarts_pawOnWheel,newDrops,tbt]=checkForSufficientChewing(finaldata.success_reachStarts_pawOnWheel,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese_pawOnWheel,true);
finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
finaldata.drop_reachStarts_pawOnWheel(newDrops==1)=1;

% Chewing power and duration thresholds
disp('Plotting which duration and power thresholds distinguish drop vs success');
if isfield(tbt,'all_reachBatch')
    tbt.currentstudysuccess=(tbt.('success_reachStarts') + tbt.('reachBatch_success_reachStarts') + tbt.('success_reachStarts_pawOnWheel') + tbt.('reachBatch_success_reachStarts_pawOnWheel')) > 0.5;
    tbt.currentstudydrop=(tbt.('drop_reachStarts') + tbt.('reachBatch_drop_reachStarts') + tbt.('drop_reachStarts_pawOnWheel') + tbt.('reachBatch_drop_reachStarts_pawOnWheel')) > 0.5;
else
    tbt.currentstudysuccess=(tbt.('success_reachStarts') + tbt.('success_reachStarts_pawOnWheel')) > 0.5;
    tbt.currentstudydrop=(tbt.('drop_reachStarts') + tbt.('drop_reachStarts_pawOnWheel')) > 0.5;
end
whichIsReachStarts_noPawOnWheel=tbt.('success_reachStarts')>0.5;
whichIsReachStarts_PawOnWheel=tbt.('success_reachStarts_pawOnWheel')>0.5;
[out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,'currentstudysuccess','currentstudydrop',overweightFP,whichIsReachStarts_noPawOnWheel,whichIsReachStarts_PawOnWheel);
chewingPower=out1.chewingPower;
chewingDuration=out1.chewingDuration;
rawIntensity=out1.rawIntensity;
s1=eval(out1.threshold);
nopawonwheel_thresh1=out1.threshold;
chewingPower=out2.chewingPower;
chewingDuration=out2.chewingDuration;
s2=eval(out2.threshold);
nopawonwheel_thresh2=out2.threshold;
theseAreSuccess=s1 & s2;
maybeDropMaybeSuccess=s1~=s2;
disp([out1.movieFrameInds' out2.movieFrameInds'])
tbt=adjustTbtUsingThresh(out1.movieFrameInds,tbt,theseAreSuccess,false,finaldata,maybeDropMaybeSuccess);

chewingPower=out1.paw_chewingPower;
chewingDuration=out1.paw_chewingDuration;
rawIntensity=out1.paw_rawIntensity;
s1=eval(out1.threshold);
chewingPower=out2.paw_chewingPower;
chewingDuration=out2.paw_chewingDuration;
s2=eval(out2.threshold);
theseAreSuccess=s1 & s2;
maybeDropMaybeSuccess=s1~=s2;
tbt=adjustTbtUsingThresh(out1.paw_movieFrameInds,tbt,theseAreSuccess,true,finaldata,maybeDropMaybeSuccess);

% if isfield(tbt,'all_reachBatch')
%     tbt.currentstudysuccess=(tbt.('success_reachStarts_pawOnWheel') + tbt.('reachBatch_success_reachStarts_pawOnWheel'))>0.5;
%     tbt.currentstudydrop=(tbt.('drop_reachStarts_pawOnWheel') + tbt.('reachBatch_drop_reachStarts_pawOnWheel'))>0.5;
% else
%     tbt.currentstudysuccess=tbt.('success_reachStarts_pawOnWheel');
%     tbt.currentstudydrop=tbt.('drop_reachStarts_pawOnWheel');
% end

% [out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,'currentstudysuccess','currentstudydrop',overweightFP,whichIsReachStarts_PawOnWheel,whichIsReachStarts_PawOnWheel);
% chewingPower=out1.chewingPower;
% chewingDuration=out1.chewingDuration;
% rawIntensity=out1.rawIntensity;
% if useThreshFromNoPawOnWheel==true
%     s1=eval(nopawonwheel_thresh1);
% else
%     s1=eval(out1.threshold);
% end
% chewingPower=out2.chewingPower;
% chewingDuration=out2.chewingDuration;
% if useThreshFromNoPawOnWheel==true
%     s2=eval(nopawonwheel_thresh2);
% else
%     s2=eval(out2.threshold);
% end
% theseAreSuccess=s1 & s2;
% disp([out1.movieFrameInds' out2.movieFrameInds'])
% tbt=adjustTbtUsingThresh(out1.movieFrameInds,tbt,theseAreSuccess,true,finaldata);

tbt=rmfield(tbt,'currentstudysuccess');
tbt=rmfield(tbt,'currentstudydrop');

end

function tbt=setAllAbove0to1(tbt,ifFieldContains)

f=fieldnames(tbt);
for i=1:length(f)
    currfield=f{i};
    if ~isempty(regexp(currfield,ifFieldContains,'once'))
        tempie=zeros(size(tbt.(currfield)));
        tempie(tbt.(currfield)>0)=1;
        tbt.(currfield)=tempie;
    end
end

end

function tbt=adjustTbtUsingThresh(movieframes,tbt,theseAreSuccess,isPawOnWheel,finaldata,maybeDropMaybeSuccess)

for i=1:length(movieframes)
    movieframe=movieframes(i);
    temp=tbt.movieframeinds;
    [~,mi]=min(abs(temp-movieframe),[],'all','omitnan','linear');
    [a,b]=ind2sub(size(temp),mi);
    if isempty(a) || isempty(b)
        error('Could not find matching movie frame in tbt IN postAnalysis_checkForChewedPellet.m, adjustTbtUsingThresh');
    end
    if theseAreSuccess(i)==true
        % is success
        % if finaldata also thinks that this is a success, make sure that
        % tbt says success
        [~,movind]=nanmin(abs(finaldata.movieframeinds-movieframe));
        if isPawOnWheel==true
            if finaldata.success_reachStarts_pawOnWheel(movind)>0.5
                tbt.success_reachStarts_pawOnWheel(a,b)=1;
                tbt.drop_reachStarts_pawOnWheel(a,b)=0;
            end
        else
            if finaldata.success_reachStarts(movind)>0.5
                tbt.success_reachStarts(a,b)=1;
                tbt.drop_reachStarts(a,b)=0;
            end
        end
    else
        % is drop
        if isPawOnWheel==true
            if tbt.drop_reachStarts_pawOnWheel(a,b)==0
                disp('switching to drop based on power and duration 2D');
            end
            tbt.drop_reachStarts_pawOnWheel(a,b)=1;
            tbt.success_reachStarts_pawOnWheel(a,b)=0;
            if maybeDropMaybeSuccess(i)==1
                tbt.maybeDrop_reachStarts_pawOnWheel(a,b)=1;
            end
        else
            if tbt.drop_reachStarts(a,b)==0
                disp('switching to drop based on power and duration 2D');
            end
            tbt.drop_reachStarts(a,b)=1;
            tbt.success_reachStarts(a,b)=0;
            if maybeDropMaybeSuccess(i)==1
                tbt.maybeDrop_reachStarts(a,b)=1;
            end
        end
    end  
end

end

function [allmi,ms]=findAllMovieInds(movieframe,allmovieframes)

allmi=[];
ms=[];
for i=1:1000
    [m,mi]=min(abs(allmovieframes-movieframe),[],'all','omitnan','linear');
    [a,b]=ind2sub(size(allmovieframes),mi);
    allmi=[allmi; [a,b]];
    ms=[ms; m];
    allmovieframes(a,b)=movieframe+1000;
    if m>0.5 && ~isempty(allmi)
        break
    end
end
allmi=allmi(ms<0.5,:);
ms=ms(ms<0.5);

end

function tbt=adjustTbtAccordingly(currReachInd,finaldata,tbt,flipped,currFlip,isPawOnWheel)

movieframe=finaldata.movieframeinds(currReachInd);
temp=tbt.movieframeinds;
[allmi,ms]=findAllMovieInds(movieframe,temp); % in case there are multiple instances of this time point at end and beginning of neighboring trials in tbt
% [~,mi]=min(abs(temp-movieframe),[],'all','omitnan','linear');
% [a,b]=ind2sub(size(temp),mi);
% if isempty(a) || isempty(b)
%     error('Could not find matching movie frame in tbt IN postAnalysis_checkForChewedPellet.m');
% end

for i=1:size(allmi,1)
    a=allmi(i,1);
    b=allmi(i,2);
    if currFlip==true && flipped(currReachInd)==true
        % should flip tbt from success to drop
        % but tbt should already be flipped
    elseif currFlip==true && flipped(currReachInd)==false
        % flip now although didn't flip before
        if isPawOnWheel==true
            tbt.drop_reachStarts_pawOnWheel(a,b)=1;
            tbt.success_reachStarts_pawOnWheel(a,b)=0;
        else
            tbt.drop_reachStarts(a,b)=1;
            tbt.success_reachStarts(a,b)=0;
        end
    elseif currFlip==false && flipped(currReachInd)==false
        % don't flip
    elseif currFlip==false && flipped(currReachInd)==true
        % previously flipped success to drop, but now flip back
        if isPawOnWheel==true
            tbt.drop_reachStarts_pawOnWheel(a,b)=0;
            tbt.success_reachStarts_pawOnWheel(a,b)=1;
        else
            tbt.drop_reachStarts(a,b)=0;
            tbt.success_reachStarts(a,b)=1;
        end
    end
end

end

function [reaches,newDrops,tbt]=checkForSufficientChewing(reaches,chewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,flipped,isPawOnWheel)

fi=find(reaches==1);
newDrops=zeros(size(reaches));
flippedBack=0;
newFlips=0;
for i=1:length(fi)
    currFlip=false;
    currReachInd=fi(i);
    % is there enough chewing within X seconds of this reach
    if currReachInd+withinXInds>length(chewing)
        chewInds=sum(chewing(currReachInd:end)>0.5);
    else
        chewInds=sum(chewing(currReachInd:currReachInd+withinXInds)>0.5);
    end
    if chewInds<minIndToPelletChew % not enough chewing to be consistent with eating pellet
        reaches(currReachInd)=0; % not a successful reach
        newDrops(currReachInd)=1; % actually a drop
        currFlip=true;
    end
    if dropIfChewingBefore==1 && newDrops(currReachInd)==0
        % was mouse chewing BEFORE reach?
        if currReachInd-priorXInds<1
            chewInds_before=sum(chewing(1:currReachInd-1)>0.5);
        else
            chewInds_before=sum(chewing(currReachInd-priorXInds:currReachInd-1)>0.5);
        end
        if chewInds_before>floor((minIndToPelletChew/withinXInds)*priorXInds)
            % mouse was chewing before reach
            % did mouse chew long enough after reach, consistent with
            % consumption of full pellet?
            if chewInds<minIndMoreStringent % not enough chewing to be consistent with eating pellet
                reaches(currReachInd)=0; % not a successful reach
                newDrops(currReachInd)=1; % actually a drop
                currFlip=true;
            end
        end
    end
    tbt=adjustTbtAccordingly(currReachInd,finaldata,tbt,flipped,currFlip,isPawOnWheel);
    if currFlip==false && flipped(currReachInd)==true
        flippedBack=flippedBack+1;
    elseif currFlip==true && flipped(currReachInd)==false
        newFlips=newFlips+1;
    end
end
disp(['flipped back ' num2str(flippedBack)]);
disp(['new flips ' num2str(newFlips)]);

end
