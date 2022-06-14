function [tbt,finaldata]=postAnalysis_checkForChewedPellet(tbt,finaldata,savehandles,zoneVals,eat,usePreviousChewingThresh)

overweightFP=6.5; 
doOr=false;
% overweightFP=1; 
% will overweight false positives if not equal to 1, FP are real drops called a success
% false positives will count overweightFP times more than true positives
% useThreshFromNoPawOnWheel=true; % if true, apply threshold from no paw on wheel reaches to paw on wheel reaches
removeZscore=true;
userSetsThresh=false;
userThresh=1.8*10^5;
if userSetsThresh==true
    ans=questdlg('User sets thresh ... continue?');
    if ~strcmp(ans,'Yes')
        return
    end
end

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
chewFrequency=settings.chew.chewFrequency;
fps=settings.movie_fps;
% Convert to inds
minIndToPelletChew=floor(minTimePelletChew/(1/fps));
withinXInds=floor(withinXSeconds/(1/fps));
priorXInds=floor(priorSeconds/(1/fps));
minIndMoreStringent=floor(minTimeMoreStringent/(1/fps));

% Plot eat zone for double check of chewing threshold
tempie=zoneVals.eatZone;
eatzone=tempie(savehandles.discardFirstNFrames+1:end);
if removeZscore==true
    params.Fs=settings.movie_fps;
    params.tapers=settings.chew.tapers;
    params.fpass=settings.chew.fpass; % in Hz
    [S,t,f]=mtspecgramc(eatzone(~isnan(eatzone)),[5 0.25],params);
    chewingpower=nanmean(S(:,f>=chewFrequency(1) & f<=chewFrequency(2)),2);
    chewingpower=(chewingpower./nanmax(chewingpower))*nanmax(eat.chewingpower);
    
    frameTimes=0:(1/settings.movie_fps):(length(eatzone(~isnan(eatzone)))-1)*(1/settings.movie_fps);
    new_chewingInFrames=mapToFrames(chewingpower,t,frameTimes);
end
figure(); plot(eatzone,'Color','k'); 
hold on; 
plot(eat.isChewing*nanmax(eatzone),'Color','g');
figure(); plot(eat.chewingInFrames,'Color','k');
hold on;
if removeZscore==true
    plot(new_chewingInFrames,'Color','r');
end
plot(eat.isChewing*nanmax(eat.chewingInFrames),'Color','g');
ylabel('chewing power');
if removeZscore==true
    useRed='Will use red. ';
else
    useRed=[];
end
if usePreviousChewingThresh==false
    newchewthresh=input(['Enter new chewing threshold or empty to keep current thresh. ' useRed]);
else
    newchewthresh=[];
end
if isempty(newchewthresh) || usePreviousChewingThresh==true
    changedChewThresh=false;
else
    changedChewThresh=true;
    disp('Changing chew threshold.');
    if removeZscore==true
        eat.chewingInFrames=new_chewingInFrames;
    end
    eat.isChewing=eat.chewingInFrames>newchewthresh;
    figure(); plot(eat.chewingInFrames,'Color','k');
    hold on;
    plot(eat.isChewing*nanmax(eat.chewingInFrames),'Color','g');
    title('With new chew thresh');
    eat.movieFrames=[1:length(eat.isChewing)]+savehandles.discardFirstNFrames;
    finaldata.isChewing=mapUsingMovieFrames(finaldata.movieframeinds,eat.isChewing,eat.movieFrames);
end

% Set up fields for storing some ambiguity about drop vs. success
tbt.maybeDrop_reachStarts_pawOnWheel=zeros(size(tbt.cue));
tbt.maybeDrop_reachStarts=zeros(size(tbt.cue));

% Fix any reach values below 1
tbt=setAllAbove0to1(tbt,'reach');

% Check that for each reach classified as a success, there is at least this
% much chewing time after the cue, 

% Take backups
if ~isfield(finaldata,'success_reachStarts_backup')
    finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
    finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
    finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
    finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
end
finaldata.flippedThese=finaldata.success_reachStarts~=finaldata.success_reachStarts_backup;
finaldata.flippedThese_pawOnWheel=finaldata.success_reachStarts_pawOnWheel~=finaldata.success_reachStarts_pawOnWheel_backup;

% If want to start by considering every potential reach a success, then
% will trim accordingly
finaldata.success_reachStarts=(finaldata.success_reachStarts + finaldata.drop_reachStarts + finaldata.success_reachStarts_backup + finaldata.drop_reachStarts_backup) > 0.5;
finaldata.success_reachStarts_pawOnWheel=(finaldata.success_reachStarts_pawOnWheel + finaldata.drop_reachStarts_pawOnWheel + finaldata.success_reachStarts_pawOnWheel_backup + finaldata.drop_reachStarts_pawOnWheel_backup) > 0.5;

% % First for reaches where paw does not start on wheel
% finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
[finaldata.success_reachStarts,newDrops,tbt]=checkForSufficientChewing(finaldata.success_reachStarts,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese,false);
% finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
finaldata.drop_reachStarts(newDrops==1)=1;
finaldata.drop_reachStarts(newDrops==0)=0;

% For reaches where paw does start on wheel
% finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
[finaldata.success_reachStarts_pawOnWheel,newDrops,tbt]=checkForSufficientChewing(finaldata.success_reachStarts_pawOnWheel,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese_pawOnWheel,true);
% finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
finaldata.drop_reachStarts_pawOnWheel(newDrops==1)=1;
finaldata.drop_reachStarts_pawOnWheel(newDrops==0)=0;

% Chewing power and duration thresholds
disp('Plotting which duration and power thresholds distinguish drop vs success');
% tbt.currentclassifysuccess=(tbt.('success_reachStarts_backup') + tbt.('success_reachStarts_pawOnWheel_backup')) > 0.5;
% tbt.currentclassifydrop=(tbt.('drop_reachStarts_backup') + tbt.('drop_reachStarts_pawOnWheel_backup')) > 0.5;
tbt.currentclassifysuccess=(tbt.('success_reachStarts') + tbt.('success_reachStarts_pawOnWheel')) > 0.5;
tbt.currentclassifydrop=(tbt.('drop_reachStarts') + tbt.('drop_reachStarts_pawOnWheel')) > 0.5;
% if isfield(tbt,'all_reachBatch')
%     tbt.currentstudysuccess=(tbt.('success_reachStarts') + tbt.('reachBatch_success_reachStarts') + tbt.('success_reachStarts_pawOnWheel') + tbt.('reachBatch_success_reachStarts_pawOnWheel')) > 0.5;
%     tbt.currentstudydrop=(tbt.('drop_reachStarts') + tbt.('reachBatch_drop_reachStarts') + tbt.('drop_reachStarts_pawOnWheel') + tbt.('reachBatch_drop_reachStarts_pawOnWheel')) > 0.5;
% else
    tbt.currentstudysuccess=(tbt.('success_reachStarts') + tbt.('success_reachStarts_pawOnWheel')) > 0.5;
    tbt.currentstudydrop=(tbt.('drop_reachStarts') + tbt.('drop_reachStarts_pawOnWheel')) > 0.5;
% end
whichIsReachStarts_noPawOnWheel=tbt.('success_reachStarts')>0.5;
whichIsReachStarts_PawOnWheel=tbt.('success_reachStarts_pawOnWheel')>0.5;
[out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,'currentstudysuccess','currentstudydrop','currentclassifysuccess','currentclassifydrop',overweightFP,whichIsReachStarts_noPawOnWheel,whichIsReachStarts_PawOnWheel);
chewingPower=out1.chewingPower(out1.isCurrReachStart==1);
chewingDuration=out1.chewingDuration(out1.isCurrReachStart==1);
rawIntensity=out1.rawIntensity(out1.isCurrReachStart==1);
if isempty(out1.predictions)
    tbt=rmfield(tbt,'currentstudysuccess');
    tbt=rmfield(tbt,'currentstudydrop');
    tbt=rmfield(tbt,'currentclassifysuccess');
    tbt=rmfield(tbt,'currentclassifydrop');
    return
end
predictions=out1.predictions(out1.isCurrReachStart==1);
currmovieFrameInds=out1.movieFrameInds(out1.isCurrReachStart==1);
s1=eval(out1.threshold);
nopawonwheel_thresh1=out1.threshold;
chewingPower=out2.chewingPower(out2.isCurrReachStart==1);
chewingDuration=out2.chewingDuration(out2.isCurrReachStart==1);
predictions=out2.predictions(out2.isCurrReachStart==1);
s2=eval(out2.threshold);
nopawonwheel_thresh2=out2.threshold;
if doOr==true
    theseAreSuccess=s1 | s2;
else
    theseAreSuccess=s1 & s2;
end
maybeDropMaybeSuccess=s1~=s2;
if userSetsThresh
    theseAreSuccess=rawIntensity > userThresh;
end
% disp([out1.movieFrameInds(out1.isCurrReachStart==1)' out2.movieFrameInds(out2.isCurrReachStart==1)'])
tbt=adjustTbtUsingThresh(currmovieFrameInds,tbt,theseAreSuccess,false,finaldata,maybeDropMaybeSuccess);

chewingPower=out1.chewingPower(out1.isCurrReachPaw==1);
chewingDuration=out1.chewingDuration(out1.isCurrReachPaw==1);
rawIntensity=out1.rawIntensity(out1.isCurrReachPaw==1);
predictions=out1.predictions(out1.isCurrReachPaw==1);
currmovieFrameInds=out1.movieFrameInds(out1.isCurrReachPaw==1);
s1=eval(out1.threshold);
chewingPower=out2.chewingPower(out2.isCurrReachPaw==1);
chewingDuration=out2.chewingDuration(out2.isCurrReachPaw==1);
predictions=out2.predictions(out2.isCurrReachPaw==1);
s2=eval(out2.threshold);
if doOr==true
    theseAreSuccess=s1 | s2;
else
    theseAreSuccess=s1 & s2;
end
maybeDropMaybeSuccess=s1~=s2;
if userSetsThresh
    theseAreSuccess=rawIntensity > userThresh;
end
tbt=adjustTbtUsingThresh(currmovieFrameInds,tbt,theseAreSuccess,true,finaldata,maybeDropMaybeSuccess);

tbt=rmfield(tbt,'currentstudysuccess');
tbt=rmfield(tbt,'currentstudydrop');
tbt=rmfield(tbt,'currentclassifysuccess');
tbt=rmfield(tbt,'currentclassifydrop');

end

function dataByFrames=mapToFrames(data,times,frameTimes)

dataByFrames=nan(size(frameTimes));

for i=1:length(times)
    [~,mi]=min(abs(times(i)-frameTimes));
    dataByFrames(mi)=data(i);
end

dataByFrames=fillInNans(dataByFrames);

end

function data=fillInNans(data)

inds=find(~isnan(data));
for i=1:length(inds)
    currind=inds(i);
    if i==1
        % fill in before
        data(1:currind-1)=data(currind);
    elseif i==length(inds)
        halfLength=floor((currind-inds(i-1))/2);
        data(inds(i-1)+1:inds(i-1)+1+halfLength)=data(inds(i-1));
        data(inds(i-1)+2+halfLength:currind-1)=data(currind);
        % fill in after
        data(currind+1:end)=data(currind);
    else
        % fill in with recent
        halfLength=floor((currind-inds(i-1))/2);
        data(inds(i-1)+1:inds(i-1)+1+halfLength)=data(inds(i-1));
        data(inds(i-1)+2+halfLength:currind-1)=data(currind);
    end
end
if any(isnan(data))
    error('Failed to replace all nans');
end     

end

function finalData=mapUsingMovieFrames(movieFramesInFinalData,dataToMap,framesInOrigData)

finalData=nan(size(movieFramesInFinalData));
for i=1:length(movieFramesInFinalData)
    [~,mi]=nanmin(abs(framesInOrigData-movieFramesInFinalData(i)));
    finalData(i)=dataToMap(mi);    
end

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
    [allmi,ms]=findAllMovieInds(movieframe,temp);
    for j=1:size(allmi,1)
        a=allmi(j,1);
        b=allmi(j,2);
        if theseAreSuccess(i)==true
            % is success
            % if finaldata also thinks that this is a success, make sure that
            % tbt says success
            [~,movind]=nanmin(abs(finaldata.movieframeinds-movieframe));
            if isPawOnWheel==true
                if finaldata.success_reachStarts_pawOnWheel(movind)>0.5
                    tbt.success_reachStarts_pawOnWheel(a,b)=1;
                    tbt.drop_reachStarts_pawOnWheel(a,b)=0;
                else
                    tbt.success_reachStarts_pawOnWheel(a,b)=1;
                    tbt.drop_reachStarts_pawOnWheel(a,b)=0;
                    tbt.maybeDrop_reachStarts_pawOnWheel(a,b)=1;
                end
                if maybeDropMaybeSuccess(i)==1
                    tbt.maybeDrop_reachStarts_pawOnWheel(a,b)=1;
                end
            else
                if finaldata.success_reachStarts(movind)>0.5
                    tbt.success_reachStarts(a,b)=1;
                    tbt.drop_reachStarts(a,b)=0;
                else
                    tbt.success_reachStarts(a,b)=1;
                    tbt.drop_reachStarts(a,b)=0;
                    tbt.maybeDrop_reachStarts(a,b)=1;
                end
                if maybeDropMaybeSuccess(i)==1
                    tbt.maybeDrop_reachStarts(a,b)=1;
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
allmi=allmi(ms<1,:);
ms=ms(ms<1);

end

function [tbt,new_success]=adjustTbtAccordingly(currReachInd,finaldata,tbt,flipped,currFlip,isPawOnWheel)

movieframe=finaldata.movieframeinds(currReachInd);
temp=tbt.movieframeinds;
[allmi,ms]=findAllMovieInds(movieframe,temp); % in case there are multiple instances of this time point at end and beginning of neighboring trials in tbt
% [~,mi]=min(abs(temp-movieframe),[],'all','omitnan','linear');
% [a,b]=ind2sub(size(temp),mi);
% if isempty(a) || isempty(b)
%     error('Could not find matching movie frame in tbt IN postAnalysis_checkForChewedPellet.m');
% end

new_success=false;
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
        % check whether is drop in tbt and make success
        if isPawOnWheel==true
            if tbt.drop_reachStarts_pawOnWheel(a,b)==1
                tbt.drop_reachStarts_pawOnWheel(a,b)=0;
                tbt.success_reachStarts_pawOnWheel(a,b)=1;
                new_success=true;
            end
        else
            if tbt.drop_reachStarts(a,b)==1
                tbt.drop_reachStarts(a,b)=0;
                tbt.success_reachStarts(a,b)=1;
                new_success=true;
            end
        end
        
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
newSuccesses=0;
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
    [tbt,new_success]=adjustTbtAccordingly(currReachInd,finaldata,tbt,flipped,currFlip,isPawOnWheel);
    if new_success==true
        newSuccesses=newSuccesses+1;
    end
    if currFlip==false && flipped(currReachInd)==true
        flippedBack=flippedBack+1;
    elseif currFlip==true && flipped(currReachInd)==false
        newFlips=newFlips+1;
    end
end
disp(['flipped back ' num2str(flippedBack)]);
disp(['new flips ' num2str(newFlips)]);
disp(['new successes ' num2str(newSuccesses)]);

end
