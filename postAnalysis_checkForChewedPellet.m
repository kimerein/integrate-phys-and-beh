function [tbt,finaldata]=postAnalysis_checkForChewedPellet(tbt,finaldata,savehandles,zoneVals,eat)

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
[out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,'success_reachStarts','drop_reachStarts');
chewingPower=out1.chewingPower;
chewingDuration=out1.chewingDuration;
s1=eval(out1.threshold);
chewingPower=out2.chewingPower;
chewingDuration=out2.chewingDuration;
s2=eval(out2.threshold);
figure();
polyg_x=[10 0 0 nanmax(out1.chewingPower)+1 nanmax(out1.chewingPower)+1 nanmax(out1.chewingPower)+1];
polyg_y=[110 160 nanmax(out2.chewingDuration)+1 nanmax(out2.chewingDuration)+1 nanmax(out2.chewingDuration)+1 110];
in=inpolygon(out1.chewingPower,out2.chewingDuration,polyg_x,polyg_y);
scatter(out1.chewingPower(in),out2.chewingDuration(in),'r+'); hold on;
plot(polyg_x,polyg_y);
scatter(out1.chewingPower(~in),out2.chewingDuration(~in),'bo');
theseAreSuccess=s1 & s2 & in;
disp([out1.movieFrameInds' out2.movieFrameInds'])
tbt=adjustTbtUsingThresh(out1.movieFrameInds,tbt,theseAreSuccess,false,finaldata);

[out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,'success_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel');
chewingPower=out1.chewingPower;
chewingDuration=out1.chewingDuration;
s1=eval(out1.threshold);
chewingPower=out2.chewingPower;
chewingDuration=out2.chewingDuration;
s2=eval(out2.threshold);
figure();
polyg_x=[10 0 0 nanmax(out1.chewingPower)+1 nanmax(out1.chewingPower)+1 nanmax(out1.chewingPower)+1];
polyg_y=[110 160 nanmax(out2.chewingDuration)+1 nanmax(out2.chewingDuration)+1 nanmax(out2.chewingDuration)+1 110];
in=inpolygon(out1.chewingPower,out2.chewingDuration,polyg_x,polyg_y);
scatter(out1.chewingPower(in),out2.chewingDuration(in),'r+'); hold on;
plot(polyg_x,polyg_y);
scatter(out1.chewingPower(~in),out2.chewingDuration(~in),'bo');
theseAreSuccess=s1 & s2 & in;
disp([out1.movieFrameInds' out2.movieFrameInds'])
tbt=adjustTbtUsingThresh(out1.movieFrameInds,tbt,theseAreSuccess,true,finaldata);

end

function tbt=adjustTbtUsingThresh(movieframes,tbt,theseAreSuccess,isPawOnWheel,finaldata)

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
        else
            if tbt.drop_reachStarts(a,b)==0
                disp('switching to drop based on power and duration 2D');
            end
            tbt.drop_reachStarts(a,b)=1;
            tbt.success_reachStarts(a,b)=0;
        end
    end  
end

end

function tbt=adjustTbtAccordingly(currReachInd,finaldata,tbt,flipped,currFlip,isPawOnWheel)

movieframe=finaldata.movieframeinds(currReachInd);
temp=tbt.movieframeinds;
[~,mi]=min(abs(temp-movieframe),[],'all','omitnan','linear');
[a,b]=ind2sub(size(temp),mi);
if isempty(a) || isempty(b)
    error('Could not find matching movie frame in tbt IN postAnalysis_checkForChewedPellet.m');
end

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
