function tbt=postAnalysis_checkForChewedPellet(tbt,finaldata)
% function finaldata=checkForChewedPellet(finaldata)

settings=autoReachAnalysisSettings();
disp(['priorToReach_chewWindow is ' num2str(settings.priorToReach_chewWindow)]);
disp(['minTimeToChew_afterReach is ' num2str(settings.minTimeToChew_afterReach)]);
disp(['withinXSeconds is ' num2str(settings.withinXSeconds)]);
disp(['minTimeToChewPellet is ' num2str(settings.minTimeToChewPellet)]);

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
finaldata.flippedThese_pawOnWheel=finaldata.success_reachStarts_pawOnWheel~=finaldata.success_reachStarts_backup_pawOnWheel;
finaldata.success_reachStarts=finaldata.success_reachStarts_backup;
finaldata.drop_reachStarts=finaldata.drop_reachStarts_backup;
finaldata.success_reachStarts_pawOnWheel=finaldata.success_reachStarts_pawOnWheel_backup;
finaldata.drop_reachStarts_pawOnWheel=finaldata.drop_reachStarts_pawOnWheel_backup;

% First for reaches where paw does not start on wheel
finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
[finaldata.success_reachStarts,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese);
finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
finaldata.drop_reachStarts(newDrops==1)=1;

% For reaches where paw does start on wheel
finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
[finaldata.success_reachStarts_pawOnWheel,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts_pawOnWheel,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,finaldata.flippedThese_pawOnWheel);
finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
finaldata.drop_reachStarts_pawOnWheel(newDrops==1)=1;

end

function adjustTbtAccordingly(currReachInd,finaldata,tbt)

movieframe=finaldata.movieframeinds(currReachInd);
temp=tbt.movieframeinds;
[a,b]=find(temp==movieframe);
if isempty(a) || isempty(b)
    error('Could not find matching movie frame in tbt IN postAnalysis_checkForChewedPellet.m');
end

end

function [reaches,newDrops]=checkForSufficientChewing(reaches,chewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent,finaldata,tbt,flipped)

fi=find(reaches==1);
newDrops=zeros(size(reaches));
flippedBack=0;
newFlip=0;
stillFlip=0;
currFlip=false;
for i=1:length(fi)
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
        adjustTbtAccordingly(currReachInd,finaldata,tbt);
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
                adjustTbtAccordingly(currReachInd,finaldata,tbt);
                currFlip=true;
            end
        end
    end
    % did analysis code originally flip this from success to drop
    if flipped(currReachInd)==1 && currFlip==true
        stillFlip=stillFlip+1;
    elseif flipped(currReachInd)==1 && currFlip==false
        flippedBack=flippedBack+1;
    elseif flipped(currReachInd)==0 && currFlip==true
        newFlip=newFlip+1;
    end
end

disp('Flipped this many back ');
disp(flippedBack);
disp('Still flip this many ');
disp(stillFlip);
disp('This many new flips ');
disp(newFlip);

end
