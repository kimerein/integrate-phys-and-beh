function recodeMissVGrab(varargin)

if isempty(varargin)
    currentVid='Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20220616\dLight_172\O2 output\VID_20130726_234313_processed_data';
    datestr='20220616';
    mousename='dLight_172';
    f_pr=regexp(currentVid,'_processed_data');
    fslash=regexp(currentVid,'\');
    aviName=currentVid(fslash(end)+1:f_pr-1);
    % placeForO2data=['Z:\MICROSCOPE\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];
    placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    alignment=[];
else
    if length(varargin)==3
        currentVid=varargin{1};
        placeForO2data=varargin{2};
        alignment=varargin{3};
    end 
end

%currentVid='Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210714\O2 output\2012-02-10 23-21-33-C_processed_data';
%datestr='20210714';
%mousename='Mar_3';
timeForDislodged=0.4;

% f_pr=regexp(currentVid,'_processed_data');
% fslash=regexp(currentVid,'\');
% aviName=currentVid(fslash(end)+1:f_pr-1);
% placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
% % placeForO2data=['Z:\MICROSCOPE\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];

load([currentVid '\tbt.mat'])
backuptbt=tbt;
load([currentVid '\final_aligned_data.mat'])
ti=mode(diff(alignment.timesfromarduino));
if ti==0
    temp=diff(alignment.timesfromarduino);
    temp(temp==0)=nan;
    ti=mode(temp);
end
ti=ti/1000;
indsForDislodged=ceil(timeForDislodged/ti);
disp(['will allow ' num2str(timeForDislodged) ' sec for mouse to dislodge pellet / ' num2str(indsForDislodged) ' indices']);
[newtbt,alignment]=postAnalysis_checkForPelletDislodged(tbt,alignment,indsForDislodged);
newtbt=addReachBatchesToSingleTbt(newtbt,'cueZone_onVoff',0.25,0,[]);
tbt=newtbt;

save([currentVid '\backuptbt.mat'],'backuptbt');
save([currentVid '\tbt.mat'],'tbt');
fid=fopen([currentVid '\fixed_miss_v_grab.txt'],'wt');
fclose(fid);

% if have already fixed success v drop, recode all reaches using success v
% drop threshold
fixDropVSuccess(currentVid,placeForO2data,alignment);

end

function [missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,whichfield,indsForDislodged)

temp=finaldata.(whichfield);
newGrabs=zeros(size(finaldata.(whichfield)));
temp(temp>0.5)=1;
misses=find(temp>0.5);
reachEnds=finaldata.reachEnds;
missType=zeros(size(misses));
movieframeinds=nan(size(misses));
for i=1:length(misses)
    % for each miss, check whether pellet dislodged within indsForDislodged
    % find closest subsequent reachEnd
    fend=find(reachEnds(misses(i):end)>0.5,1,'first')+misses(i)-1;
    movieframeinds(i)=finaldata.movieframeinds(misses(i));
    if isempty(fend)
        disp('could not find end of reach');
        fend=length(reachEnds);
    end
    if fend+indsForDislodged>length(finaldata.pelletPresent)
        if all(finaldata.pelletPresent(fend+1:end)>0.5)
            % mouse did not move pellet, so miss
            missType(i)=true;
            newGrabs(misses(i))=0;
        else
            % mouse did move pellet, so grab
            missType(i)=false;
            newGrabs(misses(i))=1;
        end
    else
        if all(finaldata.pelletPresent(fend+1:fend+indsForDislodged)>0.5)
            % mouse did not move pellet, so miss
            missType(i)=true;
            newGrabs(misses(i))=0;
        else
            % mouse did move pellet, so grab
            missType(i)=false;
            newGrabs(misses(i))=1;
        end
    end
end

end

function [tbt,finaldata]=postAnalysis_checkForPelletDislodged(tbt,finaldata,indsForDislodged)

% Fix any reach values below 1
tbt=setAllAbove0to1(tbt,'reach');

% Take backups
if ~isfield(finaldata,'success_reachStarts_backup')
    finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
    finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
    finaldata.miss_reachStarts_backup=finaldata.miss_reachStarts;
    finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
    finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
    finaldata.miss_reachStarts_pawOnWheel_backup=finaldata.miss_reachStarts_pawOnWheel;
end

% If pellet is dislodged within 10 indices of end of reach start, then
% consider this a grab rather than a miss
% get all misses
% First for reaches where paw does not start on wheel
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'miss_reachStarts',indsForDislodged);
finaldata.miss_reachStarts(newGrabs==1)=0;
finaldata.success_reachStarts(newGrabs==1)=1; % converts misses to successes
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,false,finaldata);
% For reaches where paw does start on wheel
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'miss_reachStarts_pawOnWheel',indsForDislodged);
finaldata.miss_reachStarts_pawOnWheel(newGrabs==1)=0;
finaldata.success_reachStarts_pawOnWheel(newGrabs==1)=1;
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,true,finaldata);

tbt=addReachBatchesToSingleTbt(tbt,'cueZone_onVoff',0.25,0,[]);

% double check drops and successes
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'success_reachStarts',indsForDislodged);
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,false,finaldata);
% For reaches where paw does start on wheel
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'success_reachStarts_pawOnWheel',indsForDislodged);
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,true,finaldata);
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'drop_reachStarts',indsForDislodged);
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,false,finaldata);
% For reaches where paw does start on wheel
[missType,movieframeinds,newGrabs]=whetherPelletDislodged(finaldata,'drop_reachStarts_pawOnWheel',indsForDislodged);
tbt=adjustTbtUsingThresh(movieframeinds,tbt,~missType,true,finaldata);

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

function tbt=adjustTbtUsingThresh(movieframes,tbt,theseAreSuccess,isPawOnWheel,finaldata)

for i=1:length(movieframes)
    movieframe=movieframes(i);
    temp=tbt.movieframeinds;
    [allmi,ms]=findAllMovieInds(movieframe,temp);
    for j=1:size(allmi,1)
        a=allmi(j,1);
        b=allmi(j,2);
        if theseAreSuccess(i)==true % if mouse actually grabbed pellet
            % is grab
            % make sure that tbt says success 
            if isPawOnWheel==true
                tbt.success_reachStarts_pawOnWheel(a,b)=1;
%                 tbt.drop_reachStarts_pawOnWheel(a,b)=1;
                tbt.miss_reachStarts_pawOnWheel(a,b)=0;
            else
                tbt.success_reachStarts(a,b)=1;
%                 tbt.drop_reachStarts(a,b)=1;
                tbt.miss_reachStarts(a,b)=0;
            end
        else
            % is actually a miss
            if isPawOnWheel==true
                tbt.success_reachStarts_pawOnWheel(a,b)=0;
                tbt.drop_reachStarts_pawOnWheel(a,b)=0;
                tbt.miss_reachStarts_pawOnWheel(a,b)=1;
            else
                tbt.success_reachStarts(a,b)=0;
                tbt.drop_reachStarts(a,b)=0;
                tbt.miss_reachStarts(a,b)=1;
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
