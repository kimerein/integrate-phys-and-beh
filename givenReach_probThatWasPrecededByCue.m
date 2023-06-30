function [was_reach_preceded_by_cue,trial_number,reach_time]=givenReach_probThatWasPrecededByCue(dataset, whichreachfield, cueWithinXsec, timestep, cueAtInd)

% for each reach, was it was preceded by a cue?
% also return trial number
% and when reach occurred
was_reach_preceded_by_cue=[];
trial_number=[];
reach_time=[];
reaches=dataset.(whichreachfield);
if iscell(reaches)
    reaches=reaches{1};
end
if cueWithinXsec<0
    cueWithinXind=floor(abs(cueWithinXsec)/timestep);
    cueAfterReach=true;
else
    cueWithinXind=floor(cueWithinXsec/timestep);
    cueAfterReach=false;
end
timelabelsfromcue=0:timestep:(size(reaches,2)-1)*timestep;
timelabelsfromcue=timelabelsfromcue-cueAtInd*timestep;
for i=1:size(reaches,1)
    temp=reaches(i,:);
    f=find(temp>0.5);
    if ~isempty(f)
        if cueAfterReach==true
            % of all the reaches on this trial, what fraction are within
            % cueWithinXind of cue and precede cue?
            % i.e., P(cue follows reach | reach)
            was_reach_preceded_by_cue=[was_reach_preceded_by_cue (cueAtInd-f<cueWithinXind & f-cueAtInd<0)];
        else
            % of all the reaches on this trial, what fraction are within
            % cueWithinXind of cue and follow cue?
            % i.e., P(cue precedes reach | reach)
            was_reach_preceded_by_cue=[was_reach_preceded_by_cue (f-cueAtInd<cueWithinXind & f-cueAtInd>0)];
        end
        trial_number=[trial_number ones(1,length(f))*i];
        reach_time=[reach_time timelabelsfromcue(f)];
    end
end

end