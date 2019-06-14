function boutDurations=getBoutDuration(tbt,whichBout,startInd,takeTrials)

bouts=tbt.(whichBout);

if ~isempty(takeTrials)
    bouts=bouts(takeTrials==1,:);
end

boutDurations=nan(1,size(bouts,1));
bouts(:,1:startInd-1)=0;

for i=1:size(bouts,1)
    % find first bout
    f=find(bouts(i,:)>0.5,1,'first');
    if isempty(f)
        continue
    end
    % get bout duration
    f2=find(bouts(i,f:end)<0.5,1,'first');
    if ~isempty(f2)
        boutDurations(i)=(f2+f-1)-f;
    else
        boutDurations(i)=size(bouts,2)-f;
    end
end
boutDurations=boutDurations.*mode(diff(nanmean(tbt.times,1)));