function behavior_tbt=getChewingEnds(behavior_tbt)

behavior_tbt.chewingEnds=nan(size(behavior_tbt.cue));
temp=behavior_tbt.isChewing;
for i=1:size(temp,1)
    currchew=temp(i,:);
    behavior_tbt.chewingEnds(i,2:end)=diff(currchew)<-0.5;
end