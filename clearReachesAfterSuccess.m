function alltbt=clearReachesAfterSuccess(alltbt)

% get minimum time for mouse to chew and consume full pellet
settings=autoReachAnalysisSettings();
chewTime=settings.chew.minTimeToChewPellet;
timeBin=mode(nanmean(alltbt.times(:,2:end),1)-nanmean(alltbt.times(:,1:end-1),1));
chewInds=floor(chewTime/timeBin);

for i=1:size(alltbt.all_reachBatch,1)
    temp_allReaches=alltbt.all_reachBatch(i,:);
    temp_chewing=alltbt.isChewing(i,:);
    temp_success=alltbt.reachBatch_success_reachStarts(i,:);
    f=find(temp_success==1);
    temp_timesAfterSuccess=zeros(size(temp_success));
    
    % check for reaches within X secs of successful reach, while mouse chewing
    for j=1:length(f)
        temp_timesAfterSuccess(f(j):f(j)+chewInds-1)=1;
    end        
    temp_timesAfterSuccess=temp_timesAfterSuccess(1:length(temp_success));
    chewReaches=temp_allReaches>0.05 & temp_chewing>0.05 & temp_timesAfterSuccess>0.05;
    temp_nonChewReaches=temp_allReaches;
    temp_nonChewReaches(chewReaches==1)=0;
    alltbt.all_reachBatch(i,:)=temp_nonChewReaches;
end
    
    