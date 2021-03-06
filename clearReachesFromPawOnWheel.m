function alltbt=clearReachesFromPawOnWheel(alltbt)

alltbt.reachBatch_all_pawOnWheel=alltbt.reachBatch_success_reachStarts_pawOnWheel + alltbt.reachBatch_drop_reachStarts_pawOnWheel + alltbt.reachBatch_miss_reachStarts_pawOnWheel;

for i=1:size(alltbt.all_reachBatch,1)
    temp_allReaches=alltbt.all_reachBatch(i,:);
    temp_pawOnWheel=alltbt.pawOnWheel(i,:);
    
    % shift reaches while paw on wheel
    temp_nonPawReaches=temp_allReaches;
    temp_pawReaches=temp_allReaches;
    temp_pawReaches(temp_pawOnWheel<0.05)=0;
    temp_nonPawReaches(temp_pawOnWheel>0.05)=0;
    alltbt.all_reachBatch(i,:)=temp_nonPawReaches;
    alltbt.reachBatch_all_pawOnWheel(i,:)=alltbt.reachBatch_all_pawOnWheel(i,:) + temp_pawReaches;
    
    % shift successes while paw on wheel
    temp_nonPawReaches=alltbt.reachBatch_success_reachStarts(i,:);
    temp_pawReaches=alltbt.reachBatch_success_reachStarts(i,:);
    temp_pawReaches(temp_pawOnWheel<0.05)=0;
    temp_nonPawReaches(temp_pawOnWheel>0.05)=0;
    alltbt.reachBatch_success_reachStarts(i,:)=temp_nonPawReaches;
    alltbt.reachBatch_success_reachStarts_pawOnWheel(i,:)=alltbt.reachBatch_success_reachStarts_pawOnWheel(i,:) + temp_pawReaches;
    
    % shift drops while paw on wheel
    temp_nonPawReaches=alltbt.reachBatch_drop_reachStarts(i,:);
    temp_pawReaches=alltbt.reachBatch_drop_reachStarts(i,:);
    temp_pawReaches(temp_pawOnWheel<0.05)=0;
    temp_nonPawReaches(temp_pawOnWheel>0.05)=0;
    alltbt.reachBatch_drop_reachStarts(i,:)=temp_nonPawReaches;
    alltbt.reachBatch_drop_reachStarts_pawOnWheel(i,:)=alltbt.reachBatch_drop_reachStarts_pawOnWheel(i,:) + temp_pawReaches;
    
    % shift misses while paw on wheel
    temp_nonPawReaches=alltbt.reachBatch_miss_reachStarts(i,:);
    temp_pawReaches=alltbt.reachBatch_miss_reachStarts(i,:);
    temp_pawReaches(temp_pawOnWheel<0.05)=0;
    temp_nonPawReaches(temp_pawOnWheel>0.05)=0;
    alltbt.reachBatch_miss_reachStarts(i,:)=temp_nonPawReaches;
    alltbt.reachBatch_miss_reachStarts_pawOnWheel(i,:)=alltbt.reachBatch_miss_reachStarts_pawOnWheel(i,:) + temp_pawReaches;    
end

alltbt.reachBatch_all_pawOnWheel(alltbt.reachBatch_all_pawOnWheel>0.5)=1;
alltbt.reachBatch_success_reachStarts_pawOnWheel(alltbt.reachBatch_success_reachStarts_pawOnWheel>0.5)=1;
alltbt.reachBatch_drop_reachStarts_pawOnWheel(alltbt.reachBatch_drop_reachStarts_pawOnWheel>0.5)=1;
alltbt.reachBatch_miss_reachStarts_pawOnWheel(alltbt.reachBatch_miss_reachStarts_pawOnWheel>0.5)=1;