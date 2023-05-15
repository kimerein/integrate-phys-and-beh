function alignment=putReachBatchesBackIntoContinuousAlignment(tbt,alignment)

batchfields={'reachBatch_miss_reachStarts','reachBatch_success_reachStarts_pawOnWheel','reachBatch_drop_reachStarts_pawOnWheel','reachBatch_miss_reachStarts_pawOnWheel',...
             'reachBatch_success_reachStarts','reachBatch_drop_reachStarts','all_reachBatch'};
frames=tbt.movieframeinds;
for i=1:length(batchfields)
    curr=tbt.(batchfields{i});
    for j=1:size(curr,1)
        f=find(curr(j,:)>0.5);
        for k=1:length(f)
            alignment=putBack(alignment,frames(f(k)),batchfields{i});
        end
    end
end

end

function alignment=putBack(alignment,frameNumber,whichField)

if ~isfield(alignment,whichField)
    alignment.(whichField)=zeros(size(alignment.cue));
end
temp=alignment.(whichField);
[~,mi]=nanmin(abs(alignment.movieframeinds-frameNumber));
temp(mi)=1;
alignment.(whichField)=temp;

end