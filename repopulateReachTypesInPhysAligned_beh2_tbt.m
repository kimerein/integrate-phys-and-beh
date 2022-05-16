function beh2_tbt=repopulateReachTypesInPhysAligned_beh2_tbt(fixedtbt,beh2_tbt)

reachTypeFieldNames={'success_reachStarts','drop_reachStarts','miss_reachStarts','success_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel','miss_reachStarts_pawOnWheel'...
    'reachBatch_miss_reachStarts','reachBatch_success_reachStarts_pawOnWheel','reachBatch_drop_reachStarts_pawOnWheel','reachBatch_miss_reachStarts_pawOnWheel','reachBatch_success_reachStarts',...
    'reachBatch_drop_reachStarts','all_reachBatch'};

[fixedtbt,beh2_tbt]=commonReference_usingDistractor(fixedtbt,beh2_tbt,'movie_distractor');

% reference_into field refers to which row of other behavior tbt

% repopulate all reach types
for j=1:length(reachTypeFieldNames)
    currfield=reachTypeFieldNames{j};
    temp=beh2_tbt.(currfield);
    fixedtemp=fixedtbt.(currfield);
    % for each row of fixedtbt, repopulate rows of beh2_tbt
    for i=1:size(fixedtbt.reference_into_beh2trialinds,1)
        rowInOtherBehTbt=fixedtbt.reference_into_beh2trialinds(i,1);
        if isnan(rowInOtherBehTbt)
            continue
        end
        temp(rowInOtherBehTbt,:)=fixedtemp(i,:);
    end
    beh2_tbt.(currfield)=temp;
end

end