function alltbt=fixReachesStuckAtOne(alltbt)

fixTheseFields={'all_reachBatch','pelletmissingreach_reachStarts','reachBatch_all_pawOnWheel',...
                'reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts'...
                'reachStarts','reachStarts_pelletPresent'};

% mouse cannot possibly reach more than 50 times in 1 trial
% this would be a reach rate of over 5 Hz, which is impossible
for i=1:length(fixTheseFields)
    currfield=fixTheseFields{i};
    temp=alltbt.(currfield);
    temp(sum(temp,2,'omitnan')>50,:)=zeros(size(temp(sum(temp,2,'omitnan')>50,:))); % just discard these trials, they are very infrequent
    alltbt.(currfield)=temp;
end

end