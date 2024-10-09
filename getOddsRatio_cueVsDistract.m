function [alltbt,metadata,trialTypes,distract_tbt,trialTypes_distract,metadata_distract]=getOddsRatio_cueVsDistract(alltbt,metadata,trialTypes,distract_tbt,trialTypes_distract,metadata_distract)

% for each training day (unique sessid), fit logistic regression to see
% effect of cue vs. distract on reaching
u=unique(metadata.sessid);
metadata.odds_ratio=nan(size(metadata.hit_rates));
metadata_distract.odds_ratio=nan(size(metadata_distract.hit_rates));
% get hit rates after cue, hit rates after distractor
for i=1:length(u)
    curr_sessid=u(i);
    cue_hit_rate=mode(metadata.hit_rates(metadata.sessid==curr_sessid));
    distract_hit_rate=mode(metadata_distract.hit_rates(metadata_distract.sessid==curr_sessid));
    % get # trials
    ntrialscue=nansum(metadata.sessid==curr_sessid);
    ntrialsdistract=nansum(metadata_distract.sessid==curr_sessid);
    % reconstruct data from summary stats: hit rates and n trials
    hit_occurred=[ones(round(cue_hit_rate.*ntrialscue),1); zeros(round((1-cue_hit_rate).*ntrialscue),1); ...
                  ones(round(distract_hit_rate.*ntrialsdistract),1); zeros(round((1-distract_hit_rate).*ntrialsdistract),1)];
    cue_type=[zeros(size([ones(round(cue_hit_rate.*ntrialscue),1); zeros(round((1-cue_hit_rate).*ntrialscue),1)])); ...
              ones(size([ones(round(distract_hit_rate.*ntrialsdistract),1); zeros(round((1-distract_hit_rate).*ntrialsdistract),1)]))];
    [b,dev,stats]=glmfit(cue_type,hit_occurred,'binomial','link','logit');
    odds_ratio=exp(-b(2)); % cue relative to distractor
    metadata.odds_ratio(metadata.sessid==curr_sessid)=odds_ratio;
    metadata_distract.odds_ratio(metadata_distract.sessid==curr_sessid)=odds_ratio;
end
alltbt.odds_ratio=metadata.odds_ratio;
trialTypes.odds_ratio=metadata.odds_ratio;
distract_tbt.odds_ratio=metadata_distract.odds_ratio;
trialTypes_distract.odds_ratio=metadata_distract.odds_ratio;

end