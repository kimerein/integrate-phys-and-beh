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
    cue_FA_rate=mode(metadata.FA_rates_preCue(metadata.sessid==curr_sessid));
    distract_FA_rate=mode(metadata_distract.FA_rates_preCue(metadata_distract.sessid==curr_sessid));
    % get # trials
    ntrialscue=sum(metadata.sessid==curr_sessid,[],'omitnan');
    ntrialsdistract=sum(metadata_distract.sessid==curr_sessid,[],'omitnan');
    % calculate the number of successes
    %x1=cue_hit_rate*ntrialscue;
    %x2=distract_hit_rate*ntrialsdistract;
    % calculate odds for each cue
    odds1=cue_hit_rate/(1 - cue_hit_rate);
    odds2=distract_hit_rate/(1 - distract_hit_rate);
    % calculate the odds ratio
    odds_ratio=odds1/odds2;
    
%     % reconstruct data from summary stats: hit rates and n trials
%     hit_occurred=[ones(ceil(cue_hit_rate.*ntrialscue),1); zeros(ceil((1-cue_hit_rate).*ntrialscue),1); ...
%                   ones(ceil(distract_hit_rate.*ntrialsdistract),1); zeros(ceil((1-distract_hit_rate).*ntrialsdistract),1)];
%     FA_occurred=[ones(ceil(cue_FA_rate.*ntrialscue),1); zeros(ceil((1-cue_FA_rate).*ntrialscue),1); ...
%                   ones(ceil(distract_FA_rate.*ntrialsdistract),1); zeros(ceil((1-distract_FA_rate).*ntrialsdistract),1)];
%     cue_type=[zeros(size([ones(ceil(cue_hit_rate.*ntrialscue),1); zeros(ceil((1-cue_hit_rate).*ntrialscue),1)])); ...
%               ones(size([ones(ceil(distract_hit_rate.*ntrialsdistract),1); zeros(ceil((1-distract_hit_rate).*ntrialsdistract),1)]))];
%     cue_type_FA=[zeros(size([ones(ceil(cue_FA_rate.*ntrialscue),1); zeros(ceil((1-cue_FA_rate).*ntrialscue),1)])); ...
%               ones(size([ones(ceil(distract_FA_rate.*ntrialsdistract),1); zeros(ceil((1-distract_FA_rate).*ntrialsdistract),1)]))];
%     [b,dev,stats]=glmfit(cue_type,hit_occurred,'binomial','link','logit');
%     odds_ratio=exp(-b(2)); % cue relative to distractor
%     [b_FA,dev_FA,stats_FA]=glmfit(cue_type_FA,FA_occurred,'binomial','link','logit');
%     FA_odds_ratio=exp(-b_FA(2)); % cue relative to distractor
    metadata.odds_ratio(metadata.sessid==curr_sessid)=odds_ratio;
    metadata_distract.odds_ratio(metadata_distract.sessid==curr_sessid)=odds_ratio;
end
alltbt.odds_ratio=metadata.odds_ratio;
trialTypes.odds_ratio=metadata.odds_ratio;
distract_tbt.odds_ratio=metadata_distract.odds_ratio;
trialTypes_distract.odds_ratio=metadata_distract.odds_ratio;

end