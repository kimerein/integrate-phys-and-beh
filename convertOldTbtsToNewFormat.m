function tbt=convertOldTbtsToNewFormat(tbt)

% tbt.optoOn=zeros(size(tbt.times));
% tbt.optoZone=zeros(size(tbt.times));
% tbt.cueZone_onVoff=tbt.cue;
% tbt.times_wrt_trial_start=tbt.times-repmat(min(tbt.times,[],2,'omitnan'),1,size(tbt.times,2));
% tbt.reachBatch_miss_reachStarts=tbt.miss_reachStarts;
% tbt.reachBatch_success_reachStarts_pawOnWheel=tbt.success_reachStarts_pawOnWheel;
% tbt.reachBatch_drop_reachStarts_pawOnWheel=tbt.drop_reachStarts_pawOnWheel;
% tbt.reachBatch_miss_reachStarts_pawOnWheel=tbt.miss_reachStarts_pawOnWheel;
% tbt.reachBatch_success_reachStarts=tbt.success_reachStarts;
% tbt.reachBatch_drop_reachStarts=tbt.drop_reachStarts;
% tbt.all_reachBatch=tbt.reachStarts;


figure(); plot(nanmean(tbt.times_tbt,1));
pause;
f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    if size(temp,2)<486
        tbt.(f{i})=[temp nan(size(temp,1),486-size(temp,2))];
    else
        tbt.(f{i})=temp(:,1:486);
    end
end

tbt.movie_distractor=tbt.distractor_tbt;
tbt.pelletLoaded=tbt.pelletLoaded_tbt;
tbt.timesfromarduino=tbt.times_tbt;
tbt.optoOn=zeros(size(tbt.times_tbt));
tbt.optoZone=zeros(size(tbt.times_tbt));
tbt.isChewing=tbt.eating_tbt;
tbt.movieframeinds=tbt.movieframeinds_tbt;
tbt.cueZone_onVoff=tbt.cue_tbt;
tbt.reachStarts=tbt.reachStarts_tbt;
tbt.pelletPresented=tbt.pelletPresented_tbt;
tbt.reachStarts_pelletPresent=tbt.reach_pelletPresent_tbt;
tbt.success_reachStarts=tbt.success_tbt;
tbt.drop_reachStarts=tbt.drop_tbt;
tbt.miss_reachStarts=tbt.miss_tbt;
tbt.pelletmissingreach_reachStarts=tbt.reach_wout_pellet_tbt;
tbt.eating=tbt.eating_tbt;
tbt.pawOnWheel=tbt.paw_from_wheel_tbt;
tbt.success_reachStarts_pawOnWheel=tbt.success_pawOnWheel_tbt;
tbt.drop_reachStarts_pawOnWheel=tbt.drop_pawOnWheel_tbt;
tbt.miss_reachStarts_pawOnWheel=tbt.miss_pawOnWheel_tbt;
tbt.reach_ongoing=tbt.reach_ongoing_tbt;
tbt.times=tbt.times_tbt;
tbt.times_wrt_trial_start=tbt.times_tbt-repmat(min(tbt.times_tbt,[],2,'omitnan'),1,size(tbt.times_tbt,2));
tbt.reachBatch_miss_reachStarts=tbt.miss_tbt;
tbt.reachBatch_success_reachStarts_pawOnWheel=tbt.success_pawOnWheel_tbt;
tbt.reachBatch_drop_reachStarts_pawOnWheel=tbt.drop_pawOnWheel_tbt;
tbt.reachBatch_miss_reachStarts_pawOnWheel=tbt.miss_pawOnWheel_tbt;
tbt.reachBatch_success_reachStarts=tbt.success_tbt;
tbt.reachBatch_drop_reachStarts=tbt.drop_tbt;
tbt.all_reachBatch=tbt.reachStarts_tbt;