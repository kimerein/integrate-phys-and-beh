function behaviorControl_wrapper(data_loc_array)

for i=1:size(data_loc_array,1)
    if strcmp(data_loc_array{i,3},'no_spikes')
        continue
    end
    % load behavior tbt
    load([data_loc_array{355,6} sep 'beh2_tbt.mat']);
    beh2_tbt=getChewingEnds(beh2_tbt);
    % get number of confirmatory reaches

end

end

function getConfirmReaches(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'all_reachBatch','cueZone_onVoff','first',reachWindow,[]);
% get # of reaches in behWindow, which is in seconds relative to peak of
% alignmentCompanion
nansum(dataout.y(any(~isnan(alignmentCompanion.y),2),beh:));

end