function behaviorControl_wrapper(data_loc_array)

for i=1:size(data_loc_array,1)
    if strcmp(data_loc_array{i,3},'no_spikes')
        continue
    end
    % load behavior tbt
    load([data_loc_array{355,6} sep 'beh2_tbt.mat']);
    beh2_tbt=getChewingEnds(beh2_tbt);
    % get number of confirmatory reaches
    % get duration of chewing
end

end

function chewendtimes=getChewEnd(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'chewingEnds','cueZone_onVoff','first',reachWindow,[]);
close all;
% get distribution of times when chewing ends
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
startAtInd=find(dataout.x>=eventTime+behWindow(1),1,'first');
temp=dataout.y(any(~isnan(alignmentCompanion.y),2),startAtInd:find(dataout.x<=eventTime+behWindow(2),1,'last'));
chewendtimes=nan(size(temp,1),1);
for i=1:size(temp,1)
    f=find(temp(i,:)>0.5,1,'first');
    chewendtimes(i)=dataout.x(startAtInd-1+f)-eventTime;
end

end

function reaches=getConfirmReaches(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'all_reachBatch','cueZone_onVoff','first',reachWindow,[]);
close all;
% get # of reaches in behWindow, which is in seconds relative to peak of
% alignmentCompanion
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
reaches=nansum(dataout.y(any(~isnan(alignmentCompanion.y),2),find(dataout.x>=eventTime+behWindow(1),1,'first'):find(dataout.x<=eventTime+behWindow(2),1,'last')));

end