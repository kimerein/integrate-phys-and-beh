function checkWaveformsDuringOpto(phys_tbt,spikes)

goodUnits=[];
f=fieldnames(phys_tbt);
for i=1:length(f)
    if ~isempty(regexp(f{i},'unit')) 
        if strcmp(f{i},'unitTimes') || strcmp(f{i},'unitsum')
            continue
        end
        if ~isempty(regexp(f{i},'pval')) || ~isempty(regexp(f{i},'avAlignedToOpto'))
            continue
        end
        [rstart,rend]=regexp(f{i},'unit');
        s=f{i};
        goodUnits=[goodUnits str2num(s(rend+1:end))];  
    end
end

for i=1:size(phys_tbt.opto,1)
    if phys_tbt.hasOpto(i)==1
        % get times in phys when opto was on
        firstopto=find(phys_tbt.opto(i,:)>0.01,1,'first');
        if isempty(firstopto)
            continue
        end
        firstafteropto=find(phys_tbt.opto(i,firstopto:end)<0.01,1,'first');
        if isempty(firstafteropto)
            continue
        else
            firstafteropto=firstafteropto+firstopto-1;
        end
        timewindow=phys_tbt.phys_timepoints(i,[firstopto firstafteropto-1]); 
        timewindow_noopto=phys_tbt.phys_timepoints(i,[1 firstopto-1]);
        % get av waveform of spikes occurring in this time range
        for j=1:length(goodUnits)
            optowvfm=nanmean(spikes.waveforms(spikes.unwrapped_times>timewindow(1) & spikes.unwrapped_times<timewindow(2) & spikes.assigns==goodUnits(j),:,:),1);
            nooptowvfm=nanmean(spikes.waveforms(spikes.unwrapped_times>timewindow_noopto(1) & spikes.unwrapped_times<timewindow_noopto(2) & spikes.assigns==goodUnits(j),:,:),1);
        end
        optowvfm=reshape(optowvfm,size(spikes.waveforms,2)*size(spikes.waveforms,3),1);
        nooptowvfm=reshape(nooptowvfm,size(spikes.waveforms,2)*size(spikes.waveforms,3),1);
    else
        
    end
end

end