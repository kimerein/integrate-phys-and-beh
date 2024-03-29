function phys_tbt=checkWaveformsDuringOpto(phys_tbt,spikes)

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

optowvfms=zeros(length(goodUnits),size(spikes.waveforms,2)*size(spikes.waveforms,3));
optowvfms_count=zeros(length(goodUnits),1);
nooptowvfms=zeros(length(goodUnits),size(spikes.waveforms,2)*size(spikes.waveforms,3));
nooptowvfms_count=zeros(length(goodUnits),1);
firstafteropto=[];
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
        timewindow_noopto=phys_tbt.phys_timepoints(i,[find(~isnan(phys_tbt.phys_timepoints(i,:)),1,'first') firstopto-1]);
        % get av waveform of spikes occurring in this time range
        for j=1:length(goodUnits)
            if nansum(spikes.unwrapped_times>timewindow(1) & spikes.unwrapped_times<timewindow(2) & spikes.assigns==goodUnits(j))>0
                optowvfm=nanmean(spikes.waveforms(spikes.unwrapped_times>timewindow(1) & spikes.unwrapped_times<timewindow(2) & spikes.assigns==goodUnits(j),:,:),1);
                optowvfm=reshape(optowvfm,1,size(spikes.waveforms,2)*size(spikes.waveforms,3));
                optowvfms(j,:)=addm(optowvfms(j,:),optowvfm);
                optowvfms_count(j)=optowvfms_count(j)+1;
            end    
            if nansum(spikes.unwrapped_times>timewindow_noopto(1) & spikes.unwrapped_times<timewindow_noopto(2) & spikes.assigns==goodUnits(j))>0
                nooptowvfm=nanmean(spikes.waveforms(spikes.unwrapped_times>timewindow_noopto(1) & spikes.unwrapped_times<timewindow_noopto(2) & spikes.assigns==goodUnits(j),:,:),1);
                nooptowvfm=reshape(nooptowvfm,1,size(spikes.waveforms,2)*size(spikes.waveforms,3));
                nooptowvfms(j,:)=addm(nooptowvfms(j,:),nooptowvfm);
                nooptowvfms_count(j)=nooptowvfms_count(j)+1;
            end
        end
    else
        if isempty(firstafteropto)
            timewindow_noopto=phys_tbt.phys_timepoints(i,[find(~isnan(phys_tbt.phys_timepoints(i,:)),1,'first') find(~isnan(phys_tbt.phys_timepoints(i,:)),1,'last')]);
        else
            timewindow_noopto=phys_tbt.phys_timepoints(i,[find(~isnan(phys_tbt.phys_timepoints(i,:)),1,'first') firstafteropto]);
        end
        % get av waveform of spikes occurring in this time range
        if length(timewindow_noopto)<2
            % just don't add this waveform
            continue
        end
        for j=1:length(goodUnits)
            if nansum(spikes.unwrapped_times>timewindow_noopto(1) & spikes.unwrapped_times<timewindow_noopto(2) & spikes.assigns==goodUnits(j))>0
                nooptowvfm=nanmean(spikes.waveforms(spikes.unwrapped_times>timewindow_noopto(1) & spikes.unwrapped_times<timewindow_noopto(2) & spikes.assigns==goodUnits(j),:,:),1);
                nooptowvfm=reshape(nooptowvfm,1,size(spikes.waveforms,2)*size(spikes.waveforms,3));
                nooptowvfms(j,:)=addm(nooptowvfms(j,:),nooptowvfm);
                nooptowvfms_count(j)=nooptowvfms_count(j)+1;
            end
        end
    end
end
for j=1:length(goodUnits)
    optowvfms(j,:)=optowvfms(j,:)./optowvfms_count(j);
    nooptowvfms(j,:)=nooptowvfms(j,:)./nooptowvfms_count(j);
    phys_tbt.(['unit' num2str(goodUnits(j)) '_nooptowvfm'])=nooptowvfms(j,:);
    phys_tbt.(['unit' num2str(goodUnits(j)) '_optowvfm'])=optowvfms(j,:);
end

end

function C=addm(A,B)

tmp=cat(3,A,B);
C=nansum(tmp,3);

end