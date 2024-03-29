function [halfwidth, peakToTrough, amp, avWaveforms]=getSUWaveformFeatures(spikes,Fs)

% Finds features of waveform on channel where spike was biggest

a=unique(spikes.assigns);

avWaveforms=zeros(length(a),size(spikes.waveforms,2),size(spikes.waveforms,3));
for j=1:length(a)
    avWaveforms(j,:,:)=mean(spikes.waveforms(spikes.assigns==a(j),:,:),1);
    % set first element of average waveform to 0
    avWaveforms(j,:,:)=avWaveforms(j,:,:)-repmat(avWaveforms(j,1,:),[1 size(avWaveforms,2) 1]);
end
% find ch where spike was biggest
evChsForAssigns=nan(1,length(a));
for i=1:length(a)
    minAcrossTime=min(avWaveforms(i,:,:),[],2,'omitnan'); % extracellular recording so negative deflection
    [~,evChsForAssigns(i)]=min(minAcrossTime,[],3,'omitnan');
end
    
unitAvs=zeros(length(a),1);
peakToTrough=zeros(length(a),1);
amp=zeros(length(a),1);
for j=1:length(a)
    currEventCh=evChsForAssigns(j);
    shift1=double(-avWaveforms(j,:,currEventCh))-min(double(-avWaveforms(j,:,currEventCh)));
    [peak,peakInd]=findpeaks(shift1,'SORTSTR','descend','NPEAKS',1);
    halfAmp=(peak-shift1(1))/2;
    amp(j)=peak-shift1(1);
    shift2=shift1-halfAmp;
    point1=find(shift2(peakInd:-1:1)<0,1,'first');
    if isempty(point1)
        point1=1;
    end
    point1=peakInd-point1+1;
    point2=find(shift2(peakInd:end)<0,1,'first');
    if isempty(point2)
        point2=length(shift2(peakInd:end));
    end
    point2=peakInd+point2-1;
    unitAvs(j)=point2-point1;
    
    [peak,peakInd]=findpeaks(double(-avWaveforms(j,:,currEventCh)),'SORTSTR','descend','NPEAKS',1);
    if length(avWaveforms(j,:,currEventCh))-peakInd<2
        troughInd=length(avWaveforms(j,:,currEventCh));
    else
        [trough,troughInd]=findpeaks(double(avWaveforms(j,peakInd:end,currEventCh)),'SORTSTR','descend','NPEAKS',1);
    end
    if isempty(troughInd)
        troughInd=length(avWaveforms(j,:,currEventCh));
    else
        troughInd(1)=troughInd(1)+peakInd-1;
    end
    peakToTrough(j)=(troughInd(1)-peakInd(1))/Fs;   
end
unitAvs=unitAvs/32000;
halfwidth=unitAvs;