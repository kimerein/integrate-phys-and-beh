function spikes=setOtherUnitsTo_notGoodUnits(spikes,goodUnits)

temp=spikes.labels(:,2);
temp(~ismember(spikes.labels(:,1),goodUnits))=1;
spikes.labels(:,2)=temp;
