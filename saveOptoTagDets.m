function saveOptoTagDets(dir,physiology_tbt,currAssign,onCh)

optoTag.pval1=physiology_tbt.(['unit' num2str(currAssign) '_pval1']);
optoTag.pval2=physiology_tbt.(['unit' num2str(currAssign) '_pval2']);
optoTag.avAlignedToOpto=physiology_tbt.(['unit' num2str(currAssign) '_avAlignedToOpto']);
optoTag.optoOnInUnitTimes=physiology_tbt.optoOnInUnitTimes;
optoTag.unitTimes=physiology_tbt.unitTimes;
% see spikes opto aligned for waveforms inside and out of opto
save([dir sep 'unit' num2str(currAssign) 'onCh' num2str(onCh) '_optoTag'],'optoTag');

end