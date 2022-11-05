function saveUnitDets(dir,physiology_tbt,currAssign,onCh)

unitdets.dontUseTrials=physiology_tbt.(['unit' num2str(currAssign) '_dontUseTrials']);
unitdets.meanFR=physiology_tbt.(['unit' num2str(currAssign) '_meanFR']);
unitdets.depth=physiology_tbt.(['unit' num2str(currAssign) '_depth']);
unitdets.inStructure=physiology_tbt.(['unit' num2str(currAssign) '_inStructure']);
unitdets.halfwidth=physiology_tbt.(['unit' num2str(currAssign) '_halfwidth']);
unitdets.peakToTrough=physiology_tbt.(['unit' num2str(currAssign) '_peakToTrough']);
unitdets.avWaveforms=physiology_tbt.(['unit' num2str(currAssign) '_avWaveforms']);
unitdets.amp=physiology_tbt.(['unit' num2str(currAssign) '_amp']);
unitdets.isFS=physiology_tbt.(['unit' num2str(currAssign) '_isFS']);
unitdets.isTAN=physiology_tbt.(['unit' num2str(currAssign) '_isTAN']);
unitdets.isSPN=physiology_tbt.(['unit' num2str(currAssign) '_isSPN']);
unitdets.isLowFRThin=physiology_tbt.(['unit' num2str(currAssign) '_isLowFRThin']);
save([dir sep 'unit' num2str(currAssign) 'onCh' num2str(onCh) '_unitdets'],'unitdets');

end