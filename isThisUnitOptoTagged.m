function addTag=isThisUnitOptoTagged(optoAligned_phys_tbt,currAssign,ifTaggedWhichTag)

ds_for_sustained=10; % downsamp bins

pval1=optoAligned_phys_tbt.(['unit' num2str(currAssign) '_pval1']);
pval2=optoAligned_phys_tbt.(['unit' num2str(currAssign) '_pval2']);
avAlignedToOpto=optoAligned_phys_tbt.(['unit' num2str(currAssign) '_avAlignedToOpto']);
optoOnInUnitTimes=optoAligned_phys_tbt.('optoOnInUnitTimes');

switch ifTaggedWhichTag
    case 'none'
        addTag='';
    case 'D1'
        if earlyIncrease(pval1) && earlyIncrease(pval2) && sustainedIncrease(avAlignedToOpto,optoOnInUnitTimes,2,ds_for_sustained)
            addTag='D1tagged';
        else
            addTag='';
        end   
    case 'A2a'
        if earlyIncrease(pval1) && earlyIncrease(pval2) && sustainedIncrease(avAlignedToOpto,optoOnInUnitTimes,2,ds_for_sustained)
            addTag='A2atagged';
        else
            addTag='';
        end
    case 'Nkx2.1'
        if earlyIncrease(pval1) && earlyIncrease(pval2) && sustainedIncrease(avAlignedToOpto,optoOnInUnitTimes,3,ds_for_sustained)
            addTag='Nkxtagged';
        else
            addTag='';
        end
    otherwise
        error('do not recognize ifTaggedWhichTag in isThisUnitOptoTagged.m');
end

suppressFigs=false;
if suppressFigs==false
    if ~strcmp(addTag,'')
        figure(); 
        plot(avAlignedToOpto,'Color','k'); hold on;
        plot(optoOnInUnitTimes,'Color','r');
        pause;
    end
end

end

function out=sustainedIncrease(avAlignedToOpto,optoOnInUnitTimes,nTimesStdevAtBaseline,ds)

avAlignedToOpto=downSampAv(avAlignedToOpto,ds);
optoOnInUnitTimes=downSampAv(optoOnInUnitTimes,ds);

meanAtBaseline=mean(avAlignedToOpto(optoOnInUnitTimes==0));
stdAtBaseline=std(avAlignedToOpto(optoOnInUnitTimes==0));
meanDuringOpto=mean(avAlignedToOpto(optoOnInUnitTimes==1));
out=meanDuringOpto>meanAtBaseline+(stdAtBaseline*nTimesStdevAtBaseline);

end

function out=earlyIncrease(pval1)

% assuming binsize in alignToOpto is binsize=20 in ms
out=any(pval1(1:3)<0.05);

end