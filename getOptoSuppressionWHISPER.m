function [confr,optfr,supp,unitName,con,opt,tim,optoOn]=getOptoSuppressionWHISPER(optoAligned_phys_tbt,unitName,ds,takeOnlyOptoStartsInRange,physiology_tbt)

if ~isempty(physiology_tbt)
    f=nan(size(physiology_tbt.opto,1),1);
    for i=1:size(physiology_tbt.opto,1)
        isf=find(physiology_tbt.opto(i,:)>0.1,1,'first');
        if ~isempty(isf)
            f(i)=isf;
        end
    end
end

tim=nanmean(optoAligned_phys_tbt.unitTimes_wrt_trial_start,1);
fr=optoAligned_phys_tbt.(unitName);
if ~isempty(takeOnlyOptoStartsInRange)
    con=nanmean(fr(optoAligned_phys_tbt.hasOpto==0,:),1);
    opt=nanmean(fr(optoAligned_phys_tbt.hasOpto'==1 & f>takeOnlyOptoStartsInRange(1) & f<takeOnlyOptoStartsInRange(2),:),1);
else
    con=nanmean(fr(optoAligned_phys_tbt.hasOpto==0,:),1);
    opt=nanmean(fr(optoAligned_phys_tbt.hasOpto==1,:),1);
end

figure(); 
plot(nanmean(optoAligned_phys_tbt.cue_times,1),nanmean(optoAligned_phys_tbt.cue,1),'Color','b');
hold on;
plot(nanmean(optoAligned_phys_tbt.cue_times,1),nanmean(optoAligned_phys_tbt.opto,1),'Color','r');
title('Cue and opto');

% find suppression when opto on
confr=nanmean(con(optoAligned_phys_tbt.optoOnInUnitTimes==1));
optfr=nanmean(opt(optoAligned_phys_tbt.optoOnInUnitTimes==1));
optoOn=optoAligned_phys_tbt.optoOnInUnitTimes;
supp=optfr/confr;

figure();
plot(downSampAv(tim,ds),downSampAv(con,ds),'Color','k');
hold on;
plot(downSampAv(tim,ds),downSampAv(opt,ds),'Color','r');
plot(tim,optoAligned_phys_tbt.optoOnInUnitTimes,'Color','m');
title(num2str(supp));