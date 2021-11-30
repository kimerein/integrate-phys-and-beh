function [physiology_tbt,tbt_unit_phys,strunits]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,rmExistingUnitFields)

% binsize in ms

% rmExistingUnitFields=false;

maxTrialDuration=nanmax(physiology_tbt.cuetimes_wrt_trial_start(1:end));

temp=spikes.labels(:,1); 
useAssigns=unique(temp(ismember(spikes.labels(:,2),goodUnitLabel)));

strunits={};
if rmExistingUnitFields==true
    f=fieldnames(physiology_tbt);
    for i=1:length(f)
        if ~isempty(regexp(f{i},'unit','ONCE'))
            physiology_tbt=rmfield(physiology_tbt,f{i});
        end
    end
else
    % deconflict any unit assigns
    f=fieldnames(physiology_tbt);
    for i=1:length(f)
        fi=f{i};
        fexp=regexp(fi,'unit');
        if ~isempty(fexp)
            % check whether spikes have same unit assigns
            funderscore=regexp(fi,'_');
            if ~isempty(funderscore)
                existingAssign=str2double(fi(fexp+4:funderscore-1));
                existingInd=str2double(fi(funderscore+1:end));
            else
                existingAssign=str2double(fi(fexp+4:end));
                existingInd=0;
            end
            if any(ismember(useAssigns,existingAssign))
                % this assigns already exists in physiology_tbt
                strunits{useAssigns==existingAssign}=['unit' num2str(useAssigns(i)) '_' num2str(existingInd+1)];
            end
        end
    end
end

figure();
tbt_unit_phys=cell(1,length(useAssigns));
for i=1:length(useAssigns)
    if isempty(strunits)
        strunit=['unit' num2str(useAssigns(i))];
    else
        if isempty(strunit{i})
            strunit=['unit' num2str(useAssigns(i))];
            strunits{i}=strunit;
        else
            strunit=strunits{i};
        end
    end
    [n,c,edges,x,y,ns]=psth_wStd_trialByTrial(filtspikes(spikes,0,'assigns',useAssigns(i)),binsize,bsmooth,maxTrialDuration,[],[]); % y is smoothed if 3rd arg is true
    % nan out everything after trial ends
    for j=1:size(ns,1)-1
        triallength=physiology_tbt.phys_timepoints(j+1,1)-physiology_tbt.phys_timepoints(j,1);
        ns(j,c>triallength)=nan;
    end
    physiology_tbt.(strunit)=ns;
    tbt_unit_phys{i}=ns;
    if i==1
        if isfield(physiology_tbt,'unitsum')
            sumacrossunits=physiology_tbt.unitsum;
        else
            sumacrossunits=zeros(size(ns));
        end
    end
    sumacrossunits=sumacrossunits+ns;
    plot(x,y);
    hold all;
end
physiology_tbt.unitTimes=repmat(c,size(ns,1),1);
plot(nanmean(physiology_tbt.cuetimes_wrt_trial_start,1),nanmean(physiology_tbt.cue,1)*10,'Color','b');


if isempty(strunits) % didn't reassign any unit names
    for i=1:length(useAssigns)
        strunits{i}=['unit' num2str(useAssigns(i))];
    end
end

% make psth that is sum across units
physiology_tbt.unitsum=sumacrossunits;

end