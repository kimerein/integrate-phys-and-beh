function [dontUseTrials,averageFiringRate]=inferUnitStability(f,physiology_tbt,dsinds,percentThresh,timeStretchThresh,plotInference)

if ischar(f)
    % this is a figure
    a=openfig(f,'invisible');
    %'/Users/kim/Desktop/Example data/20210621 Mar_6/SU aligned to behavior/unit201onCh11__QC.fig'
    % load firing rate and waveform amplitude across experiment from Kim's SU QC figure
    dataObjs = findobj(a,'-property','YData','Type','Line');
    % works by assuming firing rate is the only right-aligned plot
    for i=1:length(dataObjs)
        if strcmp(dataObjs(i).Parent.YAxisLocation,'right')
                fr_time=dataObjs(i).XData;
                fr_y=dataObjs(i).YData; % firing rate across expt
                wv_time=dataObjs(i).XData;
                wv_y=dataObjs(i).YData; % waveform amplitude across expt
                break
        end
    end
else
    % pass in data directly to save time opening fig
    if isempty(f)
        dontUseTrials=[];
        averageFiringRate=0;
        return
    end
    fr_time=f.fr_time;
    fr_y=f.fr_y; % firing rate across expt
    wv_time=f.wv_time;
    wv_y=f.wv_y; % waveform amplitude across experiment
end

ds=downSampAv(fr_y,dsinds); ds_time=downSampAv(fr_time,dsinds); 
normed_ds=ds./max(ds,[],2,'omitnan');
% are there any stretches of normed_ds less than percentThresh/100
% AND also longer in time than timeStretchThresh
lowFR=normed_ds<(percentThresh/100);
testedSoFar=zeros(size(lowFR));
dropTimes=zeros(size(lowFR));
timestretchinds=floor(timeStretchThresh/(ds_time(2)-ds_time(1)));
for i=1:length(lowFR)
    fl=find(lowFR==true,1,'first');
    % is this stretch long enough and low firing rate
    % continue until firing rate high again OR END
    fnexthigh=find(lowFR(fl+timestretchinds:end)==false,1,'first');
    if isempty(fnexthigh)
        currin=length(lowFR);
    else
        currin=fl+timestretchinds-1+fnexthigh-1;
    end
    testIndsForLow=fl:fl+timestretchinds-1;
    if fl+timestretchinds-1>length(lowFR)
        testIndsForLow=fl:length(lowFR);
    end
    if all(lowFR(testIndsForLow)==true)
        dropTimes(fl:currin)=1;
        testedSoFar(fl:currin)=1;
    else
        testedSoFar(fl)=1;
    end
    lowFR(testedSoFar==1)=false;
end
timestep=ds_time(2)-ds_time(1);
        
if plotInference==true
    figure();
    scatter(ds_time,dropTimes);
    hold on;
    plot(ds_time,normed_ds);
end

averageFiringRate=mean(ds(dropTimes==0),'all','omitnan');

% find trials that correspond to these time periods
dontUseTrials=zeros(size(physiology_tbt.phys_timepoints,1),1);
for i=1:size(physiology_tbt.phys_timepoints,1)
    temp=physiology_tbt.phys_timepoints(i,:);
    temp=temp(~isnan(temp) & temp~=0);
    isInRange_tempBegin=isWithinTimeRange(temp(1),ds_time(dropTimes==1),timestep);
    isInRange_tempEnd=isWithinTimeRange(temp(end),ds_time(dropTimes==1),timestep);
    if isInRange_tempBegin==1 && isInRange_tempEnd==1
        % drop this trial
        dontUseTrials(i)=1;
    end
end

end

function isInRange=isWithinTimeRange(thisTime,dropTimes,timestep)

isInRange=false;
for i=1:length(dropTimes)
    currRange=[dropTimes(i)-timestep/2 dropTimes(i)+timestep/2];
    if thisTime>=currRange(1) && thisTime<=currRange(2)
        isInRange=true;
        break
    end
end

end