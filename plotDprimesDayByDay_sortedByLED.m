function [temp_LED1,temp_LED2]=plotDprimesDayByDay_sortedByLED(alltbt,trialTypes,metadata,ledField,ledVals,dayField,reachName,cueName)

% ledField is the field in trialTypes
% ledVals is a 2-element vector of values in ledField to compare in plot
% dayField is field in metadata to use as day-by-day

if length(ledVals)~=2
    error('Expected ledVals to be a 2-element vector');
end

settings=reachExpt_analysis_settings();
settings.reachAfterCueWindow_start=0;
settings.reachAfterCueWindow_end=2.5;
% filter settings
tbt_filter.sortField='mouseid';
tbt_filter.clock_progress=false;

u_mice=unique(metadata.mouseid);
nth_sess_for_each_mouse=cell(1,length(u_mice));
for i=1:length(u_mice)
    currmouse=u_mice(i);
    nth_sess_for_each_mouse{i}=unique(metadata.nth_session(metadata.mouseid==currmouse));
end

% for each session for each mouse, get dprime with and without LED
dprimes_for_each_day_per_mouse_LED1=cell(1,length(u_mice));
dprimes_for_each_day_per_mouse_LED2=cell(1,length(u_mice));
maxDays=0;
for i=1:length(u_mice)
    currmouse=u_mice(i);
    % filter alltbt by mouse and by LED value
    tbt_filter.range_values=[currmouse-0.5 currmouse+0.5];
    [temp_alltbt,temp_trialTypes,temp_metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    % first LED val
    [temp_alltbt_LED1,temp_trialTypes_LED1,temp_metadata_LED1]=filtTbt(temp_alltbt,temp_trialTypes,ledField,[ledVals(1)-0.0001 ledVals(1)+0.0001],temp_metadata,tbt_filter.clock_progress);
    [dprimes_for_each_day_per_mouse_LED1{i}]=get_dprimes(temp_alltbt_LED1,temp_trialTypes_LED1,temp_metadata_LED1,reachName,cueName,settings,dayField);
    % second LED val
    [temp_alltbt_LED2,temp_trialTypes_LED2,temp_metadata_LED2]=filtTbt(temp_alltbt,temp_trialTypes,ledField,[ledVals(2)-0.0001 ledVals(2)+0.0001],temp_metadata,tbt_filter.clock_progress);
    [dprimes_for_each_day_per_mouse_LED2{i}]=get_dprimes(temp_alltbt_LED2,temp_trialTypes_LED2,temp_metadata_LED2,reachName,cueName,settings,dayField);
    if length(dprimes_for_each_day_per_mouse_LED1{i})>maxDays
        maxDays=length(dprimes_for_each_day_per_mouse_LED1{i});
    end
end

temp_LED1=nan(maxDays,length(dprimes_for_each_day_per_mouse_LED1));
temp_LED2=nan(maxDays,length(dprimes_for_each_day_per_mouse_LED2));
for i=1:maxDays
    for j=1:length(dprimes_for_each_day_per_mouse_LED1)
        if i>length(dprimes_for_each_day_per_mouse_LED1{j})
            continue
        else
            teeemp=dprimes_for_each_day_per_mouse_LED1{j};
            temp_LED1(i,j)=teeemp(i);
        end
        if i>length(dprimes_for_each_day_per_mouse_LED2{j})
            continue
        else
            teeemp=dprimes_for_each_day_per_mouse_LED2{j};
            temp_LED2(i,j)=teeemp(i);
        end
    end
end
figure();
plot(nanmean(temp_LED1,2),'Color','k');
hold on; 
plot(nanmean(temp_LED2,2),'Color','r');



end

function [dprimes,metadata,alltbt,out]=get_dprimes(alltbt,out,metadata,reachName,cueName,settings,dayField)

[dprimes]=get_dprime_per_session(alltbt,out,metadata,reachName,cueName,settings,dayField);

metadata.dprimes=nan(size(metadata.(dayField)));
u=unique(metadata.(dayField));
for i=1:length(u)
    curru=u(i);
    metadata.dprimes(ismember(metadata.(dayField),curru))=dprimes(i);
end

alltbt.dprimes=metadata.dprimes;
out.dprimes=metadata.dprimes;

end

function [dprimes,hit_rates,FA_rates,out]=get_dprime_per_session(tbt,out,metadata,whichReach,nameOfCue,settings,dayField)

if isempty(settings)
    settings=reachExpt_analysis_settings();
end
hitWindow_start=settings.reachAfterCueWindow_start; % wrt cue onset
hitWindow_end=settings.reachAfterCueWindow_end; % wrt cue onset
FAWindow_start=settings.preCueWindow_start; % wrt trial onset
FAWindow_end=settings.preCueWindow_end; % wrt trial onset

% calculate hit rates per session

% Convert time window wrt cue onset into indices into data
cueInd=find(nanmean(tbt.(nameOfCue),1)>settings.lowThresh,1,'first');
startInds=floor(abs(hitWindow_start)/mode(diff(nanmean(tbt.times,1))));
if hitWindow_start<0
    startInds=-startInds;
end
endInds=floor(abs(hitWindow_end)/mode(diff(nanmean(tbt.times,1))));
if hitWindow_end<0
    endInds=-endInds;
end
startInds=cueInd+startInds;
endInds=cueInd+endInds;
if startInds<1
    startInds=1;
end
if endInds>size(tbt.times,2)
    endInds=size(tbt.times,2);
end

temp=tbt.(whichReach);
hits=any(temp(:,startInds:endInds),2);
% throw out all nan trials
allnantrials=all(isnan(tbt.times),2);
hit_rates=trials_per_session(metadata,hits==1,allnantrials==0,dayField);

% calculate false alarm rate

% Convert time window wrt trial onset into indices into data
if ~iscell(FAWindow_start)
    useInds=zeros(1,size(temp,2));
    startInds=floor(FAWindow_start/mode(diff(nanmean(tbt.times,1))));
    endInds=floor(FAWindow_end/mode(diff(nanmean(tbt.times,1))));
    if startInds<1
        startInds=1;
    end
    if endInds>size(tbt.times,2)
        endInds=size(tbt.times,2);
    end
    useInds(startInds:endInds)=1;
else
    useInds=zeros(1,size(temp,2));
    for i=1:length(FAWindow_start)
        currStretch=FAWindow_start{i};
        startInds=floor(currStretch(1)/mode(diff(nanmean(tbt.times,1))));
        endInds=floor(currStretch(2)/mode(diff(nanmean(tbt.times,1))));
        if startInds<1
            startInds=1;
        end
        if endInds>size(tbt.times,2)
            endInds=size(tbt.times,2);
        end
        useInds(startInds:endInds)=1;
    end
end

FAs=any(temp(:,useInds==1),2);
FA_rates=trials_per_session(metadata,FAs==1,allnantrials==0,dayField);

dprimes=dprime(hit_rates,FA_rates);

% dprime for session, organize per trial
out.dprime=nan(size(out.led,1),1);
for i=1:size(out.led,1)
    temp=metadata.(dayField);
    out.dprime(i)=dprimes(temp(i));
end

end

function fractionTrialsInGroup=trials_per_session(metadata,in_group_trials,total_trials,dayField)

[sessid,sessStartInd]=unique(metadata.(dayField));
fractionTrialsInGroup=nan(1,length(sessid));
sessStartInd=[sessStartInd; length(metadata.(dayField))+1];
for i=1:length(sessid)
    curr_sessStartInd=sessStartInd(i);
    sessinds=curr_sessStartInd:sessStartInd(i+1)-1;
    fractionTrialsInGroup(i)=sum(in_group_trials(sessinds)==1)/sum(total_trials(sessinds)==1);
end

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end