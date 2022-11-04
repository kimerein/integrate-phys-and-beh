function plotBehaviorSessEvents_wrapper(alltbt,metadata,trialTypes,sessid_rangevalues)





sess_filt.alltbt=alltbt; sess_filt.metadata=metadata; sess_filt.trialTypes=trialTypes;
tbt_filter.sortField='sessid'; tbt_filter.range_values=sessid_rangevalues; % sess 45
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
[sess_filt.alltbt,sess_filt.trialTypes,sess_filt.metadata]=filtTbt(sess_filt.alltbt,sess_filt.trialTypes,tbt_filter.sortField,tbt_filter.range_values,sess_filt.metadata,tbt_filter.clock_progress);
sess_filt.alltbt.timesfromarduino=sess_filt.alltbt.times;
sess_filt.settings=plotCueTriggered_settings();
sess_filt.settings.plotfields={'cueZone_onVoff','optoZone','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts'};
sess_filt.settings.plotevents=sess_filt.settings.plotfields;
sess_filt.settings.eventOutlines={'b','m',[0 0.7500 0],'r','c',[0.8 0.8 0.8]};
figure(); plot(sess_filt.alltbt.optoZone(1:20,:)');
temp=input('Thresh for optoZone: ',"s");
if strcmp(temp,'optoOn')
    temp=0.5;
    sess_filt.settings.plotfields{2}='optoOn';
    sess_filt.settings.plotevents{2}=sess_filt.settings.plotfields{2};
else
    temp=eval(temp);
end
sess_filt.settings.eventThresh={[0.5],[temp],[0.5],[0.5],[0.5],[0.5]};
sess_filt.settings.eventColors={'b','none',[0 0.7500 0],'r','c',[0.8 0.8 0.8]};
sess_filt.settings.firstN={'all',[5],'all','all','all','all'};
sess_filt.settings.histoplotfields={'cueZone_onVoff','all_reachBatch'};
sess_filt.settings.shading_type=[];
plotBehavior(sess_filt.alltbt,'cueZone_onVoff',false,1:size(sess_filt.alltbt.cueZone_onVoff,1),sess_filt.settings);