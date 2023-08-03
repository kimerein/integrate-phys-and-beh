function iterateThroughDprimeChunksForMouse(alltbt,metadata,trialTypes,whichMouse,dprime_bins,saveIntoWhichFolder,saveName)

alltbt.mouseid=metadata.mouseid;
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set

metadata.takemice=ismember(metadata.mouseid,[whichMouse]);
alltbt.takemice=metadata.takemice; trialTypes.takemice=metadata.takemice;
tbt_filter.sortField='takemice';
tbt_filter.range_values=[0.5 1.5]; 
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name; temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';  temp=temp(~isspace(temp)); tbt_filter.name=temp; tbt_filter.clock_progress=true;
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

for i=1:length(dprime_bins)-1
    try
        currdprimerange=dprime_bins(i:i+1);
        metadata=backup.metadata; alltbt=backup.alltbt; trialTypes=backup.trialTypes;
        tbt_filter.sortField='dprimes';
        tbt_filter.range_values=currdprimerange;
        [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

        test.nInSequence=[2];
        trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
        trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
        test.templateSequence2_cond=eval(trial1);
        test.templateSequence2_end=eval(trial2);
        test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
        test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
        saveDir2=[saveDir '\' test.event_name];
        mkdir(saveDir2);
        fakeCueInd=50; skipCorrected=true;
        [dataset,~]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
        plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
        xlim([3.255-1 3.255+3]);

        prompt = {'Reach window:','Change window:'};
        dlgtitle = 'Input';
        answer = inputdlg(prompt,dlgtitle);
        if isempty(dataset.realDistributions.templateSequence1_cond)
            continue
        end

        temp1=answer{1}; r=regexp(temp1,'\.'); temp1name=[temp1(1:r-1) 'poi' temp1(r+1:end)];
        temp2=answer{2}; r=regexp(temp2,'\.'); temp2name=[temp2(1:r-1) 'poi' temp2(r+1:end)];

        scriptToMakeOutcomeFigure4(alltbt,trialTypes,metadata,str2num(answer{1}),str2num(answer{2}),saveDir,saveIntoWhichFolder,tbt_filter,[temp2name temp1name saveName 'dpbin' num2str(i) 'simpleMouse' num2str(whichMouse)]);
    catch
    end
end