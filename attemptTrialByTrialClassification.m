function attemptTrialByTrialClassification(dd,success_Response,failure_Response,response_to_plot1,response_to_plot2,timeWindow,squishToSess,useFirstMapping,chooseShuffle,loadIn)

% timeWindow is in seconds wrt peak of aligncomp

if isempty(loadIn)
    useTensorLabels=true;

    a=load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\unitbyunit_names_to_match_cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_all_glm_coef_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_metrics_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\unitnames_glm.mat');
    indexGLMcellsIntoUnitNames=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names);
    if useTensorLabels==true
        disp('Using Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\rank 2\idx.mat');
        load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\rank 2\idx.mat');
        a.cued_success_Response.idx=idx;
    else
        load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\consensus_idx_from_glm_when_normByGLMcoefIntegral.mat');
        idx=nan(size(a.cued_success_Response.unitbyunit_x,1),1);
        idx(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)))=idx_from_glm(~isnan(indexGLMcellsIntoUnitNames));
        a.cued_success_Response.idx=idx;
    end

    whichGLMinds=[21:26];
    a.cued_success_Response=addMetricsToResponse(a.cued_success_Response,py_metrics,py_all_glm_coef,indexGLMcellsIntoUnitNames,whichGLMinds);
    r{1}=a.cued_success_Response;

    if isempty(success_Response)
        % choose type of response to plot
        response_to_plot=response_to_plot1; %'cued_success'; %'all_success'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
        % [inStructure isFS isTAN isSPN isLowFRThin]
        plotUnitCriteria=[1 0 0 1 0]; % -100 is a wildcard, else 0 (false) and 1 (true)
        getCriteriaForUnitsToPlot(plotUnitCriteria);
        % read in some units
        dd_more=cell(1,length(dd));
        for i=1:length(dd)
            dd_more{i}=[dd{i} sep response_to_plot];
        end
        whichUnitsToGrab='_'; % '_' for all units, or can be something like 'D1tagged'
        settings=settingsForStriatumUnitPlots;
        settings.maxUnitsPerSess=30;
        settings.keepAllSingleTrials=true;
        success_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settings,[]);
    end
    if isempty(failure_Response)
        % choose type of response to plot
        response_to_plot=response_to_plot2; %'uncued_success'; %'all_drop'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
        % read in some units
        dd_more=cell(1,length(dd));
        for i=1:length(dd)
            dd_more{i}=[dd{i} sep response_to_plot];
        end
        failure_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settings,[]);
    end

    r{2}=success_Response;
    r{3}=failure_Response;
    r=matchAllUnits(r);
    success_Response=removeUnitFromResponse(success_Response,r{1}.excluded==1);
    failure_Response=removeUnitFromResponse(failure_Response,r{1}.excluded==1);
    r{2}=[];
    r{3}=[];
    idx=r{1}.idx;

    % make aligncomp peaks the same
    ti=nanmean(success_Response.aligncomp_x,1);
    aligncompy=nanmean(success_Response.aligncomp_y,1);
    [~,apeak]=nanmax(aligncompy); peakAt=ti(apeak);
    success_Response.aligncomp_x=success_Response.aligncomp_x-peakAt;
    success_Response.unitbyunit_x=success_Response.unitbyunit_x-peakAt;
    ti=nanmean(failure_Response.aligncomp_x,1);
    aligncompy=nanmean(failure_Response.aligncomp_y,1);
    [~,apeak]=nanmax(aligncompy); peakAt=ti(apeak);
    failure_Response.aligncomp_x=failure_Response.aligncomp_x-peakAt;
    failure_Response.unitbyunit_x=failure_Response.unitbyunit_x-peakAt;

    % get probability of response in time window for each unit
    [~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
    temp=nanmean(success_Response.aligncomp_x,1);
    alignSuccessTime=temp(alignPeakInd);
    [~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
    temp=nanmean(failure_Response.aligncomp_x,1);
    alignFailureTime=temp(alignPeakInd);
    [~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
    successRange=[startAt endAt];
    [~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
    failureRange=[startAt endAt];

    unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
    fromWhichUnit_success=success_Response.fromWhichUnit;
    fromWhichTrial_success=success_Response.fromWhichTrial;
    fromWhichSess_success=success_Response.fromWhichSess_forTrials;
    unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
    fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
    fromWhichTrial_success=fromWhichTrial_success(success_Response.isEventInThisTrial==1);
    fromWhichSess_success=fromWhichSess_success(success_Response.isEventInThisTrial==1);
    unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
    fromWhichUnit_failure=failure_Response.fromWhichUnit;
    fromWhichTrial_failure=failure_Response.fromWhichTrial;
    fromWhichSess_failure=failure_Response.fromWhichSess_forTrials;
    unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
    fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
    fromWhichTrial_failure=fromWhichTrial_failure(failure_Response.isEventInThisTrial==1);
    fromWhichSess_failure=fromWhichSess_failure(failure_Response.isEventInThisTrial==1);

    % Make unit trial IDs for trials from different sessions
    trialoffset=nanmax([fromWhichTrial_success; fromWhichTrial_failure])+1;
    sessids=unique([fromWhichSess_success; fromWhichSess_failure]);
    fromWhichTrialID_success=nan(size(fromWhichTrial_success));
    fromWhichTrialID_failure=nan(size(fromWhichTrial_failure));
    for i=1:length(unique(sessids))
        fromWhichTrialID_success(fromWhichSess_success==sessids(i))=fromWhichTrial_success(fromWhichSess_success==sessids(i))+(i-1)*trialoffset;
        fromWhichTrialID_failure(fromWhichSess_failure==sessids(i))=fromWhichTrial_failure(fromWhichSess_failure==sessids(i))+(i-1)*trialoffset;
    end

    units=unique(success_Response.fromWhichUnit);

    % cue coef
    cuecoef=r{1}.glmcoef_index21+r{1}.glmcoef_index22+r{1}.glmcoef_index23+r{1}.glmcoef_index24+r{1}.glmcoef_index25+r{1}.glmcoef_index26;
    cuecoef_thresh_above=prctile(cuecoef,70);
    cuecoef_thresh_below=prctile(cuecoef,30);

    % For each trial, put together units belonging to idx==1 or idx==2
    isidx1unit_cue=units(idx==1 & cuecoef>cuecoef_thresh_above);
    isidx2unit_cue=units(idx==2 & cuecoef>cuecoef_thresh_above);
    isidx1unit_uncue=units(idx==1 & cuecoef<=cuecoef_thresh_below);
    isidx2unit_uncue=units(idx==2 & cuecoef<=cuecoef_thresh_below);

    % shuffle
    idxShuffle=idx(randperm(length(idx)));
    isidx1unit_cueSHUFFLE=units(idxShuffle==1 & cuecoef>cuecoef_thresh_above);
    isidx2unit_cueSHUFFLE=units(idxShuffle==2 & cuecoef>cuecoef_thresh_above);
    isidx1unit_uncueSHUFFLE=units(idxShuffle==1 & cuecoef<=cuecoef_thresh_below);
    isidx2unit_uncueSHUFFLE=units(idxShuffle==2 & cuecoef<=cuecoef_thresh_below);

    % shuffle2
    cuecoefShuffle=cuecoef(randperm(length(cuecoef)));
    isidx1unit_cueSHUFFLE2=units(idx==1 & cuecoefShuffle>cuecoef_thresh_above);
    isidx2unit_cueSHUFFLE2=units(idx==2 & cuecoefShuffle>cuecoef_thresh_above);
    isidx1unit_uncueSHUFFLE2=units(idx==1 & cuecoefShuffle<=cuecoef_thresh_below);
    isidx2unit_uncueSHUFFLE2=units(idx==2 & cuecoefShuffle<=cuecoef_thresh_below);

    % shuffle3
    isidx1unit_cueSHUFFLE3=units(idxShuffle==1 & cuecoefShuffle>cuecoef_thresh_above);
    isidx2unit_cueSHUFFLE3=units(idxShuffle==2 & cuecoefShuffle>cuecoef_thresh_above);
    isidx1unit_uncueSHUFFLE3=units(idxShuffle==1 & cuecoefShuffle<=cuecoef_thresh_below);
    isidx2unit_uncueSHUFFLE3=units(idxShuffle==2 & cuecoefShuffle<=cuecoef_thresh_below);

    uniqueTrialIDs=unique(fromWhichTrialID_success);
    idx1_fr_success_cue=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_cue=nan(length(uniqueTrialIDs),1);
    idx1_fr_success_uncue=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_uncue=nan(length(uniqueTrialIDs),1);
    idx1_n_success_cue=nan(length(uniqueTrialIDs),1);
    idx2_n_success_cue=nan(length(uniqueTrialIDs),1);
    idx1_n_success_uncue=nan(length(uniqueTrialIDs),1);
    idx2_n_success_uncue=nan(length(uniqueTrialIDs),1);
    idx1_fr_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_fr_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_n_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_n_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_n_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_n_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_sess_success_cue=nan(length(uniqueTrialIDs),1);
    idx2_sess_success_cue=nan(length(uniqueTrialIDs),1);
    idx1_sess_success_uncue=nan(length(uniqueTrialIDs),1);
    idx2_sess_success_uncue=nan(length(uniqueTrialIDs),1);

    idx1_fr_success_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_fr_success_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_n_success_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_n_success_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_n_success_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_n_success_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);

    idx1_fr_success_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_fr_success_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_fr_success_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_n_success_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_n_success_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_n_success_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_n_success_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    for i=1:length(uniqueTrialIDs)
        idx1_fr_success_cue(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_cue(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx1_fr_success_uncue(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_uncue(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));

        idx1_n_success_cue(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_cue(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx1_n_success_uncue(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_uncue(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));

        idx1_sess_success_cue(i)=nanmean(fromWhichSess_success(ismember(fromWhichUnit_success,isidx1unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_sess_success_cue(i)=nanmean(fromWhichSess_success(ismember(fromWhichUnit_success,isidx2unit_cue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx1_sess_success_uncue(i)=nanmean(fromWhichSess_success(ismember(fromWhichUnit_success,isidx1unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_sess_success_uncue(i)=nanmean(fromWhichSess_success(ismember(fromWhichUnit_success,isidx2unit_uncue) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));

        idx1_fr_success_cueSHUFFLE(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_cueSHUFFLE(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx1_fr_success_uncueSHUFFLE(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_uncueSHUFFLE(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));

        idx1_n_success_cueSHUFFLE(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_cueSHUFFLE(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx1_n_success_uncueSHUFFLE(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_uncueSHUFFLE(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));

        idx1_fr_success_cueSHUFFLE2(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_cueSHUFFLE2(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx1_fr_success_uncueSHUFFLE2(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_uncueSHUFFLE2(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));

        idx1_n_success_cueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_cueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx1_n_success_uncueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_uncueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));

        idx1_fr_success_cueSHUFFLE3(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_cueSHUFFLE3(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx1_fr_success_uncueSHUFFLE3(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));
        idx2_fr_success_uncueSHUFFLE3(i)=nanmean(unitfr_success(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i))));

        idx1_n_success_cueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_cueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_cueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_cueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx1_n_success_uncueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_success,isidx1unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
        idx2_n_success_uncueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_success,isidx2unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_success,uniqueTrialIDs(i)));
    end
    uniqueTrialIDs=unique(fromWhichTrialID_failure);
    idx1_fr_failure_cue=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_cue=nan(length(uniqueTrialIDs),1);
    idx1_fr_failure_uncue=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_uncue=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_cue=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_cue=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_uncue=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_uncue=nan(length(uniqueTrialIDs),1);
    idx1_fr_failure_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_fr_failure_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
    idx1_sess_failure_cue=nan(length(uniqueTrialIDs),1);
    idx2_sess_failure_cue=nan(length(uniqueTrialIDs),1);
    idx1_sess_failure_uncue=nan(length(uniqueTrialIDs),1);
    idx2_sess_failure_uncue=nan(length(uniqueTrialIDs),1);

    idx1_fr_failure_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_fr_failure_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_cueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_uncueSHUFFLE2=nan(length(uniqueTrialIDs),1);

    idx1_fr_failure_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_fr_failure_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_fr_failure_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_cueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx1_n_failure_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    idx2_n_failure_uncueSHUFFLE3=nan(length(uniqueTrialIDs),1);
    for i=1:length(uniqueTrialIDs)
        idx1_fr_failure_cue(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_cue(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx1_fr_failure_uncue(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_uncue(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));

        idx1_n_failure_cue(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_cue(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx1_n_failure_uncue(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_uncue(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));

        idx1_sess_failure_cue(i)=nanmean(fromWhichSess_failure(ismember(fromWhichUnit_failure,isidx1unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_sess_failure_cue(i)=nanmean(fromWhichSess_failure(ismember(fromWhichUnit_failure,isidx2unit_cue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx1_sess_failure_uncue(i)=nanmean(fromWhichSess_failure(ismember(fromWhichUnit_failure,isidx1unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_sess_failure_uncue(i)=nanmean(fromWhichSess_failure(ismember(fromWhichUnit_failure,isidx2unit_uncue) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));

        idx1_fr_failure_cueSHUFFLE(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_cueSHUFFLE(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx1_fr_failure_uncueSHUFFLE(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_uncueSHUFFLE(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));

        idx1_n_failure_cueSHUFFLE(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_cueSHUFFLE(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx1_n_failure_uncueSHUFFLE(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_uncueSHUFFLE(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));

        idx1_fr_failure_cueSHUFFLE2(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_cueSHUFFLE2(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx1_fr_failure_uncueSHUFFLE2(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_uncueSHUFFLE2(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));

        idx1_n_failure_cueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_cueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx1_n_failure_uncueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_uncueSHUFFLE2(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE2) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));

        idx1_fr_failure_cueSHUFFLE3(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_cueSHUFFLE3(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx1_fr_failure_uncueSHUFFLE3(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));
        idx2_fr_failure_uncueSHUFFLE3(i)=nanmean(unitfr_failure(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i))));

        idx1_n_failure_cueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_cueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_cueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_cueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx1_n_failure_uncueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_failure,isidx1unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
        idx2_n_failure_uncueSHUFFLE3(i)=nansum(ismember(fromWhichUnit_failure,isidx2unit_uncueSHUFFLE3) & ismember(fromWhichTrialID_failure,uniqueTrialIDs(i)));
    end

    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_fr_success_cue.mat','idx2_fr_success_cue','idx2_fr_success_cueSHUFFLE','idx2_fr_success_cueSHUFFLE2','idx2_fr_success_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_fr_success_cue.mat','idx1_fr_success_cue','idx1_fr_success_cueSHUFFLE','idx1_fr_success_cueSHUFFLE2','idx1_fr_success_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_fr_success_uncue.mat','idx2_fr_success_uncue','idx2_fr_success_uncueSHUFFLE','idx2_fr_success_uncueSHUFFLE2','idx2_fr_success_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_fr_success_uncue.mat','idx1_fr_success_uncue','idx1_fr_success_uncueSHUFFLE','idx1_fr_success_uncueSHUFFLE2','idx1_fr_success_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_n_success_cue.mat','idx2_n_success_cue','idx2_n_success_cueSHUFFLE','idx2_n_success_cueSHUFFLE2','idx2_n_success_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_n_success_cue.mat','idx1_n_success_cue','idx1_n_success_cueSHUFFLE','idx1_n_success_cueSHUFFLE2','idx1_n_success_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_n_success_uncue.mat','idx2_n_success_uncue','idx2_n_success_uncueSHUFFLE','idx2_n_success_uncueSHUFFLE2','idx2_n_success_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_n_success_uncue.mat','idx1_n_success_uncue','idx1_n_success_uncueSHUFFLE','idx1_n_success_uncueSHUFFLE2','idx1_n_success_uncueSHUFFLE3');
    
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_sess_success_cue.mat','idx1_sess_success_cue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_sess_success_cue.mat','idx2_sess_success_cue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_sess_success_uncue.mat','idx1_sess_success_uncue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_sess_success_uncue.mat','idx2_sess_success_uncue');

    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_sess_failure_cue.mat','idx1_sess_failure_cue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_sess_failure_cue.mat','idx2_sess_failure_cue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_sess_failure_uncue.mat','idx1_sess_failure_uncue');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_sess_failure_uncue.mat','idx2_sess_failure_uncue');

    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_fr_failure_cue.mat','idx2_fr_failure_cue','idx2_fr_failure_cueSHUFFLE','idx2_fr_failure_cueSHUFFLE2','idx2_fr_failure_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_fr_failure_cue.mat','idx1_fr_failure_cue','idx1_fr_failure_cueSHUFFLE','idx1_fr_failure_cueSHUFFLE2','idx1_fr_failure_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_fr_failure_uncue.mat','idx2_fr_failure_uncue','idx2_fr_failure_uncueSHUFFLE','idx2_fr_failure_uncueSHUFFLE2','idx2_fr_failure_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_fr_failure_uncue.mat','idx1_fr_failure_uncue','idx1_fr_failure_uncueSHUFFLE','idx1_fr_failure_uncueSHUFFLE2','idx1_fr_failure_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_n_failure_cue.mat','idx2_n_failure_cue','idx2_n_failure_cueSHUFFLE','idx2_n_failure_cueSHUFFLE2','idx2_n_failure_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_n_failure_cue.mat','idx1_n_failure_cue','idx1_n_failure_cueSHUFFLE','idx1_n_failure_cueSHUFFLE2','idx1_n_failure_cueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx2_n_failure_uncue.mat','idx2_n_failure_uncue','idx2_n_failure_uncueSHUFFLE','idx2_n_failure_uncueSHUFFLE2','idx2_n_failure_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\idx1_n_failure_uncue.mat','idx1_n_failure_uncue','idx1_n_failure_uncueSHUFFLE','idx1_n_failure_uncueSHUFFLE2','idx1_n_failure_uncueSHUFFLE3');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\successRange.mat','successRange');
    save('C:\Users\sabatini\Documents\trialbytrial classification\before exclude\failureRange.mat','failureRange');

else
    tes=load([loadIn '\idx1_fr_failure_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_fr_failure_cue',f)
        idx1_fr_failure_cue=tes.idx1_fr_success_cue;
        switch chooseShuffle
            case 1
                idx1_fr_failure_cueSHUFFLE=tes.idx1_fr_success_cueSHUFFLE;
            case 2
                idx1_fr_failure_cueSHUFFLE=tes.idx1_fr_success_cueSHUFFLE2;
            case 3
                idx1_fr_failure_cueSHUFFLE=tes.idx1_fr_success_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_fr_failure_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_fr_failure_cueSHUFFLE=idx1_fr_failure_cueSHUFFLE2;
            case 3
                idx1_fr_failure_cueSHUFFLE=idx1_fr_failure_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_fr_failure_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_fr_failure_uncue',f)
        idx1_fr_failure_uncue=tes.idx1_fr_success_uncue;
        switch chooseShuffle
            case 1
                idx1_fr_failure_uncueSHUFFLE=tes.idx1_fr_success_uncueSHUFFLE;
            case 2
                idx1_fr_failure_uncueSHUFFLE=tes.idx1_fr_success_uncueSHUFFLE2;
            case 3
                idx1_fr_failure_uncueSHUFFLE=tes.idx1_fr_success_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_fr_failure_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_fr_failure_uncueSHUFFLE=idx1_fr_failure_uncueSHUFFLE2;
            case 3
                idx1_fr_failure_uncueSHUFFLE=idx1_fr_failure_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_fr_success_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_fr_success_cue',f)
        idx1_fr_success_cue=tes.idx1_fr_failure_cue;
        switch chooseShuffle
            case 1
                idx1_fr_success_cueSHUFFLE=tes.idx1_fr_failure_cueSHUFFLE;
            case 2
                idx1_fr_success_cueSHUFFLE=tes.idx1_fr_failure_cueSHUFFLE2;
            case 3
                idx1_fr_success_cueSHUFFLE=tes.idx1_fr_failure_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_fr_success_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_fr_success_cueSHUFFLE=idx1_fr_success_cueSHUFFLE2;
            case 3
                idx1_fr_success_cueSHUFFLE=idx1_fr_success_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_fr_success_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_fr_success_uncue',f)
        idx1_fr_success_uncue=tes.idx1_fr_failure_uncue;
        switch chooseShuffle
            case 1
                idx1_fr_success_uncueSHUFFLE=tes.idx1_fr_failure_uncueSHUFFLE;
            case 2
                idx1_fr_success_uncueSHUFFLE=tes.idx1_fr_failure_uncueSHUFFLE2;
            case 3
                idx1_fr_success_uncueSHUFFLE=tes.idx1_fr_failure_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_fr_success_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_fr_success_uncueSHUFFLE=idx1_fr_success_uncueSHUFFLE2;
            case 3
                idx1_fr_success_uncueSHUFFLE=idx1_fr_success_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_n_failure_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_n_failure_cue',f)
        idx1_n_failure_cue=tes.idx1_n_success_cue;
        switch chooseShuffle
            case 1
                idx1_n_failure_cueSHUFFLE=tes.idx1_n_success_cueSHUFFLE;
            case 2
                idx1_n_failure_cueSHUFFLE=tes.idx1_n_success_cueSHUFFLE2;
            case 3
                idx1_n_failure_cueSHUFFLE=tes.idx1_n_success_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_n_failure_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_n_failure_cueSHUFFLE=idx1_n_failure_cueSHUFFLE2;
            case 3
                idx1_n_failure_cueSHUFFLE=idx1_n_failure_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_n_failure_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_n_failure_uncue',f)
        idx1_n_failure_uncue=tes.idx1_n_success_uncue;
        switch chooseShuffle
            case 1
                idx1_n_failure_uncueSHUFFLE=tes.idx1_n_success_uncueSHUFFLE;
            case 2
                idx1_n_failure_uncueSHUFFLE=tes.idx1_n_success_uncueSHUFFLE2;
            case 3
                idx1_n_failure_uncueSHUFFLE=tes.idx1_n_success_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_n_failure_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_n_failure_uncueSHUFFLE=idx1_n_failure_uncueSHUFFLE2;
            case 3
                idx1_n_failure_uncueSHUFFLE=idx1_n_failure_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_n_success_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_n_success_cue',f)
        idx1_n_success_cue=tes.idx1_n_failure_cue;
        switch chooseShuffle
            case 1
                idx1_n_success_cueSHUFFLE=tes.idx1_n_failure_cueSHUFFLE;
            case 2
                idx1_n_success_cueSHUFFLE=tes.idx1_n_failure_cueSHUFFLE2;
            case 3
                idx1_n_success_cueSHUFFLE=tes.idx1_n_failure_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_n_success_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_n_success_cueSHUFFLE=idx1_n_success_cueSHUFFLE2;
            case 3
                idx1_n_success_cueSHUFFLE=idx1_n_success_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx1_n_success_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx1_n_success_uncue',f)
        idx1_n_success_uncue=tes.idx1_n_failure_uncue;
        switch chooseShuffle
            case 1
                idx1_n_success_uncueSHUFFLE=tes.idx1_n_failure_uncueSHUFFLE;
            case 2
                idx1_n_success_uncueSHUFFLE=tes.idx1_n_failure_uncueSHUFFLE2;
            case 3
                idx1_n_success_uncueSHUFFLE=tes.idx1_n_failure_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx1_n_success_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx1_n_success_uncueSHUFFLE=idx1_n_success_uncueSHUFFLE2;
            case 3
                idx1_n_success_uncueSHUFFLE=idx1_n_success_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_fr_failure_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_fr_failure_cue',f)
        idx2_fr_failure_cue=tes.idx2_fr_success_cue;
        switch chooseShuffle
            case 1
                idx2_fr_failure_cueSHUFFLE=tes.idx2_fr_success_cueSHUFFLE;
            case 2
                idx2_fr_failure_cueSHUFFLE=tes.idx2_fr_success_cueSHUFFLE2;
            case 3
                idx2_fr_failure_cueSHUFFLE=tes.idx2_fr_success_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_fr_failure_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_fr_failure_cueSHUFFLE=idx2_fr_failure_cueSHUFFLE2;
            case 3
                idx2_fr_failure_cueSHUFFLE=idx2_fr_failure_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_fr_failure_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_fr_failure_uncue',f)
        idx2_fr_failure_uncue=tes.idx2_fr_success_uncue;
        switch chooseShuffle
            case 1
                idx2_fr_failure_uncueSHUFFLE=tes.idx2_fr_success_uncueSHUFFLE;
            case 2
                idx2_fr_failure_uncueSHUFFLE=tes.idx2_fr_success_uncueSHUFFLE2;
            case 3
                idx2_fr_failure_uncueSHUFFLE=tes.idx2_fr_success_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_fr_failure_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_fr_failure_uncueSHUFFLE=idx2_fr_failure_uncueSHUFFLE2;
            case 3
                idx2_fr_failure_uncueSHUFFLE=idx2_fr_failure_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_fr_success_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_fr_success_cue',f)
        idx2_fr_success_cue=tes.idx2_fr_failure_cue;
        switch chooseShuffle
            case 1
                idx2_fr_success_cueSHUFFLE=tes.idx2_fr_failure_cueSHUFFLE;
            case 2
                idx2_fr_success_cueSHUFFLE=tes.idx2_fr_failure_cueSHUFFLE2;
            case 3
                idx2_fr_success_cueSHUFFLE=tes.idx2_fr_failure_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_fr_success_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_fr_success_cueSHUFFLE=idx2_fr_success_cueSHUFFLE2;
            case 3
                idx2_fr_success_cueSHUFFLE=idx2_fr_success_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_fr_success_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_fr_success_uncue',f)
        idx2_fr_success_uncue=tes.idx2_fr_failure_uncue;
        switch chooseShuffle
            case 1
                idx2_fr_success_uncueSHUFFLE=tes.idx2_fr_failure_uncueSHUFFLE;
            case 2
                idx2_fr_success_uncueSHUFFLE=tes.idx2_fr_failure_uncueSHUFFLE2;
            case 3
                idx2_fr_success_uncueSHUFFLE=tes.idx2_fr_failure_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_fr_success_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_fr_success_uncueSHUFFLE=idx2_fr_success_uncueSHUFFLE2;
            case 3
                idx2_fr_success_uncueSHUFFLE=idx2_fr_success_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_n_failure_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_n_failure_cue',f)
        idx2_n_failure_cue=tes.idx2_n_success_cue;
        switch chooseShuffle
            case 1
                idx2_n_failure_cueSHUFFLE=tes.idx2_n_success_cueSHUFFLE;
            case 2
                idx2_n_failure_cueSHUFFLE=tes.idx2_n_success_cueSHUFFLE2;
            case 3
                idx2_n_failure_cueSHUFFLE=tes.idx2_n_success_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_n_failure_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_n_failure_cueSHUFFLE=idx2_n_failure_cueSHUFFLE2;
            case 3
                idx2_n_failure_cueSHUFFLE=idx2_n_failure_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_n_failure_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_n_failure_uncue',f)
        idx2_n_failure_uncue=tes.idx2_n_success_uncue;
        switch chooseShuffle
            case 1
                idx2_n_failure_uncueSHUFFLE=tes.idx2_n_success_uncueSHUFFLE;
            case 2
                idx2_n_failure_uncueSHUFFLE=tes.idx2_n_success_uncueSHUFFLE2;
            case 3
                idx2_n_failure_uncueSHUFFLE=tes.idx2_n_success_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_n_failure_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_n_failure_uncueSHUFFLE=idx2_n_failure_uncueSHUFFLE2;
            case 3
                idx2_n_failure_uncueSHUFFLE=idx2_n_failure_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_n_success_cue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_n_success_cue',f)
        idx2_n_success_cue=tes.idx2_n_failure_cue;
        switch chooseShuffle
            case 1
                idx2_n_success_cueSHUFFLE=tes.idx2_n_failure_cueSHUFFLE;
            case 2
                idx2_n_success_cueSHUFFLE=tes.idx2_n_failure_cueSHUFFLE2;
            case 3
                idx2_n_success_cueSHUFFLE=tes.idx2_n_failure_cueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_n_success_cue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_n_success_cueSHUFFLE=idx2_n_success_cueSHUFFLE2;
            case 3
                idx2_n_success_cueSHUFFLE=idx2_n_success_cueSHUFFLE3;
        end
    end
    tes=load([loadIn '\idx2_n_success_uncue.mat']);
    f=fieldnames(tes); 
    if ~ismember('idx2_n_success_uncue',f)
        idx2_n_success_uncue=tes.idx2_n_failure_uncue;
        switch chooseShuffle
            case 1
                idx2_n_success_uncueSHUFFLE=tes.idx2_n_failure_uncueSHUFFLE;
            case 2
                idx2_n_success_uncueSHUFFLE=tes.idx2_n_failure_uncueSHUFFLE2;
            case 3
                idx2_n_success_uncueSHUFFLE=tes.idx2_n_failure_uncueSHUFFLE3;
        end
        
    else
        load([loadIn '\idx2_n_success_uncue.mat']);
        switch chooseShuffle
            case 1
            case 2
                idx2_n_success_uncueSHUFFLE=idx2_n_success_uncueSHUFFLE2;
            case 3
                idx2_n_success_uncueSHUFFLE=idx2_n_success_uncueSHUFFLE3;
        end
    end
    tes=load([loadIn '\successRange.mat']);
    f=fieldnames(tes); 
    if ~ismember('successRange',f)
        successRange=tes.failureRange;
    else
        load([loadIn '\successRange.mat']);
    end
    tes=load([loadIn '\failureRange.mat']);
    f=fieldnames(tes); 
    if ~ismember('failureRange',f)
        failureRange=tes.successRange;
    else
        load([loadIn '\failureRange.mat']);
    end
end

% Squish to session
if squishToSess==true
    % average together all trials from same session, keeping unit types separate
    [idx1_fr_success_cue,idx2_fr_success_cue,idx1_fr_success_uncue,idx2_fr_success_uncue,idx1_fr_success_cueSHUFFLE,idx2_fr_success_cueSHUFFLE,idx1_fr_success_uncueSHUFFLE,idx2_fr_success_uncueSHUFFLE,idx1_n_success_cue,idx2_n_success_cue,idx1_n_success_uncue,idx2_n_success_uncue,idx1_n_success_cueSHUFFLE,idx2_n_success_cueSHUFFLE,idx1_n_success_uncueSHUFFLE,idx2_n_success_uncueSHUFFLE]...
        =squishToSessions(idx1_sess_success_cue,idx2_sess_success_cue,idx1_sess_success_uncue,idx2_sess_success_uncue,idx1_fr_success_cue,idx2_fr_success_cue,idx1_fr_success_uncue,idx2_fr_success_uncue,idx1_fr_success_cueSHUFFLE,idx2_fr_success_cueSHUFFLE,idx1_fr_success_uncueSHUFFLE,idx2_fr_success_uncueSHUFFLE,...
            idx1_n_success_cue,idx2_n_success_cue,idx1_n_success_uncue,idx2_n_success_uncue,idx1_n_success_cueSHUFFLE,idx2_n_success_cueSHUFFLE,idx1_n_success_uncueSHUFFLE,idx2_n_success_uncueSHUFFLE);

    [idx1_fr_failure_cue,idx2_fr_failure_cue,idx1_fr_failure_uncue,idx2_fr_failure_uncue,idx1_fr_failure_cueSHUFFLE,idx2_fr_failure_cueSHUFFLE,idx1_fr_failure_uncueSHUFFLE,idx2_fr_failure_uncueSHUFFLE,idx1_n_failure_cue,idx2_n_failure_cue,idx1_n_failure_uncue,idx2_n_failure_uncue,idx1_n_failure_cueSHUFFLE,idx2_n_failure_cueSHUFFLE,idx1_n_failure_uncueSHUFFLE,idx2_n_failure_uncueSHUFFLE]...
        =squishToSessions(idx1_sess_failure_cue,idx2_sess_failure_cue,idx1_sess_failure_uncue,idx2_sess_failure_uncue,idx1_fr_failure_cue,idx2_fr_failure_cue,idx1_fr_failure_uncue,idx2_fr_failure_uncue,idx1_fr_failure_cueSHUFFLE,idx2_fr_failure_cueSHUFFLE,idx1_fr_failure_uncueSHUFFLE,idx2_fr_failure_uncueSHUFFLE,...
            idx1_n_failure_cue,idx2_n_failure_cue,idx1_n_failure_uncue,idx2_n_failure_uncue,idx1_n_failure_cueSHUFFLE,idx2_n_failure_cueSHUFFLE,idx1_n_failure_uncueSHUFFLE,idx2_n_failure_uncueSHUFFLE);
end

% Enough units
nunitthresh=0;
enough_success=idx1_n_success_cue>nunitthresh & idx2_n_success_cue>nunitthresh & idx1_n_success_uncue>nunitthresh & idx2_n_success_uncue>nunitthresh;
enough_failure=idx1_n_failure_cue>nunitthresh & idx2_n_failure_cue>nunitthresh & idx1_n_failure_uncue>nunitthresh & idx2_n_failure_uncue>nunitthresh;
% Enough spikes
spikethresh=0;
enoughspikes_success=idx1_fr_success_cue>=spikethresh & idx2_fr_success_cue>=spikethresh & idx1_fr_success_uncue>=spikethresh & idx2_fr_success_uncue>=spikethresh;
enoughspikes_failure=idx1_fr_failure_cue>=spikethresh & idx2_fr_failure_cue>=spikethresh & idx1_fr_failure_uncue>=spikethresh & idx2_fr_failure_uncue>=spikethresh;
enough_success=enough_success & enoughspikes_success;
enough_failure=enough_failure & enoughspikes_failure;
idx1_fr_success_cue=idx1_fr_success_cue(enough_success);
idx2_fr_success_cue=idx2_fr_success_cue(enough_success);
idx1_fr_success_uncue=idx1_fr_success_uncue(enough_success);
idx2_fr_success_uncue=idx2_fr_success_uncue(enough_success);
idx1_fr_failure_cue=idx1_fr_failure_cue(enough_failure);
idx2_fr_failure_cue=idx2_fr_failure_cue(enough_failure);
idx1_fr_failure_uncue=idx1_fr_failure_uncue(enough_failure);
idx2_fr_failure_uncue=idx2_fr_failure_uncue(enough_failure);
% SHUFFLE enough units
enough_successSHUFFLE=idx1_n_success_cueSHUFFLE>nunitthresh & idx2_n_success_cueSHUFFLE>nunitthresh & idx1_n_success_uncueSHUFFLE>nunitthresh & idx2_n_success_uncueSHUFFLE>nunitthresh;
enough_failureSHUFFLE=idx1_n_failure_cueSHUFFLE>nunitthresh & idx2_n_failure_cueSHUFFLE>nunitthresh & idx1_n_failure_uncueSHUFFLE>nunitthresh & idx2_n_failure_uncueSHUFFLE>nunitthresh;
% Enough spikes
enoughspikes_successSHUFFLE=idx1_fr_success_cueSHUFFLE>=spikethresh & idx2_fr_success_cueSHUFFLE>=spikethresh & idx1_fr_success_uncueSHUFFLE>=spikethresh & idx2_fr_success_uncueSHUFFLE>=spikethresh;
enoughspikes_failureSHUFFLE=idx1_fr_failure_cueSHUFFLE>=spikethresh & idx2_fr_failure_cueSHUFFLE>=spikethresh & idx1_fr_failure_uncueSHUFFLE>=spikethresh & idx2_fr_failure_uncueSHUFFLE>=spikethresh;
enough_successSHUFFLE=enough_successSHUFFLE & enoughspikes_successSHUFFLE;
enough_failureSHUFFLE=enough_failureSHUFFLE & enoughspikes_failureSHUFFLE;
idx1_fr_success_cueSHUFFLE=idx1_fr_success_cueSHUFFLE(enough_successSHUFFLE);
idx2_fr_success_cueSHUFFLE=idx2_fr_success_cueSHUFFLE(enough_successSHUFFLE);
idx1_fr_success_uncueSHUFFLE=idx1_fr_success_uncueSHUFFLE(enough_successSHUFFLE);
idx2_fr_success_uncueSHUFFLE=idx2_fr_success_uncueSHUFFLE(enough_successSHUFFLE);
idx1_fr_failure_cueSHUFFLE=idx1_fr_failure_cueSHUFFLE(enough_failureSHUFFLE);
idx2_fr_failure_cueSHUFFLE=idx2_fr_failure_cueSHUFFLE(enough_failureSHUFFLE);
idx1_fr_failure_uncueSHUFFLE=idx1_fr_failure_uncueSHUFFLE(enough_failureSHUFFLE);
idx2_fr_failure_uncueSHUFFLE=idx2_fr_failure_uncueSHUFFLE(enough_failureSHUFFLE);

% Y AXIS
% Decode "cued vs uncued" if passed in, e.g., cued_failure and
% uncued_failure
if useFirstMapping==true
    [n,x]=histcounts((idx2_fr_success_cue-idx1_fr_success_cue)-(idx2_fr_success_uncue-idx1_fr_success_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
    [n,x]=histcounts((idx2_fr_failure_cue-idx1_fr_failure_cue)-(idx2_fr_failure_uncue-idx1_fr_failure_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');
else
    [n,x]=histcounts((idx1_fr_success_uncue - idx2_fr_success_uncue) - (idx1_fr_success_cue - idx2_fr_success_cue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
    [n,x]=histcounts((idx1_fr_failure_uncue - idx2_fr_failure_uncue) - (idx1_fr_failure_cue - idx2_fr_failure_cue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');
end

% X AXIS
% Decode "success vs failure"
if useFirstMapping==true
    [n,x]=histcounts((idx2_fr_success_cue+idx2_fr_success_uncue)-(idx1_fr_success_cue+idx1_fr_success_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
    [n,x]=histcounts((idx2_fr_failure_cue+idx2_fr_failure_uncue)-(idx1_fr_failure_cue+idx1_fr_failure_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');
else
    [n,x]=histcounts((idx2_fr_success_cue - idx2_fr_success_uncue) + (idx1_fr_success_cue - idx1_fr_success_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
    [n,x]=histcounts((idx2_fr_failure_cue - idx2_fr_failure_uncue) + (idx1_fr_failure_cue - idx1_fr_failure_uncue),-2000-40:80:2000+40);
    [n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');
end

% dprime
if useFirstMapping==true
    temp1=(idx2_fr_success_cue-idx1_fr_success_cue)-(idx2_fr_success_uncue-idx1_fr_success_uncue);
    temp2=(idx2_fr_failure_cue-idx1_fr_failure_cue)-(idx2_fr_failure_uncue-idx1_fr_failure_uncue);
else
    temp1=(idx1_fr_success_uncue - idx2_fr_success_uncue) - (idx1_fr_success_cue - idx2_fr_success_cue);
    temp2=(idx1_fr_failure_uncue - idx2_fr_failure_uncue) - (idx1_fr_failure_cue - idx2_fr_failure_cue);
end
dp=(nanmean(temp1)-nanmean(temp2))./sqrt(nanstd(temp1,[],1).^2+nanstd(temp2,[],1).^2); disp(dp);
if useFirstMapping==true
    temp3=(idx2_fr_success_cue+idx2_fr_success_uncue)-(idx1_fr_success_cue+idx1_fr_success_uncue);
    temp4=(idx2_fr_failure_cue+idx2_fr_failure_uncue)-(idx1_fr_failure_cue+idx1_fr_failure_uncue);
else
    temp3=(idx2_fr_success_cue - idx2_fr_success_uncue) + (idx1_fr_success_cue - idx1_fr_success_uncue);
    temp4=(idx2_fr_failure_cue - idx2_fr_failure_uncue) + (idx1_fr_failure_cue - idx1_fr_failure_uncue);
end
dp=(nanmean(temp3)-nanmean(temp4))./sqrt(nanstd(temp3,[],1).^2+nanstd(temp4,[],1).^2); disp(dp);

% scatter
figure(); s=scatter(temp3./length(successRange(1):successRange(2)),temp1./length(successRange(1):successRange(2)),60,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]);
vals_x_axis=temp3./length(failureRange(1):failureRange(2)); vals_y_axis=temp1./length(failureRange(1):failureRange(2));
% [~,trmx]=rmoutliers(vals_x_axis,"ThresholdFactor",3); [~,trmy]=rmoutliers(vals_y_axis,"ThresholdFactor",2);
% vals_x_axis=vals_x_axis(~trmx & ~trmy); vals_y_axis=vals_y_axis(~trmx & ~trmy);
save('C:\Users\sabatini\Documents\confidence ellipse\vals_x_axis1.mat','vals_x_axis');
save('C:\Users\sabatini\Documents\confidence ellipse\vals_y_axis1.mat','vals_y_axis');
line([nanmean(temp3./length(successRange(1):successRange(2)))-nanstd(temp3./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp3))) nanmean(temp3./length(successRange(1):successRange(2)))+nanstd(temp3./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp3)))],...
     [nanmean(temp1./length(successRange(1):successRange(2))) nanmean(temp1./length(successRange(1):successRange(2)))],'Color','k');
line([nanmean(temp3./length(successRange(1):successRange(2))) nanmean(temp3./length(successRange(1):successRange(2)))],...
     [nanmean(temp1./length(successRange(1):successRange(2)))-nanstd(temp1./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp1))) nanmean(temp1./length(successRange(1):successRange(2)))+nanstd(temp1./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp1)))],'Color','k');

figure(); s=scatter(temp4./length(failureRange(1):failureRange(2)),temp2./length(failureRange(1):failureRange(2)),60,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]);
vals_x_axis=temp4./length(failureRange(1):failureRange(2)); vals_y_axis=temp2./length(failureRange(1):failureRange(2));
% [~,trmx]=rmoutliers(vals_x_axis,"ThresholdFactor",3); [~,trmy]=rmoutliers(vals_y_axis,"ThresholdFactor",2);
% vals_x_axis=vals_x_axis(~trmx & ~trmy); vals_y_axis=vals_y_axis(~trmx & ~trmy);
save('C:\Users\sabatini\Documents\confidence ellipse\vals_x_axis2.mat','vals_x_axis');
save('C:\Users\sabatini\Documents\confidence ellipse\vals_y_axis2.mat','vals_y_axis');
line([nanmean(temp4./length(successRange(1):successRange(2)))-nanstd(temp4./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp4))) nanmean(temp4./length(successRange(1):successRange(2)))+nanstd(temp4./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp4)))],...
     [nanmean(temp2./length(successRange(1):successRange(2))) nanmean(temp2./length(successRange(1):successRange(2)))],'Color','r');
line([nanmean(temp4./length(successRange(1):successRange(2))) nanmean(temp4./length(successRange(1):successRange(2)))],...
     [nanmean(temp2./length(successRange(1):successRange(2)))-nanstd(temp2./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp2))) nanmean(temp2./length(successRange(1):successRange(2)))+nanstd(temp2./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp2)))],'Color','r');

% SHUFFLE dprime
if useFirstMapping==true
    temp1=(idx2_fr_success_cueSHUFFLE-idx1_fr_success_cueSHUFFLE)-(idx2_fr_success_uncueSHUFFLE-idx1_fr_success_uncueSHUFFLE);
    temp2=(idx2_fr_failure_cueSHUFFLE-idx1_fr_failure_cueSHUFFLE)-(idx2_fr_failure_uncueSHUFFLE-idx1_fr_failure_uncueSHUFFLE);
else
    temp1=(idx1_fr_success_uncueSHUFFLE - idx2_fr_success_uncueSHUFFLE) - (idx1_fr_success_cueSHUFFLE - idx2_fr_success_cueSHUFFLE);
    temp2=(idx1_fr_failure_uncueSHUFFLE - idx2_fr_failure_uncueSHUFFLE) - (idx1_fr_failure_cueSHUFFLE - idx2_fr_failure_cueSHUFFLE);
end
dp=(nanmean(temp1)-nanmean(temp2))./sqrt(nanstd(temp1,[],1).^2+nanstd(temp2,[],1).^2); disp(dp);
if useFirstMapping==true
    temp3=(idx2_fr_success_cueSHUFFLE+idx2_fr_success_uncueSHUFFLE)-(idx1_fr_success_cueSHUFFLE+idx1_fr_success_uncueSHUFFLE);
    temp4=(idx2_fr_failure_cueSHUFFLE+idx2_fr_failure_uncueSHUFFLE)-(idx1_fr_failure_cueSHUFFLE+idx1_fr_failure_uncueSHUFFLE);
else
    temp3=(idx2_fr_success_cueSHUFFLE - idx2_fr_success_uncueSHUFFLE) + (idx1_fr_success_cueSHUFFLE - idx1_fr_success_uncueSHUFFLE);
    temp4=(idx2_fr_failure_cueSHUFFLE - idx2_fr_failure_uncueSHUFFLE) + (idx1_fr_failure_cueSHUFFLE - idx1_fr_failure_uncueSHUFFLE);
end
dp=(nanmean(temp3)-nanmean(temp4))./sqrt(nanstd(temp3,[],1).^2+nanstd(temp4,[],1).^2); disp(dp);
% SHUFFLE scatter
figure(); s=scatter(temp3./length(successRange(1):successRange(2)),temp1./length(successRange(1):successRange(2)),60,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]); title('SHUFFLE');
line([nanmean(temp3./length(successRange(1):successRange(2)))-nanstd(temp3./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp3))) nanmean(temp3./length(successRange(1):successRange(2)))+nanstd(temp3./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp3)))],...
     [nanmean(temp1./length(successRange(1):successRange(2))) nanmean(temp1./length(successRange(1):successRange(2)))],'Color','k');
line([nanmean(temp3./length(successRange(1):successRange(2))) nanmean(temp3./length(successRange(1):successRange(2)))],...
     [nanmean(temp1./length(successRange(1):successRange(2)))-nanstd(temp1./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp1))) nanmean(temp1./length(successRange(1):successRange(2)))+nanstd(temp1./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp1)))],'Color','k');
figure(); s=scatter(temp4./length(failureRange(1):failureRange(2)),temp2./length(failureRange(1):failureRange(2)),60,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]);
line([nanmean(temp4./length(successRange(1):successRange(2)))-nanstd(temp4./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp4))) nanmean(temp4./length(successRange(1):successRange(2)))+nanstd(temp4./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp4)))],...
     [nanmean(temp2./length(successRange(1):successRange(2))) nanmean(temp2./length(successRange(1):successRange(2)))],'Color','r');
line([nanmean(temp4./length(successRange(1):successRange(2))) nanmean(temp4./length(successRange(1):successRange(2)))],...
     [nanmean(temp2./length(successRange(1):successRange(2)))-nanstd(temp2./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp2))) nanmean(temp2./length(successRange(1):successRange(2)))+nanstd(temp2./length(successRange(1):successRange(2)),[],1)./sqrt(nansum(~isnan(temp2)))],'Color','r');
title('SHUFFLE');

% vals_x_axis=temp4; vals_y_axis=temp2;
%%% [~,trmx]=rmoutliers(vals_x_axis); [~,trmy]=rmoutliers(vals_y_axis);
%%% vals_x_axis=vals_x_axis(~trmx & ~trmy); vals_y_axis=vals_y_axis(~trmx & ~trmy);

% Try all decodes
currSigns=[1 1 1 1;...
    1 1 1 -1;...
    1 1 -1 1;...
    1 1 -1 -1;...
    1 -1 1 1;...
    1 -1 1 -1;...
    1 -1 -1 1;...
    1 -1 -1 -1;...
    -1 1 1 1;...
    -1 1 1 -1;...
    -1 1 -1 1;...
    -1 1 -1 -1;...
    -1 -1 1 1;...
    -1 -1 1 -1;...
    -1 -1 -1 1;...
    -1 -1 -1 -1];
dps=nan(size(currSigns,1),2);
if useFirstMapping==true
    close all;
    
    for i=1:size(currSigns,1)
        [dp_x,dp_y]=tryOtherDecodes(idx2_fr_success_cue,idx1_fr_success_cue,idx2_fr_success_uncue,idx1_fr_success_uncue,idx2_fr_failure_cue,idx1_fr_failure_cue,idx2_fr_failure_uncue,idx1_fr_failure_uncue,currSigns(i,:),successRange,failureRange);
        dps(i,1)=dp_x; dps(i,2)=dp_y;
    end
    figure();
    imagesc([-1 -1; dps; 1 1]); % just for colorbar range
end

% save
save('C:\Users\sabatini\Documents\trialbytrial classification\idx2_fr_success_cue.mat','idx2_fr_success_cue','idx2_fr_success_cueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx1_fr_success_cue.mat','idx1_fr_success_cue','idx1_fr_success_cueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx2_fr_success_uncue.mat','idx2_fr_success_uncue','idx2_fr_success_uncueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx1_fr_success_uncue.mat','idx1_fr_success_uncue','idx1_fr_success_uncueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx2_fr_failure_cue.mat','idx2_fr_failure_cue','idx2_fr_failure_cueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx1_fr_failure_cue.mat','idx1_fr_failure_cue','idx1_fr_failure_cueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx2_fr_failure_uncue.mat','idx2_fr_failure_uncue','idx2_fr_failure_uncueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\idx1_fr_failure_uncue.mat','idx1_fr_failure_uncue','idx1_fr_failure_uncueSHUFFLE');
save('C:\Users\sabatini\Documents\trialbytrial classification\currSigns.mat','currSigns');
save('C:\Users\sabatini\Documents\trialbytrial classification\dps.mat','dps');

end

function [idx1_fr_success_cue,idx2_fr_success_cue,idx1_fr_success_uncue,idx2_fr_success_uncue,idx1_fr_success_cueSHUFFLE,idx2_fr_success_cueSHUFFLE,idx1_fr_success_uncueSHUFFLE,idx2_fr_success_uncueSHUFFLE,idx1_n_success_cue,idx2_n_success_cue,idx1_n_success_uncue,idx2_n_success_uncue,idx1_n_success_cueSHUFFLE,idx2_n_success_cueSHUFFLE,idx1_n_success_uncueSHUFFLE,idx2_n_success_uncueSHUFFLE]=squishToSessions(idx1_sess_success_cue,idx2_sess_success_cue,idx1_sess_success_uncue,idx2_sess_success_uncue,idx1_fr_success_cueINPUT,idx2_fr_success_cueINPUT,idx1_fr_success_uncueINPUT,idx2_fr_success_uncueINPUT,idx1_fr_success_cueSHUFFLEINPUT,idx2_fr_success_cueSHUFFLEINPUT,idx1_fr_success_uncueSHUFFLEINPUT,idx2_fr_success_uncueSHUFFLEINPUT,...
    idx1_n_success_cueINPUT,idx2_n_success_cueINPUT,idx1_n_success_uncueINPUT,idx2_n_success_uncueINPUT,idx1_n_success_cueINPUTSHUFFLE,idx2_n_success_cueINPUTSHUFFLE,idx1_n_success_uncueINPUTSHUFFLE,idx2_n_success_uncueINPUTSHUFFLE)

uniqueTrialIDs=unique([idx1_sess_success_cue; idx2_sess_success_cue; idx1_sess_success_uncue; idx2_sess_success_uncue]);
idx1_fr_success_cue=nan(length(uniqueTrialIDs),1);
idx2_fr_success_cue=nan(length(uniqueTrialIDs),1);
idx1_fr_success_uncue=nan(length(uniqueTrialIDs),1);
idx2_fr_success_uncue=nan(length(uniqueTrialIDs),1);

idx1_n_success_cue=nan(length(uniqueTrialIDs),1);
idx2_n_success_cue=nan(length(uniqueTrialIDs),1);
idx1_n_success_uncue=nan(length(uniqueTrialIDs),1);
idx2_n_success_uncue=nan(length(uniqueTrialIDs),1);

idx1_n_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx2_n_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx1_n_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx2_n_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);

idx1_fr_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx2_fr_success_cueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx1_fr_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);
idx2_fr_success_uncueSHUFFLE=nan(length(uniqueTrialIDs),1);

for i=1:length(usess)
    idx1_fr_success_cue(i)=nanmean(idx1_fr_success_cueINPUT(ismember(idx1_sess_success_cue,uniqueTrialIDs(i))));
    idx2_fr_success_cue(i)=nanmean(idx2_fr_success_cueINPUT(ismember(idx2_sess_success_cue,uniqueTrialIDs(i))));
    idx1_fr_success_uncue(i)=nanmean(idx1_fr_success_uncueINPUT(ismember(idx1_sess_success_uncue,uniqueTrialIDs(i))));
    idx2_fr_success_uncue(i)=nanmean(idx2_fr_success_uncueINPUT(ismember(idx2_sess_success_uncue,uniqueTrialIDs(i))));

    idx1_n_success_cue(i)=nanmean(idx1_n_success_cueINPUT(ismember(idx1_sess_success_cue,uniqueTrialIDs(i))));
    idx2_n_success_cue(i)=nanmean(idx2_n_success_cueINPUT(ismember(idx2_sess_success_cue,uniqueTrialIDs(i))));
    idx1_n_success_uncue(i)=nanmean(idx1_n_success_uncueINPUT(ismember(idx1_sess_success_uncue,uniqueTrialIDs(i))));
    idx2_n_success_uncue(i)=nanmean(idx2_n_success_uncueINPUT(ismember(idx2_sess_success_uncue,uniqueTrialIDs(i))));

    idx1_fr_success_cueSHUFFLE(i)=nanmean(idx1_fr_success_cueSHUFFLEINPUT(ismember(idx1_sess_success_cue,uniqueTrialIDs(i))));
    idx2_fr_success_cueSHUFFLE(i)=nanmean(idx2_fr_success_cueSHUFFLEINPUT(ismember(idx2_sess_success_cue,uniqueTrialIDs(i))));
    idx1_fr_success_uncueSHUFFLE(i)=nanmean(idx1_fr_success_uncueSHUFFLEINPUT(ismember(idx1_sess_success_uncue,uniqueTrialIDs(i))));
    idx2_fr_success_uncueSHUFFLE(i)=nanmean(idx2_fr_success_uncueSHUFFLEINPUT(ismember(idx2_sess_success_uncue,uniqueTrialIDs(i))));

    idx1_n_success_cueSHUFFLE(i)=nanmean(idx1_n_success_cueINPUTSHUFFLE(ismember(idx1_sess_success_cue,uniqueTrialIDs(i))));
    idx2_n_success_cueSHUFFLE(i)=nanmean(idx2_n_success_cueINPUTSHUFFLE(ismember(idx2_sess_success_cue,uniqueTrialIDs(i))));
    idx1_n_success_uncueSHUFFLE(i)=nanmean(idx1_n_success_uncueINPUTSHUFFLE(ismember(idx1_sess_success_uncue,uniqueTrialIDs(i))));
    idx2_n_success_uncueSHUFFLE(i)=nanmean(idx2_n_success_uncueINPUTSHUFFLE(ismember(idx2_sess_success_uncue,uniqueTrialIDs(i))));
end

end

function [dp_xaxis,dp_yaxis]=tryOtherDecodes(idx2_fr_success_cue,idx1_fr_success_cue,idx2_fr_success_uncue,idx1_fr_success_uncue,idx2_fr_failure_cue,idx1_fr_failure_cue,idx2_fr_failure_uncue,idx1_fr_failure_uncue,currSigns,successRange,failureRange)

idx2_fr_success_cue=idx2_fr_success_cue.*currSigns(1);
idx1_fr_success_cue=idx1_fr_success_cue.*currSigns(2);
idx2_fr_success_uncue=idx2_fr_success_uncue.*currSigns(3);
idx1_fr_success_uncue=idx1_fr_success_uncue.*currSigns(4);

idx2_fr_failure_cue=idx2_fr_failure_cue.*currSigns(1);
idx1_fr_failure_cue=idx1_fr_failure_cue.*currSigns(2);
idx2_fr_failure_uncue=idx2_fr_failure_uncue.*currSigns(3);
idx1_fr_failure_uncue=idx1_fr_failure_uncue.*currSigns(4);

% Decode "cued vs uncued" if passed in, e.g., cued_failure and
% uncued_failure
[n,x]=histcounts((idx2_fr_success_cue-idx1_fr_success_cue)-(idx2_fr_success_uncue-idx1_fr_success_uncue),-2000-40:80:2000+40);
[n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
[n,x]=histcounts((idx2_fr_failure_cue-idx1_fr_failure_cue)-(idx2_fr_failure_uncue-idx1_fr_failure_uncue),-2000-40:80:2000+40);
[n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');

% Decode "success vs failure"
[n,x]=histcounts((idx2_fr_success_cue+idx2_fr_success_uncue)-(idx1_fr_success_cue+idx1_fr_success_uncue),-2000-40:80:2000+40);
[n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color','k');
[n,x]=histcounts((idx2_fr_failure_cue+idx2_fr_failure_uncue)-(idx1_fr_failure_cue+idx1_fr_failure_uncue),-2000-40:80:2000+40);
[n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color','r');

% dprime
temp1=(idx2_fr_success_cue-idx1_fr_success_cue)-(idx2_fr_success_uncue-idx1_fr_success_uncue);
temp2=(idx2_fr_failure_cue-idx1_fr_failure_cue)-(idx2_fr_failure_uncue-idx1_fr_failure_uncue);
dp_yaxis=(nanmean(temp1)-nanmean(temp2))./sqrt(nanstd(temp1,[],1).^2+nanstd(temp2,[],1).^2); disp(dp_yaxis);
temp3=(idx2_fr_success_cue+idx2_fr_success_uncue)-(idx1_fr_success_cue+idx1_fr_success_uncue);
temp4=(idx2_fr_failure_cue+idx2_fr_failure_uncue)-(idx1_fr_failure_cue+idx1_fr_failure_uncue);
dp_xaxis=(nanmean(temp3)-nanmean(temp4))./sqrt(nanstd(temp3,[],1).^2+nanstd(temp4,[],1).^2); disp(dp_xaxis);

% scatter
figure(); s=scatter(temp3./length(successRange(1):successRange(2)),temp1./length(successRange(1):successRange(2)),60,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]);
figure(); s=scatter(temp4./length(failureRange(1):failureRange(2)),temp2./length(failureRange(1):failureRange(2)),60,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.4); hold on;
line([-5 5],[-5 5]); line([5 -5],[-5 5]);

% pause;
close all;

end

function [dp_all,dpfr_all,isSig_all,pval_all,p_success_unitbyunit_all,p_failure_unitbyunit_all,fr_success_unitbyunit_all,fr_failure_unitbyunit_all,timebinMeans,frsd_success_unitbyunit_all,frsd_failure_unitbyunit_all]=withinSessDprime(success_Response,failure_Response,timeBins)

for i=1:size(timeBins,1)
    [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=dprimeInTimeSessBySess(success_Response,failure_Response,timeBins(i,:));
    if i==1
        dp_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        dpfr_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        isSig_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        pval_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        p_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        p_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        fr_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        fr_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        frsd_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        frsd_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        timebinMeans=nan(size(timeBins,1),1);
    end
    dp_all(i,:,:)=dp;
    dpfr_all(i,:,:)=dpfr;
    if ~isempty(isSig)
        isSig_all(i,:,:)=isSig;
        pval_all(i,:,:)=pval;
    end
    p_success_unitbyunit_all(i,:,:)=p_success_unitbyunit;
    p_failure_unitbyunit_all(i,:,:)=p_failure_unitbyunit;
    fr_success_unitbyunit_all(i,:,:)=fr_success_unitbyunit;
    fr_failure_unitbyunit_all(i,:,:)=fr_failure_unitbyunit;
    frsd_success_unitbyunit_all(i,:,:)=frsd_success_unitbyunit;
    frsd_failure_unitbyunit_all(i,:,:)=frsd_failure_unitbyunit;
    timebinMeans(i)=mean(timeBins(i,:));
end

end

function [dp_all,dpfr_all,isSig_all,pval_all,p_success_unitbyunit_all,p_failure_unitbyunit_all,fr_success_unitbyunit_all,fr_failure_unitbyunit_all,timebinMeans]=dprimesOverTime(success_Response,failure_Response,timeBins)

for i=1:size(timeBins,1)
    [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit]=dprimeInTime(success_Response,failure_Response,timeBins(i,:));
    if i==1
        dp_all=nan(size(timeBins,1),length(dp));
        dpfr_all=nan(size(timeBins,1),length(dp));
        isSig_all=nan(size(timeBins,1),length(dp));
        pval_all=nan(size(timeBins,1),length(dp));
        p_success_unitbyunit_all=nan(size(timeBins,1),length(dp));
        p_failure_unitbyunit_all=nan(size(timeBins,1),length(dp));
        fr_success_unitbyunit_all=nan(size(timeBins,1),length(dp));
        fr_failure_unitbyunit_all=nan(size(timeBins,1),length(dp));
        timebinMeans=nan(size(timeBins,1),1);
    end
    dp_all(i,:)=dp;
    dpfr_all(i,:)=dpfr;
    if ~isempty(isSig)
        isSig_all(i,:)=isSig;
        pval_all(i,:)=pval;
    end
    p_success_unitbyunit_all(i,:)=p_success_unitbyunit;
    p_failure_unitbyunit_all(i,:)=p_failure_unitbyunit;
    fr_success_unitbyunit_all(i,:)=fr_success_unitbyunit;
    fr_failure_unitbyunit_all(i,:)=fr_failure_unitbyunit;
    timebinMeans(i)=mean(timeBins(i,:));
end

end

function [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=dprimeInTimeSessBySess(success_Response,failure_Response,timeWindow)

isSig=[];
pval=[];
nBoots=100;

% get probability of response in time window for each unit
[~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
temp=nanmean(success_Response.aligncomp_x,1);
alignSuccessTime=temp(alignPeakInd);
[~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
temp=nanmean(failure_Response.aligncomp_x,1);
alignFailureTime=temp(alignPeakInd);
[~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
successRange=[startAt endAt];
[~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
failureRange=[startAt endAt];

unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
fromWhichSess_success=success_Response.fromWhichSess_forTrials;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
fromWhichSess_success=fromWhichSess_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
fromWhichSess_failure=failure_Response.fromWhichSess_forTrials;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
fromWhichSess_failure=fromWhichSess_failure(failure_Response.isEventInThisTrial==1);

units=unique(success_Response.fromWhichUnit);
% [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure);
% [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
% maxdp=3;
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
% dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;

[fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure);
dpfr=dprime_from_FR(fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit);

end

function [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit]=dprimeInTime(success_Response,failure_Response,timeWindow)

isSig=[];
pval=[];
nBoots=100;

% get probability of response in time window for each unit
[~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
temp=nanmean(success_Response.aligncomp_x,1);
alignSuccessTime=temp(alignPeakInd);
[~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
temp=nanmean(failure_Response.aligncomp_x,1);
alignFailureTime=temp(alignPeakInd);
[~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
successRange=[startAt endAt];
[~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
failureRange=[startAt endAt];

unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);

units=unique(success_Response.fromWhichUnit);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);

% [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
% maxdp=3;
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
% dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;

[fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
dpfr=dprime_from_FR(fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit);

end

function dp=dprime_from_FR(data1_me,data2_me,data1_sd,data2_sd)

% RMS sd discriminability index
% suboptimal but simple
dp=(data1_me-data2_me)./sqrt(data1_sd.^2+data2_sd.^2);

end

function [isSig,pval,real_success_per_unit,real_failure_per_unit]=getSignificantUnits_Differences(units,success_Response,failure_Response,nBoots,timeWindow1,timeWindow2,vsOrJustTimeWindow2)

boot_success_per_unit=nan(length(units),nBoots);
boot_failure_per_unit=nan(length(units),nBoots);
% get real prob values
timeWindow=timeWindow1; %[-2 0];
[~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
temp=nanmean(success_Response.aligncomp_x,1);
alignSuccessTime=temp(alignPeakInd);
[~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
temp=nanmean(failure_Response.aligncomp_x,1);
alignFailureTime=temp(alignPeakInd);
[~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
successRange=[startAt endAt];
[~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
failureRange=[startAt endAt];
unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
units=unique(success_Response.fromWhichUnit);
[p_success_unitbyunit_BEFORE,p_failure_unitbyunit_BEFORE]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
timeWindow=timeWindow2; %[1 3];
[~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
temp=nanmean(success_Response.aligncomp_x,1);
alignSuccessTime=temp(alignPeakInd);
[~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
temp=nanmean(failure_Response.aligncomp_x,1);
alignFailureTime=temp(alignPeakInd);
[~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
successRange=[startAt endAt];
[~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
failureRange=[startAt endAt];
unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
units=unique(success_Response.fromWhichUnit);
[p_success_unitbyunit_AFTER,p_failure_unitbyunit_AFTER]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
real_success_per_unit=p_success_unitbyunit_BEFORE-p_success_unitbyunit_AFTER;
real_failure_per_unit=p_failure_unitbyunit_BEFORE-p_failure_unitbyunit_AFTER;
for j=1:nBoots
    % shuffle and bootstrap
    timeWindow=[-2 0];
    [~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
    temp=nanmean(success_Response.aligncomp_x,1);
    alignSuccessTime=temp(alignPeakInd);
    [~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
    temp=nanmean(failure_Response.aligncomp_x,1);
    alignFailureTime=temp(alignPeakInd);
    [~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
    successRange=[startAt endAt];
    [~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
    failureRange=[startAt endAt];
    unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
    fromWhichUnit_success=success_Response.fromWhichUnit;
    unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
    fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
    unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
    fromWhichUnit_failure=failure_Response.fromWhichUnit;
    unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
    fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
    % shuffle
    randie=cell(length(units),1);
    for i=1:length(units)
        currunit=units(i);
        resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
        randie{i}=randperm(length(resp));
        shuffle_resp=resp(randie{i});
        howmanysuccesses=nansum(fromWhichUnit_success==currunit);
        unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
        unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
    end
    % get prob
    [p_success_unitbyunit_BEFORE,p_failure_unitbyunit_BEFORE]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    % new time window
    timeWindow=[1 3];
    [~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
    successRange=[startAt endAt];
    [~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
    failureRange=[startAt endAt];
    unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
    unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
    unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
    unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
    % shuffle
    for i=1:length(units)
        currunit=units(i);
        resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
        shuffle_resp=resp(randie{i});
        howmanysuccesses=nansum(fromWhichUnit_success==currunit);
        unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
        unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
    end
    % get prob
    [p_success_unitbyunit_AFTER,p_failure_unitbyunit_AFTER]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    switch vsOrJustTimeWindow2
        case 'versus'
            boot_success_per_unit(:,j)=p_success_unitbyunit_BEFORE-p_success_unitbyunit_AFTER;
            boot_failure_per_unit(:,j)=p_failure_unitbyunit_BEFORE-p_failure_unitbyunit_AFTER;
        case 'justTimeWindow2'
            boot_success_per_unit(:,j)=p_success_unitbyunit_AFTER;
            boot_failure_per_unit(:,j)=p_failure_unitbyunit_AFTER;
    end
end
% get significance for each unit
% is significant if real difference between success and failure is greater
% than 97.5th percentile of shuffled differences OR less than 2.5th
% percentile of shuffled differences
greaterthanwhichprctile=nan(length(units),1);
lessthanwhichprctile=nan(length(units),1);
pval=nan(length(units),1);
getptiles=1:100;
for i=1:length(units)
    ptiles=prctile(boot_success_per_unit(i,:)-boot_failure_per_unit(i,:),getptiles);
    temp=getptiles(find(ptiles<(real_success_per_unit(i)-real_failure_per_unit(i)),1,'last'));
    if ~isempty(temp)
        greaterthanwhichprctile(i)=temp;
    else
        greaterthanwhichprctile(i)=0;
    end
    temp=getptiles(find(ptiles>(real_success_per_unit(i)-real_failure_per_unit(i)),1,'first'));
    if ~isempty(temp)
        lessthanwhichprctile(i)=temp;
    else
        lessthanwhichprctile(i)=100;
    end
    if greaterthanwhichprctile(i)>50
        if greaterthanwhichprctile(i)==100
            pval(i)=1/100;
            continue
        end
        % above median
        pval(i)=((100-(greaterthanwhichprctile(i)+1))/100);
    else
        if lessthanwhichprctile(i)==0
            pval(i)=1/100;
            continue
        end
        % below median
        pval(i)=(((lessthanwhichprctile(i)-1)-0)/100);
    end
end
isSig=pval<0.025;

end

function [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots)

boot_success_per_unit=nan(length(units),nBoots);
boot_failure_per_unit=nan(length(units),nBoots);
% get real prob values
[real_success_per_unit,real_failure_per_unit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
for j=1:nBoots
    % shuffle and bootstrap
    [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    [boot_success_per_unit(:,j),boot_failure_per_unit(:,j)]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
end
% get significance for each unit
% is significant if real difference between success and failure is greater
% than 97.5th percentile of shuffled differences OR less than 2.5th
% percentile of shuffled differences
greaterthanwhichprctile=nan(length(units),1);
lessthanwhichprctile=nan(length(units),1);
pval=nan(length(units),1);
getptiles=1:100;
for i=1:length(units)
    ptiles=prctile(boot_success_per_unit(i,:)-boot_failure_per_unit(i,:),getptiles);
    temp=getptiles(find(ptiles<(real_success_per_unit(i)-real_failure_per_unit(i)),1,'last'));
    if ~isempty(temp)
        greaterthanwhichprctile(i)=temp;
    else
        greaterthanwhichprctile(i)=0;
    end
    temp=getptiles(find(ptiles>(real_success_per_unit(i)-real_failure_per_unit(i)),1,'first'));
    if ~isempty(temp)
        lessthanwhichprctile(i)=temp;
    else
        lessthanwhichprctile(i)=100;
    end
    if greaterthanwhichprctile(i)>50
        if greaterthanwhichprctile(i)==100
            pval(i)=1/100;
            continue
        end
        % above median
        pval(i)=((100-(greaterthanwhichprctile(i)+1))/100);
    else
        if lessthanwhichprctile(i)==0
            pval(i)=1/100;
            continue
        end
        % below median
        pval(i)=(((lessthanwhichprctile(i)-1)-0)/100);
    end
end
isSig=pval<0.025;

end

function [out2point5,out97point5,outbinmids]=bootstrapPerctiles(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots)

p_bins=0:0.01:1.001;
p_bins(end)=1.001;
perc2point5=nan(length(p_bins)-1,nBoots);
perc97point5=nan(length(p_bins)-1,nBoots);
binmids=nan(length(p_bins)-1,nBoots);
for j=1:nBoots
    [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    % get 95% CI for bootstrap (i.e., shuffle)
    for i=1:length(p_bins)-1
        successbin=[p_bins(i) p_bins(i+1)];
        binmids(i,j)=mean(successbin);
        P=prctile(p_failure_unitbyunit(p_success_unitbyunit>=successbin(1) & p_success_unitbyunit<successbin(2)),[2.5 97.5]);
        perc2point5(i,j)=P(1);
        perc97point5(i,j)=P(2);
    end
end
out2point5=mean(perc2point5,2,'omitnan');
out97point5=mean(perc97point5,2,'omitnan');
outbinmids=mean(binmids,2,'omitnan');

end

function [fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

fr_success_unitbyunit=nan(length(units),1);
fr_failure_unitbyunit=nan(length(units),1);
frsd_success_unitbyunit=nan(length(units),1);
frsd_failure_unitbyunit=nan(length(units),1);
disp([num2str(length(units)) ' units']);
for i=1:length(units)
    fr_success_unitbyunit(i)=nanmean(unitfr_success(fromWhichUnit_success==units(i)));
    fr_failure_unitbyunit(i)=nanmean(unitfr_failure(fromWhichUnit_failure==units(i)));
    frsd_success_unitbyunit(i)=nanstd(unitfr_success(fromWhichUnit_success==units(i)));
    frsd_failure_unitbyunit(i)=nanstd(unitfr_failure(fromWhichUnit_failure==units(i)));
end

end

function [fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure)

uSess=unique([fromWhichSess_success; fromWhichSess_failure]);
fr_success_unitbyunit=nan(length(units),1);
fr_failure_unitbyunit=nan(length(units),1);
frsd_success_unitbyunit=nan(length(units),1);
frsd_failure_unitbyunit=nan(length(units),1);
for j=1:length(uSess)
    currSess=uSess(j);
    for i=1:length(units)
        fr_success_unitbyunit(i,j)=nanmean(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess));
        fr_failure_unitbyunit(i,j)=nanmean(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess));
        frsd_success_unitbyunit(i,j)=nanstd(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess));
        frsd_failure_unitbyunit(i,j)=nanstd(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess));
    end
end

end

function [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure)

uSess=unique([fromWhichSess_success; fromWhichSess_failure]);
p_success_unitbyunit=nan(length(units),length(uSess));
p_failure_unitbyunit=nan(length(units),length(uSess));
for j=1:length(uSess)
    currSess=uSess(j);
    for i=1:length(units)
        p_success_unitbyunit(i,j)=nansum(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess)>0.5)./nansum(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess);
        p_failure_unitbyunit(i,j)=nansum(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess)>0.5)./nansum(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess);
    end
end

end

function [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

p_success_unitbyunit=nan(length(units),1);
p_failure_unitbyunit=nan(length(units),1);
disp([num2str(length(units)) ' units']);
for i=1:length(units)
    p_success_unitbyunit(i)=nansum(unitfr_success(fromWhichUnit_success==units(i))>0.5)./nansum(fromWhichUnit_success==units(i));
    p_failure_unitbyunit(i)=nansum(unitfr_failure(fromWhichUnit_failure==units(i))>0.5)./nansum(fromWhichUnit_failure==units(i));
end

end

function [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

for i=1:length(units)
    currunit=units(i);
    resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
    shuffle_resp=resp(randperm(length(resp)));
    howmanysuccesses=nansum(fromWhichUnit_success==currunit);
    unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
    unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
end

end



