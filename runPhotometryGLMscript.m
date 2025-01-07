function runPhotometryGLMscript(data_loc_array)

for i=125:253
% for i=1:253
    % Get photometry signals, pDMSt and NAcc
    switch data_loc_array{i,14}
        case 'green'
            nowproc='red';
        case 'red'
            nowproc='green';
        otherwise
            continue
    end
    if ~exist(data_loc_array{i,6},'dir')
        continue
    end
    if ~exist([data_loc_array{i,15} 'NAcc'],'dir')
        mkdir([data_loc_array{i,15} 'NAcc']);
        a=load([data_loc_array{i,6} sep 'behavior_tbt.mat']);
        behavior_tbt=a.behavior_tbt;
        a=load([data_loc_array{i,6} sep 'photometry_tbt.mat']);
        photometry_tbt=a.photometry_tbt;
        saveBehaviorAlignmentsPhotometry(photometry_tbt,behavior_tbt,[data_loc_array{i,15} 'NAcc'],'',nowproc);
    end 
    
    % Check that alignComp will not give trouble
    % pDMSt
    a=load([data_loc_array{i,15} sep 'cue' sep 'ch__cueAligned']);
    alignComp=a.alignComp;
    dataout=a.dataout;
    f=find(~isnan(alignComp.y),1,'last');
    if alignComp.y(f)>0
        f2=find(fliplr(alignComp.y(1:f))<=0,1,'first');
        alignComp.y(f-f2:end)=0;
        save([data_loc_array{i,15} sep 'cue' sep 'ch__cueAligned'],'alignComp','dataout');
    end
    % NAcc
    a=load([data_loc_array{i,15} 'NAcc' sep 'cue' sep 'ch__cueAligned']);
    alignComp=a.alignComp;
    dataout=a.dataout;
    f=find(~isnan(alignComp.y),1,'last');
    if alignComp.y(f)>0
        f2=find(fliplr(alignComp.y(1:f))<=0,1,'first');
        alignComp.y(f-f2:end)=0;
        save([data_loc_array{i,15} 'NAcc' sep 'cue' sep 'ch__cueAligned.mat'],'alignComp','dataout');
    end

    % Set up pDMSt GLM
    if ~exist([data_loc_array{i,15} sep 'forPhotoglm'],'dir')
        GLM_forPhoto_analysis(i,data_loc_array,10,'','');
        % Add reach behavior event
        a=load([data_loc_array{i,15} sep 'forPhotoglm' sep 'behEvents.mat']);
        behEvents=a.behEvents;
        behEvents=[behEvents(1:9,:); behEvents(4,:)+behEvents(5,:)+behEvents(6,:)>0.5; behEvents(10,:)];
        save([data_loc_array{i,15} sep 'forPhotoglm' sep 'behEvents.mat'],'behEvents');
    end
    
    % Set up NAcc GLM
    if ~exist([data_loc_array{i,15} sep 'forPhotoglmNAcc'],'dir')
        GLM_forPhoto_analysis(i,data_loc_array,10,'NAcc','NAcc');
        % Add reach behavior event
        a=load([data_loc_array{i,15} sep 'forPhotoglmNAcc' sep 'behEvents.mat']);
        behEvents=a.behEvents;
        behEvents=[behEvents(1:9,:); behEvents(4,:)+behEvents(5,:)+behEvents(6,:)>0.5; behEvents(10,:)];
        save([data_loc_array{i,15} sep 'forPhotoglmNAcc' sep 'behEvents.mat'],'behEvents');
    end

    close all;
end

end