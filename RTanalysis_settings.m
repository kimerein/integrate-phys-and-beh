function outSettings=RTanalysis_settings(varargin)

persistent settings

displ=false;
if ~isempty(varargin)
    if strcmp(varargin{1},'display settings');
        displ=true;
    end
    if strcmp(varargin{2},'clear')
        settings=[];
    end
end

if isempty(settings)
    % default settings
    settings.returnPreCueRTs=0; % also return "reaction times" given that time moves backwards from cue
    settings.divideByBaseRate=0; % note that reaction times pick first reach after cue -- higher baseline reaching rate will mean more RTs near 0 ... do you want to normalize reaction time distribution by base rate? if yes, then 1
    settings.subtractBaselineReaching=0; % if is 1, will subtract off the baseline reaching rate, effectively zeroing out time bins where distribution falls at or below baseline
    settings.throwOutAllAfter1stBaseline=0; % if 1, will throw out all reaction times after reaction time distribution first returns to baseline reach rate
    settings.noRTlessThan=0; % no reaction times can be less than this (in seconds); if 0, this is the classic definition of reaction time
    settings.excludePawOnWheelDuringCue=0; % 1 if want to exclude trials where mouse reached before/during cue, else 0
    settings.longRT_ifNoReach=1; % if mouse does not reach in this trial, fill in reaction time as longer than trial length IF this is set to 1, else throw out trial
    settings.reverse_trial_order=0; % 1 if want to flip trial order (fake time goes backwards)
    settings.longITIthresh=20; % in sec, trials lasting longer than this will be considered long ITI trials
    settings.multipleOptoTimes=0; % if there were multiple different opto timings
end

if displ==true
    % Order in fldname must match order in prompt
    fldname={'returnPreCueRTs','divideByBaseRate','subtractBaselineReaching','throwOutAllAfter1stBaseline','noRTlessThan','excludePawOnWheelDuringCue','longRT_ifNoReach','reverse_trial_order','longITIthresh','multipleOptoTimes'};
    prompt={'1 if want to return negative reaction times (before cue):','1 if want to normalize RT distribution by base reaching rate:',...
            '1 if want to subtract off baseline reaching rate:','1 if want to throw out all RTs after reach rate returns to baseline rate:',...
            'no reaction time can be less than this (sec), wrt cue onset:','1 if want to exclude trials in which mouse paw was already on wheel at cue onset:',...
            '1 if want to fill in RT longer than trial length if mouse did not reach in trial:','1 if want to reverse trial order (for control if time goes backwards):'...
            'threshold separating long ITI trials from short ITI trials (sec):','1 if there were multiple different opto timings:'};
    dlgtitle='Check or modify RT analysis-specific settings';
    dims=[1 35];
    definput={num2str(settings.(fldname{1})),num2str(settings.(fldname{2})),num2str(settings.(fldname{3})),num2str(settings.(fldname{4})),num2str(settings.(fldname{5})),num2str(settings.(fldname{6})),num2str(settings.(fldname{7})),num2str(settings.(fldname{8})),num2str(settings.(fldname{9})),num2str(settings.(fldname{10}))};
    answer=inputdlg(prompt,dlgtitle,dims,definput);
    for i=1:length(answer)
        if ischar(settings.(fldname{i}))
            if ~strcmp(settings.(fldname{i}),answer{i})
                settings.(fldname{i})=answer{i};
            end
        else
            % all others are numbers
            if str2num(answer{i})~=settings.(fldname{i})
                settings.(fldname{i})=str2num(answer{i});
            end
        end
    end
end

outSettings=settings;