function outSettings=reachExpt_analysis_settings(varargin)
% first argument is whether to display settings, second argument overwrites
% any default settings

% persistent settings
settings=[];

displ=false;
if ~isempty(varargin)
    if strcmp(varargin{1},'display settings')
        displ=true;
    end
end

if isempty(settings)
    % default settings
    settings.nameOfCue='cueZone_onVoff'; % name of cue in tbt data structure
    settings.lowThresh=0.05; % threshold for detecting events in trial-by-trial behavior data structure
    settings.durationOfWheelTurn=1.029; % duration of pellet presenter wheel turn, in seconds
    settings.wheelStopsThisManySecsBeforeCue=0.1; % time window between pellet presenter wheel stopping and cue onset, in seconds
    settings.cueDuration=0.25; % duration of cue, in seconds
    settings.maxTrialDuration=35; % maximum duration of a trial, in seconds
    settings.timeSlop=0.1; % max possible error in timing due to frame rate in movie and/or alignment approach, in seconds
    settings.reachAfterCueWindow_start=0; % define start of time window in which to say mouse "reached after cue", seconds from cue onset
    settings.reachAfterCueWindow_end=0+1.5; % define end of time window in which to say mouse "reached after cue", seconds from cue onset
    settings.preCueWindow_start=0; % define start of time window from trial onset, in seconds
    settings.preCueWindow_end=1.5; % define end of time window from trial onset, in seconds
    settings.check_for_human=1; % 1 if want to exclude data without humanChecked.txt
    settings.discardPreemptive=1; % 1 if want to discard data sets with preemptive reaching
    settings.doRealign=1; % in tbt, cue can be any time point when cue is on; 1 if want to realign data based on cue starts (i.e., align all cue onsets)
    settings.useOptoZone=0; % 1 if want to use manually defined optoZone in video instead of Arduino-based optoOn
    settings.maxDelayUntilOpto=13.25; % in seconds, max time from trial onset until opto turns on
    settings.isOrchestra=1; % will suppress figures if running on server
    settings.putCueAtInd=94; % if, when doing cue realignment, want to put the cue onset (i.e., cueZone_onVoff max) at a specific index; only used if doRealign==1
    % if putCueAtInd is empty, will just use first trial's cue position as
    % the default
    settings.tryForFiles={'optoOnHere','nth_session','optoThresh','preemptCue','dateFromTextFile'}; % look for these files in each directory
end

if ~isempty(varargin) && length(varargin)>1
    if isstruct(varargin{2})
        set=varargin{2};
        % passing in some settings
        f=fieldnames(set);
        for i=1:length(f)
            settings.(f{i})=set.(f{i});
            disp([f{i} ' is ' ]);
            disp(settings.(f{i}));
        end
    end
end

if displ==true
    % Order in fldname must match order in prompt
    fldname={'nameOfCue','lowThresh','durationOfWheelTurn','wheelStopsThisManySecsBeforeCue','cueDuration','maxTrialDuration','timeSlop','reachAfterCueWindow_start','reachAfterCueWindow_end',...
             'preCueWindow_start','preCueWindow_end','check_for_human','discardPreemptive','doRealign','useOptoZone','maxDelayUntilOpto','putCueAtInd'};
    prompt={'name of cue:','tbt thresh:','duration of wheel turn (s):','time between wheel stop and cue (s):','cue duration (s):','max trial duration (s):',...
            'time uncertainty (s):','"reached after cue" window s from cue on START:','"reached after cue" window s from cue on END:',...
            'pre-cue window from trial onset START (s):','pre-cue window from trial onset END (s):','1 if need human check to include data:','1 if discard days with preemptive reaching:','1 if realign all cues to cue onset:','1 if want to use movie opto zone instead of opto from Arduino:','max delay til opto (s):',...
            'if doRealign is true, will put cue onsets at this index from trial start (empty for default):'};
    dlgtitle='Check or modify experiment-specific settings';
    dims=[1 70];
    definput={settings.(fldname{1}),num2str(settings.(fldname{2})),num2str(settings.(fldname{3})),num2str(settings.(fldname{4})),num2str(settings.(fldname{5})),num2str(settings.(fldname{6})),...
              num2str(settings.(fldname{7})),num2str(settings.(fldname{8})),num2str(settings.(fldname{9})),...
              num2str(settings.(fldname{10})),num2str(settings.(fldname{11})),num2str(settings.(fldname{12})),num2str(settings.(fldname{13})),num2str(settings.(fldname{14})),num2str(settings.(fldname{15})),num2str(settings.(fldname{16})),num2str(settings.(fldname{17}))};
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

if displ==true
    strtoshow=[];
    for i=1:length(settings.tryForFiles)
        strtoshow=[strtoshow ' ' settings.tryForFiles{i}];
    end
    answer=questdlg(['Will look for these files in each data directory: ' strtoshow]);
    switch answer
        case 'Yes'
        case 'No'
            disp('Re-run with correct settings. Can modify in reachExpt_analysis_settings.m.');
    end
end

outSettings=settings;