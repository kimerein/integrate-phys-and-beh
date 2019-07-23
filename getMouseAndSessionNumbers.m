function getMouseAndSessionNumbers(exptDir,mouseDatabase)

% get mouse ID database
if isempty(mouseDatabase)
    mouseDatabase=makeMouseDatabase;
elseif isstr(mouseDatabase)
    a=load(mouseDatabase);
    mouseDatabase=a.mouseDatabase;
else
    % passed in mouseDatabase
end
    
% this should be a structure with the following fields
% 'mouse_name'          mouse's name
% 'mouse_id'            mouse's ID number
% 'description'         description of mouse genotype, injections and
%                       experimental goals
% 'trainers'            who trained this mouse
% 'cue start date'      first date of cue training on rig
% 'sac date'            date mouse was sac'ced
% 'data_directories'    a list of all directories containing behavior data
%                       related to this mouse
% 'low_speed_videos'    a list of all low-speed videos (.AVI) with behavior
%                       data for this mouse

% use this list of data_directories and low_speed_videos associated with
% each mouse to figure out mouse's id for analyzed data

ls=dir(expt_dir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if ~isempty(regexp(thisname,'processed_data')) && thisisdir==1
        a=load([expt_dir '\' thisname '\tbt.mat']);
        tbt=a.tbt;
        interptimes=0:0.035:17;
        figure(); 
        plot(interptimes,nanmean(tbt.reachStarts,1),'Color','k'); 
        hold on; 
        plot(interptimes,nanmean(tbt.cueZone_onVoff,1),'Color','b');
        th=questdlg('Pre-emptive reaching?','Does mouse reach consistently before cue?');
        if isempty(th)
            break
        elseif strcmp(th,'Yes')
            preemptCue=true;
        elseif strcmp(th,'No')
            preemptCue=false;
        elseif strcmp(th,'Cancel')
            break
        end
        save([expt_dir '\' thisname '\preemptCue.mat'],'preemptCue');
    end
end