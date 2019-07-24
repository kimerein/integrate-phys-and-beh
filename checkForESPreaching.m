function checkForESPreaching(expt_dir)

overwriteExisting=false;

ls=dir(expt_dir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if ~isempty(regexp(thisname,'processed_data','ONCE')) && thisisdir==1
        if exist([expt_dir '\' thisname '\preemptCue.mat'],'file')
            if overwriteExisting==false
                continue
            end
        end         
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
        close all;
    end
end