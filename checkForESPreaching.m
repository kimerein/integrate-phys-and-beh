function checkForESPreaching(expt_dir,within_dir)

overwriteExisting=false;
isExternalCue=false;

if ~iscell(expt_dir)
    ls=dir(expt_dir);
else
    ls=expt_dir;
end
if ~iscell(expt_dir)
    doingdir=true;
else
    doingdir=false;
end
for i=1:length(ls)
    if doingdir==true
        thisname=ls(i).name;
        thisisdir=ls(i).isdir;
    else
        currdir=ls{i};
        temp=regexp(currdir,'\');
        thisname=currdir(temp(end)+1:end);
        thisisdir=isempty(regexp(thisname,'\.','ONCE'));
        expt_dir=currdir(1:temp(end));
    end
    if ~isempty(regexp(thisname,'processed_data','ONCE')) && thisisdir==1
        if exist([expt_dir '\' thisname '\preemptCue.mat'],'file')
            if overwriteExisting==false
                continue
            end
        end         
        a=load([expt_dir '\' thisname '\tbt.mat']);
        tbt=a.tbt;
        interptimes=0:0.035:17;
        if isExternalCue==true
            cue_tbt=tbt;
            tbt=realignToDistractor(tbt,[],[],false);
            wasDistract=any(tbt.movie_distractor(:,92:end)>0.05,2);
            % only consider pre-cue period, blind viewer to post-cue period
            figure(); plot(nanmean(cue_tbt.all_reachBatch(:,1:85),1),'Color','r'); hold on;
            plot(nanmean(tbt.all_reachBatch(wasDistract==1,1:85),1),'Color','b');
        else
            figure();
            plot(interptimes,nanmean(tbt.reachStarts,1),'Color','k');
            hold on;
            plot(interptimes,nanmean(tbt.cueZone_onVoff,1),'Color','b');
        end
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