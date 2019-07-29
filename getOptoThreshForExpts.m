function getOptoThreshForExpts(expt_dir)

overwriteExistingThresh=false;

ls=dir(expt_dir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if ~isempty(regexp(thisname,'processed_data','once')) && thisisdir==1
        if exist([expt_dir '\' thisname '\optoThresh.mat'],'file') && exist([expt_dir '\' thisname '\optoOnHere.mat'],'file')
            if overwriteExistingThresh==false
                continue
            end
        end                
        a=load([expt_dir '\' thisname '\final_aligned_data.mat']);
        alignment=a.alignment;
        figure(); 
        plot(alignment.optoZone);
        th=input('Opto thresh for this? (Enter -10 if no opto here.) ');
        optoThresh=th;
        hold on;
        line([0 length(alignment.optoZone)],[optoThresh optoThresh],'Color','r');
        if optoThresh==-10
            optoOnHere=0;
        else
            optoOnHere=1;
        end
        save([expt_dir '\' thisname '\optoThresh.mat'],'optoThresh');
        save([expt_dir '\' thisname '\optoOnHere.mat'],'optoOnHere');
        pause;
        close all;
    end
end