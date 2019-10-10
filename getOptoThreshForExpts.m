function getOptoThreshForExpts(expt_dir,optoFromArduino)

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
        if optoFromArduino==true
            a=load([expt_dir '\' thisname '\tbt.mat']);
            tbt=a.tbt;
        end
        figure(); 
        if optoFromArduino==false
            plot(alignment.optoZone);
        else
            plot(alignment.optoOn);
        end
        th=input('Opto thresh for this? (Enter -10 if no opto here.) ');
        optoThresh=th;
        hold on;
        if optoFromArduino==false
            line([0 length(alignment.optoZone)],[optoThresh optoThresh],'Color','r');
        else
            line([0 length(alignment.optoOn)],[optoThresh optoThresh],'Color','r');
        end
        if optoThresh==-10
            optoOnHere=0;
        else
            optoOnHere=1;
        end
        save([expt_dir '\' thisname '\optoThresh.mat'],'optoThresh');
        save([expt_dir '\' thisname '\optoOnHere.mat'],'optoOnHere');
        if optoFromArduino==true
            alignment.optoZone=alignment.optoOn;
            save([expt_dir '\' thisname '\final_aligned_data.mat'],'alignment')
            tbt.optoZone=tbt.optoOn;
            save([expt_dir '\' thisname '\tbt.mat'],'tbt')
        end
        pause;
        close all;
    end
end