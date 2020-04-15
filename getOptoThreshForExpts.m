function getOptoThreshForExpts(expt_dir,optoFromArduino,arduinoOptoForWhichMice)

overwriteExistingThresh=false;

currOptoFromArduino=optoFromArduino;

if ~iscell(expt_dir)
    ls=dir(expt_dir);
else
    ls=expt_dir;
end
for i=1:length(ls)
    currOptoFromArduino=optoFromArduino;
    if ~iscell(expt_dir)
        thisname=ls(i).name;
        thisisdir=ls(i).isdir;
    else
        currdir=ls{i};
        temp=regexp(currdir,'\');
        thisname=currdir(temp(end)+1:end);
        thisisdir=isempty(regexp(thisname,'\.','ONCE'));
    end
    if ~isempty(regexp(thisname,'processed_data','once')) && thisisdir==1
        if exist([expt_dir '\' thisname '\optoThresh.mat'],'file') && exist([expt_dir '\' thisname '\optoOnHere.mat'],'file')
            if overwriteExistingThresh==false
                continue
            end
        end                
        a=load([expt_dir '\' thisname '\final_aligned_data.mat']);
        alignment=a.alignment;
        if ~isempty(arduinoOptoForWhichMice)
            a=load([expt_dir '\' thisname '\mouse_id.mat']);
            mouse_id=a.mouse_id;
            if ismember(mouse_id,arduinoOptoForWhichMice)
                currOptoFromArduino=true;
            end
        end
        if currOptoFromArduino==true
            a=load([expt_dir '\' thisname '\tbt.mat']);
            tbt=a.tbt;
        end
        figure(); 
        if currOptoFromArduino==false
            plot(alignment.optoZone);
        else
            plot(alignment.optoOn);
        end
        th=input([thisname ' Opto thresh for this? (Enter -10 if no opto here.) ']);
        optoThresh=th;
        hold on;
        if currOptoFromArduino==false
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
        if currOptoFromArduino==true
            alignment.optoZone=alignment.optoOn;
            save([expt_dir '\' thisname '\final_aligned_data.mat'],'alignment')
            tbt.optoZone=tbt.optoOn;
            save([expt_dir '\' thisname '\tbt.mat'],'tbt')
        end
        pause;
        close all;
    end
end