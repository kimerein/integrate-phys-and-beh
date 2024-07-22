function getOptoThreshForExpts(expt_dir,optoFromArduino,arduinoOptoForWhichMice)

overwriteExistingThresh=true;
check_for_human=1;

currOptoFromArduino=optoFromArduino;

if ~iscell(expt_dir)
    ls=dir(expt_dir);
else
    ls=expt_dir;
end
back_expt_dir=expt_dir;
for i=1:length(ls)
    currOptoFromArduino=optoFromArduino;
    if ~iscell(back_expt_dir)
        thisname=ls(i).name;
        thisisdir=ls(i).isdir;
    else
        currdir=ls{i};
        temp=regexp(currdir,'\');
        thisname=currdir(temp(end)+1:end);
        thisisdir=isempty(regexp(thisname,'\.','ONCE'));
        rar=regexp(currdir,'\');
        expt_dir=currdir(1:rar(end)-1);
    end
    if ~isempty(regexp(thisname,'processed_data','once')) && thisisdir==1
        if exist([expt_dir '\' thisname '\optoThresh.mat'],'file') && exist([expt_dir '\' thisname '\optoOnHere.mat'],'file')
            if overwriteExistingThresh==false
                continue
            end
        end    
        if check_for_human==1
            if exist([expt_dir '\' thisname '\humanchecked.txt'], 'file')==2 || exist([expt_dir '\' thisname '\humanchecked_afterResampleFix.txt'], 'file')==2
            else
%                 disp(['Not including ' expt_dir '\' thisname]);
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
        th=input([expt_dir '\' thisname ' Opto thresh for this? (Enter -10 if no opto here.) ']);
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