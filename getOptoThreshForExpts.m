function getOptoThreshForExpts(expt_dir)

ls=dir(expt_dir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if ~isempty(regexp(thisname,'processed_data')) && thisisdir==1
        a=load([expt_dir '\' thisname '\final_aligned_data.mat']);
        alignment=a.alignment;
        figure(); 
        plot(alignment.optoZone);
        th=input('Opto thresh for this? ');
        optoThresh=th;
        save([expt_dir '\' thisname '\optoThresh.mat'],'optoThresh');
    end
end