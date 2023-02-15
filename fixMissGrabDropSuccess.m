function fixMissGrabDropSuccess(datadir)

ls_base=dir(datadir);
for mousen=3:length(ls_base)
    if ~isfolder([ls_base(mousen).folder sep ls_base(mousen).name])
        continue
    end
    ls=dir([ls_base(mousen).folder sep ls_base(mousen).name]);
    mousename=ls_base(mousen).name;
    disp(['Processing ' mousename]);
    for i=3:length(ls)
        if isfolder([ls(i).folder sep ls(i).name])
            proc_data_folder=[ls(i).folder sep ls(i).name];
            f_pr=regexp(proc_data_folder,'_processed_data');
            fslash=regexp(proc_data_folder,'\');
            aviName=proc_data_folder(fslash(end)+1:f_pr-1);
            datestr=getDatestr([ls(i).folder sep ls(i).name]);
            placeForO2data=getPlaceOfO2data(mousename,aviName,datestr);
            recodeMissVGrab(proc_data_folder,placeForO2data,[]);
            fixDropVSuccess(proc_data_folder,placeForO2data,[]);
        end
    end
end
 
end

function getDatestr(foldername)

lsfilesinfolder=dir(foldername);
for i=3:length(lsfilesinfolder)
    if ~isempty(regexp(lsfilesinfolder(i).name,'.txt', 'once'))
        if isfile(lsfilesinfolder(i).name)
            if ~isempty(regexp(lsfilesinfolder(i).name,'2016', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2017', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2018', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2019', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2020', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2021', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2022', 'once')) || ...
               ~isempty(regexp(lsfilesinfolder(i).name,'2023', 'once'))
                    
            end
        end

    end
end

end

function placeForO2data=getPlaceOfO2data(mousename,aviName,datestr)

switch mousename
    case 'mitch_right'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
    case 'mitch_mismatch'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
    case 'mitch_thin'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
    case 'mitch_wide'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
    case 'Oct_B'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Oct_D'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Sep_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Feb_1'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Feb_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Nov_grey'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
end

end