function fixMissGrabDropSuccess(datadir)

usePrevChewingThresh=true;

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
%             try
                proc_data_folder=[ls(i).folder sep ls(i).name];
                f_pr=regexp(proc_data_folder,'_processed_data');
                fslash=regexp(proc_data_folder,'\');
                aviName=proc_data_folder(fslash(end)+1:f_pr-1);
                datestr=getDatestr([ls(i).folder sep ls(i).name]);
                placeForO2data=getPlaceOfO2data(mousename,aviName,datestr);
                if ~alreadyDidMissVGrab([ls(i).folder sep ls(i).name])
                    recodeMissVGrab(proc_data_folder,placeForO2data,[],usePrevChewingThresh);
                    % Note that recodeMissVGrab calls fixDropVSuccess at end
                    close all;
                else
                    disp(['Already did miss v grab for ' ls(i).folder sep ls(i).name ' so skipping']);
                    if ~alreadyDidSuccessVDrop([ls(i).folder sep ls(i).name])
                        fixDropVSuccess(proc_data_folder,placeForO2data,[],usePrevChewingThresh);
                        close all;
                    else
                        disp(['Already did success v drop for ' ls(i).folder sep ls(i).name ' so skipping']);
                    end
                end
%             catch
%                 disp(['Error in ' ls(i).folder sep ls(i).name]);
%             end
        end
    end
end
 
end

function didit=alreadyDidMissVGrab(foldername)

lsfilesinfolder=dir(foldername);
didit=false;
for i=3:length(lsfilesinfolder)
    if ~isempty(regexp(lsfilesinfolder(i).name,'.txt', 'once'))
         if ~isempty(regexp(lsfilesinfolder(i).name,'fixed_miss_v_grab', 'once'))
                    didit=true;
                    return        
         end

    end
end

end

function didit=alreadyDidSuccessVDrop(foldername)

lsfilesinfolder=dir(foldername);
didit=false;
for i=3:length(lsfilesinfolder)
    if ~isempty(regexp(lsfilesinfolder(i).name,'.txt', 'once'))
         if ~isempty(regexp(lsfilesinfolder(i).name,'fixed_tbt_success_v_drop', 'once'))
                    didit=true;
                    return        
         end

    end
end

end

function justdate=getDatestr(foldername)

lsfilesinfolder=dir(foldername);
justdate=[];
for i=3:length(lsfilesinfolder)
    if ~isempty(regexp(lsfilesinfolder(i).name,'.txt', 'once'))
         if ~isempty(regexp(lsfilesinfolder(i).name,'2016', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2017', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2018', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2019', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2020', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2021', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2022', 'once')) || ...
            ~isempty(regexp(lsfilesinfolder(i).name,'2023', 'once'))
                    fullname=lsfilesinfolder(i).name;
                    uptoext=regexp(fullname,'.txt');
                    justdate=fullname(1:uptoext-1);
                    return        
         end

    end
end

end

function placeForO2data=getPlaceOfO2data(mousename,aviName,datestr)

switch mousename
    case 'bi_agouti'
        placeForO2data=['Z:\MICROSCOPE\Kim\JuliaG\Behavior Expts\' mousename '\' datestr '\O2 output\' aviName];
    case 'bi_both'
        placeForO2data=['Z:\MICROSCOPE\Kim\JuliaG\Behavior Expts\' mousename '\' datestr '\O2 output\' aviName];
    case 'mitch_right'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\' aviName];
%         placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\' aviName];
    case 'mitch_mismatch'
%         placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\analysis\' aviName];
    case 'mitch_thin'
%         placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\' aviName];
    case 'mitch_wide'
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
%         placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\' aviName];
    case 'mitch_none'
%         placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\O2 output\' aviName];
        placeForO2data=['Z:\MICROSCOPE\Kim\Mitchell\Behavior\' mousename '\' datestr '\' aviName];
    case 'ssm12'
        placeForO2data=['Z:\MICROSCOPE\Senmiao\Behavior\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'ssm13'
        placeForO2data=['Z:\MICROSCOPE\Senmiao\Behavior\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'ssm14'
        placeForO2data=['Z:\MICROSCOPE\Senmiao\Behavior\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Oct_B'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'July_3'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'May'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Nov_dark'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Nov_ON'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Oct_A'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Sep_3'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Oct_E'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Oct_D'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Sep_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName]; 
    case 'Feb_1'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Feb_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Feb_3'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Feb_4'
%         placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\' aviName];
    case 'Nov_grey'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Dec_d'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Jan_1'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'July_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Nov_2'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Nov_stripe'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case 'Oct_0'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    case '3F_white'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
%         placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\low speed\O2 output\' aviName];
%         placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\O2 output\' aviName];
%         placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\low speed\' aviName];
    case '2F_white'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\O2 output\' aviName];
    case '3F_spot'
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\' mousename '\' datestr '\O2 output\' aviName];
    otherwise
        placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
end

end