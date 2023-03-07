function [unitbyunit_names,excluded,D1orD2taggingExpt,D1taggedCells,A2ataggedCells,fromWhichSess,fromWhichUnit]=alignToReadInUnitNames(datadir,onlyTakeTheseUnits,settings)

% settings=settingsForStriatumUnitPlots;
maxUnitsPerSess=settings.maxUnitsPerSess;
if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
unitbyunit_names={};
excluded=zeros(maxUnitsPerSess*length(dd),1);
D1orD2taggingExpt=nan(length(dd),1); % will be 1 for D1, 2 for A2a tagging session
D1taggedCells=zeros(maxUnitsPerSess*length(dd),1);
A2ataggedCells=zeros(maxUnitsPerSess*length(dd),1);
fromWhichSess=nan(maxUnitsPerSess*length(dd),1);
fromWhichUnit=nan(maxUnitsPerSess*length(dd),1);
excluded_count=1;
units_count=1;
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    disp(['reading in from ' datadir]);
    ls=dir(datadir);
    for i=3:length(ls)
        a=[];
        if contains(ls(i).name,'trainingSet') || contains(ls(i).name,'testSet') % ignore these files
            continue
        end
        
        if ismac()
            a=load([ls(i).folder '/' ls(i).name]);
        else
            a=load([ls(i).folder '\' ls(i).name]);
        end
        
        if ~isempty(regexp(ls(i).name, 'D1tagged'))
            D1orD2taggingExpt(j)=1;
            D1taggedCells(excluded_count)=1;
        elseif ~isempty(regexp(ls(i).name, 'A2atagged'))
            D1orD2taggingExpt(j)=2;
            A2ataggedCells(excluded_count)=1;
        end

        if ~isempty(onlyTakeTheseUnits)
            if isempty(regexp(ls(i).name, onlyTakeTheseUnits))
                % doesn't contain string, continue and skip this unit
                excluded(excluded_count)=1;
                excluded_count=excluded_count+1;
                continue
            end
        end

        if ~isfield(settings,'isPhotometry')
            settings.isPhotometry=false;
        end

        if settings.isPhotometry==false
            [unitTest,dontUseTrials]=doUnitTest(ls(i).folder, ls(i).name);
            if unitTest
            else
                if settings.isPhotometry==true
                else
                    excluded(excluded_count)=1;
                    excluded_count=excluded_count+1;
                    continue % failed unit test
                end
            end
        end

        if isempty(a.dataout)
            disp([ls(i).folder '\' ls(i).name ' is empty ... skipping']);
            excluded(excluded_count)=1;
            excluded_count=excluded_count+1;
            continue
        end

        if settings.isPhotometry==true
            dontUseTrials=zeros(size(a.dataout.y,1),1);
        end

        if length(dontUseTrials)~=size(a.dataout.y,1)
            % as a default, take all trials
            % but I might want to modify this behavior later
            % sometimes dontUseTrials is empty if figure opening failed,
            % or other
            disp('dontUseTrials length does not match size of a.dataout.y in line 91 of alignToCompanion.m, but taking all trials');
        end
        excluded_count=excluded_count+1;

%         if excludeHigherFR
%             if nanmean(nanmean(a.dataout.y(:,1:beforeOptoBaseline),1),2)>excludeAboveFR
%                 disp(['excluding unit ' ls(i).name ' based on high FR']);
%                 excluded=[excluded 1];
%                 continue
%             else
%                 excluded=[excluded 0];
%             end
%         else
%             if nanmean(nanmean(a.dataout.y(:,1:beforeOptoBaseline),1),2)>excludeAboveFR
%                 excluded=[excluded 1];
%             else
%                 excluded=[excluded 0];
%             end
%         end
        
        s=ls(i).name;
        r=regexp(s,'_');
        unitbyunit_names{units_count}=s(1:r-1);
        fromWhichSess(units_count)=j;
        fromWhichUnit(units_count)=units_count;

%         disp(['Added ' ls(i).name]);
        disp(['Added ' s(1:r-1)]);
        units_count=units_count+1;
    end
end
excluded_count=excluded_count-1;
units_count=units_count-1;
excluded=excluded(1:excluded_count);
D1taggedCells=D1taggedCells(1:excluded_count);
A2ataggedCells=A2ataggedCells(1:excluded_count);