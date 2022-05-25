function [cueCD,uncueCD,lateUncueCD]=alignToCompanion(datadir,getCDs,cueCD,uncueCD,lateUncueCD)

excludeHigherFR=true;
% getCDs=true;

ls=dir(datadir);
unitbyunit_x=[];
unitbyunit_y=[];
aligncomp_x=[];
aligncomp_y=[];
padsize=1000;
testForAlignment=false;
unitbaseline=150;
for i=3:length(ls)
    a=load([ls(i).folder '\' ls(i).name]);
    
    if excludeHigherFR
        if nanmean(nanmean(a.dataout.y,1),2)>4
            disp('excluding unit based on high FR');
            continue
        end
    end
    
    % realign first!
    timestep_for_aligncomp=mode(diff(a.alignComp.x));
    timestep_for_unit=mode(diff(a.dataout.x));
    [~,ma]=nanmax(nanmean(a.alignComp.y,1));
    [~,mi]=nanmin(abs(a.dataout.x-a.alignComp.x(ma)));
    if testForAlignment==true
        a.dataout.y(:,mi)=a.dataout.y(:,mi)+100;
    end
    % ensure that each unit data has the same length baseline before mi
    if mi<unitbaseline
        a.dataout.x=[nan(1,unitbaseline-mi) a.dataout.x];
        a.dataout.y=[nan(size(a.dataout.y,1),unitbaseline-mi) a.dataout.y];
    elseif unitbaseline<mi
        a.dataout.x=a.dataout.x(mi-unitbaseline:end);
        a.dataout.y=a.dataout.y(:,mi-unitbaseline:end);
    end
    if ~isempty(unitbyunit_x)
        upTo2=size(aligncomp_x,2);
    else
        upTo2=length(a.alignComp.x)+(padsize-ma);
    end
    if upTo2-(padsize-ma)>length(a.alignComp.x)
        a.alignComp.x=[a.alignComp.x nan(1,upTo2-(padsize-ma)-length(a.alignComp.x))];
        a.alignComp.y=[a.alignComp.y nan(size(a.alignComp.y,1),upTo2-(padsize-ma)-size(a.alignComp.y,2))];
    end
    aligncomp_x=[aligncomp_x; [nan(1,padsize-ma) a.alignComp.x(1:upTo2-(padsize-ma))]];
    aligncomp_y=[aligncomp_y; [nan(1,padsize-ma) nanmean(a.alignComp.y(:,1:upTo2-(padsize-ma)),1)]];
    if ~isempty(unitbyunit_x)
        upTo=size(unitbyunit_x,2);
    else
        upTo=length(a.dataout.x);
    end
    if upTo>size(a.dataout.x,2)
        % truncate
        unitbyunit_x=unitbyunit_x(:,1:size(a.dataout.x,2));
        unitbyunit_y=unitbyunit_y(:,1:size(a.dataout.x,2));
        upTo=size(unitbyunit_x,2);
    end
    unitbyunit_x=[unitbyunit_x; [a.dataout.x(1:upTo)]]; 
    unitbyunit_y=[unitbyunit_y; [nanmean(a.dataout.y(:,1:upTo),1)]];

    disp(['Added ' ls(i).name]);
end

ds=3;
if ds~=1
    unitbyunit_x=downSampMatrix(unitbyunit_x,ds);
    unitbyunit_y=downSampMatrix(unitbyunit_y,ds);
    aligncomp_x=downSampMatrix(aligncomp_x,ds);
    aligncomp_y=downSampMatrix(aligncomp_y,ds);
end

figure();
plot(nanmean(unitbyunit_x,1),unitbyunit_y');
hold on;
plot(nanmean(aligncomp_x,1),aligncomp_y','Color','b');

figure(); 
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1),'Color','k');
hold on;
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)-nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)+nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');

plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');
hold on;
plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)-nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');
plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)+nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');

if getCDs==true
    cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[0 1.5]);
    uncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-1.5 0]);
    lateUncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[4 9.5]);
else
    cueCD=[];
    uncueCD=[];
    lateUncueCD=[];
end

plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD);

end

function CD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,windowWrtAlignCompMax)

% Find cue coding dimension
% Find max of alignment companion
[~,ma]=nanmax(nanmean(aligncomp_y,1));
aligncomp_times=nanmean(aligncomp_x,1);
subtractOff=aligncomp_times(ma);
aligncomp_times=aligncomp_times-subtractOff;
unit_times=nanmean(unitbyunit_x,1);
unit_times=unit_times-subtractOff;
[~,beginIndsForCD]=nanmin(abs(unit_times-windowWrtAlignCompMax(1)));
[~,endIndsForCD]=nanmin(abs(unit_times-windowWrtAlignCompMax(2)));
indsForCD=beginIndsForCD:endIndsForCD;
disp(['Calculating CD for time window ' num2str(unit_times(indsForCD(1))) ' to ' num2str(unit_times(indsForCD(end)))]);

% Find population vector during window
CD=nanmean(unitbyunit_y(:,indsForCD),2);

end

function plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD)

[~,ma]=nanmax(nanmean(aligncomp_y,1));
aligncomp_times=nanmean(aligncomp_x,1);
subtractOff=aligncomp_times(ma);
aligncomp_times=aligncomp_times-subtractOff;
unit_times=nanmean(unitbyunit_x,1);
unit_times=unit_times-subtractOff;

evolveProjectOntoCue=nan(1,size(unitbyunit_y,2));
evolveProjectOntoUncue=nan(1,size(unitbyunit_y,2));
for i=1:size(unitbyunit_y,2)
    evolveProjectOntoCue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),cueCD);
    evolveProjectOntoUncue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),uncueCD);
end

cmap=colormap('cool');
stepintocmap=ceil(size(cmap,1)/nansum(~isnan(evolveProjectOntoCue)));
indsintocmap=1:stepintocmap:nansum(~isnan(evolveProjectOntoCue));

figure();
k=1;
for i=1:length(evolveProjectOntoCue)
    if isnan(evolveProjectOntoCue(i)) | isnan(evolveProjectOntoUncue(i))
        continue
    end
    scatter(evolveProjectOntoUncue(i),evolveProjectOntoCue(i),[],cmap(indsintocmap(k),:));
    k=k+1;
    hold on;
end

end

function proj=projectTimePointOntoCD(activity,CD)

% Project activity onto CD
% Find unit vector of CD
unit_CD=CD./norm(CD);
disp(norm(unit_CD));
proj=dot(activity,unit_CD);

end