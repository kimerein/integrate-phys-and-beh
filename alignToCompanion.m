function [cueCD,uncueCD,lateUncueCD]=alignToCompanion(datadir,getCDs,cueCD,uncueCD,lateUncueCD)

excludeHigherFR=false;
cutAtTime=3; % stop plotting this many seconds after max of alignment companion
ds=6; % downsample bin size
onlyTakeTheseUnits=''; % if is not empty, will only take units with this string in filename, e.g., 'D1tagged'
% or make cutAtTime empty to plot all time points
% getCDs=true;

unitbyunit_x=[];
unitbyunit_y=[];
aligncomp_x=[];
aligncomp_y=[];
padsize=1000;
testForAlignment=false;
unitbaseline=300; %150;
if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    ls=dir(datadir);
    for i=3:length(ls)
        a=load([ls(i).folder '\' ls(i).name]);
        
        if excludeHigherFR
            if nanmean(nanmean(a.dataout.y,1),2)>4
                disp(['excluding unit ' ls(i).name ' based on high FR']);
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
end

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
%     % for cue
    cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 1.25]); % cue duration is 250 ms
%     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 0.25]); % cue duration is 250 ms
    uncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-1.75 -0.25]);
    lateUncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[4 9.5]);
    % for reach
%     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-1 0]); % around reach
% %     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 0.25]); % cue duration is 250 ms
%     uncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[3 9]); % after reach
%     lateUncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[4 9.5]);
else
    % use CDs passed in
end

% If want to cut off the plotting at a particular time after the
% aligncomp_y max, use cutAtTime
% else
% cutAtTime=[];
plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD,cutAtTime);

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

function plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD,cutAtTime)

[~,ma]=nanmax(nanmean(aligncomp_y,1));
aligncomp_times=nanmean(aligncomp_x,1);
subtractOff=aligncomp_times(ma);
aligncomp_times=aligncomp_times-subtractOff;
unit_times=nanmean(unitbyunit_x,1);
unit_times=unit_times-subtractOff;

if ~isempty(cutAtTime)
    % find ind into aligncomp for cutoff
    [~,cutoffind]=nanmin(abs(aligncomp_times-cutAtTime));
    % find ind into unitbyunit for cutoff
    [~,cutoffind2]=nanmin(abs(unit_times-cutAtTime));
    aligncomp_x=aligncomp_x(:,1:cutoffind);
    aligncomp_times=aligncomp_times(1:cutoffind);
    aligncomp_y=aligncomp_y(:,1:cutoffind);
    unitbyunit_x=unitbyunit_x(:,1:cutoffind2);
    unit_times=unit_times(1:cutoffind2);
    unitbyunit_y=unitbyunit_y(:,1:cutoffind2);    
end

evolveProjectOntoCue=nan(1,size(unitbyunit_y,2));
evolveProjectOntoUncue=nan(1,size(unitbyunit_y,2));
for i=1:size(unitbyunit_y,2)
    evolveProjectOntoCue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),cueCD);
    evolveProjectOntoUncue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),uncueCD);
end

cmap=colormap('cool');
stepintocmap=size(cmap,1)/nansum(~isnan(evolveProjectOntoCue));
indsintocmap=0:stepintocmap:(nansum(~isnan(evolveProjectOntoCue))-1)*stepintocmap;
indsintocmap(1)=0.0001;

figure();
k=1;
for i=1:length(evolveProjectOntoCue)
    if isnan(evolveProjectOntoCue(i)) | isnan(evolveProjectOntoUncue(i))
        continue
    end
    scatter(evolveProjectOntoUncue(i),evolveProjectOntoCue(i),[],cmap(ceil(indsintocmap(k)),:));
    k=k+1;
    hold on;
end
daspect([1 1 1]);
xlabel('Projection onto uncued CD');
ylabel('Projection onto cued CD');
linemi=min([evolveProjectOntoUncue evolveProjectOntoCue],[],2,'omitnan');
linema=max([evolveProjectOntoUncue evolveProjectOntoCue],[],2,'omitnan');
line([linemi linema],[linemi linema],'Color',[0.8 0.8 0.8]);

figure(); plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');
hold on; 
k=1;
for i=1:length(evolveProjectOntoCue)
    if isnan(evolveProjectOntoCue(i)) | isnan(evolveProjectOntoUncue(i))
        continue
    end
    scatter(nanmean(unitbyunit_x(:,i),1),nanmean(unitbyunit_y(:,i),1),[],cmap(ceil(indsintocmap(k)),:));
    k=k+1;
    hold on;
end
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1),'Color','k');

end

function proj=projectTimePointOntoCD(activity,CD)

% Project activity onto CD
% Find unit vector of CD
unit_CD=CD./norm(CD);
proj=dot(activity,unit_CD);

end