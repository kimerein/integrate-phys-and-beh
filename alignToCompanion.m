function alignToCompanion(datadir)

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
end

ds=1;
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