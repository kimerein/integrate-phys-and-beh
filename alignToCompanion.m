function alignToCompanion(datadir)

ls=dir(datadir);
unitbyunit_x=[];
unitbyunit_y=[];
aligncomp_x=[];
aligncomp_y=[];
padsize=1000;
for i=3:length(ls)
    a=load([ls(i).folder '\' ls(i).name]);
    % realign first!
    [~,ma]=nanmax(nanmean(a.alignComp.y,1));
    if ~isempty(unitbyunit_x)
        upTo=size(unitbyunit_x,2);
        upTo2=size(aligncomp_x,2);
    else
        upTo=length(a.dataout.x)+(padsize-ma);
        upTo2=length(a.alignComp.x)+(padsize-ma);
    end
    unitbyunit_x=[unitbyunit_x; [nan(1,padsize-ma) a.dataout.x(1:upTo-(padsize-ma))]];
    unitbyunit_y=[unitbyunit_y; [nan(1,padsize-ma) nanmean(a.dataout.y(:,1:upTo-(padsize-ma)),1)]];
    aligncomp_x=[aligncomp_x; [nan(1,padsize-ma) a.alignComp.x(1:upTo2-(padsize-ma))]];
    aligncomp_y=[aligncomp_y; [nan(1,padsize-ma) nanmean(a.alignComp.y(:,1:upTo2-(padsize-ma)),1)]];
end

ds=5;
if ds~=1
    unitbyunit_x=downSampMatrix(unitbyunit_x,ds);
    unitbyunit_y=downSampMatrix(unitbyunit_y,ds);
    aligncomp_x=downSampMatrix(aligncomp_x,ds);
    aligncomp_y=downSampMatrix(aligncomp_y,ds);
end

figure();
plot(nanmean(unitbyunit_x,1),unitbyunit_y');
hold on;
plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');

figure();
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1),'Color','k');
hold on;
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)-nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)+nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');

plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');
hold on;
plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)-nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');
plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)+nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');