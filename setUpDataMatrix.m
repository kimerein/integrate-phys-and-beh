function dataMatrix=setUpDataMatrix(data1,data2,data3,data4,data5,takePointsBeforeZero,takePointsAfterZero)

% set up as units X time X conditions
u1=getUnitbyUnitResponses(data1,takePointsBeforeZero,takePointsAfterZero);
u2=getUnitbyUnitResponses(data2,takePointsBeforeZero,takePointsAfterZero);
u3=getUnitbyUnitResponses(data3,takePointsBeforeZero,takePointsAfterZero);
u4=getUnitbyUnitResponses(data4,takePointsBeforeZero,takePointsAfterZero);
u5=getUnitbyUnitResponses(data5,takePointsBeforeZero,takePointsAfterZero);

% data matrix
dataMatrix=cat(3,u1,u2,u3,u4,u5);

end

function out=getUnitbyUnitResponses(data,takePointsBeforeZero,takePointsAfterZero)

% find alignment companion peak
aligncomp=mean(data.aligncomp_y,1,'omitnan');
aligncomp_x=mean(data.aligncomp_x,1,'omitnan');
[~,ma]=max(aligncomp,[],2,'omitnan');
maxTime=aligncomp_x(ma);
[~,mi]=min(abs(mean(data.unitbyunit_x,1,'omitnan')-maxTime));
if mi-takePointsBeforeZero<1
    amountToPadFront=nan(size(data.unitbyunit_y,1),abs(mi-takePointsBeforeZero)+1);
    startInd=1;
else
    amountToPadFront=nan(size(data.unitbyunit_y,1),0);
    startInd=mi-takePointsBeforeZero;
end
if mi+takePointsAfterZero-1>size(data.unitbyunit_y,2)
    amountToPadBack=nan(size(data.unitbyunit_y,1),mi+takePointsAfterZero-1-size(data.unitbyunit_y,2));
    endInd=size(data.unitbyunit_y,2);
else
    amountToPadBack=nan(size(data.unitbyunit_y,1),0);
    endInd=mi+takePointsAfterZero-1;
end
out=[amountToPadFront data.unitbyunit_y(:,startInd:endInd) amountToPadBack];

end