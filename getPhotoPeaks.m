function photolocs=getPhotoPeaks(Response,peakOrDip,inRange)

% inRange is in seconds wrt peak of aligncomp_y
[~,ma]=max(mean(Response.aligncomp_y,1,'omitnan'),[],2,'omitnan');
times=mean(Response.aligncomp_x,1);
timeOfAlign=times(ma);
disp(['Aligned at ' num2str(timeOfAlign) 'seconds']);
startRange=timeOfAlign+inRange(1);
endRange=timeOfAlign+inRange(2);
phototimes=mean(Response.unitbyunit_x,1);
[~,startInd]=nanmin(abs(phototimes-startRange));
[~,endInd]=nanmin(abs(phototimes-endRange));
points=nan(size(Response.unitbyunit_y,1),1);
for i=1:size(Response.unitbyunit_y,1)
    R=Response.unitbyunit_y(i,:);
    switch peakOrDip
        case 'peak'
            [~,ma]=max(R(startInd:endInd),[],'all','omitnan');
            ma=ma+startInd-1;
        case 'dip'
            [~,ma]=min(R(startInd:endInd),[],'all','omitnan');
            ma=ma+startInd-1;
        otherwise
            error('Do not recognize argument passed in as peakOrDip in getPhotoPeaks.m');
    end
    points(i)=ma;
end
photolocs=phototimes(points);