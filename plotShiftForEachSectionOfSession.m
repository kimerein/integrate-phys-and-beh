function [uncued,cued]=plotShiftForEachSectionOfSession(data1,data2,doPlot)

if length(data1.cued)>length(data2.cued)
    % interp
    orig_length_data2=length(data2.cued);
    data2.cued=interp1(linspace(1,length(data1.cued),orig_length_data2),data2.cued,1:length(data1.cued));
    data2.uncued=interp1(linspace(1,length(data1.cued),orig_length_data2),data2.uncued,1:length(data1.cued));
elseif length(data2.cued)>length(data1.cued)
    % interp
    orig_length_data1=length(data1.cued);
    data1.cued=interp1(linspace(1,length(data2.cued),orig_length_data1),data1.cued,1:length(data2.cued));
    data1.uncued=interp1(linspace(1,length(data2.cued),orig_length_data1),data1.uncued,1:length(data2.cued));
end

if doPlot==true
    figure();
    quiver(0,0,nanmean(data1.uncued-data2.uncued),nanmean(data1.cued-data2.cued),'Color','k');
end
uncued=nanmean(data1.uncued-data2.uncued);
cued=nanmean(data1.cued-data2.cued);

end