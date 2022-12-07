function GLM_analysis(whichSess,dd)
% wrapper for first pass at GLM, calls B's poissModel

if length(whichSess)>1
    dd=dd(whichSess);
end

whichUnitsToGrab='_'; plotUnitCriteria=[-100 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
setForUn=settingsForStriatumUnitPlots;
if setForUn.keepAllSingleTrials~=true
    error('need trial by trial data for GLM analysis');
end
response_to_plot='cue';
if length(whichSess)>1
    dd_more=cell(1,length(dd));
    dd_m=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
        dd_m{i}=[dd{i}];
    end
    ResponseCued=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
    grabOtherBehaviorEvents(dd_m);
else
    ResponseCued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
    grabOtherBehaviorEvents(dd{whichSess});
end
ResponseCued.unitbyunit_x=downSampMatrix(ResponseCued.unitbyunit_x,downSampBy);
ResponseCued.unitbyunit_y=downSampMatrix(ResponseCued.unitbyunit_y,downSampBy);
ResponseCued=makeUnitsUnique(ResponseCued);




% temp=Response.unitbyunit_y; temp(isnan(temp))=0;
% behEvents=zeros(size(temp));
% behEvents(:,50)=1;
% poissModel(temp,0:0.04:(size(temp,2)-1)*0.04,behEvents);

end

function grabOtherBehaviorEvents(datadir)

if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    disp(['reading in beh events from ' datadir]);
    ls=dir(datadir);
%     for i=3:length(ls)
%         a=[];
%         if ismac()
%             a=load([ls(i).folder '/' ls(i).name]);
%         else
%             a=load([ls(i).folder '\' ls(i).name]);
%         end

end
end