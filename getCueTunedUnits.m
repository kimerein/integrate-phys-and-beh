function [cuetuningindex,isCueTuned]=getCueTunedUnits(cueR,uncuedReachR,method,meanormax,ds,cuebaserange,cueonrange,reachbaserange,reachonrange)

% ds=25;
smoo=[]; %100; %80; %42; % NOTHING ABOVE 100 BCZ WILL SMEAR IN OPTO RESPONSE
noGaussSmooth=true;
Zscore=false;
% cuebaserange=[5 9]; %[-3 -2]; %[-0.75 -0.54];
% cueonrange=[-0.375 1.25]; %[-0.375 0.375]; %[-0.25 0.8];
% reachbaserange=[-3 -2];
% reachonrange=[-2 0];
% reachbaserange=cuebaserange;
% reachonrange=cueonrange;

% match which units and get responses unit by unit
r=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cueR,uncuedReachR,'meanAcrossUnits',ds,true);
cueresp=r.response1; uncuedreachresp=r.response2;
if Zscore==true
    sd=std([cueresp.unitbyunit(:,cueresp.t>=4 & cueresp.t<9) uncuedreachresp.unitbyunit(:,uncuedreachresp.t>=4 & uncuedreachresp.t<9)],[],2,'omitnan');
    sd(sd<1/30)=1/30;
    cueresp.unitbyunit=cueresp.unitbyunit./repmat(sd,1,size(cueresp.unitbyunit,2));
    cueresp.unitbyunit(isinf(cueresp.unitbyunit))=nan;

    uncuedreachresp.unitbyunit=uncuedreachresp.unitbyunit./repmat(sd,1,size(uncuedreachresp.unitbyunit,2));
    uncuedreachresp.unitbyunit(isinf(uncuedreachresp.unitbyunit))=nan;
end
if ~isempty(smoo)
    for i=1:size(cueresp.unitbyunit,1)
        temp=cueresp.unitbyunit(i,:);
        if noGaussSmooth==true
            cueresp.unitbyunit(i,:)=smooth(temp,smoo);
        else
            cueresp.unitbyunit(i,:)=smoothdata(temp,'gaussian',smoo);
        end
    end
    for i=1:size(uncuedreachresp.unitbyunit,1)
        temp=uncuedreachresp.unitbyunit(i,:);
         if noGaussSmooth==true
            uncuedreachresp.unitbyunit(i,:)=smooth(temp,smoo);
         else
            uncuedreachresp.unitbyunit(i,:)=smoothdata(temp,'gaussian',smoo);
         end
    end
end

% get cue response vs response during uncued reach
t=cueresp.unittimes;
switch meanormax
    case 'max'
        beforecue=max(cueresp.unitbyunit(:,t>cuebaserange(1) & t<cuebaserange(2)),[],2,'omitnan'); 
        duringcue=max(cueresp.unitbyunit(:,t>cueonrange(1) & t<cueonrange(2)),[],2,'omitnan'); 
    case 'mean'
        beforecue=mean(cueresp.unitbyunit(:,t>cuebaserange(1) & t<cuebaserange(2)),2,'omitnan'); 
        duringcue=mean(cueresp.unitbyunit(:,t>cueonrange(1) & t<cueonrange(2)),2,'omitnan'); 
end
cuer=duringcue-beforecue;

t=cueresp.unittimes;
switch meanormax
    case 'max'
        beforereach=max(uncuedreachresp.unitbyunit(:,t>reachbaserange(1) & t<reachbaserange(2)),[],2,'omitnan');
        duringreach=max(uncuedreachresp.unitbyunit(:,t>reachonrange(1) & t<reachonrange(2)),[],2,'omitnan');
    case 'mean'
        beforereach=mean(uncuedreachresp.unitbyunit(:,t>reachbaserange(1) & t<reachbaserange(2)),2,'omitnan');
        duringreach=mean(uncuedreachresp.unitbyunit(:,t>reachonrange(1) & t<reachonrange(2)),2,'omitnan');
end
reachr=duringreach-beforereach;

% define cue tuned units using cue tuning index
switch method
    case 'justcue_v_justuncue'
        cuetuningindex=duringcue-duringreach;
    case 'vs_uncued_reach'
        cuer(abs(cuer)<0.0001)=0;
        reachr(abs(reachr)<0.0001)=0;
        cuetuningindex=(cuer-reachr)./(cuer+reachr);
    case'vs_uncued_reach_no_index'
        cuer(abs(cuer)<0.0001)=0;
        reachr(abs(reachr)<0.0001)=0;
        cuetuningindex=(cuer-reachr);
    case 'cue_vs_baseline'
        cuetuningindex=(duringcue-beforecue)./(duringcue+beforecue);
        cuetuningindex(cuetuningindex>1)=1;
        cuetuningindex(cuetuningindex<-1)=-1;
        fill_insufficient_data_w_nans=false;
        if fill_insufficient_data_w_nans
            cuetuningindex(cuetuningindex==1)=nan;
            cuetuningindex(cuetuningindex==-1)=nan;
        end
    case 'cue_vs_baseline_no_index'
        cuetuningindex=(duringcue-beforecue);
    case 'justcue'
        cuetuningindex=duringcue;
end
isCueTuned=cuetuningindex>0;

end