function combineFigure4DataAcrossMice(dd,useOptoForThisGp,whichtoplot_eventCond,timeStep,cueind,doCDF,doLEDcdf,calcCued,atleast_n_trials,ds)

% combineFigure4DataAcrossMice(dd([1:9 11:16]),useOptoForThisGp([1:9 11:16]),'cued success_seqLength3_win2and3_minus2to0_optoDurforcuepoi4then1NOOPREwin3',0.035,94,false,false,false,200,1);

rr_noLED_tri1_uncued=[];
rr_noLED_tri1_cued=[];
rr_noLED_trinext_uncued=[];
rr_noLED_trinext_cued=[];

rr_LED_tri1_uncued=[];
rr_LED_tri1_cued=[];
rr_LED_trinext_uncued=[];
rr_LED_trinext_cued=[];

if iscell(whichtoplot_eventCond)
    % more than one event condition to combine
    thesetoplot=cell(1,length(dd)*length(whichtoplot_eventCond));
    newdd=cell(1,length(dd)*length(whichtoplot_eventCond));
    newUseOpto=nan(1,length(dd)*length(whichtoplot_eventCond));
    k=1;
    for i=1:length(dd)
        for j=1:length(whichtoplot_eventCond)
            thesetoplot{k}=whichtoplot_eventCond{j};
            newdd{k}=dd{i};
            newUseOpto(k)=useOptoForThisGp(i);
            k=k+1;
        end
    end
    dd=newdd;
    useOptoForThisGp=newUseOpto;
end

returnThis_control=cell(length(dd));
returnThis_LED=cell(length(dd));

forcdf_rawreach_trial1=[];
forcdf_rawreach_trial2=[];

s_trial1_alltrials_uncued=[];
s_trial1_alltrials_cued=[];
s_alltrials_uncued=[];
s_alltrials_cued=[];
sLED_trial1_alltrials_uncued=[];
sLED_trial1_alltrials_cued=[];
sLED_alltrials_uncued=[];
sLED_alltrials_cued=[];

m_trial1_alltrials_uncued=[];
m_trial1_alltrials_cued=[];
m_alltrials_uncued=[];
m_alltrials_cued=[];
mLED_trial1_alltrials_uncued=[];
mLED_trial1_alltrials_cued=[];
mLED_alltrials_uncued=[];
mLED_alltrials_cued=[];
returnControlPDF=[];
returnLEDPDF=[];
k=1;
for i=1:length(dd)
    currdir=dd{i};
    if doCDF==false
        % Get mouse by mouse
        % and sess by sess
        % rows are sessid, column 1 is sessid, column2 is mouseid
        a=load(fullfile(currdir,'sessIDandMouseID_NOOPRE.mat'));
        sessIDandMouseID=a.sessIDandMouseID;

        if iscell(whichtoplot_eventCond)
            % more than one event condition to combine
            whichtoplot=thesetoplot{k};
            k=k+1;
        else
            whichtoplot=whichtoplot_eventCond;
        end

        if ~exist(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'),'file')
            continue
        end
        
        a=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'));
        returnThis_control{i}=a.returnThis;

        a=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThis.mat'));
        returnThis_LED{i}=a.returnThis;

        % combine PDFs
        if i==2
            [~,~,returnControlPDF]=combineReachPlotDatasets(returnThis_control{1},[1 2],returnThis_control{2},[1 2],'k'); close all;
            if ~all(isnan(returnThis_LED{i}.data1_mean{1}))
                [~,~,returnLEDPDF]=combineReachPlotDatasets(returnThis_LED{1},[1 2],returnThis_LED{2},[1 2],'k'); close all;
            end
        elseif i>2
            [~,~,returnControlPDF]=combineReachPlotDatasets(returnControlPDF,[1 2],returnThis_control{i},[1 2],'k'); close all;
            if ~all(isnan(returnThis_LED{i}.data1_mean{1}))
                [~,~,returnLEDPDF]=combineReachPlotDatasets(returnLEDPDF,[1 2],returnThis_LED{i},[1 2],'k'); close all;
            end
        end

        a=load(fullfile(currdir,whichtoplot,'returnout.mat'));
        rout=a.returnout;
        if isempty(rout.reachrates_noLED)
            % nothing to add
            continue
        end
        

        temp=rout.reachrates_noLED.trial1_alltrials_uncued; temp=temp(:).';
        rr_noLED_tri1_uncued=[rr_noLED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_noLED.trial1_alltrials_cued; temp=temp(:).';
        rr_noLED_tri1_cued=[rr_noLED_tri1_cued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_uncued; temp=temp(:).';
        rr_noLED_trinext_uncued=[rr_noLED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_cued; temp=temp(:).';
        rr_noLED_trinext_cued=[rr_noLED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_noLED_tri1_uncued);
        rr_noLED_tri1_uncued=rr_noLED_tri1_uncued(~ina);
        rr_noLED_tri1_cued=rr_noLED_tri1_cued(~ina);
        rr_noLED_trinext_uncued=rr_noLED_trinext_uncued(~ina);
        rr_noLED_trinext_cued=rr_noLED_trinext_cued(~ina);

        if isempty(rout.reachrates_LED) || useOptoForThisGp(i)==0
            % nothing to add
            [sbys,mbym]=getMbyM_SbyS(sessIDandMouseID,rout,atleast_n_trials);
            s_trial1_alltrials_uncued=[s_trial1_alltrials_uncued sbys.reachrates_noLED.trial1_alltrials_uncued'];
            s_trial1_alltrials_cued=[s_trial1_alltrials_cued sbys.reachrates_noLED.trial1_alltrials_cued'];
            s_alltrials_uncued=[s_alltrials_uncued sbys.reachrates_noLED.alltrials_uncued'];
            s_alltrials_cued=[s_alltrials_cued sbys.reachrates_noLED.alltrials_cued'];

            m_trial1_alltrials_uncued=[m_trial1_alltrials_uncued mbym.reachrates_noLED.trial1_alltrials_uncued];
            m_trial1_alltrials_cued=[m_trial1_alltrials_cued mbym.reachrates_noLED.trial1_alltrials_cued];
            m_alltrials_uncued=[m_alltrials_uncued mbym.reachrates_noLED.alltrials_uncued];
            m_alltrials_cued=[m_alltrials_cued mbym.reachrates_noLED.alltrials_cued];

            sLED_trial1_alltrials_uncued=[sLED_trial1_alltrials_uncued sbys.reachrates_LED.trial1_alltrials_uncued'];
            sLED_trial1_alltrials_cued=[sLED_trial1_alltrials_cued sbys.reachrates_LED.trial1_alltrials_cued'];
            sLED_alltrials_uncued=[sLED_alltrials_uncued sbys.reachrates_LED.trial1_alltrials_uncued'];
            sLED_alltrials_cued=[sLED_alltrials_cued sbys.reachrates_LED.trial1_alltrials_cued'];

            mLED_trial1_alltrials_uncued=[mLED_trial1_alltrials_uncued mbym.reachrates_LED.trial1_alltrials_uncued];
            mLED_trial1_alltrials_cued=[mLED_trial1_alltrials_cued mbym.reachrates_LED.trial1_alltrials_cued];
            mLED_alltrials_uncued=[mLED_alltrials_uncued mbym.reachrates_LED.trial1_alltrials_uncued];
            mLED_alltrials_cued=[mLED_alltrials_cued mbym.reachrates_LED.trial1_alltrials_cued];
            continue
        end
        temp=rout.reachrates_LED.trial1_alltrials_uncued; temp=temp(:).';
        rr_LED_tri1_uncued=[rr_LED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_LED.trial1_alltrials_cued; temp=temp(:).';
        rr_LED_tri1_cued=[rr_LED_tri1_cued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_uncued; temp=temp(:).';
        rr_LED_trinext_uncued=[rr_LED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_cued; temp=temp(:).';
        rr_LED_trinext_cued=[rr_LED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_LED_tri1_uncued);
        rr_LED_tri1_uncued=rr_LED_tri1_uncued(~ina);
        rr_LED_tri1_cued=rr_LED_tri1_cued(~ina);
        rr_LED_trinext_uncued=rr_LED_trinext_uncued(~ina);
        rr_LED_trinext_cued=rr_LED_trinext_cued(~ina);

        [sbys,mbym]=getMbyM_SbyS(sessIDandMouseID,rout,atleast_n_trials);
        s_trial1_alltrials_uncued=[s_trial1_alltrials_uncued sbys.reachrates_noLED.trial1_alltrials_uncued'];
        s_trial1_alltrials_cued=[s_trial1_alltrials_cued sbys.reachrates_noLED.trial1_alltrials_cued'];
        s_alltrials_uncued=[s_alltrials_uncued sbys.reachrates_noLED.alltrials_uncued'];
        s_alltrials_cued=[s_alltrials_cued sbys.reachrates_noLED.alltrials_cued'];

        m_trial1_alltrials_uncued=[m_trial1_alltrials_uncued mbym.reachrates_noLED.trial1_alltrials_uncued];
        m_trial1_alltrials_cued=[m_trial1_alltrials_cued mbym.reachrates_noLED.trial1_alltrials_cued];
        m_alltrials_uncued=[m_alltrials_uncued mbym.reachrates_noLED.alltrials_uncued];
        m_alltrials_cued=[m_alltrials_cued mbym.reachrates_noLED.alltrials_cued];

        sLED_trial1_alltrials_uncued=[sLED_trial1_alltrials_uncued sbys.reachrates_LED.trial1_alltrials_uncued'];
        sLED_trial1_alltrials_cued=[sLED_trial1_alltrials_cued sbys.reachrates_LED.trial1_alltrials_cued'];
        sLED_alltrials_uncued=[sLED_alltrials_uncued sbys.reachrates_LED.alltrials_uncued'];
        sLED_alltrials_cued=[sLED_alltrials_cued sbys.reachrates_LED.alltrials_cued'];

        mLED_trial1_alltrials_uncued=[mLED_trial1_alltrials_uncued mbym.reachrates_LED.trial1_alltrials_uncued];
        mLED_trial1_alltrials_cued=[mLED_trial1_alltrials_cued mbym.reachrates_LED.trial1_alltrials_cued];
        mLED_alltrials_uncued=[mLED_alltrials_uncued mbym.reachrates_LED.alltrials_uncued];
        mLED_alltrials_cued=[mLED_alltrials_cued mbym.reachrates_LED.alltrials_cued];
    else
        if iscell(whichtoplot_eventCond)
            % more than one event condition to combine
            whichtoplot=thesetoplot{k};
            k=k+1;
        else
            whichtoplot=whichtoplot_eventCond;
        end

        % do CDF only
        if doLEDcdf==true
            a=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThisCDF.mat'));
            rthisCDF=a.returnThisCDF;
            forcdf_rawreach_trial1=[forcdf_rawreach_trial1; rthisCDF.trial1_rawReachMatrix];
            forcdf_rawreach_trial2=[forcdf_rawreach_trial2; rthisCDF.trial2_rawReachMatrix];
        else
            a=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThisCDF.mat'));
            rthisCDF=a.returnThisCDF;
            forcdf_rawreach_trial1=[forcdf_rawreach_trial1; rthisCDF.trial1_rawReachMatrix];
            forcdf_rawreach_trial2=[forcdf_rawreach_trial2; rthisCDF.trial2_rawReachMatrix];
        end
    end
end

if doCDF==true
    returnThisCDF.trial1_rawReachMatrix=forcdf_rawreach_trial1;
    returnThisCDF.trial2_rawReachMatrix=forcdf_rawreach_trial2;

    if isempty(ds)
        ds=1;
    end
    % Reach rate over trials
    [n,x]=cityscape_hist(mean(returnThisCDF.trial1_rawReachMatrix,1,'omitnan'),0:timeStep:(size(returnThisCDF.trial1_rawReachMatrix,2)-1)*(timeStep));
    figure(); plot(x,n,'Color','k');
    [n2,x2]=cityscape_hist(mean(returnThisCDF.trial2_rawReachMatrix,1,'omitnan'),0:timeStep:(size(returnThisCDF.trial2_rawReachMatrix,2)-1)*(timeStep));
    hold on; plot(x2,n2,'Color','r');
    title('Reach rate average over trials');

    % Reach PDF averaged over trials
    temp=downSampMatrix(returnThisCDF.trial1_rawReachMatrix,ds); temp=temp./repmat(nansum(temp,2),1,size(temp,2)); % make each row its own pdf
    ds_times=0:timeStep*ds:(size(temp,2)-1)*(timeStep*ds);
    [n,x]=cityscape_hist(mean(temp,1,'omitnan'),ds_times);
    figure(); plot(x,n,'Color','k');
    temp=downSampMatrix(returnThisCDF.trial2_rawReachMatrix,ds); temp=temp./repmat(nansum(temp,2),1,size(temp,2)); % make each row its own pdf
    [n2,x2]=cityscape_hist(mean(temp,1,'omitnan'),ds_times);
    hold on; plot(x2,n2,'Color','r');
    ylabel('Average of each trial as PDF');

    % plot reach probability per bin
    temp=downSampMatrix(returnThisCDF.trial1_rawReachMatrix,ds); temp(temp>0)=1;
    ds_times=0:timeStep*ds:(size(temp,2)-1)*(timeStep*ds);
%     returnThisCDF.trial1_rawReachMatrix=temp;
    returnThisCDF.trial1_rawReachMatrix=downSampMatrix(returnThisCDF.trial1_rawReachMatrix,ds);
    [n,x]=cityscape_hist(mean(temp,1,'omitnan'),ds_times);
    figure(); plot(x,n,'Color','k');
    temp=downSampMatrix(returnThisCDF.trial2_rawReachMatrix,ds); temp(temp>0)=1;
%     returnThisCDF.trial2_rawReachMatrix=temp;
    returnThisCDF.trial2_rawReachMatrix=downSampMatrix(returnThisCDF.trial2_rawReachMatrix,ds);
    [n2,x2]=cityscape_hist(mean(temp,1,'omitnan'),ds_times);
    hold on; plot(x2,n2,'Color','r');
    ylabel('Probability that mouse reaches over trials');

    stopPDFat=9.5;
    if ~isempty(stopPDFat)
        figure(); plot(x,n./nansum(n(x<stopPDFat)),'Color','k'); hold on;
        plot(x2,n2./nansum(n2(x2<stopPDFat)),'Color','r');
        xlim([min(x,[],'all','omitnan') stopPDFat]);
        ylabel('Reach PDF');
    else
        figure(); plot(x,n./nansum(n),'Color','k'); hold on;
        plot(x2,n2./nansum(n2),'Color','r');
        ylabel('Reach PDF');
    end


    getMeanAndBootstrapForCDF(timeStep*ds,returnThisCDF,floor(cueind/ds));
else
    % plot pdfs
    [returnControlPDF.data1_se{1},returnControlPDF.data2_se{1},returnControlPDF.data1_mean{1},returnControlPDF.data2_mean{1},returnControlPDF.time_for_x]=downsamplePDF(returnControlPDF.data1_se{1},returnControlPDF.data2_se{1},returnControlPDF.data1_mean{1},returnControlPDF.data2_mean{1},returnControlPDF.time_for_x,ds);
    figure();
    plotPDF(returnControlPDF.data1_mean{1},returnControlPDF.data1_se{1},returnControlPDF.time_for_x,returnControlPDF.n,'k'); hold on;
    plotPDF(returnControlPDF.data2_mean{1},returnControlPDF.data2_se{1},returnControlPDF.time_for_x,returnControlPDF.n,'b');
    title('Control trial sequence');
    if ~isempty(returnLEDPDF)
        [returnLEDPDF.data1_se{1},returnLEDPDF.data2_se{1},returnLEDPDF.data1_mean{1},returnLEDPDF.data2_mean{1},returnLEDPDF.time_for_x]=downsamplePDF(returnLEDPDF.data1_se{1},returnLEDPDF.data2_se{1},returnLEDPDF.data1_mean{1},returnLEDPDF.data2_mean{1},returnLEDPDF.time_for_x,ds);
        figure();
        plotPDF(returnLEDPDF.data1_mean{1},returnLEDPDF.data1_se{1},returnLEDPDF.time_for_x,returnLEDPDF.n,'r'); hold on;
        plotPDF(returnLEDPDF.data2_mean{1},returnLEDPDF.data2_se{1},returnLEDPDF.time_for_x,returnLEDPDF.n,'b');
        title('LED trial sequence');
    end


    % do bootstrapped scatter
    figure();
    [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_noLED_tri1_uncued,rr_noLED_tri1_cued,'k','k',false,calcCued);
    plotMeAndSe(rr_noLED_tri1_uncued,rr_noLED_tri1_cued,'k',2,false,calcCued);
    [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_noLED_trinext_uncued,rr_noLED_trinext_cued,[0.5 0.5 0.5],[0.5 0.5 0.5],false,calcCued);
    plotMeAndSe(rr_noLED_trinext_uncued,rr_noLED_trinext_cued,[0.5 0.5 0.5],2,false,calcCued);

    figure();
    [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_LED_tri1_uncued,rr_LED_tri1_cued,'k','r',false,calcCued);
    plotMeAndSe(rr_LED_tri1_uncued,rr_LED_tri1_cued,'k',2,false,calcCued);
    [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_LED_trinext_uncued,rr_LED_trinext_cued,'r','r',false,calcCued);
    plotMeAndSe(rr_LED_trinext_uncued,rr_LED_trinext_cued,'r',2,false,calcCued);

    figure();
    plotByBy(s_trial1_alltrials_uncued,s_trial1_alltrials_cued,s_alltrials_uncued,s_alltrials_cued,'k',calcCued,false);
    plotByBy(sLED_trial1_alltrials_uncued,sLED_trial1_alltrials_cued,sLED_alltrials_uncued,sLED_alltrials_cued,'r',calcCued,false);
    title('Sess by sess');

    figure();
    plotMeAndSe(s_alltrials_uncued-s_trial1_alltrials_uncued,s_alltrials_cued-s_trial1_alltrials_cued,'k',2,false,calcCued);
    hold on;
    plotMeAndSe(sLED_alltrials_uncued-sLED_trial1_alltrials_uncued,sLED_alltrials_cued-sLED_trial1_alltrials_cued,'r',2,false,calcCued);
    title('Sess by sess mean and se');


    figure();
    plotByBy(m_trial1_alltrials_uncued,m_trial1_alltrials_cued,m_alltrials_uncued,m_alltrials_cued,'k',calcCued,true);
    plotByBy(mLED_trial1_alltrials_uncued,mLED_trial1_alltrials_cued,mLED_alltrials_uncued,mLED_alltrials_cued,'r',calcCued,true);
    title('Mouse by mouse');

    stats_compareChangeInRate(m_alltrials_uncued-m_trial1_alltrials_uncued,m_alltrials_cued-m_trial1_alltrials_cued,mLED_alltrials_uncued-mLED_trial1_alltrials_uncued,mLED_alltrials_cued-mLED_trial1_alltrials_cued);
    disp('AND NOW SESS BY SESS');
    stats_compareChangeInRate(s_alltrials_uncued-s_trial1_alltrials_uncued,s_alltrials_cued-s_trial1_alltrials_cued,sLED_alltrials_uncued-sLED_trial1_alltrials_uncued,sLED_alltrials_cued-sLED_trial1_alltrials_cued);
end

end

function [data1_se,data2_se,data1_mean,data2_mean,timeBins]=downsamplePDF(data1_se,data2_se,data1_mean,data2_mean,timeBins,ds)

data1_se=data1_se.*sqrt(ds/ds^2);
data2_se=data2_se.*sqrt(ds/ds^2);
data1_se=downSampAv(data1_se,ds);
data2_se=downSampAv(data2_se,ds);

data1_mean=downSampAv(data1_mean,ds);
data2_mean=downSampAv(data2_mean,ds);
timeBins=downSampAv(timeBins,ds);

end

function plotPDF(mecombo,combo_se,time_for_x,ntrials,linecol)

[n,x]=cityscape_hist(mecombo,time_for_x);
plot(x,n,'Color',linecol);
hold on;
[n,x]=cityscape_hist(mecombo-combo_se,time_for_x);
plot(x,n,'Color',linecol);
[n,x]=cityscape_hist(mecombo+combo_se,time_for_x);
plot(x,n,'Color',linecol);
xlabel('Time (sec)'); ylabel(['Reach rate (reaches per sec) ' num2str(ntrials) ' trials']);

end

function stats_compareChangeInRate(uncued_change,cued_change,uncued_changeLED,cued_changeLED)

ina=isnan(uncued_change) | isnan(uncued_changeLED);
uncued_change=uncued_change(~ina);
uncued_changeLED=uncued_changeLED(~ina);
ina=isnan(cued_change) | isnan(cued_changeLED);
cued_change=cued_change(~ina);
cued_changeLED=cued_changeLED(~ina);

if isempty(cued_changeLED)
    return
end

p=signrank(uncued_change,uncued_changeLED);
disp(['p of signrank for uncued change with and without LED: ' num2str(p)]);
p=signrank(cued_change,cued_changeLED);
disp(['p of signrank for cued change with and without LED: ' num2str(p)]);

end

function [sbys,mbym]=getMbyM_SbyS(sessIDandMouseID,rout,atleast_n_trials)

sbys.reachrates_noLED.trial1_alltrials_uncued=mean(rout.reachrates_noLED.trial1_alltrials_uncued,2,'omitnan');
sbys.reachrates_noLED.trial1_alltrials_cued=mean(rout.reachrates_noLED.trial1_alltrials_cued,2,'omitnan');
sbys.reachrates_noLED.alltrials_uncued=mean(rout.reachrates_noLED.alltrials_uncued,2,'omitnan');
sbys.reachrates_noLED.alltrials_cued=mean(rout.reachrates_noLED.alltrials_cued,2,'omitnan');
ns=sum(~isnan(rout.reachrates_noLED.trial1_alltrials_uncued),2,'omitnan');

mbym.reachrates_noLED.trial1_alltrials_uncued=getmbym(sessIDandMouseID,rout.reachrates_noLED.trial1_alltrials_uncued,atleast_n_trials);
mbym.reachrates_noLED.trial1_alltrials_cued=getmbym(sessIDandMouseID,rout.reachrates_noLED.trial1_alltrials_cued,atleast_n_trials);
mbym.reachrates_noLED.alltrials_uncued=getmbym(sessIDandMouseID,rout.reachrates_noLED.alltrials_uncued,atleast_n_trials);
mbym.reachrates_noLED.alltrials_cued=getmbym(sessIDandMouseID,rout.reachrates_noLED.alltrials_cued,atleast_n_trials);

sbys.reachrates_noLED.trial1_alltrials_uncued(ns<atleast_n_trials)=nan;
sbys.reachrates_noLED.trial1_alltrials_cued(ns<atleast_n_trials)=nan;
sbys.reachrates_noLED.alltrials_uncued(ns<atleast_n_trials)=nan;
sbys.reachrates_noLED.alltrials_cued(ns<atleast_n_trials)=nan;

if isempty(rout.reachrates_LED)
    sbys.reachrates_LED.trial1_alltrials_uncued=nan(size(sbys.reachrates_noLED.trial1_alltrials_uncued));
    sbys.reachrates_LED.trial1_alltrials_cued=nan(size(sbys.reachrates_noLED.trial1_alltrials_cued));
    sbys.reachrates_LED.alltrials_uncued=nan(size(sbys.reachrates_noLED.alltrials_uncued));
    sbys.reachrates_LED.alltrials_cued=nan(size(sbys.reachrates_noLED.alltrials_cued));

    mbym.reachrates_LED.trial1_alltrials_uncued=nan(size(mbym.reachrates_noLED.trial1_alltrials_uncued));
    mbym.reachrates_LED.trial1_alltrials_cued=nan(size(mbym.reachrates_noLED.trial1_alltrials_cued));
    mbym.reachrates_LED.alltrials_uncued=nan(size(mbym.reachrates_noLED.alltrials_uncued));
    mbym.reachrates_LED.alltrials_cued=nan(size(mbym.reachrates_noLED.alltrials_cued));
    return
end

sbys.reachrates_LED.trial1_alltrials_uncued=mean(rout.reachrates_LED.trial1_alltrials_uncued,2,'omitnan');
sbys.reachrates_LED.trial1_alltrials_cued=mean(rout.reachrates_LED.trial1_alltrials_cued,2,'omitnan');
sbys.reachrates_LED.alltrials_uncued=mean(rout.reachrates_LED.alltrials_uncued,2,'omitnan');
sbys.reachrates_LED.alltrials_cued=mean(rout.reachrates_LED.alltrials_cued,2,'omitnan');
ns=sum(~isnan(rout.reachrates_LED.trial1_alltrials_uncued),2,'omitnan');

mbym.reachrates_LED.trial1_alltrials_uncued=getmbym(sessIDandMouseID,rout.reachrates_LED.trial1_alltrials_uncued,atleast_n_trials);
mbym.reachrates_LED.trial1_alltrials_cued=getmbym(sessIDandMouseID,rout.reachrates_LED.trial1_alltrials_cued,atleast_n_trials);
mbym.reachrates_LED.alltrials_uncued=getmbym(sessIDandMouseID,rout.reachrates_LED.alltrials_uncued,atleast_n_trials);
mbym.reachrates_LED.alltrials_cued=getmbym(sessIDandMouseID,rout.reachrates_LED.alltrials_cued,atleast_n_trials);

sbys.reachrates_LED.trial1_alltrials_uncued(ns<atleast_n_trials)=nan;
sbys.reachrates_LED.trial1_alltrials_cued(ns<atleast_n_trials)=nan;
sbys.reachrates_LED.alltrials_uncued(ns<atleast_n_trials)=nan;
sbys.reachrates_LED.alltrials_cued(ns<atleast_n_trials)=nan;

end

function plotByBy(trial1_uncued,trial1_cued,trialn_uncued,trialn_cued,col,calcCued,plotLines)

if calcCued==true
    scatter(trialn_uncued-trial1_uncued,(trialn_cued-trialn_uncued)-(trial1_cued-trial1_uncued),[],col); hold on;
    scatter(mean(trialn_uncued-trial1_uncued,'all','omitnan'),mean((trialn_cued-trialn_uncued)-(trial1_cued-trial1_uncued),'all','omitnan'),80,col,'filled'); 
else
    scatter(trialn_uncued-trial1_uncued,trialn_cued-trial1_cued,[],col); hold on;
    scatter(mean(trialn_uncued-trial1_uncued,'all','omitnan'),mean(trialn_cued-trial1_cued,'all','omitnan'),80,col,'filled');
end
if plotLines==true
    for i=1:length(trialn_uncued)
        if calcCued==true
            line([0 trialn_uncued(i)-trial1_uncued(i)],[0 (trialn_cued(i)-trialn_uncued(i))-(trial1_cued(i)-trial1_uncued(i))],'Color',col); hold on;
        else
            line([0 trialn_uncued(i)-trial1_uncued(i)],[0 trialn_cued(i)-trial1_cued(i)],'Color',col); hold on;
        end
    end
    if calcCued==true
        line([0 mean(trialn_uncued-trial1_uncued,'all','omitnan')],[0 mean((trialn_cued-trialn_uncued)-(trial1_cued-trial1_uncued),'all','omitnan')],'Color',col,'LineWidth',4); hold on;
    else
        line([0 mean(trialn_uncued-trial1_uncued,'all','omitnan')],[0 mean(trialn_cued-trial1_cued,'all','omitnan')],'Color',col,'LineWidth',4); hold on;
    end
end

end

function rrout=getmbym(sessIDandMouseID,rr,atleast_n_trials)

u=unique(sessIDandMouseID(:,2));
rrout=nan(1,length(u));
for i=1:length(u)
    currmouse=u(i);
    temp=rr(sessIDandMouseID(:,2)==currmouse,:);
    % if fewer than atleast_n_trials not nan trials, drop this mouse
    nonan=nansum(~isnan(temp(1:end)));
    if nonan<atleast_n_trials
        continue
    end
    rrout(i)=mean(temp,'all','omitnan');
end

end

function [x,y]=readInFigData(figname)

% Open the .fig file
fig = openfig(figname,'reuse'); 

% Get the handles of all the child objects
axes_handles = get(fig,'Children');

li=get(axes_handles(1),'Children');
[x,y]=getquivers(li);

end

function [xdata,ydata]=getquivers(li)

% this format from kim's figs
% first quiver is average
xdata=[]; ydata=[];
for i=1:length(li)
    if ~isa(li(i), 'matlab.graphics.chart.primitive.Quiver')
        continue
    end
    if li(i).LineWidth==4
        % skip it, was av
        continue
    end
    xdata=[xdata li(i).UData];
    ydata=[ydata li(i).VData];
end

end

function plotMeAndSe(data1,data2,c,linewidth,suppressOutput,calcCued)
% make inputs vectors if they are not
data1=data1(1:end);
data2=data2(1:end);
if suppressOutput==false
    if calcCued==true
        % cued window is binomial bcz window small
        % use normal approximation
        % variance=np(1-p)
        p=nanmean(data2>0.01);
        n=nansum(~isnan(data2));
        cuedVar=p*(1-p); % var of av rather than sum
        cuedSD=sqrt(cuedVar);
        % scale factor to change back to rate
        scalefac=mean(data2,'all','omitnan')/p;
        cuedRateSD=cuedSD*scalefac;
        cuedRateSE=cuedRateSD./sqrt(n);
        line([nanmean(data1)-nanstd(data1,[],2)./sqrt(nansum(~isnan(data1))) nanmean(data1)+nanstd(data1,[],2)./sqrt(nansum(~isnan(data1)))],[nanmean(data2) nanmean(data2)],'Color',c,'LineWidth',linewidth);
        hold on;
        line([nanmean(data1) nanmean(data1)],[nanmean(data2)-cuedRateSE nanmean(data2)+cuedRateSE],'Color',c,'LineWidth',linewidth);
    else
        line([nanmean(data1)-nanstd(data1,[],2)./sqrt(nansum(~isnan(data1))) nanmean(data1)+nanstd(data1,[],2)./sqrt(nansum(~isnan(data1)))],[nanmean(data2) nanmean(data2)],'Color',c,'LineWidth',linewidth);
        hold on;
        line([nanmean(data1) nanmean(data1)],[nanmean(data2)-nanstd(data2,[],2)./sqrt(nansum(~isnan(data2))) nanmean(data2)+nanstd(data2,[],2)./sqrt(nansum(~isnan(data2)))],'Color',c,'LineWidth',linewidth);
    end
end

end

function getMeanAndBootstrapForCDF(timeStep,returnThisCDF,f)

startAtCue=false;
cuelengthadjust=0.25; % duration of cue
subtractPreCue=false;
startAtPrecue=true;
preCueWindow=[-2 -1]-cuelengthadjust;
cutCDFat=9-cuelengthadjust; % cut cdf at this time

maxTrile=cutCDFat; % make this empty if want whole trial length

% timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(returnThisCDF.trial1_rawReachMatrix,2)-1)*timeStep;
% find cue ind
% [~,f]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
if startAtPrecue==true
    cueTime=timeBinsForReaching(f);
    precuecutoff=cueTime+preCueWindow(1);

    temp=returnThisCDF.trial1_rawReachMatrix;
    temp(:,timeBinsForReaching<precuecutoff)=0;
    returnThisCDF.trial1_rawReachMatrix=temp;
    
    temp=returnThisCDF.trial2_rawReachMatrix;
    temp(:,timeBinsForReaching<precuecutoff)=0;
    returnThisCDF.trial2_rawReachMatrix=temp;
end
if startAtCue==true
    cueTime=timeBinsForReaching(f);
    preCueWindow=[cueTime+preCueWindow(1) cueTime+preCueWindow(2)];
else
    % make preCueWindow wrt cueTime=0
    cueTime=timeBinsForReaching(f);
    preCueWindow=[cueTime+preCueWindow(1) cueTime+preCueWindow(2)];
    cueTime=0;
end
[d1,d2]=plotCDF_rawReaches(returnThisCDF.trial1_rawReachMatrix,returnThisCDF.trial2_rawReachMatrix,timeBinsForReaching,cueTime,['CDF Raw Reaches: trial 1 (black) vs later (red)'],subtractPreCue,preCueWindow,maxTrile);

% Bootstrap CDFs
dat2=returnThisCDF.trial1_rawReachMatrix;
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs=nan(nRuns,length(timeBinsForReaching));
bootPDFs=nan(nRuns,length(timeBinsForReaching));
whichTrialsTaken=cell(1,nRuns);
for i=1:nRuns
    takeTheseForBoot=randi(size(dat2,1),1,takeIndsForBootstrap); % with replacement
    whichTrialsTaken{i}=takeTheseForBoot;
    sub_dat2=dat2(takeTheseForBoot,:);
    sub_dat2=sum(sub_dat2,1,'omitnan');
    if isempty(maxTrile)
    else
        sub_dat2(timeBinsForReaching>maxTrile)=0;
    end
    cond2_cdf=accumulateDistribution(sub_dat2);
    cond2_cdf=cond2_cdf./nanmax(cond2_cdf);
    bootCDFs(i,:)=cond2_cdf;
    bootPDFs(i,:)=mean(sub_dat2,1,'omitnan')./sum(mean(sub_dat2,1,'omitnan'),'all','omitnan');
end
% Show bootstrapped 95% CI
sorted_bootCDFs=nan(size(bootCDFs));
fifthPerc=nan(1,size(bootCDFs,2));
ninetyfifthPerc=nan(1,size(bootCDFs,2));
for i=1:size(bootCDFs,2)
    sorted_bootCDFs(:,i)=sort(bootCDFs(:,i));
    fifthPerc(i)=prctile(sorted_bootCDFs(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootCDFs(:,i),95);
end
plot(timeBinsForReaching,fifthPerc,'Color','k'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','k');

dat2=returnThisCDF.trial2_rawReachMatrix;
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs2nd=nan(nRuns,length(timeBinsForReaching));
bootPDFs2nd=nan(nRuns,length(timeBinsForReaching));
for i=1:nRuns
    % takeTheseForBoot=randi(size(dat2,1),1,takeIndsForBootstrap); % with replacement
    takeTheseForBoot=whichTrialsTaken{i};
    sub_dat2=dat2(takeTheseForBoot,:);
    sub_dat2=sum(sub_dat2,1,'omitnan');
    if isempty(maxTrile)
    else
        sub_dat2(timeBinsForReaching>maxTrile)=0;
    end
    cond2_cdf=accumulateDistribution(sub_dat2);
    cond2_cdf=cond2_cdf./nanmax(cond2_cdf);
    bootCDFs2nd(i,:)=cond2_cdf;
    bootPDFs2nd(i,:)=mean(sub_dat2,1,'omitnan')./sum(mean(sub_dat2,1,'omitnan'),'all','omitnan');
end
% Show bootstrapped 95% CI
sorted_bootCDFs=nan(size(bootCDFs2nd));
fifthPerc=nan(1,size(bootCDFs2nd,2));
ninetyfifthPerc=nan(1,size(bootCDFs2nd,2));
for i=1:size(bootCDFs2nd,2)
    sorted_bootCDFs(:,i)=sort(bootCDFs2nd(:,i));
    fifthPerc(i)=prctile(sorted_bootCDFs(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootCDFs(:,i),95);
end
plot(timeBinsForReaching,fifthPerc,'Color','r'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','r');

% Now get bootstrapped earth mover's function
earthmovers=bootCDFs2nd-bootCDFs;
% Show bootstrapped 95% CI
sorted_bootEM=nan(size(earthmovers));
fifthPerc=nan(1,size(earthmovers,2));
ninetyfifthPerc=nan(1,size(earthmovers,2));
for i=1:size(earthmovers,2)
    sorted_bootEM(:,i)=sort(earthmovers(:,i));
    fifthPerc(i)=prctile(sorted_bootEM(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootEM(:,i),95);
end
figure();
plot(timeBinsForReaching,nanmean(earthmovers,1),'Color','r'); hold on;
plot(timeBinsForReaching,fifthPerc,'Color','r'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','r');

% Bootstrapped difference in PDFs
pdfdiff=bootPDFs2nd-bootPDFs;
% Show bootstrapped 95% CI
sorted_bootPD=nan(size(pdfdiff));
fifthPerc=nan(1,size(pdfdiff,2));
ninetyfifthPerc=nan(1,size(pdfdiff,2));
for i=1:size(pdfdiff,2)
    sorted_bootPD(:,i)=sort(pdfdiff(:,i));
    fifthPerc(i)=prctile(sorted_bootPD(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootPD(:,i),95);
end
figure();
plot(timeBinsForReaching,nanmean(pdfdiff,1),'Color','r'); hold on;
plot(timeBinsForReaching,fifthPerc,'Color','r'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','r');
title('Diff between PDFs from bootstrap');

end


function [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(approach2_alltrials_uncued,approach2_alltrials_cued,colorForBootstrapPoints,scatterPointsEdgeColor,suppressOutput,calcCued)

stillsuppressbootstrap=false;

if isempty(approach2_alltrials_uncued) || isempty(approach2_alltrials_cued)
    uncued_mean_out=[];
    cued_mean_out=[];
    bootMeans=[];
    return
end

% can't directly subtract uncued from cued on a trial by trial
% basis bcz cued is basically a binomial, either 1 or 0
% need to get average rate before doing subtraction
% calculate true cued, because reaching in cued window is actually
% the sum of uncued and cued

% altogether_prob_cued=nanmean(approach2_alltrials_cued,2);
% altogether_prob_uncued=nanmean(approach2_alltrials_uncued,2);
altogether_prob_cued=approach2_alltrials_cued(1:end); % better to bootstrap across trials, not sessions, because in some sessions, mouse drops a lot
altogether_prob_uncued=approach2_alltrials_uncued(1:end);
takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
% nan trials are from cases where a session had fewer trials, just nanned
% to fill in matrix
% disp(['dropping this many trials because of nan ' num2str(nansum(~takeTrials))]);
altogether_prob_cued=altogether_prob_cued(takeTrials==1);
altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
% Show bootstrapped 95% CI
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
nRuns=100;
bootMeans=nan(2,nRuns);
for i=1:nRuns
    takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
    sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
    sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
    % can do subtraction here after calculating rate
    if calcCued==true
        bootMeans(1,i)=nanmean(sub_prob_uncued);
        bootMeans(2,i)=nanmean(sub_prob_cued)-nanmean(sub_prob_uncued);
    else
        bootMeans(1,i)=nanmean(sub_prob_uncued);
        bootMeans(2,i)=nanmean(sub_prob_cued);
    end
end
if suppressOutput==false && stillsuppressbootstrap==false
    s=scatter(bootMeans(1,:),bootMeans(2,:),20,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints); hold on;
    s.AlphaData = 0.5*ones(1,size(bootMeans,2));
    s.MarkerFaceAlpha = 'flat';
    if calcCued==true
        scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued)-nanmean(altogether_prob_uncued),50,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints);
    else
        scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints);
    end
end
if calcCued==true
    uncued_mean_out=nanmean(altogether_prob_uncued);
    cued_mean_out=nanmean(altogether_prob_cued)-nanmean(altogether_prob_uncued);
else
    uncued_mean_out=nanmean(altogether_prob_uncued);
    cued_mean_out=nanmean(altogether_prob_cued);
end

end

function [out1,out2]=plotCDF_rawReaches(data1,data2,timesteps,cueTime,tit,subtractPreCue,preCueWindow,maxTrialLength)

% make raw reaching data a timeseries (i.e., a pdf)
% then simply accumulate distribution, selecting only time points after the
% cue

% maxTrialLength=[]; %9.5; % in sec, or empty if want to include all reaches

data1=sum(data1,1,'omitnan');
data2=sum(data2,1,'omitnan');

if subtractPreCue==true
    data1=data1-mean(data1(timesteps>=preCueWindow(1) & timesteps<=preCueWindow(2)),2,'omitnan');
    data2=data2-mean(data2(timesteps>=preCueWindow(1) & timesteps<=preCueWindow(2)),2,'omitnan');
    data1(data1<0)=0;
    data2(data2<0)=0;
end
if ~isempty(maxTrialLength)
    data1(timesteps>maxTrialLength)=0;
    data2(timesteps>maxTrialLength)=0;
end

suppPlots=false;

if suppPlots==false
    cond1_cdf=accumulateDistribution(data1(timesteps>cueTime));
    figure();
    plot(timesteps(timesteps>cueTime),cond1_cdf./nanmax(cond1_cdf),'Color','k');
    xlabel('CDF');
    ylabel('Count');
    title(tit);
    out1=cond1_cdf./nanmax(cond1_cdf);
    
    hold on;
    cond2_cdf=accumulateDistribution(data2(timesteps>cueTime));
    plot(timesteps(timesteps>cueTime),cond2_cdf./nanmax(cond2_cdf),'Color','r');
    out2=cond2_cdf./nanmax(cond2_cdf);
end

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=sum(data(1:i),'all','omitnan');
end

end