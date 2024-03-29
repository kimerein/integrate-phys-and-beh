function tuningOutput=plotUnitSummariesAfterTCAlabels(varargin)

tuningOutput=[];
if length(varargin)==9
    groupLabelsFromTCA=varargin{1};
    cuez=varargin{2};
    cued_success_Response=varargin{3};
    cued_failure_Response=varargin{4};
    uncued_success_Response=varargin{5};
    uncued_failure_Response=varargin{6};
    isSig=varargin{7};
    doingCued=varargin{8};
    justAvsOrTuning=varargin{9};
    useTheseCuezBins=[];
elseif length(varargin)==10
    groupLabelsFromTCA=varargin{1};
    cuez=varargin{2};
    cued_success_Response=varargin{3};
    cued_failure_Response=varargin{4};
    uncued_success_Response=varargin{5};
    uncued_failure_Response=varargin{6};
    isSig=varargin{7};
    doingCued=varargin{8};
    justAvsOrTuning=varargin{9};
    useTheseCuezBins=varargin{10};
end

switch justAvsOrTuning
    case 'justAvs'
        basesubtract=false;
        individBase=true;
        basetimewindow=[-3.5 -2.5]; %[9 12.5]; %[4 9];

        plotAll=false;
        Zscore=false;
        minmaxnorm=false;
        normToBase=false; 
        chopOutliers=true;
        smoo=1; %30; % only used for tuning %30; %6; %smoo=3; %smoo=42;
        smoothBeforeResids=true; 
        smooBef=100; %15; %30; %83;
        getResiduals=false; % but need this to get rid of mid-range
        ds=1; %1;
        removeInsufficientBaseline=true; % will nan out units that don't have at least X seconds of baseline before aligncomp max
        atLeastXBaseline=3; %0.75; % in sec
        removeInsufficientPostBaseline=true;
        atLeastXAfterBaseline=10;
        addbeginnings=false;
        dropTrialBeginnings=false;

        dropTrialBeginnings=true;
        moreThanXSecBeforeAlign=1;
        addbeginnings=true;
        atLeastXBaseline=0;
    case 'tuning'
        % for cue tuned plots
        % doingCued='uncuedOverCued'; % 'cued' or 'uncued' or 'cuedOverUncued' or 'uncuedOverCued'
        basesubtract=false;
        individBase=false;
        basetimewindow=[9 16]; %[-5 -4]; %[9 12.5]; %[4 9];

        plotAll=false;
        Zscore=false;
        normToBase=false; % doesn't work well, too much noise
        minmaxnorm=false;
        smoo=1; %6; %smoo=3; %smoo=42;
        chopOutliers=true;
        smoothBeforeResids=true;
        smooBef=100; %150; %10; %30;
        getResiduals=false; % but need this to get rid of mid-range
        dsForCuez=1;
        removeInsufficientBaseline=true; % will nan out units that don't have at least X seconds of baseline before aligncomp max
        atLeastXBaseline=3; % in sec
        removeInsufficientPostBaseline=true;
        atLeastXAfterBaseline=10;
        addbeginnings=false;
        dropTrialBeginnings=false;
end

% plot all SU
% although ugly, the raw raw data actually shows effects (maybe for
% supplement)
% cuezbins=prctile(cuez,0:5:100); 
% cuezbins=prctile(cuez,[0:10:90 92 94 96 97 98 99 100]); cuezbins(1)=cuezbins(1)-0.0001; cuezbins(end)=cuezbins(end)+0.0001;

% temp=prctile(cuez(groupLabelsFromTCA==1),[0:25:75 80 85 90 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==2),[0:25:75 80 85 90 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==1),[0 12.5 22 30 45 72 85 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==2),[0 12.5 22 30 45 72 85 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;

%     temp=prctile(cuez(groupLabelsFromTCA==1),[0 39 50 72 85 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%     temp=prctile(cuez(groupLabelsFromTCA==2),[0 39 50 72 85 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2

switch doingCued
    case 'justAvs'
    case 'cued'
%         smooBef=10; %100;
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 39 50 72 85 90 94 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 39 50 72 85 90 94 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 50 72 78 85 90 94 98 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 50 72 78 85 90 94 98 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 25 40 55 70:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 25 40 55 70:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
        temp=prctile(cuez(groupLabelsFromTCA==1),[0 25 40 55:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
        temp=prctile(cuez(groupLabelsFromTCA==2),[0 25 40 55:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
    case 'uncued'
        temp=prctile(cuez(groupLabelsFromTCA==1),[0 42 50 72 83 88 96 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 42th prctile is 0 cuez for grp 1
        temp=prctile(cuez(groupLabelsFromTCA==2),[0 42 50 72 77 85 98 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
    case 'cuedOverUncued'
        basesubtract=false; % [0 4 10 15 50 85 96 100] [0 10 20 50 70 80 96 100] [0 10 20 50 70 80 96 100]
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 10 20 50 60 70 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 10 20 50 60 70 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
        temp=prctile(cuez(groupLabelsFromTCA==1),[0 15 30 45 55 65 80 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
        temp=prctile(cuez(groupLabelsFromTCA==2),[0 15 30 45 55 65 80 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 15 30 45 55 65 80 94 98 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 39th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 15 30 45 55 65 80 94 98 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
    case 'uncuedOverCued'
%         basesubtract=false;
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 10 20 50 60 70 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 10 20 50 60 70 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; 
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 6 12 50 88 94 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; temp=sort(unique(temp)); cuezbins{1}=temp;
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 6 12 50 88 94 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; temp=sort(unique(temp)); cuezbins{2}=temp; 
        temp=prctile(cuez(groupLabelsFromTCA==1),[0 50 100]); temp(2)=0; temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; temp=sort(unique(temp)); cuezbins{1}=temp;
        temp=prctile(cuez(groupLabelsFromTCA==2),[0 50 100]); temp(2)=0; temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; temp=sort(unique(temp)); cuezbins{2}=temp; 
        if ~isempty(useTheseCuezBins)
            cuezbins=useTheseCuezBins;
        end
%         cuezbins{1}=[-0.5 0.5 1.5];
%         cuezbins{2}=[-0.5 0.5 1.5];

%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 45 55 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; temp=sort(unique(temp)); cuezbins{2}=temp; 
%         temp=prctile(cuez(groupLabelsFromTCA==1),[0 17 22 40 60 82 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp; % 42th prctile is 0 cuez for grp 1
%         temp=prctile(cuez(groupLabelsFromTCA==2),[0 17 22 40 60 82 90 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp; % 28th prctile is 0 cuez for grp 2
end
% temp=prctile(cuez(groupLabelsFromTCA==1),[1:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==2),[1:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;

if isempty(isSig)
    isSig=ones(size(groupLabelsFromTCA));
end
groupLabelsFromTCA(isSig~=1)=-10; % omit these in plotting

% else plotBinRange
% cuezbins=prctile(cuez,[0 20 40 60 70 80 82 84 86 88 90 91 92 93 95 97 100]);
% plotAll=false;
% Zscore=false; % no because miss cued success activity for grp 2
% smoo=6;
% getResiduals=true; % but need this to get rid of mid-range

if dropTrialBeginnings==true
    cued_success_Response=dropTrialBeginningsFunc(cued_success_Response,moreThanXSecBeforeAlign);
    cued_failure_Response=dropTrialBeginningsFunc(cued_failure_Response,moreThanXSecBeforeAlign);
    uncued_success_Response=dropTrialBeginningsFunc(uncued_success_Response,moreThanXSecBeforeAlign);
    uncued_failure_Response=dropTrialBeginningsFunc(uncued_failure_Response,moreThanXSecBeforeAlign);
end

if addbeginnings==true
    cued_success_Response=addLastTrialToNextBeginning(cued_success_Response);
    cued_failure_Response=addLastTrialToNextBeginning(cued_failure_Response);
    uncued_success_Response=addLastTrialToNextBeginning(uncued_success_Response);
    uncued_failure_Response=addLastTrialToNextBeginning(uncued_failure_Response);
end

if removeInsufficientBaseline==true
   % find units with insufficient pre-event baseline
   [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=removeInsuffBase(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,atLeastXBaseline);
end

if removeInsufficientPostBaseline==true
   [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=removeInsufPostBase(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,atLeastXAfterBaseline);
end

if basesubtract==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=baseSubResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,basetimewindow,individBase);
end

if smoothBeforeResids==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=smoothResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,smooBef);
end

if chopOutliers==true
    aboveThisFR=100;
    [cued_success_Response,~,uncued_success_Response,~]=chopOuts(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,aboveThisFR,groupLabelsFromTCA==1);
    [~,cued_failure_Response,~,uncued_failure_Response]=chopOuts(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,aboveThisFR,groupLabelsFromTCA==2);
end

if getResiduals==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=getTrialTypeSpecificResiduals(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end

if Zscore==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=ZscoreResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
elseif minmaxnorm==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=maxNorm(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end

if normToBase==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=normByBaseWindow(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,basetimewindow);
end

switch justAvsOrTuning
    case 'justAvs'
        % ds=6;
        cued_success_Response.idx=groupLabelsFromTCA; exclu=cued_success_Response.excluded; gpLab=1;
        sub1=subResponse(cued_success_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(cued_success_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
        outCuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedSucc,'k','r'); title('Cued success');
        cued_failure_Response.idx=groupLabelsFromTCA; exclu=cued_failure_Response.excluded; gpLab=1;
        sub1=subResponse(cued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(cued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
        outCuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedFail,'k','r'); title('Cued failure');
        uncued_success_Response.idx=groupLabelsFromTCA; exclu=uncued_success_Response.excluded; gpLab=1;
        sub1=subResponse(uncued_success_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(uncued_success_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
        outUncuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedSucc,'k','r'); title('Uncued success');
        uncued_failure_Response.idx=groupLabelsFromTCA; exclu=uncued_failure_Response.excluded; gpLab=1;
        sub1=subResponse(uncued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(uncued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
        outUncuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedFail,'k','r'); title('Uncued failure');

        r.response1=outCuedSucc.response1; r.response2=outCuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Cued grp 1');
        r.response1=outCuedSucc.response2; r.response2=outCuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Cued grp 2');
        r.response1=outUncuedSucc.response1; r.response2=outUncuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Uncued grp 1');
        r.response1=outUncuedSucc.response2; r.response2=outUncuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Uncued grp 2');
        
        makeDirectionQuiverPlot(outCuedSucc,outCuedFail,outUncuedSucc,outUncuedFail,'response2',[1 5]);

    case 'tuning'
        % failure_off={'unit91onCh1_A2atagged','unit97onCh28_A2atagged','unit98onCh31_A2atagged','unit99onCh28_A2atagged','unit147onCh22_A2atagged','unit159onCh27_A2atagged','unit160onCh27_A2atagged','unit162onCh27_A2atagged','unit163onCh27_A2atagged','unit208onCh23_A2atagged','unit208onCh25_A2atagged','unit209onCh21_A2atagged','unit209onCh23_A2atagged','unit215onCh27_A2atagged','unit217onCh27_A2atagged','unit227onCh30_A2atagged','unit233onCh30_A2atagged'};
        % plotSU_contextAndOutcome('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_1\20210803\SU aligned to behavior',failure_off);

        % cuezbins=-2:0.5:3;
        % cuezbins(1)=-2.0001; cuezbins(end)=3.0001;
        % cuezbins=prctile(cuez,0:10:100);
        % cuezbins=prctile(cuez,[0 10 20 30 40 50 60 70 75 80 85 87.5 90 92.5 95 97.5 100]);
        % cuezbins=prctile(cuez,[0 20 40 60 70 80 82 84 86 88 90 91 92 93 95 97 100]);
        basesubtract=false;
        basetimewindow=[9 12.5];
        [grp1_succ,grp2_succ]=plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,'cued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
        [grp1_fail,grp2_fail]=plotByCuez(cued_failure_Response,cuez,groupLabelsFromTCA,'cued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
        [grp1_succ_uncue,grp2_succ_uncue]=plotByCuez(uncued_success_Response,cuez,groupLabelsFromTCA,'uncued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
        [grp1_fail_uncue,grp2_fail_uncue]=plotByCuez(uncued_failure_Response,cuez,groupLabelsFromTCA,'uncued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
        tuningOutput.grp1_succ=grp1_succ;
        tuningOutput.grp2_succ=grp2_succ;
        tuningOutput.grp1_fail=grp1_fail;
        tuningOutput.grp2_fail=grp2_fail;
        tuningOutput.grp1_succ_uncue=grp1_succ_uncue;
        tuningOutput.grp2_succ_uncue=grp2_succ_uncue;
        tuningOutput.grp1_fail_uncue=grp1_fail_uncue;
        tuningOutput.grp2_fail_uncue=grp2_fail_uncue;
        % plotDiffOfBycuez(grp1_succ,grp1_fail,[]);
        % plotDiffOfBycuez(grp1_succ,grp1_fail,[-1.5 0]);
        % plotDiffOfBycuez(grp1_succ,grp1_fail,[0 2.1]);
        % violinPlots(grp1_fail_uncue,[1 4]);
        % violinPlots(grp1_succ_uncue,[1 4]);
        % violinPlots(grp2_fail_uncue,[1 4]);

%         load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\dprimes_single_trials\cued failure vs uncued failure\pvals_for_curr_cued_success_Response.mat')
%         nonsigcuez=cuez; nonsigcuez(pvals_for_curr_cued_success_Response>0.05)=nan;
%         [grp1_succsig,grp2_succsig]=plotByCuez(cued_success_Response,nonsigcuez,groupLabelsFromTCA,'cued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
%         [grp1_failsig,grp2_failsig]=plotByCuez(cued_failure_Response,nonsigcuez,groupLabelsFromTCA,'cued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
%         [grp1_succ_uncuesig,grp2_succ_uncuesig]=plotByCuez(uncued_success_Response,nonsigcuez,groupLabelsFromTCA,'uncued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
%         [grp1_fail_uncuesig,grp2_fail_uncuesig]=plotByCuez(uncued_failure_Response,nonsigcuez,groupLabelsFromTCA,'uncued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);

%         disp('doing grp1succ vs grp1_succ_uncue');
%         plotOutsOverlayed(grp1_succ,grp1_succ_uncue); pause;
%         disp('doing grp1fail vs grp1_fail_uncue');
%         plotOutsOverlayed(grp1_fail,grp1_fail_uncue); pause;
%         disp('doing grp2succ vs grp2_succ_uncue');
%         plotOutsOverlayed(grp2_succ,grp2_succ_uncue); pause;
%         disp('doing grp2fail vs grp2_fail_uncue');
%         plotOutsOverlayed(grp2_fail,grp2_fail_uncue);
end

end

function makeDirectionQuiverPlot(outCuedSucc,outCuedFail,outUncuedSucc,outUncuedFail,whichResp,plotWindow)

if ~strcmp(whichResp,'response1')
    [~,mi]=nanmin(abs(outCuedFail.response2.unittimes-plotWindow(1)));
    plotInds(1)=mi;
    [~,mi]=nanmin(abs(outCuedFail.response2.unittimes-plotWindow(2)));
    plotInds(2)=mi;
    dir1_av=nanmean(outCuedSucc.response2.me(plotInds(1):plotInds(2)));
    dir1_ubyu=nanmean(outCuedSucc.response2.unitbyunit(:,plotInds(1):plotInds(2)),2);
    dir2_av=nanmean(outUncuedSucc.response2.me(plotInds(1):plotInds(2)));
    dir2_ubyu=nanmean(outUncuedSucc.response2.unitbyunit(:,plotInds(1):plotInds(2)),2);
    dir3_av=nanmean(outCuedFail.response2.me(plotInds(1):plotInds(2)));
    dir3_ubyu=nanmean(outCuedFail.response2.unitbyunit(:,plotInds(1):plotInds(2)),2);
    dir4_av=nanmean(outUncuedFail.response2.me(plotInds(1):plotInds(2)));
    dir4_ubyu=nanmean(outUncuedFail.response2.unitbyunit(:,plotInds(1):plotInds(2)),2);

    figure();
    quiver(0,0,0,dir1_av,'Color',[0.5 0.5 0.5]); hold on;
    quiver(0,0,dir2_av,0,'Color',[0.5 0.5 0.5]);
    quiver(0,0,-dir3_av,0,'Color',[0.5 0.5 0.5]);
    quiver(0,0,0,-dir4_av,'Color',[0.5 0.5 0.5]);
    quiver(0,0,[dir2_av-dir3_av],[dir1_av-dir4_av],'Color','k');
    for i=1:size(dir1_ubyu)
        scatter([dir2_ubyu(i,:)-dir3_ubyu(i,:)],[dir1_ubyu(i,:)-dir4_ubyu(i,:)],[],'k');
    end
    xlabel('Uncued success positive, cued failure negative');
    ylabel('Cued success positive, uncued failure negative');
    return
end
[~,mi]=nanmin(abs(outCuedFail.response1.unittimes-plotWindow(1)));
plotInds(1)=mi;
[~,mi]=nanmin(abs(outCuedFail.response1.unittimes-plotWindow(2)));
plotInds(2)=mi;
dir1_av=nanmean(outCuedSucc.response1.me(plotInds(1):plotInds(2)));
dir1_ubyu=nanmean(outCuedSucc.response1.unitbyunit(:,plotInds(1):plotInds(2)),2);
dir2_av=nanmean(outUncuedSucc.response1.me(plotInds(1):plotInds(2)));
dir2_ubyu=nanmean(outUncuedSucc.response1.unitbyunit(:,plotInds(1):plotInds(2)),2);
dir3_av=nanmean(outCuedFail.response1.me(plotInds(1):plotInds(2)));
dir3_ubyu=nanmean(outCuedFail.response1.unitbyunit(:,plotInds(1):plotInds(2)),2);
dir4_av=nanmean(outUncuedFail.response1.me(plotInds(1):plotInds(2)));
dir4_ubyu=nanmean(outUncuedFail.response1.unitbyunit(:,plotInds(1):plotInds(2)),2);

figure();
quiver(0,0,0,dir1_av,'Color',[0.5 0.5 0.5]); hold on;
quiver(0,0,dir2_av,0,'Color',[0.5 0.5 0.5]);
quiver(0,0,-dir3_av,0,'Color',[0.5 0.5 0.5]);
quiver(0,0,0,-dir4_av,'Color',[0.5 0.5 0.5]);
quiver(0,0,[dir2_av-dir3_av],[dir1_av-dir4_av],'Color','k');
for i=1:size(dir1_ubyu)
    scatter([dir2_ubyu(i,:)-dir3_ubyu(i,:)],[dir1_ubyu(i,:)-dir4_ubyu(i,:)],[],'k');
end
xlabel('Uncued success positive, cued failure negative'); 
ylabel('Cued success positive, uncued failure negative'); 

% figure();
% for i=1:size(dir1_ubyu)
%     s=scatter([dir2_ubyu(i,:)-dir3_ubyu(i,:)],[dir1_ubyu(i,:)-dir4_ubyu(i,:)],[],'filled','k'); hold on;
%     s.AlphaData=0.2;
%     s.MarkerFaceAlpha='flat';
% end

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=removeInsufPostBase(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,atLeastXBaseline)

shouldExclude=zeros(size(cued_success_Response.unitbyunit_y,1),1);

currR=cued_success_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
[~,mi]=nanmin(abs(time-timeOfMaxAlignComp));
for i=1:size(currR.unitbyunit_y,1)
    lastNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'last');
    lastNonNanTime=time(lastNonNan);
    if lastNonNanTime-timeOfMaxAlignComp<atLeastXBaseline
        shouldExclude(i)=1;
    end
    % if isnan at timeOfMaxAlignComp
    if isnan(currR.unitbyunit_y(i,mi))
        shouldExclude(i)=1;
    end
end

currR=cued_failure_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
[~,mi]=nanmin(abs(time-timeOfMaxAlignComp));
for i=1:size(currR.unitbyunit_y,1)
    lastNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'last');
    lastNonNanTime=time(lastNonNan);
    if lastNonNanTime-timeOfMaxAlignComp<atLeastXBaseline
        shouldExclude(i)=1;
    end
    % if isnan at timeOfMaxAlignComp
    if isnan(currR.unitbyunit_y(i,mi))
        shouldExclude(i)=1;
    end
end

currR=uncued_success_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
[~,mi]=nanmin(abs(time-timeOfMaxAlignComp));
for i=1:size(currR.unitbyunit_y,1)
    lastNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'last');
    lastNonNanTime=time(lastNonNan);
    if lastNonNanTime-timeOfMaxAlignComp<atLeastXBaseline
        shouldExclude(i)=1;
    end
    % if isnan at timeOfMaxAlignComp
    if isnan(currR.unitbyunit_y(i,mi))
        shouldExclude(i)=1;
    end
end

currR=uncued_failure_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
[~,mi]=nanmin(abs(time-timeOfMaxAlignComp));
for i=1:size(currR.unitbyunit_y,1)
    lastNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'last');
    lastNonNanTime=time(lastNonNan);
    if lastNonNanTime-timeOfMaxAlignComp<atLeastXBaseline
        shouldExclude(i)=1;
    end
    % if isnan at timeOfMaxAlignComp
    if isnan(currR.unitbyunit_y(i,mi))
        shouldExclude(i)=1;
    end
end

% remove units with insufficient baseline from all responses
for i=1:length(shouldExclude)
    if shouldExclude(i)==1
        cued_success_Response.unitbyunit_y(i,:)=nan(size(cued_success_Response.unitbyunit_y(i,:)));
        cued_failure_Response.unitbyunit_y(i,:)=nan(size(cued_failure_Response.unitbyunit_y(i,:)));
        uncued_success_Response.unitbyunit_y(i,:)=nan(size(uncued_success_Response.unitbyunit_y(i,:)));
        uncued_failure_Response.unitbyunit_y(i,:)=nan(size(uncued_failure_Response.unitbyunit_y(i,:)));
    end
end

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=removeInsuffBase(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,atLeastXBaseline)

shouldExclude=zeros(size(cued_success_Response.unitbyunit_y,1),1);

currR=cued_success_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
for i=1:size(currR.unitbyunit_y,1)
    firstNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'first');
    firstNonNanTime=time(firstNonNan);
    if timeOfMaxAlignComp-firstNonNanTime<atLeastXBaseline
        shouldExclude(i)=1;
    end
end

currR=cued_failure_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
for i=1:size(currR.unitbyunit_y,1)
    firstNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'first');
    firstNonNanTime=time(firstNonNan);
    if timeOfMaxAlignComp-firstNonNanTime<atLeastXBaseline
        shouldExclude(i)=1;
    end
end

currR=uncued_success_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
for i=1:size(currR.unitbyunit_y,1)
    firstNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'first');
    firstNonNanTime=time(firstNonNan);
    if timeOfMaxAlignComp-firstNonNanTime<atLeastXBaseline
        shouldExclude(i)=1;
    end
end

currR=uncued_failure_Response;
time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
for i=1:size(currR.unitbyunit_y,1)
    firstNonNan=find(~isnan(currR.unitbyunit_y(i,:)),1,'first');
    firstNonNanTime=time(firstNonNan);
    if timeOfMaxAlignComp-firstNonNanTime<atLeastXBaseline
        shouldExclude(i)=1;
    end
end

% remove units with insufficient baseline from all responses
for i=1:length(shouldExclude)
    if shouldExclude(i)==1
        cued_success_Response.unitbyunit_y(i,:)=nan(size(cued_success_Response.unitbyunit_y(i,:)));
        cued_failure_Response.unitbyunit_y(i,:)=nan(size(cued_failure_Response.unitbyunit_y(i,:)));
        uncued_success_Response.unitbyunit_y(i,:)=nan(size(uncued_success_Response.unitbyunit_y(i,:)));
        uncued_failure_Response.unitbyunit_y(i,:)=nan(size(uncued_failure_Response.unitbyunit_y(i,:)));
    end
end

end

function plotOutsOverlayed(out1,out2)

baseSubDiffers=false;
forvio_timewindow=[2 4.5]; %[2 5];
smoothbeforediff=false;
smoothbeforebin=50;

if smoothbeforediff==true
    for i=1:length(out1.allunits)
        data1=out1.allunits{i};
        data2=out2.allunits{i};
        disp('smoo data1');
        for j=1:size(data1,1)
            if mod(j,100)==0
                disp(['j: ' num2str(j)]);
            end
            data1(j,:)=smooth(data1(j,:),smoothbeforebin,'moving');
        end
        disp('smoo data2');
        for j=1:size(data2,1)
            if mod(j,100)==0
                disp(['j: ' num2str(j)]);
            end
            data2(j,:)=smooth(data2(j,:),smoothbeforebin,'moving');
        end
        out1.allunits{i}=data1;
        out2.allunits{i}=data2;
    end
end

cmap=getCmapWithRed(1:length(out2.allunits)+1); hold on;
differs=cell(1,length(out1.allunits));
differstimes=cell(1,length(out1.allunits));
all_differs=[];
all_cuez=[];
whichgp=[];
for i=1:length(out1.allunits)
%     figure();
    data1=out1.allunits{i};
    data2=out2.allunits{i};
%     plot(out1.time{i},nanmean(data1,1),'Color',cmap(i,:));
%     hold on;
%     plot(out2.time{i},nanmean(data2,1),'Color',cmap(i,:));
%     scatter(out2.time{i},nanmean(data2,1),4,cmap(i,:));
    si=min(size(data1,2),size(data2,2));
%     differs{i}=nanmean(data1(:,1:si),1)-nanmean(data2(:,1:si),1);
    differs{i}=data1(:,1:si)-data2(:,1:si); temp=differs{i};
    t1=out1.time{i}; 
    all_differs=[all_differs; nanmean(temp(:,t1(1:si)>=forvio_timewindow(1) & t1(1:si)<forvio_timewindow(2)),2)];
    all_cuez=[all_cuez; out1.cuez{i}];
    whichgp=[whichgp; ones(size(out1.cuez{i})).*i];
    if baseSubDiffers==true
        temp=differs{i};
%         base=nanmean(nanmean(temp(:,out1.time{i}>-5 & out1.time{i}<-3),2),1);
%         base=nanmean(nanmean(temp(:,out1.time{i}>-5 & out1.time{i}<-3),2),1);
        base=nanmean(nanmean(temp(:,out1.time{i}>7 & out1.time{i}<12),2),1);
        differs{i}=differs{i}-base;
%         base=nanmean(temp(:,out1.time{i}>-5 & out1.time{i}<-3),2);
%         differs{i}=differs{i}-repmat(base,1,size(differs{i},2));
    end
    differstimes{i}=t1(1:si);
end

figure(); scatter(all_cuez,all_differs);
disp(['THIS IS SIGNRANK: ' num2str(signrank(all_differs))]);
try
    disp(['THIS IS SIGNRANK set 1: ' num2str(signrank(all_differs(whichgp==1)))]);
    disp(['THIS IS SIGNRANK set 2: ' num2str(signrank(all_differs(whichgp==2)))]);
catch
end

% temp=differs{1};
% base=nanmean(temp(:,differstimes{1}>-5 & differstimes{1}<-3),2);
% figure(); plot(differstimes{1},(temp-repmat(base,1,size(temp,2)))','Color',cmap(1,:));

figure();
smoo=10;
% forvio_timewindow=[smoo*0.06 4];
% forvio_timewindow=[1 4];
forvio=cell(1,2);
k=1;
alldatas=[];
alldatas_labels=[];
excludeForVio=zeros(length(out1.allunits),1);
for i=1:length(out1.allunits)
    if ~isempty(smoo)
        toplot=smooth(nanmean(differs{i},1),smoo);
        se_plus_toplot=smooth(nanmean(differs{i},1)+nanstd(differs{i},[],1)./sqrt(size(differs{i},1)),smoo);
        se_minus_toplot=smooth(nanmean(differs{i},1)-nanstd(differs{i},[],1)./sqrt(size(differs{i},1)),smoo);
    end
    plot(differstimes{i},toplot,'Color',cmap(i,:)); hold on;
    plot(differstimes{i},se_plus_toplot,'Color',cmap(i,:));
    plot(differstimes{i},se_minus_toplot,'Color',cmap(i,:));
%     if i==1 || i==length(out1.allunits)
%         temp=differs{i};
%         forvio{k}=nanmean(temp(:,differstimes{i}>=forvio_timewindow(1) & differstimes{i}<forvio_timewindow(2)),2);
%         k=k+1;
%     end
    temp=differs{i};
    forvio{i}=nanmean(temp(:,differstimes{i}>=forvio_timewindow(1) & differstimes{i}<forvio_timewindow(2)),2);
    if isempty(forvio{i})
        excludeForVio(i)=1;
    end
    alldatas=[alldatas; forvio{i}];
    alldatas_labels=[alldatas_labels; i.*ones(size(forvio{i}))];
end

% figure();
% violin(forvio,'facecolor',cmap([1 end],:),'medc',[],'facealpha',1,'bw',0.15);
% p=ranksum(forvio{1},forvio{2});
% disp(['pval from ranksum is ' num2str(p)]);

figure();
violin(forvio(excludeForVio==0),'facecolor',cmap(excludeForVio==0,:),'medc',[],'facealpha',1); %,'bw',0.4);
try
    [p,tbl,stats]=anova1(alldatas,alldatas_labels);
    results=multcompare(stats);
    results_tbl = array2table(results,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
end

end

function violinPlots(out,inTimeWindow)

figure();
dataforviolin=cell(1,length(out.allunits));
alldatas=[];
alldatas_labels=[];
for i=1:length(out.allunits)
    data=out.allunits{i};
    dataforviolin{i}=nanmean(data(:,out.time{i}>=inTimeWindow(1) & out.time{i}<=inTimeWindow(2)),2);
    alldatas=[alldatas; dataforviolin{i}];
    alldatas_labels=[alldatas_labels; i.*ones(size(dataforviolin{i}))];
end
cmap=getCmapWithRed(1:length(out.allunits)+1); hold on;
% violin(dataforviolin,'facecolor',cmap,'edgecolor','none','bw',0.1,'mc','k','medc','r-.')
violin(dataforviolin,'facecolor',cmap,'edgecolor','none','medc',[],'facealpha',1);
ylabel('Firing rate','FontSize',14);
% hold on;
% plot all points
% for i=1:length(out.allunits)
%     scatter(i+(rand(size(dataforviolin{i}))./8)-0.125/2,dataforviolin{i},2,cmap(i,:),'filled');
% end

[p,tbl,stats]=anova1(alldatas,alldatas_labels);
results=multcompare(stats);
results_tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

end

function currR=dropTrialBeginningsFunc(currR,moreThanXSecBeforeAlign)

time=nanmean(currR.unitbyunit_x,1);
temp=nanmean(currR.aligncomp_x,1);
[~,ma]=nanmax(nanmean(currR.aligncomp_y,1));
timeOfMaxAlignComp=temp(ma);
xsecbefore=timeOfMaxAlignComp-moreThanXSecBeforeAlign;
f=find(time<xsecbefore,1,'last');
faligncomp=find(nanmean(currR.aligncomp_x,1)<xsecbefore,1,'last');
temp=currR.unitbyunit_y;
currR.unitbyunit_y=temp(:,f+1:end);
temp=currR.unitbyunit_x;
currR.unitbyunit_x=temp(:,f+1:end);
temp=currR.aligncomp_x;
currR.aligncomp_x=temp(:,faligncomp+1:end);
temp=currR.aligncomp_y;
currR.aligncomp_y=temp(:,faligncomp+1:end);

end

function resp=addLastTrialToNextBeginning(resp)

trialEnds=[13.5 18];
times=nanmean(resp.unitbyunit_x,1);
batch1begins=resp.unitbyunit_y(:,times>trialEnds(1)-6.5 & times<trialEnds(1));
batch2begins=resp.unitbyunit_y(:,times>trialEnds(2)-6.5 & times<trialEnds(2));
% batch2begins=[];
togbatch=cat(3,batch1begins,batch2begins);
togbatch=reshape(nanmean(togbatch,3),size(batch1begins,1),size(batch1begins,2));
timesbatch=times(times>trialEnds(1)-6.5 & times<trialEnds(1));
timesbatch=timesbatch-nanmax(timesbatch);
timesbatch=timesbatch+nanmin(times);
resp.unitbyunit_x=[repmat(timesbatch,size(resp.unitbyunit_x,1),1) resp.unitbyunit_x];
resp.unitbyunit_y=[togbatch resp.unitbyunit_y];

end

function plotDiffOfBycuez(data1,data2,zeroAtWindow)

figure();
cmap=getCmapWithRed(1:length(data1.bycuez)+1);
if length(data1.bycuez{1})>length(data2.bycuez{1})
    for i=1:length(data1.bycuez)
        temp=data1.bycuez{i};
        data1.bycuez{i}=temp(1:length(data2.bycuez{1}));
        temp=data1.time{i};
        data1.time{i}=temp(1:length(data2.bycuez{1}));
    end
elseif length(data1.bycuez{1})<length(data2.bycuez{1})
    for i=1:length(data2.bycuez)
        temp=data2.bycuez{i};
        data2.bycuez{i}=temp(1:length(data1.bycuez{1}));
        temp=data2.time{i};
        data2.time{i}=temp(1:length(data1.bycuez{1}));
    end
end
for i=1:length(data1.bycuez)
    if ~isempty(zeroAtWindow)
        temp=data2.bycuez{i}-data1.bycuez{i};
        base=nanmean(temp(data1.time{i}>zeroAtWindow(1) & data1.time{i}<zeroAtWindow(2)));
        plot(data1.time{i},temp-base,'Color',cmap(i,:)); hold on;
    else
        plot(data1.time{i},data2.bycuez{i}-data1.bycuez{i},'Color',cmap(i,:)); hold on;
    end
end

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=normByBaseWindow(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,baseWindow)

[cuetime,t]=getCueTimeForThisResponse(cued_success_Response);
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y./repmat(mean(cued_success_Response.unitbyunit_y(:,t>baseWindow(1) & t<baseWindow(2)),2,'omitnan'),1,size(cued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(cued_failure_Response);
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y./repmat(mean(cued_failure_Response.unitbyunit_y(:,t>baseWindow(1) & t<baseWindow(2)),2,'omitnan'),1,size(cued_failure_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_success_Response);
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y./repmat(mean(uncued_success_Response.unitbyunit_y(:,t>baseWindow(1) & t<baseWindow(2)),2,'omitnan'),1,size(uncued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_failure_Response);
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y./repmat(mean(uncued_failure_Response.unitbyunit_y(:,t>baseWindow(1) & t<baseWindow(2)),2,'omitnan'),1,size(uncued_failure_Response.unitbyunit_y,2));
      
end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=maxNorm(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

minmaxwindow=[-2 5]; % in sec wrt cue
[cuetime,t]=getCueTimeForThisResponse(cued_success_Response);
cued_success_Response.unitbyunit_y(cued_success_Response.unitbyunit_y<0.001)=0;
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y./repmat(max(cued_success_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(cued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(cued_failure_Response);
cued_failure_Response.unitbyunit_y(cued_failure_Response.unitbyunit_y<0.001)=0;
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y./repmat(max(cued_failure_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(cued_failure_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_success_Response);
uncued_success_Response.unitbyunit_y(uncued_success_Response.unitbyunit_y<0.001)=0;
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y./repmat(max(uncued_success_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(uncued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_failure_Response);
uncued_failure_Response.unitbyunit_y(uncued_failure_Response.unitbyunit_y<0.001)=0;
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y./repmat(max(uncued_failure_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(uncued_failure_Response.unitbyunit_y,2));
      
end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=chopOuts(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,aboveThisFR,justForThese)

backup=cued_success_Response.unitbyunit_y; cued_success_Response.unitbyunit_y(cued_success_Response.unitbyunit_y>aboveThisFR)=aboveThisFR; cued_success_Response.unitbyunit_y(~justForThese,:)=backup(~justForThese,:);
backup=cued_failure_Response.unitbyunit_y; cued_failure_Response.unitbyunit_y(cued_failure_Response.unitbyunit_y>aboveThisFR)=aboveThisFR; cued_failure_Response.unitbyunit_y(~justForThese,:)=backup(~justForThese,:);
backup=uncued_success_Response.unitbyunit_y; uncued_success_Response.unitbyunit_y(uncued_success_Response.unitbyunit_y>aboveThisFR)=aboveThisFR; uncued_success_Response.unitbyunit_y(~justForThese,:)=backup(~justForThese,:);
backup=uncued_failure_Response.unitbyunit_y; uncued_failure_Response.unitbyunit_y(uncued_failure_Response.unitbyunit_y>aboveThisFR)=aboveThisFR; uncued_failure_Response.unitbyunit_y(~justForThese,:)=backup(~justForThese,:);

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=smoothResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,smoo)

doGauss=false;
if doGauss==true
    temp=cued_success_Response.unitbyunit_y;
    for i=1:size(temp,1)
        temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
    end
    cued_success_Response.unitbyunit_y=temp;

    temp=cued_failure_Response.unitbyunit_y;
    for i=1:size(temp,1)
        temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
    end
    cued_failure_Response.unitbyunit_y=temp;

    temp=uncued_success_Response.unitbyunit_y;
    for i=1:size(temp,1)
        temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
    end
    uncued_success_Response.unitbyunit_y=temp;

    temp=uncued_failure_Response.unitbyunit_y;
    for i=1:size(temp,1)
        temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
    end
    uncued_failure_Response.unitbyunit_y=temp;
else
    temp=cued_success_Response.unitbyunit_y;
    for i=1:size(temp,1)
        tempie=smooth(temp(i,:)',smoo,'moving')';
        f=find(~isnan(tempie),1,'last');
        tempie(f-smoo:end)=nan;
        f=find(~isnan(tempie),1,'first');
        tempie(1:f+smoo)=nan;
        temp(i,:)=tempie;
    end
    cued_success_Response.unitbyunit_y=temp;
    disp('smooth 1 of 4 done');

    temp=cued_failure_Response.unitbyunit_y;
    for i=1:size(temp,1)
        tempie=smooth(temp(i,:)',smoo,'moving')';
        f=find(~isnan(tempie),1,'last');
        tempie(f-smoo:end)=nan;
        f=find(~isnan(tempie),1,'first');
        tempie(1:f+smoo)=nan;
        temp(i,:)=tempie;
    end
    cued_failure_Response.unitbyunit_y=temp;
    disp('smooth 2 of 4 done');

    temp=uncued_success_Response.unitbyunit_y;
    for i=1:size(temp,1)
        tempie=smooth(temp(i,:)',smoo,'moving')';
        f=find(~isnan(tempie),1,'last');
        tempie(f-smoo:end)=nan;
        f=find(~isnan(tempie),1,'first');
        tempie(1:f+smoo)=nan;
        temp(i,:)=tempie;
    end
    uncued_success_Response.unitbyunit_y=temp;
    disp('smooth 3 of 4 done');

    temp=uncued_failure_Response.unitbyunit_y;
    for i=1:size(temp,1)
        tempie=smooth(temp(i,:)',smoo,'moving')';
        f=find(~isnan(tempie),1,'last');
        tempie(f-smoo:end)=nan;
        f=find(~isnan(tempie),1,'first');
        tempie(1:f+smoo)=nan;
        temp(i,:)=tempie;
    end
    uncued_failure_Response.unitbyunit_y=temp;
    disp('smooth 4 of 4 done');
end

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=baseSubResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,basewindow,individBase)

% individBase=false;

% basewindow=[5 9]; % in sec wrt cue
t=nanmean(cued_success_Response.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(cued_success_Response.aligncomp_y,1));
talign=nanmean(cued_success_Response.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;
t2=nanmean(cued_failure_Response.unitbyunit_x,1); t2=t2-cuetime;
t3=nanmean(uncued_success_Response.unitbyunit_x,1); t3=t3-cuetime;
t4=nanmean(uncued_failure_Response.unitbyunit_x,1); t4=t4-cuetime;
if individBase==false
    base=mean([cued_success_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)) cued_failure_Response.unitbyunit_y(:,t2>basewindow(1) & t2<basewindow(2)) uncued_success_Response.unitbyunit_y(:,t3>basewindow(1) & t3<basewindow(2)) uncued_failure_Response.unitbyunit_y(:,t4>basewindow(1) & t4<basewindow(2))],2,'omitnan');
    cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y-repmat(base,1,size(cued_success_Response.unitbyunit_y,2));
    cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y-repmat(base,1,size(cued_failure_Response.unitbyunit_y,2));
    uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y-repmat(base,1,size(uncued_success_Response.unitbyunit_y,2));
    uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y-repmat(base,1,size(uncued_failure_Response.unitbyunit_y,2));
else
    base=mean(cued_success_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)),2,'omitnan');
    cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y-repmat(base,1,size(cued_success_Response.unitbyunit_y,2));
    base=mean(cued_failure_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)),2,'omitnan');
    cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y-repmat(base,1,size(cued_failure_Response.unitbyunit_y,2));
    base=mean(uncued_success_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)),2,'omitnan');
    uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y-repmat(base,1,size(uncued_success_Response.unitbyunit_y,2));
    base=mean(uncued_failure_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)),2,'omitnan');
    uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y-repmat(base,1,size(uncued_failure_Response.unitbyunit_y,2));
end

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=ZscoreResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

sdwindow=[4 9]; % in sec wrt cue
t=nanmean(cued_success_Response.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(cued_success_Response.aligncomp_y,1));
talign=nanmean(cued_success_Response.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;
t2=nanmean(cued_failure_Response.unitbyunit_x,1); t2=t2-cuetime;
t3=nanmean(uncued_success_Response.unitbyunit_x,1); t3=t3-cuetime;
t4=nanmean(uncued_failure_Response.unitbyunit_x,1); t4=t4-cuetime;
sd=std([cued_success_Response.unitbyunit_y(:,t>sdwindow(1) & t<sdwindow(2)) cued_failure_Response.unitbyunit_y(:,t2>sdwindow(1) & t2<sdwindow(2)) uncued_success_Response.unitbyunit_y(:,t3>sdwindow(1) & t3<sdwindow(2)) uncued_failure_Response.unitbyunit_y(:,t4>sdwindow(1) & t4<sdwindow(2))],[],2,'omitnan');
sd(sd<0.01)=0.01;
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y./repmat(sd,1,size(cued_success_Response.unitbyunit_y,2));
cued_success_Response.unitbyunit_y(cued_success_Response.unitbyunit_y>60)=60;
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y./repmat(sd,1,size(cued_failure_Response.unitbyunit_y,2));
cued_failure_Response.unitbyunit_y(cued_failure_Response.unitbyunit_y>60)=60;
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y./repmat(sd,1,size(uncued_success_Response.unitbyunit_y,2));
uncued_success_Response.unitbyunit_y(uncued_success_Response.unitbyunit_y>60)=60;
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y./repmat(sd,1,size(uncued_failure_Response.unitbyunit_y,2));
uncued_failure_Response.unitbyunit_y(uncued_failure_Response.unitbyunit_y>60)=60;

end

function [cuetime,t]=getCueTimeForThisResponse(r1)

t=nanmean(r1.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(r1.aligncomp_y,1));
talign=nanmean(r1.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;

end

function [r1,r2,r3,r4]=getTrialTypeSpecificResiduals(r1,r2,r3,r4)

smoothBeforeNonspecific=true;
smooBefore=200;

[cuetime,t]=getCueTimeForThisResponse(r1);

% assumes same binning for all
timesteps=[mode(diff(nanmean(r1.unitbyunit_x,1))) mode(diff(nanmean(r2.unitbyunit_x,1))) mode(diff(nanmean(r3.unitbyunit_x,1))) mode(diff(nanmean(r4.unitbyunit_x,1)))];
if ~all(diff(timesteps)<0.001)
    error('Binning must match in all responses');
end

resp1=r1.unitbyunit_y(:,t>=-3 & t<=9); 
[cuetime,t]=getCueTimeForThisResponse(r2); resp2=r2.unitbyunit_y(:,t>=-3 & t<=9.01);
[cuetime,t]=getCueTimeForThisResponse(r3); resp3=r3.unitbyunit_y(:,t>=-3 & t<=9.01);
[cuetime,t]=getCueTimeForThisResponse(r4); resp4=r4.unitbyunit_y(:,t>=-3 & t<=9.01);
resp2=resp2(:,1:size(resp1,2));
resp3=resp3(:,1:size(resp1,2));
resp4=resp4(:,1:size(resp1,2));
allResp=nan(size(resp1,1),size(resp1,2),4);
allResp(:,:,1)=resp1; allResp(:,:,2)=resp2; allResp(:,:,3)=resp3; allResp(:,:,4)=resp4;

forNonspec=allResp;
if smoothBeforeNonspecific==true
    for i=1:size(forNonspec,3)
        for j=1:size(forNonspec,1)
            forNonspec(j,:,i)=smoothdata(reshape(forNonspec(j,:,i),1,size(forNonspec,2)),'gaussian',smooBefore);
        end
    end
end
nonspecific=reshape(min(noNansOrInfs(forNonspec),[],3,'omitnan'),size(forNonspec,1),size(forNonspec,2));
for i=1:4
    allResp(:,:,i)=allResp(:,:,i)-nonspecific;
end

% put back
[cuetime,t]=getCueTimeForThisResponse(r1); r1.unitbyunit_y(:,t>=-3 & t<=9)=reshape(allResp(:,:,1),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r2); f=find(t>=-3,1,'first'); r2.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,2),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r3); f=find(t>=-3,1,'first'); r3.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,3),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r4); f=find(t>=-3,1,'first'); r4.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,4),size(allResp,1),size(allResp,2));

end

function data=noNansOrInfs(data)

data(isnan(data))=mean(data(isfinite(data)),'all');
data(isinf(data))=mean(data(isfinite(data)),'all');

end

function plotOverlayedResponses(out,c1,c2)

figure();
plot(out.response1.t,out.response1.me,'Color',c1); hold on; 
plot(out.response1.t,out.response1.plusSe,'Color',c1); plot(out.response1.t,out.response1.minusSe,'Color',c1);
plot(out.response2.t,out.response2.me,'Color',c2); hold on; 
plot(out.response2.t,out.response2.plusSe,'Color',c2); plot(out.response2.t,out.response2.minusSe,'Color',c2);

end

function cmap=getCmapWithRed(cuezbins)

cmap=colormap(jet(255));
cmapstep=floor(255/(length(cuezbins)-1));
cmap=cmap(end:-cmapstep:1,:);
cmap=cmap(1:length(cuezbins)-1,:);
cmap=flipud(cmap);

end

function [grp1,grp2]=plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,addToTit,ds,smoo,cmapname,cuezbins,basesubtract,basetimewindow,plotAll)

alph=0.8;
plotBackwards=false;
doMaxInstead=false;
maxTimeBin=0.1;
noGaussSmooth=true;

backup_cuezbins=cuezbins;
if iscell(backup_cuezbins)
    cuezbins=backup_cuezbins{1};
end

bycuez=cell(length(cuezbins)-1,1);
bycuez_plusse=cell(length(cuezbins)-1,1);
bycuez_minusse=cell(length(cuezbins)-1,1);
time=cell(length(cuezbins)-1,1);
gpLab=1;
for i=1:length(cuezbins)-1
    temp=cued_success_Response;
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',gpLab); f=find(exclu~=1); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=1;
    sub1.incuezrange=incuezrange(temp.idx==gpLab); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=1); sub1.excluded(f(incuezrange~=1))=1;
    if all(isnan(sub1.unitbyunit_y),'all')
        continue
    end
    if doMaxInstead==true
        out=plotVariousSUResponsesAlignedToBeh('maxAcrossUnits',sub1,ds,maxTimeBin,true);
    else
        out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    end
    if ~isempty(smoo)
        if noGaussSmooth==true
            out.me=smooth(out.me,smoo); out.plusSe=smooth(out.plusSe,smoo); out.minusSe=smooth(out.minusSe,smoo); 
        else
            out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
        end
        if plotAll==true
            for j=1:size(out.unitbyunit,1)
                if noGaussSmooth==true
                    out.unitbyunit(j,:)=smooth(out.unitbyunit(j,:),smoo);
                else
                    out.unitbyunit(j,:)=smoothdata(out.unitbyunit(j,:),'gaussian',smoo);
                end
            end
        end
    end
    grp1.allunits{i}=out.unitbyunit;
    grp1.cuez{i}=cuez((cuez>cuezbins(i) & cuez<=cuezbins(i+1)) & temp.idx==gpLab);
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    if plotAll==true
        bycuez{i}=out.unitbyunit;
        time{i}=out.unittimes;
    else
        bycuez{i}=out.me;
        bycuez_plusse{i}=out.plusSe;
        bycuez_minusse{i}=out.minusSe;
        time{i}=out.t;
    end
end
grp1.bycuez=bycuez; grp1.time=time;
figure();
if strcmp(cmapname,'jet')
%     cmap=colormap(jet(length(cuezbins)-1));
    cmap=getCmapWithRed(cuezbins);
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
if plotBackwards==true
    for i=length(cuezbins)-1:-1:1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
else
    for i=1:length(cuezbins)-1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==false
            plot(time{i},bycuez_plusse{i},'Color',cmap(i,:)); 
            plot(time{i},bycuez_minusse{i},'Color',cmap(i,:)); 
        end
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
end
title(['groupLabel is 1 ' addToTit]);
txt=['grp1: '];
for i=1:length(grp1.allunits)
    txt=[txt num2str(size(grp1.allunits{i},1)) ' units, '];
end
text(1,1,txt);

if iscell(backup_cuezbins)
    cuezbins=backup_cuezbins{2};
end

bycuez=cell(length(cuezbins)-1,1);
gpLab=2;
for i=1:length(cuezbins)-1
    temp=cued_success_Response;
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',gpLab); f=find(exclu~=1); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=1;
    sub1.incuezrange=incuezrange(temp.idx==gpLab); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=1); sub1.excluded(f(incuezrange~=1))=1;
    if all(isnan(sub1.unitbyunit_y),'all')
        continue
    end
    if doMaxInstead==true
        out=plotVariousSUResponsesAlignedToBeh('maxAcrossUnits',sub1,ds,maxTimeBin,true);
    else
        out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    end
    if ~isempty(smoo)
        if noGaussSmooth==true
            out.me=smooth(out.me,smoo); out.plusSe=smooth(out.plusSe,smoo); out.minusSe=smooth(out.minusSe,smoo); 
        else
            out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
        end
        if plotAll==true
            for j=1:size(out.unitbyunit,1)
                if noGaussSmooth==true
                    out.unitbyunit(j,:)=smooth(out.unitbyunit(j,:),smoo);
                else
                    out.unitbyunit(j,:)=smoothdata(out.unitbyunit(j,:),'gaussian',smoo);
                end
            end
        end
    end
    grp2.allunits{i}=out.unitbyunit;
    grp2.cuez{i}=cuez((cuez>cuezbins(i) & cuez<=cuezbins(i+1)) & temp.idx==gpLab);
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    if plotAll==true
        bycuez{i}=out.unitbyunit;
        time{i}=out.unittimes;
    else
        bycuez{i}=out.me;
        bycuez_plusse{i}=out.plusSe;
        bycuez_minusse{i}=out.minusSe;
        time{i}=out.t;
    end
end
grp2.bycuez=bycuez; grp2.time=time;
figure();
if strcmp(cmapname,'jet')
%     cmap=colormap(jet(length(cuezbins)-1));
    cmap=getCmapWithRed(cuezbins);
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
if plotBackwards==true
    for i=length(cuezbins)-1:-1:1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
else
    for i=1:length(cuezbins)-1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==false
            plot(time{i},bycuez_plusse{i},'Color',cmap(i,:)); 
            plot(time{i},bycuez_minusse{i},'Color',cmap(i,:)); 
        end
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
        hold on;
    end
end
title(['groupLabel is 2 ' addToTit]);

txt=['grp2: '];
for i=1:length(grp2.allunits)
    txt=[txt num2str(size(grp2.allunits{i},1)) ' units, '];
end
text(1,1,txt);

end