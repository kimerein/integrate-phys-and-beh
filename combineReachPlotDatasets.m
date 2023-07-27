function [mecombo,combo_se,returnout]=combineReachPlotDatasets(plotReturn1,whichTrial1,plotReturn2,whichTrial2,linecol)

if length(whichTrial1)>1
    % combine trial 1 to trial 1 and trial 2 to trial 2

    % Standard error of raw reaching was calculated as
    % std(raw reach rate) / sqrt(n trials)

    % 1. Get back to standard deviation
    temp1=plotReturn1.(['data' num2str(1) '_se']);
    sd1=temp1{1}.*sqrt(plotReturn1.n);
    temp2=plotReturn2.(['data' num2str(1) '_se']);
    sd2=temp2{1}.*sqrt(plotReturn2.n);

    temp1_trialn=plotReturn1.(['data' num2str(2) '_se']);
    sd1_trialn=temp1_trialn{1}.*sqrt(plotReturn1.n);
    temp2_trialn=plotReturn2.(['data' num2str(2) '_se']);
    sd2_trialn=temp2_trialn{1}.*sqrt(plotReturn2.n);

    % 2. Variance from standard deviation
    var1=sd1.^2;
    var2=sd2.^2;

    var1_trialn=sd1_trialn.^2;
    var2_trialn=sd2_trialn.^2;

    % 3. Combine variances
    varcombo=combineVar(var1,var2,plotReturn1.n,plotReturn2.n);

    varcombo_trialn=combineVar(var1_trialn,var2_trialn,plotReturn1.n,plotReturn2.n);

    % 4. Take sqrt to get combined standard deviation
    combo_sd=sqrt(varcombo);
    
    combo_sd_trialn=sqrt(varcombo_trialn);

    % 5. Divide through by sqrt(n trials) to get se
    combo_se=combo_sd./sqrt(plotReturn1.n + plotReturn2.n);

    combo_se_trialn=combo_sd_trialn./sqrt(plotReturn1.n + plotReturn2.n);

    % 5. Get trial-weighted combined mean
    temp1=plotReturn1.(['data' num2str(1) '_mean']);
    temp2=plotReturn2.(['data' num2str(1) '_mean']);
    mecombo=combineMean(temp1{1},temp2{1},plotReturn1.n,plotReturn2.n);

    temp1_trialn=plotReturn1.(['data' num2str(2) '_mean']);
    temp2_trialn=plotReturn2.(['data' num2str(2) '_mean']);
    mecombo_trialn=combineMean(temp1_trialn{1},temp2_trialn{1},plotReturn1.n,plotReturn2.n);

    % 6. Plot combo
    if length(plotReturn1.time_for_x)~=length(plotReturn2.time_for_x)
        error('Assumed times were the same for plotReturn1 and plotReturn2');
    end

    returnout.data1_mean{1}=mecombo;
    returnout.data1_se{1}=combo_se;
    returnout.data2_mean{1}=mecombo_trialn;
    returnout.data2_se{1}=combo_se_trialn;

    returnout.time_for_x=plotReturn1.time_for_x;
    returnout.n=plotReturn1.n+plotReturn2.n;
else
    if ~ismember(whichTrial1,[1 2])
        error('whichTrial1 must be 1 or 2');
    end
    if ~ismember(whichTrial2,[1 2])
        error('whichTrial2 must be 1 or 2');
    end

    % Standard error of raw reaching was calculated as
    % std(raw reach rate) / sqrt(n trials)

    % 1. Get back to standard deviation
    temp1=plotReturn1.(['data' num2str(whichTrial1) '_se']);
    sd1=temp1{1}.*sqrt(plotReturn1.n);
    temp2=plotReturn2.(['data' num2str(whichTrial2) '_se']);
    sd2=temp2{1}.*sqrt(plotReturn2.n);

    % 2. Variance from standard deviation
    var1=sd1.^2;
    var2=sd2.^2;

    % 3. Combine variances
    varcombo=combineVar(var1,var2,plotReturn1.n,plotReturn2.n);

    % 4. Take sqrt to get combined standard deviation
    combo_sd=sqrt(varcombo);

    % 5. Divide through by sqrt(n trials) to get se
    combo_se=combo_sd./sqrt(plotReturn1.n + plotReturn2.n);

    % 5. Get trial-weighted combined mean
    temp1=plotReturn1.(['data' num2str(whichTrial1) '_mean']);
    temp2=plotReturn2.(['data' num2str(whichTrial2) '_mean']);
    mecombo=combineMean(temp1{1},temp2{1},plotReturn1.n,plotReturn2.n);

    % 6. Plot combo
    if length(plotReturn1.time_for_x)~=length(plotReturn2.time_for_x)
        error('Assumed times were the same for plotReturn1 and plotReturn2');
    end
    figure();
    [n,x]=cityscape_hist(mecombo,plotReturn1.time_for_x);
    plot(x,n,'Color',linecol);
    hold on;
    [n,x]=cityscape_hist(mecombo-combo_se,plotReturn1.time_for_x);
    plot(x,n,'Color',linecol);
    [n,x]=cityscape_hist(mecombo+combo_se,plotReturn1.time_for_x);
    plot(x,n,'Color',linecol);
    xlabel('Time (sec)'); ylabel(['Reach rate (reaches per sec) ' num2str(plotReturn1.n+plotReturn2.n) ' trials']);

    returnout.data1_mean{1}=mecombo;
    returnout.data1_se{1}=combo_se;
    returnout.data2_mean{1}=mecombo;
    returnout.data2_se{1}=combo_se;
    returnout.time_for_x=plotReturn1.time_for_x;
    returnout.n=plotReturn1.n+plotReturn2.n;
end

end

function mecombo=combineMean(me1,me2,n1,n2)

mecombo=(n1.*me1 + n2.*me2) ./ (n1+n2);

end

function varcombo=combineVar(var1,var2,n1,n2)

varcombo=((n1-1).*var1 + (n2-1).*var2)./( (n1-1) + (n2-1) );

end