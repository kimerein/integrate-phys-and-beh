function returnout=downSampReachPlotDatasets(plotReturn1,whichTrial1,downSampBy,linecol)

if ~ismember(whichTrial1,[1 2])
    error('whichTrial1 must be 1 or 2');
end

% Standard error of raw reaching was calculated as
% std(raw reach rate) / sqrt(n trials)

% 0. Get down sampled mean
temp1=plotReturn1.(['data' num2str(whichTrial1) '_mean']);
me=downSampAv(temp1{1},downSampBy);

% 1. Get back to standard deviation
temp1=plotReturn1.(['data' num2str(whichTrial1) '_se']);
sd1=temp1{1}.*sqrt(plotReturn1.n);

% 2. Variance from standard deviation
var1=sd1.^2;

% 3. Down sample variance
B=downSampBy*sum()
varDS=downSampAv(var1,downSampBy);

% 4. Take sqrt to get standard deviation
ds_sd=sqrt(varDS);

% 5. Divide through by sqrt(n trials) to get se
ds_se=ds_sd./sqrt(plotReturn1.n);



% 6. Plot downsampled
plotReturn1.time_for_x=downSampAv(plotReturn1.time_for_x,downSampBy);
figure();
[n,x]=cityscape_hist(me,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
hold on;
[n,x]=cityscape_hist(me-ds_se,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
[n,x]=cityscape_hist(me+ds_se,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
xlabel('Time (sec)'); ylabel(['Reach rate (reaches per sec) ' num2str(plotReturn1.n) ' trials']);

returnout.data1_mean{1}=me;
returnout.data1_se{1}=ds_se;
returnout.data2_mean{1}=me;
returnout.data2_se{1}=ds_se;
returnout.time_for_x=plotReturn1.time_for_x;
returnout.n=plotReturn1.n;

end

function mecombo=combineMean(me1,me2,n1,n2)

mecombo=(n1.*me1 + n2.*me2) ./ (n1+n2);

end

function varcombo=combineVar(var1,var2,n1,n2)

varcombo=((n1-1).*var1 + (n2-1).*var2)./( (n1-1) + (n2-1) );

end