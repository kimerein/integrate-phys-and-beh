function trialTypes=getLongITIs(alltbt,trialTypes,RTsettings)

firstTimes=nanmin(alltbt.timesFromSessionStart,[],2);
endTimes=nanmax(alltbt.timesFromSessionStart,[],2);

trialLengths=endTimes-firstTimes;

[n,x]=hist(trialLengths,50);

figure();
plot(x,n);
ylabel('Count');
xlabel('ITI distribution (sec)');

disp(['Using ITI threshold of ' num2str(RTsettings.longITIthresh) ' sec']);

trialTypes.isLongITI=trialLengths>=RTsettings.longITIthresh;
