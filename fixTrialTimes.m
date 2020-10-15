function alltbt=fixTrialTimes(alltbt)

for i=1:size(alltbt.times,1)
    mi=nanmin(alltbt.times(i,:));
    if mi>=35
        alltbt.times(i,:)=alltbt.times(i,:)-mi;
    end
end
