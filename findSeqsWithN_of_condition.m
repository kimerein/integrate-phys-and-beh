function isSeq=findSeqsWithN_of_condition(trialTypes, whichField, N, K, greaterThan)

temp=trialTypes.(whichField);
isSeq=nan(length(temp),1);
for i=1:length(temp)
    if i+K-1>length(temp)
        break
    end
    subtemp=temp(i:i+K-1);
    if greaterThan==true
        if nansum(subtemp==1)>N
            isSeq(i)=1;
        else
            isSeq(i)=0;
        end
    else
        if nansum(subtemp==1)<N
            isSeq(i)=1;
        else
            isSeq(i)=0;
        end
    end
end

% fix sequences that span video or session breaks

end