function isSeq=findSeqsWithN_of_condition(trialTypes, whichField, N, K, greaterThan)

temp=trialTypes.(whichField);
isSeq=nan(length(temp),1);
showCheck=true;
if showCheck==true
    figure();
    counter=1;
end
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
    % fix sequences that span video or session breaks
    if length(unique(trialTypes.sessid(i:i+K-1)))>1
        isSeq(i)=0;
    end
    if showCheck==true
        if isSeq(i)==1
            scatter([i:i+K-1] - i,ones(size(i:i+K-1))*counter, [], subtemp, 'filled');
            hold all;
            counter=counter+1;
        end
    end
end

end