function whichToTake=matchUpNames(bigList, isInNames)

whichToTake=zeros(1,length(bigList));
for i=1:length(bigList)
    for j=1:length(isInNames)
        cn=isInNames{j};
        r=regexp(cn,'_');
        subname=cn(1:r(1)-1);
        if contains(bigList{i},subname)
            whichToTake(i)=1;
            break
        end
    end
end