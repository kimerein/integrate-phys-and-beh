function sub=subResponse(Response,whichFieldToTest,test)

f=fieldnames(Response);
whichToTest=Response.(whichFieldToTest);
useThese=ismember(whichToTest,test);
for i=1:length(f)
    temp=Response.(f{i});
    if size(temp,1)==length(useThese)
        sub.(f{i})=temp(useThese,:);
    else
        % don't filter
        sub.(f{i})=temp;
    end
end