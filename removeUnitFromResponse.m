function [newResponse,whichToTakeFromUnits]=removeUnitFromResponse(Response,whichRmv)

warning('Function currently only works for unit averages, not trial by trial');

f=find(Response.excluded'==0 & ismember(1:length(Response.excluded),find(whichRmv==1)));
currIncluded=find(Response.excluded'==0);
whichToTakeFromUnits=find(ismember(currIncluded,f));

newResponse=Response;
fi=fieldnames(newResponse);
for i=1:length(fi)
    temp=newResponse.(fi{i});
    if size(temp,1)==length(currIncluded)
        % is units
        newResponse.(fi{i})=temp(~ismember(1:size(temp,1),whichToTakeFromUnits),:);
    end
end
newResponse.excluded(whichRmv)=1;

end