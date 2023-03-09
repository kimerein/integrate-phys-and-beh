function r=matchAllUnits(r)

if length(r)<2
    error('Must pass in at least two responses');
end

k=1;
for i=1:length(r)
    for j=i+1:length(r)
        disp([i j]);
        out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',r{i},r{j},[],[],[]);
        if k==1
            allrmvd=out.Response1.excluded;
        else
            allrmvd=allrmvd+out.Response1.excluded;
        end
        k=k+1;
    end
end
allrmvd=allrmvd>0;

for i=1:length(r)
    r{i}=removeUnitFromResponse(r{i},allrmvd);
end

end