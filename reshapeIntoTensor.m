function [tens,labels]=reshapeIntoTensor(mergedAllData,labels,nPointsBeforeEvent,takeNPointsAfterEvent)

% inputs from LDA_analysis
totalPoints=nPointsBeforeEvent+takeNPointsAfterEvent;
nCells=size(mergedAllData,2)/totalPoints;
tens=nan(nCells,totalPoints,length(labels));
for j=1:length(labels)
    for i=1:nCells
        tens(i,:,j)=mergedAllData(j,(i-1)*totalPoints+1:i*totalPoints);
    end
end
