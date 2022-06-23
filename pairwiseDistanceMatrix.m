function arrayDist=pairwiseDistanceMatrix(dataArray,referenceArray)

% rows are 3 dimensions
% columns are different points

if size(dataArray,2)~=size(referenceArray,2)
    error('dataArray and referenceArray must have the same number of points');
end

% scale two matrices so the pairwise distances between points matches
pairDist=nan(size(dataArray,2),size(dataArray,2));
for i=1:size(dataArray,2)
    for j=1:size(dataArray,2)
        pairDist(i,j)=sqrt((dataArray(1,i)-dataArray(1,j))^2 + (dataArray(2,i)-dataArray(2,j))^2 + (dataArray(3,i)-dataArray(3,j))^2);
    end
end
dataArrayPairwiseDist=sum(sum(pairDist,2,'omitnan'),1,'omitnan');
dataArray=dataArray./dataArrayPairwiseDist;
pairDist=nan(size(referenceArray,2),size(referenceArray,2));
for i=1:size(referenceArray,2)
    for j=1:size(referenceArray,2)
        pairDist(i,j)=sqrt((referenceArray(1,i)-referenceArray(1,j))^2 + (referenceArray(2,i)-referenceArray(2,j))^2 + (referenceArray(3,i)-referenceArray(3,j))^2);
    end
end
referenceArrayPairwiseDist=sum(sum(pairDist,2,'omitnan'),1,'omitnan');
referenceArray=referenceArray./referenceArrayPairwiseDist;

arrayDist=nan(size(dataArray,2),1);
for i=1:size(dataArray,2)
    arrayDist(i)=sqrt((dataArray(1,i)-referenceArray(1,i))^2 + (dataArray(2,i)-referenceArray(2,i))^2 + (dataArray(3,i)-referenceArray(3,i))^2);
end

figure();
imagesc(arrayDist);

disp('Abs val sum of distances between arrays');
disp(sum(abs(arrayDist),1,'omitnan'));

end