function [tensor,allLabels]=catAndExpandTensor(tensor1, tensor, allLabels1, allLabels)

togetherTensor=cat(3, [tensor1; nan(size(tensor,1),size(tensor,2),size(tensor1,3))], [nan(size(tensor1,1),size(tensor1,2),size(tensor,3)); tensor]); 
tensor=togetherTensor;
allLabels=[allLabels1; allLabels];