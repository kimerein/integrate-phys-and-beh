function [tensorfortrain,allLabels]=dataAugmentationTensor(tensor,allr2scores,augmentscheme)

% starting tensor
truetensor=tensor(:,:,1:4);

switch augmentscheme
    case 1
        % drop cells according to achievable quality of GLM fit
        abovezero_tensor=truetensor;
        abovezero_tensor(allr2scores<=0,:,:)=0;

        abovepoi05_tensor=truetensor;
        abovepoi05_tensor(allr2scores<=0.05,:,:)=0;

        abovepoi1_tensor=truetensor;
        abovepoi1_tensor(allr2scores<=0.1,:,:)=0;

        abovepoi15_tensor=truetensor;
        abovepoi15_tensor(allr2scores<=0.15,:,:)=0;

        % randomly drop 30% of the cells
        nRuns=50;
        tensorwithdrops=nan(size(truetensor,1),size(truetensor,2),nRuns*4);
        for i=1:nRuns
            dropforthistens=randsample(size(truetensor,1),floor(0.33*size(truetensor,1)));
            tensorwithdrops(:,:,(i-1)*4+1:i*4)=truetensor;
            tensorwithdrops(dropforthistens,:,(i-1)*4+1:i*4)=0;
        end

        % randomly drop 80% of the cells
        nRuns=50;
        tensorwithdrops2=nan(size(truetensor,1),size(truetensor,2),nRuns*4);
        for i=1:nRuns
            dropforthistens=randsample(size(truetensor,1),floor(0.8*size(truetensor,1)));
            tensorwithdrops2(:,:,(i-1)*4+1:i*4)=truetensor;
            tensorwithdrops2(dropforthistens,:,(i-1)*4+1:i*4)=0;
        end

        % put all together
        tensorfortrain=cat(3,truetensor,abovezero_tensor);
        tensorfortrain=cat(3,tensorfortrain,abovepoi05_tensor);
        tensorfortrain=cat(3,tensorfortrain,abovepoi1_tensor);
        tensorfortrain=cat(3,tensorfortrain,abovepoi15_tensor);
        tensorfortrain=cat(3,tensorfortrain,tensorwithdrops);
        tensorfortrain=cat(3,tensorfortrain,tensorwithdrops2);

        allLabels=repmat((0:3)',size(tensorfortrain,3)/4,1);
    case 2
        % just train cued success vs uncued success
        subtensor=truetensor(:,:,[1 3]);
        % randomly drop 90% of the cells
        nRuns=100;
        tensorwithdrops2=nan(size(subtensor,1),size(subtensor,2),nRuns*2);
        for i=1:nRuns
            dropforthistens=randsample(size(subtensor,1),floor(0.9*size(subtensor,1)));
            tensorwithdrops2(:,:,(i-1)*2+1:i*2)=subtensor;
            tensorwithdrops2(dropforthistens,:,(i-1)*2+1:i*2)=0;
        end
        tensorfortrain=cat(3,truetensor,truetensor);
        tensorfortrain=cat(3,tensorfortrain,tensorwithdrops2);

        allLabels=[(0:3)'; (0:3)'; repmat(([0 2])',size(tensorwithdrops2,3)/2,1)];
    case 3
        % just train cued failure vs uncued failure
        subtensor=truetensor(:,:,[2 4]);
        % randomly drop 90% of the cells
        nRuns=100;
        tensorwithdrops2=nan(size(subtensor,1),size(subtensor,2),nRuns*2);
        for i=1:nRuns
            dropforthistens=randsample(size(subtensor,1),floor(0.9*size(subtensor,1)));
            tensorwithdrops2(:,:,(i-1)*2+1:i*2)=subtensor;
            tensorwithdrops2(dropforthistens,:,(i-1)*2+1:i*2)=0;
        end
        tensorfortrain=cat(3,truetensor,truetensor);
        tensorfortrain=cat(3,tensorfortrain,tensorwithdrops2);

        allLabels=[(0:3)'; (0:3)'; repmat(([1 3])',size(tensorwithdrops2,3)/2,1)];
end

