function data=normalizeDataMatrix(data,whichdim,bywhat)

if length(whichdim)>1
    % expand along dimension whichdim(2)
    % then max normalize along dimension whichdim(1)
    ed=flatten(data,'expand',[whichdim(2) whichdim(1)]);
    switch bywhat
        case 'max'
            ma=max(ed,[],whichdim(1),'omitnan');
        case 'sd'
            ma=std(ed,0,whichdim(1),'omitnan');
        case 'mean'
            ma=mean(ed,whichdim(1),'omitnan');
    end
    d=[1 2 3]; missingdim=d(~ismember(d,whichdim));
    if missingdim==1
        sz=[1 size(data,2) size(data,3)];
    elseif missingdim==2
        sz=[size(data,1) 1 size(data,3)];
    elseif missingdim==3
        sz=[size(data,1) size(data,2) 1];
    else
        error('principaledCA.m works only for 3D data matrices');
    end
    data=data./repmat(ma,sz);
else
    if whichdim==1
        sz=[size(data,1) 1 1];
    elseif whichdim==2
        sz=[1 size(data,2) 1];
    elseif whichdim==3
        sz=[1 1 size(data,3)];
    else
        error('principaledCA.m works only for 3D data matrices');
    end
    switch bywhat
        case 'max'
           ma=max(data,[],whichdim,'omitnan');
        case 'sd'
           ma=std(data,0,whichdim,'omitnan');
        case 'mean'
           ma=mean(data,whichdim,'omitnan');
    end
    data=data./repmat(ma,sz);
end

end

function flatData=flatten(data,method,dim)

switch method
    case 'mean'
        flatData=mean(data,dim,'omitnan');
        sz=size(flatData);
        flatData=reshape(flatData,sz(find(sz>1,1,'first')),sz(find(sz>1,1,'last')));
    case 'max'
        flatData=max(data,[],dim,'omitnan');
        sz=size(flatData);
        flatData=reshape(flatData,sz(find(sz>1,1,'first')),sz(find(sz>1,1,'last')));
    case 'expand'
        if dim(1)==1 
            s=[2 3 1];
        elseif dim(1)==2
            s=[1 3 2];
        elseif dim(1)==3
            s=[1 2 3];
        else
            error('principaledCA.m works only for 3D data matrices');
        end
        data=permute(data,s);
        flatData=[];
        for i=1:size(data,3)
            if dim(2)==1
                flatData=[flatData; data(:,:,i)];
            elseif dim(2)==2
                flatData=[flatData data(:,:,i)];
            else
                error('failed to flatten');
            end
        end        
end

end