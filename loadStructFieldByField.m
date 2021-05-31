function out=loadStructFieldByField(varargin)

out=[];

if length(varargin)==1
    datadir=varargin{1};
    fieldsToLoad=[];
elseif length(varargin)==2
    datadir=varargin{1};
    fieldsToLoad=varargin{2};
end

doingcellarray=false;

if ismac==true
    sprtr='/';
else
    sprtr='\';
end

ls=dir(datadir);
for i=3:length(ls)
    thisname=ls(i).name;
    r=regexp(thisname,'.mat');
    thisname_withoutmat=thisname(1:r-1);
    thisisdir=ls(i).isdir;
    if thisisdir==1
        continue
    end
    if ~isempty(fieldsToLoad) 
        doesContain=false;
        for j=1:length(fieldsToLoad)
            if strcmp(thisname_withoutmat,fieldsToLoad{j})
                doesContain=true;
                indIntoFields=j;
                break
            end
        end
        if doesContain==false
            continue
        end
    end
    if ~isempty(fieldsToLoad)
        disp(['Loading field ' fieldsToLoad{indIntoFields}]);
        disp(['Loading field ' thisname]);
    end
    temp=regexp(thisname,'.mat','once');
    if ~isempty(regexp(thisname,'scell','once'))
        doingcellarray=true;
        if ~isempty(temp)
            a=load([datadir sprtr thisname]);
            indtostr=regexp(thisname,'scell','once');
            cellind=str2num(thisname(indtostr+5:temp-1));
            outcell(cellind)=a.temp;
        end
    elseif ~isempty(regexp(thisname,'cell','once'))
        doingcellarray=true;
        if ~isempty(temp)
            a=load([datadir sprtr thisname]);
            indtostr=regexp(thisname,'cell','once');
            cellind=str2num(thisname(indtostr+4:temp-1));
            outcell{cellind}=a.temp;
        end        
    else
        if ~isempty(temp)
            a=load([datadir sprtr thisname]);
            if length(fieldnames(a))==0
            else
                out.(thisname(1:temp-1))=a.temp;
            end
        end
    end
end

if doingcellarray==true
    out=outcell;
end