function [indexIntoList,backwardsIndexIntoList,unitnames_glm_notInThisList]=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names)

list=unitbyunit_names.names;
indexIntoList=nan(length(unitnames_glm),1);
backwardsIndexIntoList=nan(length(list),2);
unitnames_glm_notInThisList=zeros(length(unitnames_glm),1);

currname=unitnames_glm{1};
r=regexp(currname,'_', 'once');
if isempty(r)
    notGLMnames=true;
else
    notGLMnames=false;
end

startAt=1;
counterFor_backwardsIndexIntoList=1;
for i=1:length(unitnames_glm)
    currname=unitnames_glm{i};
    if notGLMnames==false
        r=regexp(currname,'_');
        subname=currname(1:r(1)-1);
    else
        subname=currname;
    end
    for j=startAt:length(list)
        if strcmp(subname,list{j})
            indexIntoList(i)=j;
            startAt=j;
            backwardsIndexIntoList(counterFor_backwardsIndexIntoList,:)=[i j];
            counterFor_backwardsIndexIntoList=counterFor_backwardsIndexIntoList+1;
            break
        else
            if j==length(list)
                unitnames_glm_notInThisList(i)=1;
                counterFor_backwardsIndexIntoList=counterFor_backwardsIndexIntoList+1;
            end
        end
    end
end

end