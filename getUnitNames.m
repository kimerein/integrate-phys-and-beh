function names=getUnitNames(inpt)

if ~iscell(inpt)
    dd{1}=inpt;
else
    dd=inpt;
end

unit_count=1;
names={};
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    ls=dir(datadir);
    for i=3:length(ls)
        if contains(ls(i).name,'trainingSet') || contains(ls(i).name,'testSet')
            continue
        end
        unitTest=doUnitTest(ls(i).folder, ls(i).name);
        if isempty(unitTest)
            disp(['skipping ' ls(i).name]);
            continue
        end
        if unitTest
            names{unit_count}=ls(i).name;
            unit_count=unit_count+1;
        end
    end
end
unit_count=unit_count-1;

end