function alltimesfromcue=whenCuesWrtReach(data_loc_array,whichreach,onlyFirstReach)

alltimesfromcue=[];
for i=1:size(data_loc_array,1)
    if exist([data_loc_array{i,6} '\beh2_tbt.mat'],'file')
        a=load([data_loc_array{i,6} '\beh2_tbt.mat']);
        [timesfromcue,indsfromcue]=forday_whenCuesWrtReach(a.beh2_tbt,whichreach,onlyFirstReach);
        alltimesfromcue=[alltimesfromcue; timesfromcue];
    end
    if mod(i,25)==0
        disp(i);
    end
end


end

function [timesfromcue,indsfromcue]=forday_whenCuesWrtReach(beh2_tbt,whichreach,onlyFirstReach)

% for each trial, find first reach after cue
% then get time of cue wrt this reach
temp=beh2_tbt.(whichreach);

if onlyFirstReach==true
    indsfromcue=nan(size(temp,1),1);
    for i=1:size(temp,1)
        f=find(temp(i,94:end),1,'first'); % find first reach after cue
        if ~isempty(f)
            indsfromcue(i)=f-1;
        end
    end
    timesfromcue=indsfromcue.*0.035;
else
    indsfromcue=[];
    for i=1:size(temp,1)
        f=find(temp(i,94:end)); % find all reaches after cue
        if ~isempty(f)
            indsfromcue=[indsfromcue f-1];
        end
    end
    timesfromcue=indsfromcue.*0.035;
    timesfromcue=timesfromcue';
end

% figure(); histogram(timesfromcue,100);

end