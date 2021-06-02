function [alltbt,metadata,trialTypes]=turnOffLED(alltbt,metadata,trialTypes,mouseids)

for i=1:length(mouseids)
    currmouseid=mouseids(i);
    % correct alltbt
    temp=~isnan(alltbt.optoOn(metadata.mouseid==currmouseid,:));
    alltbt.optoOn(temp)=0;
    % correct metadata
    metadata.optoOnHere(metadata.mouseid==currmouseid)=0;
    % correct trialTypes
    f=fieldnames(trialTypes);
    for j=1:length(f)
        currfield=f{j};
        if ~isempty(regexp(currfield,'led'))
            if i==1
                disp(['correcting field in trialTypes: ' currfield]);
            end
            temp=trialTypes.(currfield);
            temp(metadata.mouseid==currmouseid)=0;
            trialTypes.(currfield)=temp;
        elseif ~isempty(regexp(currfield,'optoGroup'))
            if i==1
                disp(['correcting field in trialTypes: ' currfield]);
            end
            temp=trialTypes.(currfield);
            temp(metadata.mouseid==currmouseid)=nan;
            trialTypes.(currfield)=temp;
        end
    end
end

end