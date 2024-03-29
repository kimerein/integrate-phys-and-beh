function Response=excludeTooFewTrials(Response,trial_n_cutoff,plotNtrialsHistogram)

if plotNtrialsHistogram
    figure(); 
    histogram([Response.ns],50);
    title('Histogram of trial counts across units');
end

% exclude units with too few trials
if nansum(Response.excluded==0)~=length(Response.ns)
    ('excluded field does not match Response.ns');
    return
end

% trial_n_cutoff = at least this many trials, else exclude unit
newExcluded=Response.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(newExcluded)) ' units because too few trials']);
Response.excluded(Response.excluded==0)=newExcluded;
u=unique(Response.fromWhichUnit);
Response=filterResponseToOneSU(Response,u(newExcluded==0));

end

function out=filterResponseToOneSU(Response,whichUnitToUse)

f=fieldnames(Response);
fieldLikeResponseSize=size(Response.unitbyunit_y,1);
if length(Response.fromWhichUnit)==fieldLikeResponseSize
    % took all trials
    whichUnit=Response.fromWhichUnit;
elseif length(Response.fromWhichSess)==fieldLikeResponseSize
    % took unit by unit
    whichUnit=1:fieldLikeResponseSize;
else
    error('do not recognize structure of Response in plotVariousSUResponsesAlignedToBeh.m');
end

for i=1:length(f)
    temp=Response.(f{i});
    if size(temp,1)==fieldLikeResponseSize
        % filter this field
        if length(size(temp))>1
            % 2D
            out.(f{i})=temp(ismember(whichUnit,whichUnitToUse),:);
        else
            % 1D
            out.(f{i})=temp(ismember(whichUnit,whichUnitToUse));
        end
    else
        % don't filter this field
        out.(f{i})=temp;
    end
end

end
