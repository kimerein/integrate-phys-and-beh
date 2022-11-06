function excludeTooFewTrials(Response,trial_n_cutoff,plotNtrialsHistogram)

if plotNtrialsHistogram
    figure(); 
    histogram([Response.ns],50);
    title('Histogram of trial counts across units');
end

% exclude units with too few trials
if length(Response.excluded)~=length(Response.ns)
    Response.excluded=zeros(size(Responses.ns));
end

% trial_n_cutoff = at least this many trials, else exclude unit
Response.excluded=Response.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(Response.excluded==1)) ' units because too few trials']);

[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D1untagged_cueResponse,activityD1untagged]=cutExcluded(D1untagged_cueResponse,activityD1untagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);
[D2untagged_cueResponse,activityD2untagged]=cutExcluded(D2untagged_cueResponse,activityD2untagged);

end

function [D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged)

excluded1=D1tagged_cueResponse.excluded; 
excluded2=activityD1tagged.excluded;
eitherExcluded=excluded1==1 | excluded2==1;
D1tagged_cueResponse.unitbyunit_x=D1tagged_cueResponse.unitbyunit_x(eitherExcluded==0,:);
D1tagged_cueResponse.unitbyunit_y=D1tagged_cueResponse.unitbyunit_y(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_x=D1tagged_cueResponse.aligncomp_x(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_y=D1tagged_cueResponse.aligncomp_y(eitherExcluded==0,:);
D1tagged_cueResponse.excluded=D1tagged_cueResponse.excluded(eitherExcluded==0);
D1tagged_cueResponse.ns=D1tagged_cueResponse.ns(eitherExcluded==0);
D1tagged_cueResponse.fromWhichSess=D1tagged_cueResponse.fromWhichSess(eitherExcluded==0);

activityD1tagged.unitbyunit_x=activityD1tagged.unitbyunit_x(eitherExcluded==0,:);
activityD1tagged.unitbyunit_y=activityD1tagged.unitbyunit_y(eitherExcluded==0,:);
activityD1tagged.aligncomp_x=activityD1tagged.aligncomp_x(eitherExcluded==0,:);
activityD1tagged.aligncomp_y=activityD1tagged.aligncomp_y(eitherExcluded==0,:);
activityD1tagged.excluded=activityD1tagged.excluded(eitherExcluded==0);
activityD1tagged.ns=activityD1tagged.ns(eitherExcluded==0);
activityD1tagged.fromWhichSess=activityD1tagged.fromWhichSess(eitherExcluded==0);

end