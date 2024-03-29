function Response=takeOnlyUnitsFromSess(Response,takeTheseSess)

if isfield(Response,'fromWhichSess_forTrials')
    if length(Response.fromWhichSess_forTrials)==size(Response.unitbyunit_y,1)
        sessionsMatchingRows=Response.fromWhichSess_forTrials;
    else
        sessionsMatchingRows=Response.fromWhichSess;
    end
else
    sessionsMatchingRows=Response.fromWhichSess;
end

Response.unitbyunit_x=Response.unitbyunit_x(ismember(sessionsMatchingRows,takeTheseSess),:);
Response.unitbyunit_y=Response.unitbyunit_y(ismember(sessionsMatchingRows,takeTheseSess),:);
Response.aligncomp_x=Response.aligncomp_x(ismember(sessionsMatchingRows,takeTheseSess),:);
Response.aligncomp_y=Response.aligncomp_y(ismember(sessionsMatchingRows,takeTheseSess),:);
if length(Response.excluded)==length(sessionsMatchingRows)
    Response.excluded=Response.excluded(ismember(sessionsMatchingRows,takeTheseSess));
end
if length(Response.ns)==length(sessionsMatchingRows)
    Response.ns=Response.ns(ismember(sessionsMatchingRows,takeTheseSess));
end
if length(Response.fromWhichSess)==length(sessionsMatchingRows)
    Response.fromWhichSess=Response.fromWhichSess(ismember(sessionsMatchingRows,takeTheseSess));
end

end