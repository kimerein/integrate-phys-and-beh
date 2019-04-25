function zscore_data=Zscore_by_session(data,sessid)

% all sessions with nan as sessid
% Zscore together using same mean and variance 

if any(isnan(sessid))
    disp(['There are ' num2str(sum(isnan(sessid))) ' trials with nan as sessid. Grouping these all into a separate session.']);
    m=nanmax(sessid);
    sessid(isnan(sessid))=m+1;
end

% make data a row vector
if size(data,1)>1
    data=data';
end
if size(data,1)>1
    error('data must be a 1D vector');
end

u=unique(sessid);
zscore_data=nan(size(data));
for i=1:length(u)
    curr_sessid=u(i);
    me=nanmean(data(sessid==curr_sessid));
    sd=nanstd(data(sessid==curr_sessid),[],2);
    zscore_data(sessid==curr_sessid)=(data(sessid==curr_sessid)-me)./sd;
end

