function [metadata,alltbt,out]=add_dprimes_to_tbt(varargin)

if length(varargin)==3
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
elseif length(varargin)==6
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
    reachName=varargin{4};
    cueName=varargin{5};
    settings=varargin{6};
end

[dprimes]=get_dprime_per_session(alltbt,out,metadata,reachName,cueName,settings);

metadata.dprimes=nan(size(metadata.sessid));
u=unique(metadata.sessid);
for i=1:length(u)
    curru=u(i);
    metadata.dprimes(ismember(metadata.sessid,curru))=dprimes(i);
end

alltbt.dprimes=metadata.dprimes;
out.dprimes=metadata.dprimes;