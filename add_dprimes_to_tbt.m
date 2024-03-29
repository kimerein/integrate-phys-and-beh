function [metadata,alltbt,out]=add_dprimes_to_tbt(varargin)

if length(varargin)==4
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
    dprimes=varargin{4};
    disp('input 1: alltbt, input 2: trialTypes, input 3: metadata, input 4: dprimes');
elseif length(varargin)==7
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
    dprimes=varargin{4};
    reachName=varargin{5};
    cueName=varargin{6};
    settings=varargin{7};
end

if ~isempty(dprimes)
    disp('adding dprimes passed in');
    % assume dprimes are ordered by session
    u=unique(metadata.sessid);
    for i=1:length(u)
        alltbt.dprimes(metadata.sessid==u(i))=dprimes(i);
    end
    metadata.dprimes=alltbt.dprimes;
    out.dprimes=alltbt.dprimes;
else
    % Need to make sure sessid are continuous ints from 1 to number of
    % unique sessions
    % First assign unique sessids
    u=unique(metadata.mouseid);
    j=0;
    for i=1:length(u)
        metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
        j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
    end
    u=unique(metadata.sessid); backup=metadata.sessid;
    for i=1:length(u)
        metadata.sessid(backup==u(i))=i;
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
end