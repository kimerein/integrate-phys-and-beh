function depth=getUnitDepth_forWHISPER(spikeBiggestOnCh, chmap, chSpacing, ventralMostTipDepth)

% second row is Matlab index, first row is depth on probe, where 32 is most
% dorsal, 1 is most ventral
% for A1x32Edge
% was using 20 micron spacing between chs
% chDepthMapping=[1   21; ...
%                 2   24; ...
%                 3   22; ...
%                 4   23; ...
%                 5   20; etc.

f=find(chmap(:,2)==spikeBiggestOnCh);
if isempty(f)
    disp('cannot find channel in channel map in getUnitDepth_forWHISPER.m');
    return
end
depth=ventralMostTipDepth-(chSpacing*(32-f));

end