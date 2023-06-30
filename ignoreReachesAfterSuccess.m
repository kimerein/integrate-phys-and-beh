function alltbt=ignoreReachesAfterSuccess(alltbt,metadata,timeAfterSuccess)

% will drop any no pellet reaches within timeAfterSuccess seconds from
% successful reach batch, because these likely represent chewing arm movements 
% depends to a degree on version of Kim's rig
succ=alltbt.reachBatch_success_reachStarts;
nopelletreaches=alltbt.pelletmissingreach_reachStarts;
temp=alltbt.times; temp=temp'; 
indsAfterSuccess=floor(timeAfterSuccess/mode(diff(temp(1:end))));
for i=1:size(alltbt.cue,1)
    temp=succ(i,:);
    currsess=metadata.sessid(i);
    if any(temp>0.5,2)
        f=find(temp>0.5,1,'first');
        indstouse=f:f+indsAfterSuccess;
        if indstouse(end)>size(temp,2)
            indstouse=f:size(temp,2);
            % go to next trial
            if i+1>size(alltbt.cue,1)
                break
            end
            if metadata.sessid(i+1)==currsess
                intonexttrial=f+indsAfterSuccess-size(temp,2);
                nopelletreaches(i+1,1:intonexttrial)=0;
            end
        end
        nopelletreaches(i,indstouse)=0;
    end
end
alltbt.pelletmissingreach_reachStarts=nopelletreaches;

end