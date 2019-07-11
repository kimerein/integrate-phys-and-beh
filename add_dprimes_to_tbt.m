function metadata=add_dprimes_to_tbt(alltbt,out,metadata)

[dprimes]=get_dprime_per_session(alltbt,out,metadata,'all_reachBatch','cueZone_onVoff',[]);

metadata.dprimes=nan(size(metadata.sessid));
u=unique(metadata.sessid);
for i=1:length(u)
    curru=u(i);
    metadata.dprimes(ismember(metadata.sessid,curru))=dprimes(i);
end