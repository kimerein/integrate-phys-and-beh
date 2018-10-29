function out=interpLearningCurve(isreaching_out,dprimes,breakHere)

getTheseSess=1:1:40;

sess1isAfterNReaches=20; % mouse has to at least achieve this many successful reaches before start counting sessions
% start counting from the day before this day

f=find(isreaching_out.hasSuccess>=sess1isAfterNReaches,1,'first');
f=f-1;
if f<1
    f=1;
end
isreaching_out.nth_session=isreaching_out.nth_session-isreaching_out.nth_session(f)+1;

interp_dprime=interp1(isreaching_out.nth_session(1:breakHere-1),dprimes(1:breakHere-1),getTheseSess);
out.dprime=interp_dprime;
out.nth_session=getTheseSess;

if breakHere>length(dprimes)
    return
end

isreaching_out.nth_session(breakHere:end)=isreaching_out.nth_session(breakHere:end)-isreaching_out.nth_session(breakHere)+1;
interp_dprime=interp1(isreaching_out.nth_session(breakHere:end),dprimes(breakHere:end),getTheseSess);
out.dprime_2=interp_dprime;
out.nth_session_2=getTheseSess;