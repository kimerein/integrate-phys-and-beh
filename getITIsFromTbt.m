function out=getITIsFromTbt(tbt,out)

nWheelTurns=nan(size(tbt.pelletPresented,1),1);
for i=1:size(tbt.pelletPresented,1)
    nWheelTurns(i)=sum(tbt.pelletPresented(i,:)>0.05);
end

out.nWheelTurns=nWheelTurns;