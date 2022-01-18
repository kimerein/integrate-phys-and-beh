function plotLossOverTraining(trainingoutfile)

fileID = fopen(trainingoutfile);

maxLines=100000;
tline=0;
whichit=nan(1,maxLines);
loss=nan(1,maxLines);
for i=1:maxLines
    if tline==-1
        break
    end
    tline = fgetl(fileID);
    if tline==-1
        break
    end
    % read iteration and loss
    iind=regexp(tline,'iteration: ');
    lind=regexp(tline,' loss: ');
    lrind=regexp(tline,' lr: ');
    if isempty(iind) || isempty(lind) || isempty(lrind)
        continue
    end
    whichit(i)=str2num(tline(iind+11:lind-1));
    loss(i)=str2num(tline(lind+8:lrind-1));
end
 
fclose(fileID);

figure();
plot(whichit,loss);
xlabel('Iteration');
ylabel('Loss');

end