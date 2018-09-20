function alltbt=fixOptoZone(alltbt)

alltbt.optoZone=alltbt.optoZone-repmat(min(alltbt.optoZone,[],2),1,size(alltbt.optoZone,2));

figure();
plot(alltbt.optoZone');

x=input('Enter threshold for opto on. ');

alltbt.optoOn=alltbt.optoZone>x;