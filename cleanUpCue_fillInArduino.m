function [aligned,settings]=cleanUpCue_fillInArduino(varargin)
 
isOrchestra=0;

aligned=varargin{1};
if length(varargin)>1
    minProm=varargin{2};
    minProm2=varargin{3};
else
    if isOrchestra==1
        %     minProm=prctile(aligned.cueZone,95)-prctile(aligned.cueZone,40);
        %     minProm2=prctile(aligned.cueZone,60)-prctile(aligned.cueZone,20);
        minProm=prctile(aligned.cueZone,99)-prctile(aligned.cueZone,1);
        minProm2=prctile(aligned.cueZone,98)-prctile(aligned.cueZone,2);
    else
        minProm=50000;
        minProm2=25000;
    end
end
settings.minProm=minProm;
settings.minProm2=minProm2;

temp=aligned.cueZone-median(aligned.cueZone,'omitnan');
temp=temp./max(temp,[],2,'omitnan');
temp(temp>=0.5)=1;
temp(temp<0.5)=0;
aligned.cue=temp;

[~,fmovie,widths,proms]=findpeaks(aligned.cueZone,'MinPeakProminence',minProm);
[~,farduino]=findpeaks(aligned.cue,'MinPeakHeight',0.5);
usePeak=zeros(size(aligned.cueZone));
usePeak(isnan(aligned.cueZone))=nan;
for i=1:length(farduino)
    currArduinoPeak=farduino(i);
    [~,mi]=min(abs(currArduinoPeak-fmovie));
    for j=1:50
        temp=aligned.cueZone(fmovie(mi)-(j-1))-aligned.cueZone(fmovie(mi)-j);
        if temp>minProm2
            break
        end
    end
    usePeak(fmovie(mi)-(j-1))=1;
end
aligned.cueZone_onVoff=usePeak;

figure();
plot(aligned.cueZone,'Color','k');
hold on;
plot(aligned.cue*range(aligned.cueZone)+min(aligned.cueZone),'Color','r');
plot(aligned.cueZone_onVoff*range(aligned.cueZone)+min(aligned.cueZone),'Color','b');
legend({'original from movie','from arduino','fixed'});

end

