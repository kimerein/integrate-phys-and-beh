function [trialTypes]=getTimingOfOptoOffset(alltbt,useOptoField,trialTypes,multipleOptoTimes)

thresh=0.5;

% Get timing of opto ends
opto=alltbt.(useOptoField);
optoEnd=nan(size(opto,1),1);
for i=1:size(opto,1)
    fon=find(opto(i,:)>thresh,1,'first');
    f=fon+find(opto(i,fon:end)<thresh,1,'first')-1;
    if ~isempty(f)
        optoEnd(i)=f;
    end
end

switch useOptoField
    case 'optoOn'
        % subtract cue (Arduino cue)
        cue=alltbt.cue;
        cueStart=nan(size(cue,1),1);
        for i=1:size(cue,1)
            f=find(cue(i,:)>thresh,1,'first');
            if ~isempty(f)
                cueStart(i)=f;
            end
        end
    case 'optoZone'
        % subtract cueZone_onVoff (movie cue)
        f=find(nanmean(alltbt.cueZone_onVoff,1)>thresh,1,'first');
        cueStart=repmat(f,size(opto,1),1);
end

optoEnd=optoEnd-cueStart; % get optoEnd with respect to cueStart

timeStep=mode(diff(nanmean(alltbt.times,1)));

optoEnd=optoEnd.*timeStep; % convert to real time

% Plot
[n,x]=hist(optoEnd,floor(size(opto,1)/2));
figure();
plot(x,n);
ylabel('Count');
xlabel('sec from trial onset');
pause;
    
if multipleOptoTimes==true
    answer=inputdlg({'First threshold:','Second threshold:','Third threshold:','Fourth threshold:'});
    
    figure();
    plot(x,n,'Color','k');
    hold all;
    thresh1=str2num(answer{1});
    thresh2=str2num(answer{2});
    thresh3=str2num(answer{3});
    thresh4=str2num(answer{4});
    line([thresh1 thresh1],[0 nanmax(n)]);
    line([thresh2 thresh2],[0 nanmax(n)]);
    line([thresh3 thresh3],[0 nanmax(n)]);
    line([thresh4 thresh4],[0 nanmax(n)]);
    ylabel('Count');
    xlabel('sec from trial onset');
    title('Opto OFFSET');
    
    trialTypes.optoGroupEnd=nan(size(trialTypes.led));
    trialTypes.optoGroupEnd(optoEnd<thresh1)=1;
    trialTypes.optoGroupEnd(optoEnd>=thresh1 & optoEnd<thresh2)=2;
    trialTypes.optoGroupEnd(optoEnd>=thresh2 & optoEnd<thresh3)=3;
    trialTypes.optoGroupEnd(optoEnd>=thresh3 & optoEnd<thresh4)=4;
    trialTypes.led(optoEnd>=thresh4)=0;
else
    answer=inputdlg({'Threshold (opto cannot end after this, in sec):'});
    
    figure();
    plot(x,n,'Color','k');
    hold all;
    thresh1=str2num(answer{1});
    line([thresh1 thresh1],[0 nanmax(n)]);
    ylabel('Count');
    xlabel('sec from trial onset');
    
    trialTypes.led(optoEnd>=thresh1)=0;
end


