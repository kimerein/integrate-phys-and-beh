function [trialTypes]=getTimingOfOpto(alltbt,useOptoField,trialTypes,multipleOptoTimes)

thresh=0.5;

% Get timing of opto starts
opto=alltbt.(useOptoField);
optoStart=nan(size(opto,1),1);
for i=1:size(opto,1)
    f=find(opto(i,:)>thresh,1,'first');
    if ~isempty(f)
        optoStart(i)=f;
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

optoStart=optoStart-cueStart; % get optoStart with respect to cueStart

timeStep=mode(diff(nanmean(alltbt.times,1)));

optoStart=optoStart.*timeStep; % convert to real time

% Plot
[n,x]=hist(optoStart,floor(size(opto,1)/2));
figure();
plot(x,n);
ylabel('Count');
xlabel('sec from trial onset');
pause;
    
if multipleOptoTimes==true
    answer=inputdlg({'First threshold:','Second threshold:','Third threshold'});
    
    figure();
    plot(x,n,'Color','k');
    hold all;
    thresh1=str2num(answer{1});
    thresh2=str2num(answer{2});
    thresh3=str2num(answer{3});
    line([thresh1 thresh1],[0 nanmax(n)]);
    line([thresh2 thresh2],[0 nanmax(n)]);
    line([thresh3 thresh3],[0 nanmax(n)]);
    ylabel('Count');
    xlabel('sec from trial onset');
    
    trialTypes.optoGroup=nan(size(trialTypes.led));
    trialTypes.optoGroup(optoStart<thresh1)=1;
    trialTypes.optoGroup(optoStart>=thresh1 & optoStart<thresh2)=2;
    trialTypes.optoGroup(optoStart>=thresh2 & optoStart<thresh3)=3;
    trialTypes.led(optoStart>=thresh3)=0;
else
    answer=inputdlg({'Threshold (opto cannot start after this, in sec):'});
    
    figure();
    plot(x,n,'Color','k');
    hold all;
    thresh1=str2num(answer{1});
    line([thresh1 thresh1],[0 nanmax(n)]);
    ylabel('Count');
    xlabel('sec from trial onset');
    
    trialTypes.led(optoStart>=thresh1)=0;
end


