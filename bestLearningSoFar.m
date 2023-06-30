function bestLearningSoFar(days,include_firsthalf,include_secondhalf)

backup.include_firsthalf=include_firsthalf;
backup.include_secondhalf=include_secondhalf;

day1isfirstday=true;
if day1isfirstday==true
    % day 1 is the first day w any data to focus on early learning
    alignto=find(days==1);
    ndaysafter1=length(include_firsthalf(1,alignto:end));
    for i=1:size(include_firsthalf)
        temp=include_firsthalf(i,:);
        firstnotnan=find(~isnan(temp),1,'first');
        include_firsthalf(i,:)=nan;
        if firstnotnan>alignto
            tempie=temp(firstnotnan:end);
            include_firsthalf(i,alignto:alignto+length(tempie)-1)=tempie;
            temp=include_secondhalf(i,:);
            include_secondhalf(i,:)=nan;
            include_firsthalf(i,alignto:alignto+length(tempie)-1)=temp(firstnotnan:end);
        else
            include_firsthalf(i,alignto:end)=temp(firstnotnan:firstnotnan+ndaysafter1-1);
            temp=include_secondhalf(i,:);
            include_secondhalf(i,:)=nan;
            include_secondhalf(i,alignto:end)=temp(firstnotnan:firstnotnan+ndaysafter1-1);
        end
    end
end

for i=1:size(include_firsthalf,1)
    % mouse by mouse
    % find days when average dprime of first and second halves 
    % exceeds any previous dprime for that mouse
    avformouse=nanmean([nanmean(include_firsthalf(i,:),1); nanmean(include_secondhalf(i,:),1)],1);
    higherThanBefore=getDaysHigherThanBefore(avformouse,include_firsthalf(i,:),include_secondhalf(i,:));
    % nan out days that were not higherThanBefore
    include_firsthalf(i,higherThanBefore==0)=nan;
    include_secondhalf(i,higherThanBefore==0)=nan;
end

figure(); scatter(days-0.25,nanmean(include_firsthalf,1),[],'k'); hold on;
scatter(days+0.25,nanmean(include_secondhalf,1),[],'r');
for i=1:size(include_firsthalf,2)
    line([days(i)-0.25 days(i)+0.25],[nanmean(include_firsthalf(:,i),1) nanmean(include_secondhalf(:,i),1)],'Color','k');
end
for i=1:size(include_firsthalf,2)
    if nanmean(include_secondhalf(:,i),1)>nanmean(include_firsthalf(:,i),1)
        line([days(i) days(i)],[0.9 1],'Color','r');
        hold on;
    elseif nanmean(include_secondhalf(:,i),1)<nanmean(include_firsthalf(:,i),1)
        line([days(i) days(i)],[0.9 1],'Color','k');
    end
end

figure();
for i=1:size(include_firsthalf,2)
    if ismember(days(i),[1:39])
        line([1 2],[0 nanmean(include_secondhalf(:,i),1)-nanmean(include_firsthalf(:,i),1)],'Color','k'); hold on;
        scatter(2,nanmean(include_secondhalf(:,i),1)-nanmean(include_firsthalf(:,i),1));
    end
end

disp(nanmean(nanmean(include_firsthalf(:,ismember(days,1:39)),1)));
disp(nanmean(nanmean(include_secondhalf(:,ismember(days,1:39)),1)));
disp(signrank(nanmean(include_firsthalf(:,ismember(days,1:39)),1),nanmean(include_secondhalf(:,ismember(days,1:39)),1)));

end

function higherThanBefore=getDaysHigherThanBefore(avformouse,firsthalf,secondhalf)

useDayAv=false;
nanOutAfterExpert=true;

highestSoFar=-inf;
if nanOutAfterExpert==true
    highestSoFar=3;
end
higherThanBefore=zeros(size(avformouse));
for i=1:length(avformouse)
    if nanOutAfterExpert==true
        if avformouse(i)>highestSoFar
            higherThanBefore(i)=1;
        end
    elseif useDayAv==true
        if avformouse(i)>highestSoFar
            highestSoFar=avformouse(i);
            higherThanBefore(i)=1;
        end
    else
        if firsthalf(i)>0.75 || secondhalf(i)>0.75
            continue
        end
        if firsthalf(i)>highestSoFar || secondhalf(i)>highestSoFar
            if firsthalf(i)>highestSoFar
                highestSoFar=firsthalf(i);
                higherThanBefore(i)=1;
            elseif secondhalf(i)>highestSoFar
                highestSoFar=secondhalf(i);
                higherThanBefore(i)=1;
            end
        end
    end
end
if nanOutAfterExpert==true
    higherThanBefore=~higherThanBefore;
end

end