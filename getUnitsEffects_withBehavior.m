function [unitFx,all_frs,supp]=getUnitsEffects_withBehavior(spikes,led_window,norm_window,useTrials,binsize)

% led_window gives start and end of opto on in seconds

goodunit_label=2; % will only save data for good units
% binsize=1; % in ms

unitlist=spikes.labels;
for i=1:length(unitlist)
    unitFx(i).no_led_FR=nan;
    unitFx(i).led_FR=nan;
    unitFx(i).isGoodUnit=0;
    unitFx(i).isSuppressed=false;
    unitFx(i).response_x=nan;
    unitFx(i).response_y=nan;
    unitFx(i).response_x_led=nan;
    unitFx(i).response_y_led=nan;
    if ~isempty(norm_window)
        unitFx(i).normedResponse_x=nan;
        unitFx(i).normedResponse_y=nan;
        unitFx(i).normedResponse_x_led=nan;
        unitFx(i).normedResponse_y_led=nan;
    end
end
unitAss=unitlist(:,1);
unitLabels=unitlist(:,2);
useInds=find(unitLabels==goodunit_label);
useUnits=unitAss(useInds);

if isempty(useTrials)
    useTrials=unique(spikes.trials);
    useTrials=useTrials(~isnan(useTrials));
end

all_frs=nan(length(useInds),2);
supp=nan(length(useInds),2);
for i=1:length(useInds)
    [x,y,x_led,y_led]=plotSUresponse(spikes,useUnits(i),binsize,useTrials,0);
    unitFx(useInds(i)).response_x=x;
    unitFx(useInds(i)).response_y=y;
    unitFx(useInds(i)).response_x_led=x_led;
    unitFx(useInds(i)).response_y_led=y_led;
    % get spiking activity in time window
    unitFx(useInds(i)).no_led_FR=nanmean(y(x>=led_window(1) & x<=led_window(2)));
    unitFx(useInds(i)).led_FR=nanmean(y_led(x_led>=led_window(1) & x_led<=led_window(2)));
    unitFx(useInds(i)).isGoodUnit=1;
    all_frs(i,1)=unitFx(useInds(i)).no_led_FR;
    all_frs(i,2)=unitFx(useInds(i)).led_FR;
    unitFx(useInds(i)).isSuppressed=all_frs(i,2)<all_frs(i,1);
    supp(i)=all_frs(i,2)<all_frs(i,1);
    if ~isempty(norm_window)
        norm=nanmean([y(x>=norm_window(1) & x<=norm_window(2)) y_led(x_led>=norm_window(1) & x_led<=norm_window(2))]);
        if norm~=0
            unitFx(useInds(i)).normedResponse_x=x;
            unitFx(useInds(i)).normedResponse_y=y./norm;
            unitFx(useInds(i)).normedResponse_x_led=x_led;
            unitFx(useInds(i)).normedResponse_y_led=y_led./norm;
        end
    end
    unitFx(useInds(i)).assigns=useUnits(i);
end

figure();
scatter(all_frs(:,1),all_frs(:,2));
xlabel('Firing rate no LED (Hz)');
ylabel('Firing rate LED (Hz)');
title('Effects of LED on single units');


    