function out=plot_behAndPhys_averages(unitFx,behFx,ds_forReaches)

% Get non-norm and norm phys data, separated by whether unit was suppressed
% or facilitated

j=1;
for i=1:length(unitFx)
    if isnan(unitFx(i).no_led_FR) % not a good unit
        % skip
        continue
    end
    if unitFx(i).no_led_FR>=unitFx(i).led_FR % unit is suppressed
        supp(j)=1; % yes, suppressed
    else
        supp(j)=0; % no, facilitated
    end
    if j==1
        x=unitFx(i).response_x;
    end
    if j~=1
        if length(unitFx(i).response_y)<length(y(j-1,:))
            y(j,:)=[unitFx(i).response_y nan];
            y_led(j,:)=[unitFx(i).response_y_led nan];
            y_norm(j,:)=[unitFx(i).normedResponse_y nan];
            y_norm_led(j,:)=[unitFx(i).normedResponse_y_led nan];
        elseif length(unitFx(i).response_y)>length(y(j-1,:))
            y(j,:)=unitFx(i).response_y(1:length(y(j-1,:)));
            y_led(j,:)=unitFx(i).response_y_led(1:length(y(j-1,:)));
            y_norm(j,:)=unitFx(i).normedResponse_y(1:length(y(j-1,:)));
            y_norm_led(j,:)=unitFx(i).normedResponse_y_led(1:length(y(j-1,:)));
        else
            y(j,:)=unitFx(i).response_y;
            y_led(j,:)=unitFx(i).response_y_led;
            y_norm(j,:)=unitFx(i).normedResponse_y;
            y_norm_led(j,:)=unitFx(i).normedResponse_y_led;
        end
    else
        y(j,:)=unitFx(i).response_y;
        y_led(j,:)=unitFx(i).response_y_led;
        y_norm(j,:)=unitFx(i).normedResponse_y;
        y_norm_led(j,:)=unitFx(i).normedResponse_y_led;
    end
    j=j+1;
end

% % Plot behavior with and without opto
% figure();
% temp=downSampAv(behFx.reaches_no_led,ds_forReaches);
% ma=max(temp);
% plot(downSampAv(behFx.x,ds_forReaches),temp,'Color','k');
% hold on;
% behFx.cue=behFx.cue-min(behFx.cue);
% plot(downSampAv(behFx.x,ds_forReaches),(downSampAv(behFx.cue,ds_forReaches)./max(downSampAv(behFx.cue,ds_forReaches))).*ma,'Color','b');
% xlabel('Time (s)');
% ylabel('# reaches');
% title('Reaches without opto');
% 
% figure();
% temp=downSampAv(behFx.reaches_led,ds_forReaches);
% ma=max(temp);
% plot(downSampAv(behFx.x,ds_forReaches),temp,'Color','k');
% hold on;
% plot(downSampAv(behFx.x,ds_forReaches),(downSampAv(behFx.cue,ds_forReaches)./max(downSampAv(behFx.cue,ds_forReaches))).*ma,'Color','b');
% behFx.opto=behFx.opto-min(behFx.opto);
% plot(downSampAv(behFx.x,ds_forReaches),(downSampAv(behFx.opto,ds_forReaches)./max(downSampAv(behFx.opto,ds_forReaches))).*ma,'Color','r');
% xlabel('Time (s)');
% ylabel('# reaches');
% title('Reaches with opto');

suppis=1;
figure();
getsupp=find(supp==suppis);
for i=1:length(getsupp)
    plot(x,y(getsupp(i),:),'Color','k');
    hold on;
    plot(x,y_led(getsupp(i),:),'Color','r');
end

% Plot spikes with and without opto 
% Suppressed units
takeCued=[1 3 6 7 10 16 23 34 46 53 58];
temp=1:58;
unCued=temp(~ismember(temp,takeCued));
isSuppInds=find(supp==suppis);
useThese=isSuppInds(ismember(isSuppInds,unCued));
figure();
plot(x,nanmean(y(useThese,:),1),'Color','k');
hold on;
plot(x,nanmean(y_led(useThese,:),1),'Color','r');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title('Spikes, suppressed units');
out.x=x;
out.y=nanmean(y(supp==suppis,:),1);
out.y_led=nanmean(y_led(supp==suppis,:),1);

% Facilitated units
suppis=0;
figure();
plot(x,nanmean(y(supp==suppis,:),1),'Color','k');
hold on;
plot(x,nanmean(y_led(supp==suppis,:),1),'Color','r');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title('Spikes, facilitated units');


% Plot spikes with and without opto -- normalized
% Suppressed units
suppis=1;
figure();
% fill([x fliplr(x)],[nanmean(y_norm(supp==suppis,:),1)+nanstd(y_norm(supp==suppis,:),[],1)./sqrt(sum(supp==suppis)) fliplr(nanmean(y_norm(supp==suppis,:),1)-nanstd(y_norm(supp==suppis,:),[],1)./sqrt(sum(supp==suppis)))],[0.8 0.8 0.8]);
% hold on;
plot(x,nanmean(y_norm(supp==suppis,:),1),'Color','k');
hold on;
% fill([x fliplr(x)],[nanmean(y_norm_led(supp==suppis,:),1)+nanstd(y_norm_led(supp==suppis,:),[],1) fliplr(nanmean(y_norm_led(supp==suppis,:),1)-nanstd(y_norm_led(supp==suppis,:),[],1))],[0.95 0.6 0.6]);
plot(x,nanmean(y_norm_led(supp==suppis,:),1),'Color','r');
xlabel('Time (s)');
ylabel('Norm firing rate');
title('Spikes normalized, suppressed units');

% Facilitated units
suppis=0;
figure();
% fill([x fliplr(x)],[nanmean(y_norm(supp==suppis,:),1)+nanstd(y_norm(supp==suppis,:),[],1)./sqrt(sum(supp==suppis)) fliplr(nanmean(y_norm(supp==suppis,:),1)-nanstd(y_norm(supp==suppis,:),[],1)./sqrt(sum(supp==suppis)))],[0.8 0.8 0.8]);
% hold on;
plot(x,nanmean(y_norm(supp==suppis,:),1),'Color','k');
hold on;
% fill([x fliplr(x)],[nanmean(y_norm_led(supp==suppis,:),1)+nanstd(y_norm_led(supp==suppis,:),[],1) fliplr(nanmean(y_norm_led(supp==suppis,:),1)-nanstd(y_norm_led(supp==suppis,:),[],1))],[0.95 0.6 0.6]);
plot(x,nanmean(y_norm_led(supp==suppis,:),1),'Color','r');
xlabel('Time (s)');
ylabel('Norm firing rate');
title('Spikes normalized, facilitated units');






