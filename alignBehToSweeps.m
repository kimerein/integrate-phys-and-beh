function alignBehToSweeps(beh_sweeps,spikes_sweeps)

% Assumes led is random
% Try to align based on led
[X,Y,D]=alignsignals(beh_sweeps.led,spikes_sweeps.led); 
figure();
plot(X,'Color','r');
hold on;
plot(Y,'Color','b');
leg={'beh_led','spikes_led'};
xlabel('indices');
ylabel('led on is 1');
title('Aligning beh_led to spikes_led');
legend(leg);

figure();
plot(Y,'Color','b');
hold on;
plot(X,'Color','r');
leg={'spikes_led','beh_led'};
xlabel('indices');
ylabel('led on is 1');
title('Aligning beh_led to spikes_led');