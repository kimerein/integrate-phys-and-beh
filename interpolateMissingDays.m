function [interp_reachingout,angleout]=interpolateMissingDays(isreaching_out,first_sess)

doScatter_all_days=false;
find_closest_matched_reach_rate=true;
plot_begin_to_end_vector=true;

cmap=colormap('cool');
day_range=1:20;
cmap_range=size(cmap,1);
day_maps=nan(length(day_range),3);
for i=1:length(day_range)
    f=ceil(cmap_range/length(day_range))*(i-1)+1;
    if f>cmap_range
        f=cmap_range;
    end
    day_maps(i,:)=cmap(f,:);
end

nth_session=isreaching_out.nth_session;
cued_reach_rate=isreaching_out.cued_reach_rate;
noncued_reach_rate=isreaching_out.noncued_reach_rate;
nth_session=nth_session(find(isreaching_out.nth_session==first_sess,1,'first'):end);
cued_reach_rate=cued_reach_rate(find(isreaching_out.nth_session==first_sess,1,'first'):end);
noncued_reach_rate=noncued_reach_rate(find(isreaching_out.nth_session==first_sess,1,'first'):end);

interp_cued=interp1(nth_session,cued_reach_rate,nth_session(1):nth_session(end),'linear');
interp_noncued=interp1(nth_session,noncued_reach_rate,nth_session(1):nth_session(end),'linear');
interp_nth_session=nth_session(1):nth_session(end);
interp_nth_session=interp_nth_session-min(interp_nth_session)+1;

if find_closest_matched_reach_rate==true
    test_vals=0:0.05:1.5;
    distances=nan(1,length(test_vals));
    for i=1:length(test_vals)
        distances(i)=sqrt((interp_cued(1)-test_vals(i))^2+(interp_noncued(1)-test_vals(i))^2);
    end
    [~,mi]=min(distances);
    startPoint=[test_vals(mi) test_vals(mi)];
    interp_nth_session=[0 interp_nth_session];
    interp_cued=[startPoint(2) interp_cued];
    interp_noncued=[startPoint(1) interp_noncued];
end

interp_reachingout.interp_nth_session=interp_nth_session;
interp_reachingout.interp_cued=interp_cued;
interp_reachingout.interp_noncued=interp_noncued;

figure();
for i=1:length(interp_nth_session)-1
    %line([interp_noncued(i) interp_noncued(i+1)],[interp_cued(i) interp_cued(i+1)],'Color',day_maps(i,:));
    if i>size(day_maps,1)
        d=day_maps(end,:);
    else
        d=day_maps(i,:);
    end
    quiver(interp_noncued(i),interp_cued(i),interp_noncued(i+1)-interp_noncued(i),interp_cued(i+1)-interp_cued(i),'Color',d);
    hold on;
    if doScatter_all_days==true
        scatter(interp_noncued(i),interp_cued(i),[],d);
    end
end
if i+1>size(day_maps,1)
    d=day_maps(end,:);
else
    d=day_maps(i+1,:);
end
scatter(interp_noncued(i+1),interp_cued(i+1),[],d);
% scatter(interp_noncued(1),interp_cued(1),[],day_maps(1,:));
% scatter(interp_noncued(end),interp_cued(end),[],day_maps(end,:));
xlabel('Non-cued reaching rate (Hz)');
ylabel('Cued reaching rate (Hz)');

if plot_begin_to_end_vector==true
    quiver(interp_noncued(1),interp_cued(1),interp_noncued(end)-interp_noncued(1),interp_cued(end)-interp_cued(1),'Color','k');
    [theta,rho]=cart2pol(interp_noncued(end)-interp_noncued(1),interp_cued(end)-interp_cued(1));
    disp('theta');
    disp(rad2deg(theta));
    angleout=rad2deg(theta);
end