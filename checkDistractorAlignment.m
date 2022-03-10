function checkDistractorAlignment(data1,whichTime1,whichField1,data2,whichTime2,whichField2)

dis1=data1.(whichField1);
t1=data1.(whichTime1);
dis2=data2.(whichField2);
t2=data2.(whichTime2);

figure();
offset=0;
for i=1:size(data1.times_wrt_trial_start,1)
    plot(t1(i,:),offset+dis1(i,:),'Color','k');
    hold on;
    plot(t2(i,:),offset+dis2(i,:),'Color','r');
    offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
end
    