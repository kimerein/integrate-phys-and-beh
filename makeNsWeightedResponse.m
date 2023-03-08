function outResponse=makeNsWeightedResponse(r1,r2)

outResponse=r1;
figure(); 
plot(nanmean(r1.aligncomp_y,1),'Color','b'); 
hold on; plot(nanmean(r2.aligncomp_y,1),'Color','r');
disp('assumes both responses have some number of points taken before alignment peak');
takeThisMany=min(size(r1.unitbyunit_x,2),size(r2.unitbyunit_x,2));
sumns=r1.ns+r2.ns;
temp=cat(3,repmat(r1.ns./sumns,1,takeThisMany).*r1.unitbyunit_y(:,1:takeThisMany),repmat(r2.ns./sumns,1,takeThisMany).*r2.unitbyunit_y(:,1:takeThisMany));
outResponse.unitbyunit_y=reshape(nanmean(temp,3),size(r1.unitbyunit_y,1),takeThisMany);
outResponse.unitbyunit_x=outResponse.unitbyunit_x(:,1:takeThisMany);
outResponse.ns=sumns;
