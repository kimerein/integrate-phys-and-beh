function SVMforD1vA2a(D1mods,A2amods)

X=[[D1mods.modIndex1 D1mods.modIndex2 D1mods.modIndex3]; [A2amods.modIndex1 A2amods.modIndex2 A2amods.modIndex3]];
Y=[zeros(size(D1mods.modIndex1,1),1); ...
   ones(size(A2amods.modIndex1,1),1)];
Mdl=fitcsvm(X,Y,'KernelScale','auto','Standardize',true,'OutlierFraction',0.05);
sv=Mdl.SupportVectors;
figure;
gscatter(X(:,1),X(:,3),Y,'br','xo');
% hold on;
% plot(sv(:,1),sv(:,3),'ko','MarkerSize',10);
title('Input data View 1');
figure;
gscatter(X(:,1),X(:,2),Y,'br','xo');
% hold on;
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10);
title('Input data View 2');
figure;
gscatter(X(:,3),X(:,2),Y,'br','xo');
% hold on;
% plot(sv(:,3),sv(:,2),'ko','MarkerSize',10);
title('Input data View 3');
predictions=predict(Mdl,X);
figure;
gscatter(X(:,1),X(:,3),predictions,'rb','ox');
title('Prediction View 1');
figure;
gscatter(X(:,1),X(:,2),predictions,'rb','ox');
title('Prediction View 2');
figure;
gscatter(X(:,3),X(:,2),predictions,'rb','ox');
title('Prediction View 3');

end