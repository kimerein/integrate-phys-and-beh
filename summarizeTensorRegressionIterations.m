function summarizeTensorRegressionIterations(dir_by_rank,whatRank)

dd=dir(dir_by_rank);
dd=dd(3:end);
loss=nan(1,length(dd));
sumByFac=zeros(length(dd),whatRank);
mixedLoading=nan(length(dd),1);
for i=1:length(dd)
    if ~isempty(regexp(dd(i).name,'it'))
        a=load([dd(i).folder sep dd(i).name sep 'results.mat']);
        fac0=load([dd(i).folder sep dd(i).name sep 'factor_0.mat']);
        fac1=load([dd(i).folder sep dd(i).name sep 'factor_1.mat']);
        fac2=load([dd(i).folder sep dd(i).name sep 'factor_2.mat']);
    else
        continue
    end
    loss(i)=a.results.loss_train_final;
    % factor sums
    for j=1:size(fac0.factor,2)
        sumByFac(i,j)=sumByFac(i,j)+nansum(abs(fac0.factor(:,j)-nanmean(fac0.factor(:,j))));
    end
    for j=1:size(fac1.factor,2)
        sumByFac(i,j)=sumByFac(i,j)+nansum(abs(fac1.factor(:,j)-nanmean(fac1.factor(:,j))));
    end
    for j=1:size(fac2.factor,2)
        sumByFac(i,j)=sumByFac(i,j)+nansum(abs(fac2.factor(:,j)-nanmean(fac2.factor(:,j))));
    end
    % order factors by sumByFac
    sumByFac(i,:)=sort(sumByFac(i,:),'descend');
    % mixed loading
    mixedLoading(i)=getMixedLoading(fac0.factor);
end


figure();
lossBins=0.6:0.1:1.2;
cmap=brewermap(length(lossBins)-1,"RdYlGn");
cmap=flipud(cmap);
for i=1:length(loss)
    if loss(i)<lossBins(1)
        c=cmap(1,:);
    elseif loss(i)>lossBins(end)
        c=cmap(end,:);
    else
        currlossbin=find(lossBins>=loss(i),1,'first')-1;
        c=cmap(currlossbin,:);
    end
    scatter(1:size(sumByFac,2),sumByFac(i,:),[],c,'LineWidth',2);
    hold on;
    plot(1:size(sumByFac,2),sumByFac(i,:),'LineWidth',2,'Color',c);
end

figure();
scatter(whatRank*ones(size(mixedLoading)),mixedLoading,[],'k');
% scatter(whatRank*ones(size(mixedLoading(sumByFac(:,end)>=1))),mixedLoading(sumByFac(:,end)>=1),[],'k');
% hold on;
% scatter(whatRank*ones(size(mixedLoading(sumByFac(:,end)<1))),mixedLoading(sumByFac(:,end)<1),'filled','k');
line([whatRank-0.2 whatRank+0.2],[nanmean(mixedLoading) nanmean(mixedLoading)],'Color','k');

end

function overall_mixedLoadingPenalty=getMixedLoading(neuron_fac)

pairwiseDiffs=nan(size(neuron_fac,1),size(neuron_fac,2)*size(neuron_fac,2)-size(neuron_fac,2));
pairwiseSums=nan(size(neuron_fac,1),size(neuron_fac,2)*size(neuron_fac,2)-size(neuron_fac,2));
totalFac=nansum(nansum(neuron_fac));
for i=1:size(neuron_fac,1) % across neurons
    l=1;
    for j=1:size(neuron_fac,2)
        for k=1:size(neuron_fac,2) % across pairs of factors
            % skip if these are the same
            if j==k
                continue
            end
            pairwiseDiffs(i,l)=abs(neuron_fac(i,j))-abs(neuron_fac(i,k));
            pairwiseSums(i,l)=abs(neuron_fac(i,j))+abs(neuron_fac(i,k));
            l=l+1;
        end
    end
end
mixedLoadingPenalty=nansum(abs(pairwiseDiffs),2);
sumpairs=nansum(abs(pairwiseSums),2);
mix=nansum(mixedLoadingPenalty);
if mix<0.1
    mix=0.1;
end
overall_mixedLoadingPenalty=nansum(sumpairs)/mix;

end