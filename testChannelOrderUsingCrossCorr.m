function rs=testChannelOrderUsingCrossCorr(channel_files)

rs=nan(length(channel_files),length(channel_files));
datas=cell(1,32);
for i=1:length(channel_files)
    if i==1
        a=load(channel_files{i});
        data1=a.data;
        datas{i}=data1;
    else
        data1=datas{i};
    end
    for j=i+1:length(channel_files)
        disp(['comparing ch ' num2str(i) ' and ch ' num2str(j)]);
        if i==1
            a=load(channel_files{j});
            data2=a.data;
            datas{j}=data2;
        else
            data2=datas{j};
        end
        temp=corrcoef(data1.Values,data2.Values);
        rs(i,j)=temp(1,2);
    end
end

for i=1:length(channel_files)
    for j=1:i-1
        rs(i,j)=rs(j,i);
    end
end

[maxCorrPerCh,indMaxPerCh]=nanmax(rs,[],2);
starterChs=1:32;
starterChs=starterChs';
[maxSorted,si]=sort(maxCorrPerCh);
sortedInds=indMaxPerCh(si);
sortedStarters=starterChs(si);
figure();
imagesc(maxSorted);
disp([sortedStarters sortedInds]);


figure();
imagesc(rs);

% testOrder=[13 14 10 12 11 6 9 8 2 4 32 15 7 5 3 1 30 28 26 24 31 29 16 17 27 18 25 19 22 21 23 20]; % super to deep
testOrder=[13 14 10 12 11 9 8 6 15 2 32 4 7 5 3 1 24 30 28 26 31 29 16 17 27 18 25 19 22 21 20 23]; % super to deep

reordered_rs=rs(testOrder,:);
figure();
imagesc(reordered_rs);

end