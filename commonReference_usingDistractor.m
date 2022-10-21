function [beh1_tbt,beh2_tbt]=commonReference_usingDistractor(beh1_tbt,beh2_tbt,fieldtoalign)

% Note that I selected a subset of trials for each of beh1_tbt and
% beh2_tbt, depending on which trials were present in physiology and
% photometry
% HOWEVER, movie_distractor should be the same in both behavior_tbt's
% Thus, need to align trial indices

maxTrialsLag=20; % should not be more than this many trials missing in the middle
R_thresh=0.95; % threshold for correlation coefficient, if this is the right trial alignment
% note that beh1_tbt and beh2_tbt were originally derived from the same
% array, thus should be exactly matching when correctly aligned

data1=beh1_tbt.(fieldtoalign);
data2=beh2_tbt.(fieldtoalign);
backup_data1=data1;
backup_data2=data2;
data1_trialinds=repmat([1:size(data1,1)]',1,size(data1,2)).*ones(size(data1));
data2_trialinds=repmat([1:size(data2,1)]',1,size(data2,2)).*ones(size(data2));
backup_data1_trialinds=data1_trialinds;
backup_data2_trialinds=data2_trialinds;
data1=data1';
data2=data2';
data1_trialinds=data1_trialinds';
data2_trialinds=data2_trialinds';
data1=data1(1:end);
data2=data2(1:end);
data1_trialinds=data1_trialinds(1:end);
data2_trialinds=data2_trialinds(1:end);
data1(isnan(data1))=0;
data2(isnan(data2))=0;

% Start with a global alignment
% Pad according to global alignment, then refine
[Xa,Ya,D]=alignsignals(data1,data2);
if D>0
    linedUp_trialinds_Ya=data2_trialinds; linedUp_trialinds_Xa=[zeros(1,abs(D)) data1_trialinds];
    D_trials=linedUp_trialinds_Ya(find(linedUp_trialinds_Xa>0,1,'first'))-1;
    [backup_data1_trialinds,backup_data2_trialinds]=prepend_to_data1(backup_data1_trialinds,backup_data2_trialinds,D_trials);
    [backup_data1,backup_data2]=prepend_to_data1(backup_data1,backup_data2,D_trials);
elseif D<0
    linedUp_trialinds_Xa=data1_trialinds; linedUp_trialinds_Ya=[zeros(1,abs(D)) data2_trialinds];
    D_trials=linedUp_trialinds_Xa(find(linedUp_trialinds_Ya>0,1,'first'))-1;
    [backup_data2_trialinds,backup_data1_trialinds]=prepend_to_data1(backup_data2_trialinds,backup_data1_trialinds,D_trials);
    [backup_data2,backup_data1]=prepend_to_data1(backup_data2,backup_data1,D_trials);
else
    linedUp_trialinds_Xa=data1_trialinds; linedUp_trialinds_Ya=data2_trialinds;
end
disp('First argument is black, second is red');
figure(); plot(Xa,'Color','k'); hold on; plot(Ya,'Color','r'); title('Initial global alignment'); ylabel('Trial indices');
checkAlignment(backup_data1*1.5,backup_data2); title('After global alignment');

% But trials may have been dropped in the middle of session
% So then proceed to check for best aligned trial for each trial
% individually
correct_backup_data2_trialinds=backup_data2_trialinds;
for i=1:size(backup_data1,1)
    if i>size(backup_data1,1)
        break
    end
    if i>size(backup_data2,1)
        break
    end
    R=corrcoef(backup_data1(i,:),backup_data2(i,:));
    if isnan(R(1,2)) || R(1,2)<R_thresh
        % if R is nan or less than R_thresh, find correct trial
        test_R=nan(1,length(i-maxTrialsLag:i+maxTrialsLag));
        k=1;
        for j=i-maxTrialsLag:i+maxTrialsLag
            if j<1 || j>size(backup_data2,1)
                k=k+1;
                continue
            end
            tempd1=backup_data1(i,:);
            tempd1(isnan(tempd1))=0;
            tempd2=backup_data2(j,:);
            tempd2(isnan(tempd2))=0;
            temp=corrcoef(tempd1,tempd2);
            test_R(k)=temp(1,2);
            k=k+1;
        end
        % find best alignment
        [maxie,bestalign]=nanmax(test_R);
        if maxie>R_thresh
            correct_backup_data2_trialinds(i,:)=ones(size(correct_backup_data2_trialinds(i,:))).*(i-maxTrialsLag+(bestalign-1));
        else
            % trial likely does not exist in data2
            correct_backup_data2_trialinds(i,:)=nan;
        end
    end
end
% Adjust backup_data1 and backup_data2 according to correct trial-by-trial
% alignment
for i=1:size(backup_data1,1)
    if i>size(correct_backup_data2_trialinds,1)
        break
    end
    if isnan(mode(correct_backup_data2_trialinds(i,:)))
        backup_data2(i,:)=nan;
    else
        backup_data2(i,:)=backup_data2(mode(correct_backup_data2_trialinds(i,:)),:);
    end
end
checkAlignment(backup_data1*1.5,backup_data2); title('After trial-by-trial alignment');

beh1_tbt.reference_into_beh2trialinds=nan(size(beh1_tbt.(fieldtoalign)));
beh1_tbt.this_is_which_beh=1;
data1ti=mode(backup_data1_trialinds,2);
data2ti=mode(correct_backup_data2_trialinds,2);
for i=1:size(beh1_tbt.reference_into_beh2trialinds,1)
    f=find(data1ti==i);
    if isempty(f)
        continue
    end
    if f>length(data2ti) 
        continue
    end
    beh1_tbt.reference_into_beh2trialinds(i,:)=ones(size(beh1_tbt.reference_into_beh2trialinds(i,:)))*data2ti(f);
end
beh2_tbt.reference_into_beh1trialinds=nan(size(beh2_tbt.(fieldtoalign)));
beh2_tbt.this_is_which_beh=2;
for i=1:size(beh2_tbt.reference_into_beh1trialinds,1)
    f=find(data2ti==i);
    if isempty(f)
        continue
    end
    if f>length(data1ti)
        continue
    end
    if length(f)>1
        figure(); 
        plot(data2ti);
        disp(f);
        disambig=input(['Found two matches for row ' num2str(i) '. Which to use? ']);
        f=disambig;
    end
    beh2_tbt.reference_into_beh1trialinds(i,:)=ones(size(beh2_tbt.reference_into_beh1trialinds(i,:)))*data1ti(f);
end
% figure(); plot(1.1*behavior_tbt.movie_distractor(10,:),'Color','k'); hold on; plot(beh2_tbt.movie_distractor(beh2_tbt.reference_into_beh1trialinds(10,1),:),'Color','r');

end

function [backup_data1_trialinds,backup_data2_trialinds]=prepend_to_data1(backup_data1_trialinds,backup_data2_trialinds,D)

% always prepends to first argument

backup_data1_trialinds=[nan(abs(D),size(backup_data1_trialinds,2)); backup_data1_trialinds];
backup_data2_trialinds=[backup_data2_trialinds; nan(size(backup_data1_trialinds,1)-size(backup_data2_trialinds,1),size(backup_data2_trialinds,2))];

end

function checkAlignment(dis1,dis2)

figure();
offset=1;
dis1=dis1./max(dis1(1:end),[],2,'omitnan');
dis2=dis2./max(dis2(1:end),[],2,'omitnan');
for i=1:size(dis1,1)
    if i>size(dis1,1)
        disp('skipping end of dis1');
        break
    end
    if i>size(dis2,1)
        disp('skipping end of dis2');
        break
    end
    plot(offset+dis1(i,:),'Color','k');
    hold on;
    plot(offset+dis2(i,:),'Color','r');
    if isnan(nanmax([dis1(i,:) dis2(i,:)]))
    else
%         offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
        offset=offset+1;
    end
end

end