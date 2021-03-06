function [rt_pairs,used_as_rematch_inds]=correctInputDistribution(rt_pairs,matchTo,used_as_rematch_inds)

plotOutput=0;

if ~isempty(used_as_rematch_inds)
    if matchTo==1
        temp=rt_pairs.all_rt2(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.all_rt2=temp(used_as_rematch_inds);
        temp=rt_pairs.all_rt2_triali(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.all_rt2_triali=temp(used_as_rematch_inds);
        temp=rt_pairs.real_rt_pair2(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.real_rt_pair2=temp(used_as_rematch_inds);
        temp=rt_pairs.rt_pairs2(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.rt_pairs2=temp(used_as_rematch_inds);
        temp=rt_pairs.forPairs_sequenceMatchStarts2(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.forPairs_sequenceMatchStarts2=temp(used_as_rematch_inds);
        temp=rt_pairs.sequenceMatchStarts2(rt_pairs.sequenceMatchStarts2==1);
        rt_pairs.sequenceMatchStarts2=temp(used_as_rematch_inds);
        rt_pairs.rt_pairs2_contingent=rt_pairs.rt_pairs2(rt_pairs.real_rt_pair2==true & rt_pairs.forPairs_sequenceMatchStarts2==true);
    elseif matchTo==2
        temp=rt_pairs.all_rt1(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.all_rt1=temp(used_as_rematch_inds);
        temp=rt_pairs.all_rt1_triali(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.all_rt1_triali=temp(used_as_rematch_inds);
        temp=rt_pairs.real_rt_pair1(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.real_rt_pair1=temp(used_as_rematch_inds);
        temp=rt_pairs.rt_pairs1(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.rt_pairs1=temp(used_as_rematch_inds);
        temp=rt_pairs.forPairs_sequenceMatchStarts1(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.forPairs_sequenceMatchStarts1=temp(used_as_rematch_inds);
        temp=rt_pairs.sequenceMatchStarts1(rt_pairs.sequenceMatchStarts1==1);
        rt_pairs.sequenceMatchStarts1=temp(used_as_rematch_inds);
        rt_pairs.rt_pairs1_contingent=rt_pairs.rt_pairs1(rt_pairs.real_rt_pair1==true & rt_pairs.forPairs_sequenceMatchStarts1==true);
    end
    return
end
        
if matchTo==2
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts2;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts1;
    reactionTimes1=rt_pairs.all_rt2;
    reactionTimes2=rt_pairs.all_rt1;
    rt1_triali=rt_pairs.all_rt2_triali;
    rt2_triali=rt_pairs.all_rt1_triali;
    real_rt_pairs1=rt_pairs.real_rt_pair2;
    real_rt_pairs2=rt_pairs.real_rt_pair1;
    rt_pairs1=rt_pairs.rt_pairs2;
    rt_pairs2=rt_pairs.rt_pairs1;
    forPairs_sequenceMatch1=rt_pairs.forPairs_sequenceMatchStarts2;
    forPairs_sequenceMatch2=rt_pairs.forPairs_sequenceMatchStarts1;
elseif matchTo==1
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts1;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts2;
    reactionTimes1=rt_pairs.all_rt1;
    reactionTimes2=rt_pairs.all_rt2;
    rt1_triali=rt_pairs.all_rt1_triali;
    rt2_triali=rt_pairs.all_rt2_triali;
    real_rt_pairs1=rt_pairs.real_rt_pair1;
    real_rt_pairs2=rt_pairs.real_rt_pair2;
    rt_pairs1=rt_pairs.rt_pairs1;
    rt_pairs2=rt_pairs.rt_pairs2;
    forPairs_sequenceMatch1=rt_pairs.forPairs_sequenceMatchStarts1;
    forPairs_sequenceMatch2=rt_pairs.forPairs_sequenceMatchStarts2;
end

RT_binsize=0.033; % in seconds
maxRT=nanmax([reactionTimes1(sequenceMatchStarts1==1) reactionTimes2(sequenceMatchStarts2==1)])+RT_binsize;

if plotOutput==1
    RT_bins=0:RT_binsize:maxRT;
    [n,x]=hist(rt_pairs.all_rt1(rt_pairs.sequenceMatchStarts1==1),RT_bins);
    figure();
    plot(x,n,'Color','k');
    hold all;
    scatter(x,n,[],'k');
    hold on;
    [n,x]=hist(rt_pairs.all_rt2(rt_pairs.sequenceMatchStarts2==1),RT_bins);
    plot(x,n,'Color','r');
    title('Comparison RT distributions -- before');
end

% match reaction time distribution of templateSequence group 2 
% to reaction time distribution of templateSequence group 1
RT_bins=0:RT_binsize:maxRT;
[n,x]=hist(reactionTimes1(sequenceMatchStarts1==1),RT_bins);
% get mapping of reaction times to bins for group 2
temp=reactionTimes2(sequenceMatchStarts2==1);
mapping_to_bins=nan(1,length(temp));
for i=1:length(temp)
    if isnan(temp(i))
        continue
    end
    mapping_to_bins(i)=find(temp(i)>=RT_bins(1:end-1) & temp(i)<RT_bins(2:end));
end
% resample group 2 to match reaction time distribution of group 1
group2_for_bins=cell(1,length(n));
for i=1:length(n)
    curr_n=n(i);
    curr_group2=sum(mapping_to_bins==i);
    if curr_n>curr_group2
        % bootstrap by resampling from existing group 2
        addn=curr_n-curr_group2;
        if addn==0
            continue
        end
        temp_inds=find(mapping_to_bins==i);
        if isempty(temp_inds)
            for j=1:100
                temp_inds=find(mapping_to_bins==i-j | mapping_to_bins==i+j);
                if ~isempty(temp_inds)
                    break
                end
            end
        end
        if isempty(temp_inds)
            continue
        end
        usenew=temp_inds(randsample(length(temp_inds),addn,true));
        if size(usenew,1)>1
            usenew=usenew';
        end
        group2_for_bins{i}=[find(mapping_to_bins==i) usenew];
    elseif curr_n<curr_group2
        % resample from existing group 2
        taken=curr_n;
        if taken==0
            continue
        end
        temp_inds=find(mapping_to_bins==i);
        if isempty(temp_inds)
            for j=1:100
                temp_inds=find(mapping_to_bins==i-j | mapping_to_bins==i+j);
                if ~isempty(temp_inds)
                    break
                end
            end
        end
        if isempty(temp_inds)
            continue
        end
        usethese=temp_inds(randsample(length(temp_inds),taken,true));
        group2_for_bins{i}=[usethese];
    else
        group2_for_bins{i}=[find(mapping_to_bins==i)];
    end
end
newgroup2_reactionTimes2=[];
newgroup2_rt2_triali=[];
newgroup2_sequenceMatchStarts2=[];
newgroup2_real_rt_pairs2=[];
newgroup2_rt_pairs2=[];
newgroup2_forPairs_sequenceMatch2=[];
used_as_rematch_inds=[];
for i=1:length(group2_for_bins)
    currindstouse=group2_for_bins{i};
    if size(currindstouse,1)>1
        currindstouse=currindstouse';
    end
    used_as_rematch_inds=[used_as_rematch_inds currindstouse];
    temp=reactionTimes2(sequenceMatchStarts2==1);
    newgroup2_reactionTimes2=[newgroup2_reactionTimes2 temp(currindstouse)];
    temp=rt2_triali(sequenceMatchStarts2==1);
    newgroup2_rt2_triali=[newgroup2_rt2_triali temp(currindstouse)];
    temp=real_rt_pairs2(sequenceMatchStarts2==1);
    newgroup2_real_rt_pairs2=[newgroup2_real_rt_pairs2 temp(currindstouse)];
    temp=rt_pairs2(sequenceMatchStarts2==1);
    newgroup2_rt_pairs2=[newgroup2_rt_pairs2 temp(currindstouse)];
    temp=forPairs_sequenceMatch2(sequenceMatchStarts2==1);
    newgroup2_forPairs_sequenceMatch2=[newgroup2_forPairs_sequenceMatch2 temp(currindstouse)];
    temp=sequenceMatchStarts2(sequenceMatchStarts2==1);
    newgroup2_sequenceMatchStarts2=[newgroup2_sequenceMatchStarts2 temp(currindstouse)];
end

reactionTimes2=newgroup2_reactionTimes2;
rt2_triali=newgroup2_rt2_triali;
real_rt_pairs2=newgroup2_real_rt_pairs2;
rt_pairs2=newgroup2_rt_pairs2;
forPairs_sequenceMatch2=newgroup2_forPairs_sequenceMatch2;
sequenceMatchStarts2=newgroup2_sequenceMatchStarts2;

% Flip back
if matchTo==2
    rt_pairs.sequenceMatchStarts2=sequenceMatchStarts1;
    rt_pairs.sequenceMatchStarts1=sequenceMatchStarts2;
    rt_pairs.all_rt2=reactionTimes1;
    rt_pairs.all_rt1=reactionTimes2;
    rt_pairs.all_rt2_triali=rt1_triali;
    rt_pairs.all_rt1_triali=rt2_triali;
    rt_pairs.real_rt_pair2=real_rt_pairs1;
    rt_pairs.real_rt_pair1=real_rt_pairs2;
    rt_pairs.rt_pairs2_contingent=rt_pairs1(real_rt_pairs1==true & forPairs_sequenceMatch1==true);
    rt_pairs.rt_pairs1_contingent=rt_pairs2(real_rt_pairs2==true & forPairs_sequenceMatch2==true);
elseif matchTo==1 
    rt_pairs.sequenceMatchStarts1=sequenceMatchStarts1;
    rt_pairs.sequenceMatchStarts2=sequenceMatchStarts2;
    rt_pairs.all_rt1=reactionTimes1;
    rt_pairs.all_rt2=reactionTimes2;
    rt_pairs.all_rt1_triali=rt1_triali;
    rt_pairs.all_rt2_triali=rt2_triali;
    rt_pairs.real_rt_pair1=real_rt_pairs1;
    rt_pairs.real_rt_pair2=real_rt_pairs2;
    rt_pairs.rt_pairs1_contingent=rt_pairs1(real_rt_pairs1==true & forPairs_sequenceMatch1==true);
    rt_pairs.rt_pairs2_contingent=rt_pairs2(real_rt_pairs2==true & forPairs_sequenceMatch2==true);
end

% Plot if worked
if plotOutput==1
    RT_bins=0:RT_binsize:maxRT;
    [n,x]=hist(rt_pairs.all_rt1(rt_pairs.sequenceMatchStarts1==1),RT_bins);
    figure();
    plot(x,n,'Color','k');
    hold all;
    scatter(x,n,[],'k');
    hold on;
    [n,x]=hist(rt_pairs.all_rt2(rt_pairs.sequenceMatchStarts2==1),RT_bins);
    plot(x,n,'Color','r');
    title('Comparison RT distributions -- after');
end

