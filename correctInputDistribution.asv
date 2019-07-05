function correctInputDistribution(rt_pairs,matchTo)

if matchTo==2
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts2;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts1;
    reactionTimes1=rt_pairs.all_rt2;
    reactionTimes2=rt_pairs.all_rt1;
    real_rt_pairs1=rt_pairs.real_rt_pair2;
    real_rt_pairs2=rt_pairs.real_rt_pair1;
    rt_pairs1=rt_pairs.rt_pairs2;
    rt_pairs2=rt_pairs.rt_pairs1;
    forPairs_sequenceMatch1=rt_pairs.forPairs_sequenceMatchStarts2;
    forPairs_sequenceMatch2=rt_pairs.forPairs_sequenceMatchStarts1;
elseif matchTo==1
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts2;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts1;
    reactionTimes1=rt_pairs.all_rt2;
    reactionTimes2=rt_pairs.all_rt1;
    real_rt_pairs1=rt_pairs.real_rt_pair1;
    real_rt_pairs2=rt_pairs.real_rt_pair2;
    rt_pairs1=rt_pairs.rt_pairs1;
    rt_pairs2=rt_pairs.rt_pairs2;
    forPairs_sequenceMatch1=rt_pairs.forPairs_sequenceMatchStarts1;
    forPairs_sequenceMatch2=rt_pairs.forPairs_sequenceMatchStarts2;
end

% match reaction time distribution of templateSequence group 2 
% to reaction time distribution of templateSequence group 1
RT_binsize=0.01; % in seconds
maxRT=nanmax([reactionTimes1(sequenceMatchStarts1==1) reactionTimes2(sequenceMatchStarts2==1)])+RT_binsize;
RT_bins=0:RT_binsize:maxRT;
[n,x]=hist(reactionTimes1(sequenceMatchStarts1==1),RT_bins);
% get mapping of reaction times to bins for group 2
temp=reactionTimes2(sequenceMatchStarts2==1);
mapping_to_bins=nan(1,length(temp));
for i=1:length(temp)
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
        temp_inds=find(mapping_to_bins==i);
        if isempty(temp_inds)
            temp_inds=find(mapping_to_bins==i-1 | mapping_to_bins==i+1);
        end
        usenew=temp_inds(randsample(addn,1:length(temp_inds),true));
        group2_for_bins{i}=[find(mapping_to_bins==i) usenew];
    elseif curr_n<curr_group2
        % resample from existing group 2
        taken=curr_n;
        temp_inds=find(mapping_to_bins==i);
        if isempty(temp_inds)
            temp_inds=find(mapping_to_bins==i-1 | mapping_to_bins==i+1);
        end
        usethese=temp_inds(randsample(taken,1:length(temp_inds),true));
        group2_for_bins{i}=[usethese];
    else
        group2_for_bins{i}=[find(mapping_to_bins==i)];
    end
end
newgroup2_reactionTimes2=[];
newgroup2_sequenceMatchStarts2=[];
newgroup2_real_rt_pairs2=[];
newgroup2_rt_pairs2=[];
newgroup2_forPairs_sequenceMatch2=[];
for i=1:length(group2_for_bins)
    currindstouse=group2_for_bins{i};
    temp=reactionTimes2(sequenceMatchStarts2==1);
    newgroup2_reactionTimes2=[newgroup2_reactionTimes2 temp(currindstouse)];
    temp=reactionTimes2(sequenceMatchStarts2==1);
    newgroup2_sequenceMatchStarts2=[newgroup2_sequenceMatchStarts2 temp(currindstouse)];


    out.rt_pairs1=RT_pairs1.rt_change;
out.rt_pairs2=RT_pairs2.rt_change;
out.real_rt_pair1=RT_pairs1.real_rt_pair;
out.real_rt_pair2=RT_pairs2.real_rt_pair;
out.forPairs_sequenceMatchStarts1=forPairs_sequenceMatchStarts1;
out.forPairs_sequenceMatchStarts2=forPairs_sequenceMatchStarts2;
    
out.rt_pairs1_contingent=RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & forPairs_sequenceMatchStarts1==true);
out.rt_pairs2_contingent=RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & forPairs_sequenceMatchStarts2==true);        

% Flip back
if matchTo==2
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts2;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts1;
    reactionTimes1=rt_pairs.all_rt2;
    reactionTimes2=rt_pairs.all_rt1;
    real_rt_pairs1=rt_pairs.rt_pairs2_contingent;
    real_rt_pairs2=rt_pairs.rt_pairs1_contingent;
elseif matchTo==1
    sequenceMatchStarts1=rt_pairs.sequenceMatchStarts2;
    sequenceMatchStarts2=rt_pairs.sequenceMatchStarts1;
    reactionTimes1=rt_pairs.all_rt2;
    reactionTimes2=rt_pairs.all_rt1;
    real_rt_pairs1=rt_pairs.rt_pairs2_contingent;
    real_rt_pairs2=rt_pairs.rt_pairs1_contingent;
end