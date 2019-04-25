function probabilityOfDifferentTrialTypes(alltbt,metadata,trialTypes)

n_trials_back=3;

% exclusive groups for the sake of this analysis

% group 1A
% first reach after cue is
% success, drop, miss or no pellet
f_1A{1}='after_cue_success';
f_1A{2}='after_cue_drop';
f_1A{3}='after_cue_miss';
f_1A{4}='after_cue_no_pellet';
trialTypes=add_negate_field(trialTypes,f_1A);

% group 1B
% did mouse eventually touch pellet
f_1B{1}='touched_pellet';
f_1B{2}='not_touched_pellet';
trialTypes=add_negate_field(trialTypes,f_1B);

% group 1C
% did mouse eventually consume pellet
f_1C{1}='consumed_pellet';
f_1C{2}='not_consumed_pellet';
trialTypes=add_negate_field(trialTypes,f_1C);

% group 1D
% did mouse touch pellet in cue window
f_1D{1}='touch_in_cued_window';
f_1D{2}='not_touch_in_cued_window';
trialTypes=add_negate_field(trialTypes,f_1D);

% group 2
% something about reach timing
% paw out during wheel, cued reach or reach after cue
f_2{1}='paw_during_wheel';
f_2{2}='cued_reach';
f_2{3}='not_cued_reach';
trialTypes=add_negate_field(trialTypes,f_2);

% group 3
% DMS-tail silencing on or off
f_3{1}='led';
f_3{2}='not_led';
trialTypes=add_negate_field(trialTypes,f_3);

% group 4
% mouse still chewing at trial start or
% not chewing at trial start
f_4{1}='chewing_at_trial_start';
f_4{2}='not_chewing_at_trial_start';
trialTypes=add_negate_field(trialTypes,f_4);

% group 5
% ith trial

% rank trial type by probability of each type

% group 1A w groups 2-4
% n_groups=4*3*2*2;
% curr_f_1=f_1A;
% doAnalysis(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes,alltbt,metadata);

% group 1B w groups 2-4
% n_groups=2*3*2*2;
% curr_f_1=f_1B;
% doAnalysis(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes,alltbt,metadata);

% group 1C w groups 2-4
n_groups=2*3*2*2;
curr_f_1=f_1C;
doAnalysis(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes,alltbt,metadata);

% group 1D w groups 2-4
% n_groups=2*3*2*2;
% curr_f_1=f_1D;
% doAnalysis(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes,alltbt,metadata);

end

function doAnalysis(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes,alltbt,metadata)

[group_names,trialType_count]=countTrials(curr_f_1,f_2,f_3,f_4,n_groups,trialTypes);
[~,group_names]=rankAndDisplayByProbability(group_names,trialType_count,size(trialTypes.(curr_f_1{1}),1),1);    
% do p-value comparison for these trial types with and without led on
% PREVIOUS TRIAL
n_groups=length(curr_f_1)*length(f_2)*length(f_4);
[group_names,trialType_count]=countTrials_3fields(curr_f_1,f_2,f_4,n_groups,trialTypes);
[~,group_names]=rankAndDisplayByProbability(group_names,trialType_count,size(trialTypes.(curr_f_1{1}),1),0);  
for i=1:length(group_names)
    gpname=group_names{i};
    temp=group_names{i};
    % take apart group_names into fields
    firstAmp=regexp(temp,' & ');
    if isempty(firstAmp)
        error('group_names improperly formatted');
    end
    temp=temp(firstAmp+4:end);
    secondAmp=regexp(temp,' & ');
    if isempty(secondAmp)
        error('group_names improperly formatted');
    end
    secondAmp=secondAmp+2+firstAmp;
    templateSequence_noLED{1}=trialTypes.(gpname(1:firstAmp-1))==1 & trialTypes.(gpname(firstAmp+3:secondAmp))==1 & trialTypes.(gpname(secondAmp+4:end))==1 & trialTypes.led==0; % previous trial
    templateSequence_noLED{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0; % current trial
    templateSequence_LED{1}=trialTypes.(gpname(1:firstAmp-1))==1 & trialTypes.(gpname(firstAmp+3:secondAmp))==1 & trialTypes.(gpname(secondAmp+4:end))==1 & trialTypes.led==1; % previous trial
    templateSequence_LED{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0; % current trial
    [p_RT_pairs,p_RTs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence_noLED,templateSequence_LED);
    alpha=0.2;
    if p_RT_pairs<alpha
        disp([gpname ' RT pairs p-val is ' num2str(p_RT_pairs)]);
    end
    if p_RTs<alpha
        disp([gpname ' RTs p-val is ' num2str(p_RTs)]);
    end
end

end

function [trialType_count,group_names]=rankAndDisplayByProbability(group_names,trialType_count,totalTrials,dispStuff)

trialType_count=trialType_count./totalTrials;
[trialType_count,i]=sort(trialType_count,'descend');
group_names=group_names(i);
if dispStuff==1
    disp('Probability of each trial type');
    for i=1:length(trialType_count)
        disp([num2str(trialType_count(i)) ' ' group_names{i}]);
    end
end

end

function [group_names,trialType_count]=countTrials_3fields(f_1,f_2,f_3,n_groups,trialTypes)

group_names=cell(1,length(n_groups));
trialType_count=nan(1,n_groups);
m=1;
for i=1:length(f_1)
    curr_field_1=f_1{i};
    for j=1:length(f_2)
        curr_field_2=f_2{j};
        for k=1:length(f_3)
            curr_field_3=f_3{k};
            trialType_count(m)=nansum(trialTypes.(curr_field_1)==1 & trialTypes.(curr_field_2)==1 & trialTypes.(curr_field_3)==1);
            group_names{m}=[curr_field_1 ' & ' curr_field_2 ' & ' curr_field_3];
            m=m+1;
        end
    end
end
trialType_count(isnan(trialType_count))=0;

end

function [group_names,trialType_count]=countTrials(f_1,f_2,f_3,f_4,n_groups,trialTypes)

group_names=cell(1,length(n_groups));
trialType_count=nan(1,n_groups);
m=1;
for i=1:length(f_1)
    curr_field_1=f_1{i};
    for j=1:length(f_2)
        curr_field_2=f_2{j};
        for k=1:length(f_3)
            curr_field_3=f_3{k};
            for l=1:length(f_4)
                curr_field_4=f_4{l};
                trialType_count(m)=nansum(trialTypes.(curr_field_1)==1 & trialTypes.(curr_field_2)==1 & trialTypes.(curr_field_3)==1 & trialTypes.(curr_field_4)==1);
                group_names{m}=[curr_field_1 ' & ' curr_field_2 ' & ' curr_field_3 ' & ' curr_field_4];
                m=m+1;
            end
        end
    end
end
trialType_count(isnan(trialType_count))=0;

end

function trialTypes=add_negate_field(trialTypes,f)

for i=1:length(f)
    curr_field=f{i};
    if length(curr_field)>=4
        if strcmp(curr_field(1:4),'not_')
            % this is a negation of an existing field in trialTypes
            trialTypes.(curr_field)=trialTypes.(curr_field(5:end))==0;
        end
    end
end

end