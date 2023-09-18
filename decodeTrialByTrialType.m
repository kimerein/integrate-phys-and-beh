function decodeTrialByTrialType(tbyt_cuedsuccess_vs_cuedfailure,tbyt_uncuedsuccess_vs_uncuedfailure,idx,nBoot,nUnits,nTrials,withReplacement)

% Take a random sampling of trials from n units of idx==1, take average
% Idx 1 and idx 2 units not necessarily acquired simultaneously, so
% Take a different random sampling of trials from n units of idx==2, take average

unique_units=unique([unique(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success); unique(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure); unique(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success); unique(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure)]);
% unit numbers in unique_units correspond to units 1 to 1063 in
% cued_success_Response
if length(idx)~=1063
    error('dont recognize this idx');
end
unique_units_idx=1:1063;
gp1_unit_ids=unique_units_idx(idx==1);
gp2_unit_ids=unique_units_idx(idx==2);

% cued success sample
[trialIDsPerSample_gp1_cuedsuccess,unitsPerSample_gp1_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,withReplacement);
% cued failure sample
[trialIDsPerSample_gp1_cuedfailure,unitsPerSample_gp1_cuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,withReplacement);
% uncued success sample
[trialIDsPerSample_gp1_uncuedsuccess,unitsPerSample_gp1_uncuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,withReplacement);
% uncued failure sample
[trialIDsPerSample_gp1_uncuedfailure,unitsPerSample_gp1_uncuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,withReplacement);

% cued success sample
[trialIDsPerSample_gp2_cuedsuccess,unitsPerSample_gp2_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,withReplacement);
% cued failure sample
[trialIDsPerSample_gp2_cuedfailure,unitsPerSample_gp2_cuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,withReplacement);
% uncued success sample
[trialIDsPerSample_gp2_uncuedsuccess,unitsPerSample_gp2_uncuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,withReplacement);
% uncued failure sample
[trialIDsPerSample_gp2_uncuedfailure,unitsPerSample_gp2_uncuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,withReplacement);



temp1=nan(nBoot,nTrials); temp2=nan(nBoot,nTrials);
for i=1:nBoot
    temp1(i,:)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp1_cuedsuccess(i,:)))); 
    temp2(i,:)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp2_cuedsuccess(i,:)))); 
end
figure();
scatter(-temp1,temp2,[],'g'); hold on; scatter(-nanmean(temp1),nanmean(temp2),[],'g','filled'); cuedsuccmeanx=nanmean(temp1); cuedsuccmeany=nanmean(temp2);

temp1=nan(nBoot,nTrials); temp2=nan(nBoot,nTrials);
for i=1:nBoot
    temp1(i,:)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp1_cuedfailure(i,:)))); 
    temp2(i,:)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp2_cuedfailure(i,:)))); 
end
scatter(-temp1,temp2,[],'r'); scatter(-nanmean(temp1),nanmean(temp2),[],'r','filled'); cuedfailmeanx=nanmean(temp1); cuedfailmeany=nanmean(temp2);

temp1=nan(nBoot,nTrials); temp2=nan(nBoot,nTrials);
for i=1:nBoot
    temp1(i,:)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp1_uncuedsuccess(i,:)))); 
    temp2(i,:)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp2_uncuedsuccess(i,:)))); 
end
scatter(-temp1,temp2,[],'b'); scatter(-nanmean(temp1),nanmean(temp2),[],'b','filled'); uncuedsuccmeanx=nanmean(temp1); uncuedsuccmeany=nanmean(temp2);

temp1=nan(nBoot,nTrials); temp2=nan(nBoot,nTrials);
for i=1:nBoot
    temp1(i,:)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp1_uncuedfailure(i,:)))); 
    temp2(i,:)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp2_uncuedfailure(i,:)))); 
end
scatter(-temp1,temp2,[],'y'); scatter(-nanmean(temp1),nanmean(temp2),[],'y','filled'); uncuedfailmeanx=nanmean(temp1); uncuedfailmeany=nanmean(temp2);
scatter(-(cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');

end

function [trialIDsPerSample_gp1_cuedsuccess,unitsPerSample_gp1_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,currtbyt_trialIDs,currtbyt_units,withReplacement)

% currtbyt_trialIDs=tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success
% currtbyt_units=tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success

% Random samples of nUnits of idx==1 from nUnits random single trials
unitsPerSample_gp1_cuedsuccess=nan(nBoot,nUnits);
trialIDsPerSample_gp1_cuedsuccess=nan(nBoot,nTrials);
for i=1:nBoot
    % rand sample of units
    unitsPerSample_gp1_cuedsuccess(i,:)=gp1_unit_ids(randsample(length(gp1_unit_ids),nUnits,withReplacement));
    % for these units, take a random sample of trials
    currcondition_trialIDs_fortheseunits=currtbyt_trialIDs(ismember(currtbyt_units,unitsPerSample_gp1_cuedsuccess(i,:)));
    % ensure only 1 trial per unit
    unitIDs_for_thesetrialIDs=currtbyt_units(ismember(currtbyt_units,unitsPerSample_gp1_cuedsuccess(i,:)));
    for j=1:length(unitIDs_for_thesetrialIDs)
        currunitid=unitIDs_for_thesetrialIDs(j);
        if nansum(unitIDs_for_thesetrialIDs==currunitid)>1
            % randomly choose one of the trials
            f=find(unitIDs_for_thesetrialIDs==currunitid);
            omitthesetrials=f(randsample(length(f),length(f)-1));
            unitIDs_for_thesetrialIDs(omitthesetrials)=nan;
            currcondition_trialIDs_fortheseunits(omitthesetrials)=nan;
        end
    end
    % drop nans
    keepnotnan=~isnan(currcondition_trialIDs_fortheseunits);
    currcondition_trialIDs_fortheseunits=currcondition_trialIDs_fortheseunits(keepnotnan);
    unitIDs_for_thesetrialIDs=unitIDs_for_thesetrialIDs(keepnotnan);
    trialIDsPerSample_gp1_cuedsuccess(i,:)=currcondition_trialIDs_fortheseunits(randsample(length(currcondition_trialIDs_fortheseunits),nTrials,withReplacement));
end

end