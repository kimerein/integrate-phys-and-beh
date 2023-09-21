function decodeTrialByTrialType(tbyt_cuedsuccess_vs_cuedfailure,tbyt_uncuedsuccess_vs_uncuedfailure,idx,nBoot,nUnits,nTrials,withReplacement,addThirdAxis,nanallzeros,justBootstrapTrials,collapseWithinUnit)

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

% average all trials, all units
temp1=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,gp1_unit_ids))); 
temp2=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,gp2_unit_ids)));
figure();
scatter(temp1,temp2,[],'g'); hold on;
temp1=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,gp1_unit_ids))); 
temp2=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,gp2_unit_ids)));
scatter(temp1,temp2,[],'r'); 
temp1=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,gp1_unit_ids))); 
temp2=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,gp2_unit_ids)));
scatter(temp1,temp2,[],'b'); 
temp1=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,gp1_unit_ids))); 
temp2=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,gp2_unit_ids)));
scatter(temp1,temp2,[],'y'); 

if collapseWithinUnit==true
    % trying to figure out why unit average does not look like trial
    % average
    u=unique(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success);
    for i=1:length(u)
        curru=u(i);
        temp=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,curru)));
        tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,curru))=temp;
    end
    u=unique(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure);
    for i=1:length(u)
        curru=u(i);
        temp=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,curru)));
        tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,curru))=temp;
    end
    u=unique(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success);
    for i=1:length(u)
        curru=u(i);
        temp=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,curru)));
        tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,curru))=temp;
    end
    u=unique(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure);
    for i=1:length(u)
        curru=u(i);
        temp=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,curru)));
        tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,curru))=temp;
    end
end

% cued success sample
[trialIDsPerSample_gp1_cuedsuccess,unitsPerSample_gp1_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,withReplacement,justBootstrapTrials);
% cued failure sample
[trialIDsPerSample_gp1_cuedfailure,unitsPerSample_gp1_cuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,withReplacement,justBootstrapTrials);
% uncued success sample
[trialIDsPerSample_gp1_uncuedsuccess,unitsPerSample_gp1_uncuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,withReplacement,justBootstrapTrials);
% uncued failure sample
[trialIDsPerSample_gp1_uncuedfailure,unitsPerSample_gp1_uncuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,withReplacement,justBootstrapTrials);

% cued success sample
[trialIDsPerSample_gp2_cuedsuccess,unitsPerSample_gp2_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,withReplacement,justBootstrapTrials);
% cued failure sample
[trialIDsPerSample_gp2_cuedfailure,unitsPerSample_gp2_cuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,withReplacement,justBootstrapTrials);
% uncued success sample
[trialIDsPerSample_gp2_uncuedsuccess,unitsPerSample_gp2_uncuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,withReplacement,justBootstrapTrials);
% uncued failure sample
[trialIDsPerSample_gp2_uncuedfailure,unitsPerSample_gp2_uncuedfailure]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp2_unit_ids,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,withReplacement,justBootstrapTrials);

% nan all trials that lack any data
if nanallzeros==true
    tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(all(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success==0,2),:)=nan;
    tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(all(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure==0,2),:)=nan;
    tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(all(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success==0,2),:)=nan;
    tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(all(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure==0,2),:)=nan;
end

temp1=nan(nBoot,1); temp2=nan(nBoot,1);
for i=1:nBoot
    temp1(i)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp1_cuedsuccess(i,:)) ...
          & ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,unitsPerSample_gp1_cuedsuccess(i,:)))); 
    temp2(i)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_success(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp2_cuedsuccess(i,:)) & ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success,unitsPerSample_gp2_cuedsuccess(i,:)))); 
end
figure();
if addThirdAxis==false
    scatter(temp1,temp2,[],'g'); hold on; scatter(nanmean(temp1),nanmean(temp2),[],'g','filled'); cuedsuccmeanx=nanmean(temp1); cuedsuccmeany=nanmean(temp2);
else
    scatter3(temp1,temp2,(temp1+temp2),[],'g'); hold on; scatter3(nanmean(temp1),nanmean(temp2),nanmean(temp1+temp2),[],'g','filled');
end
disp('cued success is green');

temp1=nan(nBoot,1); temp2=nan(nBoot,1);
for i=1:nBoot
    temp1(i)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp1_cuedfailure(i,:)) & ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,unitsPerSample_gp1_cuedfailure(i,:)))); 
    temp2(i)=nanmean(tbyt_cuedsuccess_vs_cuedfailure.unitfr_failure(ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp2_cuedfailure(i,:)) & ismember(tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_failure,unitsPerSample_gp2_cuedfailure(i,:)))); 
end
if addThirdAxis==false
    scatter(temp1,temp2,[],'r'); scatter(nanmean(temp1),nanmean(temp2),[],'r','filled'); cuedfailmeanx=nanmean(temp1); cuedfailmeany=nanmean(temp2);
else
    scatter3(temp1,temp2,(temp1+temp2),[],'r'); hold on; scatter3(nanmean(temp1),nanmean(temp2),nanmean(temp1+temp2),[],'r','filled');
end
disp('cued failure is red');

temp1=nan(nBoot,1); temp2=nan(nBoot,1);
for i=1:nBoot
    temp1(i)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp1_uncuedsuccess(i,:)) & ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,unitsPerSample_gp1_uncuedsuccess(i,:)))); 
    temp2(i)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_success(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_success,trialIDsPerSample_gp2_uncuedsuccess(i,:)) & ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_success,unitsPerSample_gp2_uncuedsuccess(i,:)))); 
end
if addThirdAxis==false
    scatter(temp1,temp2,[],'b'); scatter(nanmean(temp1),nanmean(temp2),[],'b','filled'); uncuedsuccmeanx=nanmean(temp1); uncuedsuccmeany=nanmean(temp2);
else
    scatter3(temp1,temp2,(temp1+temp2),[],'b'); hold on; scatter3(nanmean(temp1),nanmean(temp2),nanmean(temp1+temp2),[],'b','filled');
end
disp('uncued success is blue');

temp1=nan(nBoot,1); temp2=nan(nBoot,1);
for i=1:nBoot
    temp1(i)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp1_uncuedfailure(i,:)) & ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,unitsPerSample_gp1_uncuedfailure(i,:)))); 
    temp2(i)=nanmean(tbyt_uncuedsuccess_vs_uncuedfailure.unitfr_failure(ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichTrialID_failure,trialIDsPerSample_gp2_uncuedfailure(i,:)) & ismember(tbyt_uncuedsuccess_vs_uncuedfailure.fromWhichUnit_failure,unitsPerSample_gp2_uncuedfailure(i,:)))); 
end
if addThirdAxis==false
    scatter(temp1,temp2,[],'y'); scatter(nanmean(temp1),nanmean(temp2),[],'y','filled'); uncuedfailmeanx=nanmean(temp1); uncuedfailmeany=nanmean(temp2);
else
    scatter3(temp1,temp2,(temp1+temp2),[],'y'); hold on; scatter3(nanmean(temp1),nanmean(temp2),nanmean(temp1+temp2),[],'y','filled');
end
disp('uncued failure is yellow');
if addThirdAxis==false
    scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
end
xlabel('Gp 1 average trial firing rate (flipped)');
ylabel('Gp 2 average trial firing rate');
if addThirdAxis==true
    zlabel('Sum of gp1 and gp2 firing rate');
end

end

function [trialIDsPerSample_gp1_cuedsuccess,unitsPerSample_gp1_cuedsuccess]=sampleUnitsAndTrials(nBoot,nUnits,nTrials,gp1_unit_ids,currtbyt_trialIDs,currtbyt_units,withReplacement,justBootstrapTrials)

% currtbyt_trialIDs=tbyt_cuedsuccess_vs_cuedfailure.fromWhichTrialID_success
% currtbyt_units=tbyt_cuedsuccess_vs_cuedfailure.fromWhichUnit_success

% Random samples of nUnits of idx==1 from nUnits random single trials
unitsPerSample_gp1_cuedsuccess=nan(nBoot,nUnits);
trialIDsPerSample_gp1_cuedsuccess=nan(nBoot,nTrials);
for i=1:nBoot
    if justBootstrapTrials==false
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
    else
        % testing
        unitsPerSample_gp1_cuedsuccess=repmat(currtbyt_units(ismember(currtbyt_units,gp1_unit_ids))',nBoot,1);
        trialIDsPerSample_gp1_cuedsuccess=repmat(currtbyt_trialIDs(ismember(currtbyt_units,gp1_unit_ids))',nBoot,1);
        

        % all trials for units of this gp
%         currTrials=currtbyt_trialIDs(ismember(currtbyt_units,gp1_unit_ids));
%         tookwhich=randsample(length(currTrials),nTrials,withReplacement);
%         trialIDsPerSample_gp1_cuedsuccess(i,:)=currTrials(tookwhich);
%         f=nan(1,length(trialIDsPerSample_gp1_cuedsuccess(i,:)));
%         for j=1:length(trialIDsPerSample_gp1_cuedsuccess(i,:))
%             fallunits=find(currtbyt_trialIDs==trialIDsPerSample_gp1_cuedsuccess(i,j)); % trial ID repeated for different units
%             allunits=currtbyt_units(fallunits);
%             % randomly choose one of the units
%             whichunit=randsample(length(fallunits),1,withReplacement);
%             unitsPerSample_gp1_cuedsuccess(i,j)=allunits(whichunit);
%         end
    end
end

end