function makeFig5figs(tensor_tomatchcuedsuccess,cued_success_Response)

excludeneuronsforallplot=921:983; % alignment weird for cued
includeneurons=~ismember(1:size(tensor_tomatchcuedsuccess,1),excludeneuronsforallplot);

% all neurons, unsorted
figure();
imagesc(tensor_tomatchcuedsuccess(includeneurons,:,1));
colormap(flipud(colormap('gray')));
title('cued success all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(includeneurons,:,2));
colormap(flipud(colormap('gray')));
title('cued failure all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(includeneurons,:,3));
colormap(flipud(colormap('gray')));
title('uncued success all neurons unsorted');
figure();
imagesc(tensor_tomatchcuedsuccess(includeneurons,:,4));
colormap(flipud(colormap('gray')));
title('uncued failure all neurons unsorted');

% orderings
cuedsuccess_minus_uncuedsuccess=medfilt1(tensor_tomatchcuedsuccess(:,:,1)-tensor_tomatchcuedsuccess(:,:,3),100,[],2);
successdiff=nanmean(cuedsuccess_minus_uncuedsuccess(:,300:700),2);
cuedsuccess_minus_uncuedsuccess(isnan(cuedsuccess_minus_uncuedsuccess))=0; % for plotting
cuedsuccess_minus_uncuedsuccess(cuedsuccess_minus_uncuedsuccess>3)=3; 
cuedsuccess_minus_uncuedsuccess(cuedsuccess_minus_uncuedsuccess<-3)=-3;
successdiff_forplot=nanmean(cuedsuccess_minus_uncuedsuccess(:,300:700),2);
cuedfailure_minus_uncuedfailure=medfilt1(tensor_tomatchcuedsuccess(:,:,2)-tensor_tomatchcuedsuccess(:,:,4),100,[],2);
failurediff=nanmean(cuedfailure_minus_uncuedfailure(:,300:700),2);
cuedfailure_minus_uncuedfailure(isnan(cuedfailure_minus_uncuedfailure))=0; % for plotting
cuedfailure_minus_uncuedfailure(cuedfailure_minus_uncuedfailure>3)=3; 
cuedfailure_minus_uncuedfailure(cuedfailure_minus_uncuedfailure<-3)=-3;
failurediff_forplot=nanmean(cuedfailure_minus_uncuedfailure(:,300:700),2);
cuedsuccess_minus_cuedfailure=medfilt1(tensor_tomatchcuedsuccess(:,:,1)-tensor_tomatchcuedsuccess(:,:,2),100,[],2);
cueddiff=nanmean(cuedsuccess_minus_cuedfailure(:,300:700),2);
cuedsuccess_minus_cuedfailure(isnan(cuedsuccess_minus_cuedfailure))=0; % for plotting
cuedsuccess_minus_cuedfailure(cuedsuccess_minus_cuedfailure>3)=3; 
cuedsuccess_minus_cuedfailure(cuedsuccess_minus_cuedfailure<-3)=-3;
cueddiff_forplot=nanmean(cuedsuccess_minus_cuedfailure(:,300:700),2);
uncuedsuccess_minus_uncuedfailure=medfilt1(tensor_tomatchcuedsuccess(:,:,3)-tensor_tomatchcuedsuccess(:,:,4),100,[],2);
uncueddiff=nanmean(uncuedsuccess_minus_uncuedfailure(:,300:700),2);
uncuedsuccess_minus_uncuedfailure(isnan(uncuedsuccess_minus_uncuedfailure))=0; % for plotting
uncuedsuccess_minus_uncuedfailure(uncuedsuccess_minus_uncuedfailure>3)=3; 
uncuedsuccess_minus_uncuedfailure(uncuedsuccess_minus_uncuedfailure<-3)=-3;
uncueddiff_forplot=nanmean(uncuedsuccess_minus_uncuedfailure(:,300:700),2);

% sorted by cued succ minus cued fail
currentcombo=cueddiff_forplot(cued_success_Response.consensus_idx==2);
f=find(cued_success_Response.consensus_idx==2);
[~,si]=sort(currentcombo,'descend');
temp=tensor_tomatchcuedsuccess(f(si),:,1);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp2 cued success sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,2);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp2 cued failure sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,3);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp2 uncued success sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,4);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp2 uncued failure sorted by descending cueddiff');
currentcombo=cueddiff_forplot(cued_success_Response.consensus_idx==1);
f=find(cued_success_Response.consensus_idx==1);
[~,si]=sort(currentcombo,'descend');
temp=tensor_tomatchcuedsuccess(f(si),:,1);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp1 cued success sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,2);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp1 cued failure sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,3);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp1 uncued success sorted by descending cueddiff');
temp=tensor_tomatchcuedsuccess(f(si),:,4);
figure(); imagesc(temp);
colormap(flipud(colormap('gray')));
title('gp1 uncued failure sorted by descending cueddiff');

% succ minus fail
% cued
currentcombo=cueddiff_forplot(cued_success_Response.consensus_idx==2);
f=find(cued_success_Response.consensus_idx==2);
[~,si]=sort(currentcombo,'descend');
temp=cuedsuccess_minus_cuedfailure(f(si),:);
figure(); imagesc(temp); colormap('winter');
title('gp2 cued success minus failure');
currentcombo=cueddiff_forplot(cued_success_Response.consensus_idx==1);
f=find(cued_success_Response.consensus_idx==1);
[~,si]=sort(currentcombo,'descend');
temp=cuedsuccess_minus_cuedfailure(f(si),:);
figure(); imagesc(temp); colormap('winter');
title('gp1 cued success minus failure');
% uncued
currentcombo=uncueddiff_forplot(cued_success_Response.consensus_idx==2);
f=find(cued_success_Response.consensus_idx==2);
[~,si]=sort(currentcombo,'descend');
temp=uncuedsuccess_minus_uncuedfailure(f(si),:);
figure(); imagesc(temp); colormap('winter');
title('gp2 uncued success minus failure');
currentcombo=uncueddiff_forplot(cued_success_Response.consensus_idx==1);
f=find(cued_success_Response.consensus_idx==1);
[~,si]=sort(currentcombo,'descend');
temp=uncuedsuccess_minus_uncuedfailure(f(si),:);
figure(); imagesc(temp); colormap('winter');
title('gp1 uncued success minus failure');

% cued minus uncued
% success
currentcombo=successdiff_forplot(cued_success_Response.consensus_idx==2);
f=find(cued_success_Response.consensus_idx==2);
[~,si]=sort(currentcombo,'descend');
temp=cuedsuccess_minus_uncuedsuccess(f(si),:);
figure(); imagesc(temp); colormap('spring');
title('gp2 success cued minus uncued');
currentcombo=successdiff_forplot(cued_success_Response.consensus_idx==1);
f=find(cued_success_Response.consensus_idx==1);
[~,si]=sort(currentcombo,'descend');
temp=cuedsuccess_minus_uncuedsuccess(f(si),:);
figure(); imagesc(temp); colormap('spring');
title('gp1 success cued minus uncued');
% failure
currentcombo=failurediff_forplot(cued_success_Response.consensus_idx==2);
f=find(cued_success_Response.consensus_idx==2);
[~,si]=sort(currentcombo,'descend');
temp=cuedfailure_minus_uncuedfailure(f(si),:);
figure(); imagesc(temp); colormap('spring');
title('gp2 failure cued minus uncued');
currentcombo=failurediff_forplot(cued_success_Response.consensus_idx==1);
f=find(cued_success_Response.consensus_idx==1);
[~,si]=sort(currentcombo,'descend');
temp=cuedfailure_minus_uncuedfailure(f(si),:);
figure(); imagesc(temp); colormap('spring');
title('gp1 failure cued minus uncued');

return
% cue and outcome interaction
% cueAndOutcomeInteraction(tensor_tomatchcuedsuccess,cued_success_Response.consensus_idx,2,-5-(0.15/2):0.15:5,2);
cueAndOutcomeInteraction(tensor_tomatchcuedsuccess,cued_success_Response.consensus_idx,3,-5-0.025:0.05:5,2);
cueAndOutcomeInteraction(tensor_tomatchcuedsuccess,cued_success_Response.consensus_idx,2,-5-0.15:0.3:5,2);
cueAndOutcomeInteraction(tensor_tomatchcuedsuccess,ones(size(cued_success_Response.consensus_idx)),2,-5-0.05:0.1:5,1);

end