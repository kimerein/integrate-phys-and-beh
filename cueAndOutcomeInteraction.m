function cueAndOutcomeInteraction(tensor,consensus_idx,ds,bins)

% ds is for 2D smoothing
% bins=-5-0.05:0.1:5;

cuedsuccess_minus_uncuedsuccess=medfilt1(tensor(:,:,1)-tensor(:,:,3),300,[],2);
successdiff=nanmean(cuedsuccess_minus_uncuedsuccess(:,300:700),2);

cuedfailure_minus_uncuedfailure=medfilt1(tensor(:,:,2)-tensor(:,:,4),300,[],2);
failurediff=nanmean(cuedfailure_minus_uncuedfailure(:,300:700),2);

cuedsuccess_minus_cuedfailure=medfilt1(tensor(:,:,1)-tensor(:,:,2),100,[],2);
cueddiff=nanmean(cuedsuccess_minus_cuedfailure(:,300:700),2);

uncuedsuccess_minus_uncuedfailure=medfilt1(tensor(:,:,3)-tensor(:,:,4),100,[],2);
uncueddiff=nanmean(uncuedsuccess_minus_uncuedfailure(:,300:700),2);

% gp 2
[n,x]=histcounts(cueddiff(consensus_idx==2),-5-0.05:0.1:5);
[n_xaxis,x_xaxis]=cityscape_hist(n,x);
[n,x]=histcounts(successdiff(consensus_idx==2),-5-0.05:0.1:5);
[n_yaxis,x_yaxis]=cityscape_hist(n,x);
% real joint distribution
[joint,zeroat]=jointDistribution(cueddiff(consensus_idx==2),successdiff(consensus_idx==2),bins,bins);
temp=successdiff(consensus_idx==2);
[fakejoint,fakezeroat]=jointDistribution(repmat(cueddiff(consensus_idx==2),5,1),[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))],bins,bins);
figure();
subplot(4,1,1); plot(x_xaxis,n_xaxis); xlabel('cued succ vs fail');
subplot(4,1,2); plot(x_yaxis,n_yaxis); xlabel('succ cued vs uncued');
% subplot(4,1,3); imagesc(conv2(joint,ones(ds,ds),'same')); xlabel('true joint');
% subplot(4,1,4); imagesc(conv2(fakejoint,ones(ds,ds),'same')); xlabel('fake joint');
temp=medfilt2(joint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(zeroat(1),zeroat(2))=nan;
subplot(4,1,3); imagesc(temp); xlabel('true joint');
temp=medfilt2(fakejoint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(fakezeroat(1),fakezeroat(2))=nan;
subplot(4,1,4); imagesc(temp); xlabel('fake joint');

end

function [joint,zeroat]=jointDistribution(x,y,xbins,ybins)

dropzero=true;

xbins(end)=xbins(end)+0.0001;
ybins(end)=ybins(end)+0.0001;
joint=nan(length(xbins)-1,length(ybins)-1);
zeroat=nan(1,2);
for i=1:length(xbins)-1
    for j=1:length(ybins)-1
        currxbin=[xbins(i) xbins(i+1)];
        currybin=[ybins(j) ybins(j+1)];
        if 0>=currxbin(1) & 0<currxbin(2) & 0>=currybin(1) & 0<currybin(2) & dropzero==true
            zeroat(1)=i; zeroat(2)=j;
            continue
        end
        joint(i,j)=nansum(x>=currxbin(1) & x<currxbin(2) & y>=currybin(1) & y<currybin(2));
    end
end

if dropzero==true
    % fill with neighbors
    joint(zeroat(1),zeroat(2))=nanmean([joint(zeroat(1)-1,zeroat(2)-1) joint(zeroat(1)+1,zeroat(2)+1) joint(zeroat(1)+1,zeroat(2)-1) joint(zeroat(1)-1,zeroat(2)+1) joint(zeroat(1)-1,zeroat(2)) joint(zeroat(1),zeroat(2)-1) joint(zeroat(1)+1,zeroat(2)) joint(zeroat(1),zeroat(2)+1)]);
end

end

