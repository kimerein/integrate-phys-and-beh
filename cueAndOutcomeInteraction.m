function cueAndOutcomeInteraction(tensor,consensus_idx,ds,bins,consensus_currgroup)

% cueAndOutcomeInteraction(tensor_tomatchcuedsuccess,cued_success_Response.consensus_idx,2,-5-(0.15/2):0.15:5,2);
% ds is for 2D smoothing
% bins=-5-0.05:0.1:5;
coloroffset=0.3;

cuedsuccess_minus_uncuedsuccess=medfilt1(tensor(:,:,1)-tensor(:,:,3),50,[],2);
successdiff=nanmean(cuedsuccess_minus_uncuedsuccess(:,300:700),2);

cuedfailure_minus_uncuedfailure=medfilt1(tensor(:,:,2)-tensor(:,:,4),50,[],2);
failurediff=nanmean(cuedfailure_minus_uncuedfailure(:,300:700),2);

cuedsuccess_minus_cuedfailure=medfilt1(tensor(:,:,1)-tensor(:,:,2),50,[],2);
cueddiff=nanmean(cuedsuccess_minus_cuedfailure(:,300:700),2);

uncuedsuccess_minus_uncuedfailure=medfilt1(tensor(:,:,3)-tensor(:,:,4),50,[],2);
uncueddiff=nanmean(uncuedsuccess_minus_uncuedfailure(:,300:700),2);

% cued and succ
[n,x]=histcounts(cueddiff(consensus_idx==consensus_currgroup),bins);
[n_xaxis,x_xaxis]=cityscape_hist(n,x);
[n,x]=histcounts(successdiff(consensus_idx==consensus_currgroup),bins);
[n_yaxis,x_yaxis]=cityscape_hist(n,x);
% real joint distribution
temp1=cueddiff(consensus_idx==consensus_currgroup); temp2=successdiff(consensus_idx==consensus_currgroup); eithernan=isnan(temp1) | isnan(temp2);
[R,p]=corrcoef(temp1(~eithernan),temp2(~eithernan));
[joint,zeroat,xfit,yfit]=jointDistribution(cueddiff(consensus_idx==consensus_currgroup),successdiff(consensus_idx==consensus_currgroup),bins,bins);
temp=successdiff(consensus_idx==consensus_currgroup);
[fakejoint,fakezeroat]=jointDistribution(repmat(cueddiff(consensus_idx==consensus_currgroup),5,1),[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))],bins,bins);
temp1=repmat(cueddiff(consensus_idx==consensus_currgroup),5,1); temp2=[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))]; eithernan=isnan(temp1) | isnan(temp2);
[Rfake,pfake]=corrcoef(temp1(~eithernan),temp2(~eithernan));
figure(); 
subplot(4,1,1); plot(x_xaxis,n_xaxis); xlabel('cued succ vs fail');
subplot(4,1,2); plot(x_yaxis,n_yaxis); xlabel('succ cued vs uncued');
% subplot(4,1,3); imagesc(conv2(joint,ones(ds,ds),'same')); xlabel('true joint');
% subplot(4,1,4); imagesc(conv2(fakejoint,ones(ds,ds),'same')); xlabel('fake joint');
temp=medfilt2(joint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(zeroat(1),zeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,3); imagesc(bins,bins,temp); xlabel('true joint'); ylabel(['R ' num2str(R(1,2)) ' p ' num2str(p(1,2))]); colormap(flipud(colormap('gray')));
hold on; plot(yfit,xfit,'Color','r');
temp=medfilt2(fakejoint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(fakezeroat(1),fakezeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,4); imagesc(bins,bins,temp); xlabel('fake joint'); ylabel(['Rfake ' num2str(Rfake(1,2)) ' pfake ' num2str(pfake(1,2))]); colormap(flipud(colormap('gray')));

% cued and fail
[n,x]=histcounts(cueddiff(consensus_idx==consensus_currgroup),bins);
[n_xaxis,x_xaxis]=cityscape_hist(n,x);
[n,x]=histcounts(failurediff(consensus_idx==consensus_currgroup),bins);
[n_yaxis,x_yaxis]=cityscape_hist(n,x);
% real joint distribution
temp1=cueddiff(consensus_idx==consensus_currgroup); temp2=failurediff(consensus_idx==consensus_currgroup); eithernan=isnan(temp1) | isnan(temp2);
[R,p]=corrcoef(temp1(~eithernan),temp2(~eithernan));
[joint,zeroat,xfit,yfit]=jointDistribution(cueddiff(consensus_idx==consensus_currgroup),failurediff(consensus_idx==consensus_currgroup),bins,bins);
temp=failurediff(consensus_idx==consensus_currgroup);
[fakejoint,fakezeroat]=jointDistribution(repmat(cueddiff(consensus_idx==consensus_currgroup),5,1),[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))],bins,bins);
temp1=repmat(cueddiff(consensus_idx==consensus_currgroup),5,1); temp2=[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))]; eithernan=isnan(temp1) | isnan(temp2);
[Rfake,pfake]=corrcoef(temp1(~eithernan),temp2(~eithernan));
figure(); 
subplot(4,1,1); plot(x_xaxis,n_xaxis); xlabel('cued succ vs fail');
subplot(4,1,2); plot(x_yaxis,n_yaxis); xlabel('fail cued vs uncued');
% subplot(4,1,3); imagesc(conv2(joint,ones(ds,ds),'same')); xlabel('true joint');
% subplot(4,1,4); imagesc(conv2(fakejoint,ones(ds,ds),'same')); xlabel('fake joint');
temp=medfilt2(joint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(zeroat(1),zeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,3); imagesc(bins,bins,temp); xlabel('true joint'); ylabel(['R ' num2str(R(1,2)) ' p ' num2str(p(1,2))]); colormap(flipud(colormap('gray')));
hold on; plot(yfit,xfit,'Color','r');
temp=medfilt2(fakejoint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(fakezeroat(1),fakezeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,4); imagesc(bins,bins,temp); xlabel('fake joint'); ylabel(['Rfake ' num2str(Rfake(1,2)) ' pfake ' num2str(pfake(1,2))]); colormap(flipud(colormap('gray')));

% uncued and fail
[n,x]=histcounts(uncueddiff(consensus_idx==consensus_currgroup),bins);
[n_xaxis,x_xaxis]=cityscape_hist(n,x);
[n,x]=histcounts(failurediff(consensus_idx==consensus_currgroup),bins);
[n_yaxis,x_yaxis]=cityscape_hist(n,x);
% real joint distribution
temp1=uncueddiff(consensus_idx==consensus_currgroup); temp2=failurediff(consensus_idx==consensus_currgroup); eithernan=isnan(temp1) | isnan(temp2);
[R,p]=corrcoef(temp1(~eithernan),temp2(~eithernan));
[joint,zeroat,xfit,yfit]=jointDistribution(uncueddiff(consensus_idx==consensus_currgroup),failurediff(consensus_idx==consensus_currgroup),bins,bins);
temp=failurediff(consensus_idx==consensus_currgroup);
[fakejoint,fakezeroat]=jointDistribution(repmat(uncueddiff(consensus_idx==consensus_currgroup),5,1),[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))],bins,bins);
temp1=repmat(uncueddiff(consensus_idx==consensus_currgroup),5,1); temp2=[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))]; eithernan=isnan(temp1) | isnan(temp2);
[Rfake,pfake]=corrcoef(temp1(~eithernan),temp2(~eithernan));
figure();
subplot(4,1,1); plot(x_xaxis,n_xaxis); xlabel('uncued succ vs fail');
subplot(4,1,2); plot(x_yaxis,n_yaxis); xlabel('fail cued vs uncued');
% subplot(4,1,3); imagesc(conv2(joint,ones(ds,ds),'same')); xlabel('true joint');
% subplot(4,1,4); imagesc(conv2(fakejoint,ones(ds,ds),'same')); xlabel('fake joint');
temp=medfilt2(joint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(zeroat(1),zeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,3); imagesc(bins,bins,temp); xlabel('true joint'); ylabel(['R ' num2str(R(1,2)) ' p ' num2str(p(1,2))]); colormap(flipud(colormap('gray')));
hold on; plot(yfit,xfit,'Color','r');
temp=medfilt2(fakejoint,[ds ds]); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(fakezeroat(1),fakezeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,4); imagesc(bins,bins,temp); xlabel('fake joint'); ylabel(['Rfake ' num2str(Rfake(1,2)) ' pfake ' num2str(pfake(1,2))]); colormap(flipud(colormap('gray')));

% uncued and succ
[n,x]=histcounts(uncueddiff(consensus_idx==consensus_currgroup),bins);
[n_xaxis,x_xaxis]=cityscape_hist(n,x);
[n,x]=histcounts(successdiff(consensus_idx==consensus_currgroup),bins);
[n_yaxis,x_yaxis]=cityscape_hist(n,x);
% real joint distribution
temp1=uncueddiff(consensus_idx==consensus_currgroup); temp2=successdiff(consensus_idx==consensus_currgroup); eithernan=isnan(temp1) | isnan(temp2);
[R,p]=corrcoef(temp1(~eithernan),temp2(~eithernan));
[joint,zeroat,xfit,yfit]=jointDistribution(uncueddiff(consensus_idx==consensus_currgroup),successdiff(consensus_idx==consensus_currgroup),bins,bins);
temp=successdiff(consensus_idx==consensus_currgroup);
[fakejoint,fakezeroat]=jointDistribution(repmat(uncueddiff(consensus_idx==consensus_currgroup),5,1),[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))],bins,bins);
temp1=repmat(uncueddiff(consensus_idx==consensus_currgroup),5,1); temp2=[temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)));temp(randperm(length(temp)))]; eithernan=isnan(temp1) | isnan(temp2);
[Rfake,pfake]=corrcoef(temp1(~eithernan),temp2(~eithernan));
figure(); 
subplot(4,1,1); plot(x_xaxis,n_xaxis); xlabel('uncued succ vs fail');
subplot(4,1,2); plot(x_yaxis,n_yaxis); xlabel('succ cued vs uncued');
% subplot(4,1,3); imagesc(conv2(joint,ones(ds,ds),'same')); xlabel('true joint');
% subplot(4,1,4); imagesc(conv2(fakejoint,ones(ds,ds),'same')); xlabel('fake joint');
temp=medfilt2(joint,[ds ds]); 
% temp=imgaussfilt(joint,ds); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(zeroat(1),zeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,3); imagesc(bins,bins,temp); xlabel('true joint'); ylabel(['R ' num2str(R(1,2)) ' p ' num2str(p(1,2))]); colormap(flipud(colormap('gray')));
hold on; plot(yfit,xfit,'Color','r');
temp=medfilt2(fakejoint,[ds ds]); 
% temp=imgaussfilt(fakejoint,ds); 
% temp=conv2(temp,ones(ds,ds),'same'); 
temp(fakezeroat(1),fakezeroat(2))=nan;
temp(temp==0)=coloroffset; temp=log10(temp+0.51);
subplot(4,1,4); imagesc(bins,bins,temp); xlabel('fake joint'); ylabel(['Rfake ' num2str(Rfake(1,2)) ' pfake ' num2str(pfake(1,2))]); colormap(flipud(colormap('gray')));

end

function [x,yfit]=bestfitline(x,y)

nonan=~isnan(x) & ~isnan(y) & ~(x==0 & y==0);
p=polyfit(x(nonan),y(nonan),1);
yfit=polyval(p,x);

end

function [joint,zeroat,xfit,yfit]=jointDistribution(x,y,xbins,ybins)

dropzero=true;

[xfit,yfit]=bestfitline(x,y);

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

