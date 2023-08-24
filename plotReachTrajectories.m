function [X,Y,Z]=plotReachTrajectories(X,Y,Z,X_from_under,reachTrajTimes,smoobin,fps)

% Explanation of scaleY
% X and Z are from the side mirror
% Y is exclusively from the under mirror
% X_from_under should match well with X but is from the under mirror
% scaleY brings the scale of the under mirror distances into alignment with
% the scale of the side mirror distances
% Y is the side to side movement of paw

Xdelta_thresh=0;
Zdelta_thresh=0; %25;
Z_diff_thresh=0.5;
nindsAboveForZ=10;
realToInd=floor(2./(1/fps));
atleastthismany_notnan_points=200;
scaleY=0.25;
% startsinrange_X=[]; 
% startsinrange_Y=[]; 
% startsinrange_Z=[]; 
startsinrange_X=[0 120]; %[0 110]; %[90 110];
startsinrange_Y=[100 700]; %[0 680]; %[200 450];
startsinrange_Z=[160 1000]; %[160 300]; %1000]; %[200 350];
startTime=2; % in seconds
flipZ=false;

% Good settings for Z:\MICROSCOPE\Kim\KER Behavior\By date\High speed\20190605\March_C
% Xdelta_thresh=0;
% Zdelta_thresh=0; %25;
% Z_diff_thresh=0.5;
% nindsAboveForZ=10;
% realToInd=floor(2./(1/fps));
% atleastthismany_notnan_points=200;
% scaleY=0.25;
% startsinrange_X=[0 110]; %[90 110];
% startsinrange_Y=[0 680]; %[200 450];
% startsinrange_Z=[160 1000]; %[200 350];
% startTime=2; % in seconds

% Toss outliers
X=toss(X); X=fillNans(X);
Y=toss(Y); Y=fillNans(Y);
Z=toss(Z); Z=fillNans(Z);
X_from_under=toss(X_from_under); X_from_under=fillNans(X_from_under);

% Smooth
if ~isempty(smoobin)
    for i=1:size(X,1)
        X(i,:)=smooth(X(i,:)',smoobin);
        Y(i,:)=smooth(Y(i,:)',smoobin);
        Z(i,:)=smooth(Z(i,:)',smoobin);
        X_from_under(i,:)=smooth(X_from_under(i,:)',smoobin);
    end
end

% [X,Y,Z,X_from_under,Z_deviates]=realignToReachStarts(X,Y,Z,X_from_under,Z_diff_thresh,nindsAboveForZ,realToInd);
% USE X, i.e., FORWARD / BACKWARD, AS CRITERION INSTEAD OF UP / DOWN
backupX=X; backupZ=Z; Z=-backupX; X=backupZ;
[X,Y,Z,X_from_under,Z_deviates]=realignToReachStarts(X,Y,Z,X_from_under,Z_diff_thresh,nindsAboveForZ,realToInd);
temp=X; X=-Z; Z=temp;

X=X(~isnan(Z_deviates),:);
Y=Y(~isnan(Z_deviates),:);
Z=Z(~isnan(Z_deviates),:);
X_from_under=X_from_under(~isnan(Z_deviates),:);

figure(); plot(nanmean(X,1)-nanmean(nanmean(X,1)),'Color','k'); hold on; plot(nanmean(X_from_under,1)-nanmean(nanmean(X_from_under,1)),'Color','b');
title('X is black, X from under is blue');

% Take only reaches that start in zone defined by "startsinrange"
if ~isempty(startsinrange_X)
    startInds=ceil(startTime/(1/fps));
    startsatperch_X=any(X(:,1:startInds)>startsinrange_X(1),2) & any(X(:,1:startInds)<startsinrange_X(2),2);
    startsatperch_Y=any(Y(:,1:startInds)>startsinrange_Y(1),2) & any(Y(:,1:startInds)<startsinrange_Y(2),2);
    startsatperch_Z=any(Z(:,1:startInds)>startsinrange_Z(1),2) & any(Z(:,1:startInds)<startsinrange_Z(2),2);
%     startsatperch_Z(any(Z>400,2))=0;
    startsatperch=startsatperch_X & startsatperch_Y & startsatperch_Z;
    X=X(startsatperch,:);
    Y=Y(startsatperch,:);
    Z=Z(startsatperch,:);
    X_from_under=X_from_under(startsatperch,:);
end

% Take only the reaches with delta X greater than Xdelta_thresh
bigenoughreach=ones(size(X,1),1);
if Xdelta_thresh>0
    bigenoughreach=(max(X,[],2,'omitnan')-min(X,[],2,'omitnan'))>Xdelta_thresh;
elseif Zdelta_thresh>0
    bigenoughreach=(max(Z,[],2,'omitnan')-min(Z,[],2,'omitnan'))>Zdelta_thresh;
end
X=X(bigenoughreach==1,:);
Y=Y(bigenoughreach==1,:);
Z=Z(bigenoughreach==1,:);
X_from_under=X_from_under(bigenoughreach==1,:);

X=fillNans(X);
Y=fillNans(Y);
Z=fillNans(Z);
X_from_under=fillNans(X_from_under);

% Enough points
notnan=~isnan(X) & ~isnan(Y) & ~isnan(Z);
enoughpoints=sum(notnan,2)>atleastthismany_notnan_points;
X=X(enoughpoints,:);
Y=Y(enoughpoints,:);
Z=Z(enoughpoints,:);
X_from_under=X_from_under(enoughpoints,:);

if isempty(X)
    disp('No reaches left');
    return
end
figure(); plot(reachTrajTimes,X'); title('X');
figure(); plot(reachTrajTimes,Y'); title('Y');
figure(); plot(reachTrajTimes,Z'); title('Z');

if flipZ==true
    backupZ=Z;
    Z=-Z;
end

figure(); 
colorsUpTo=size(X,2);
cmap=colormap(cool(colorsUpTo));
for i=1:size(X,1)
%     plot3(X(i,:),Y(i,:),Z(i,:),'Color','k'); 
%     hold all; 
%     plot3(X(i,1:100),Y(i,1:100),Z(i,1:100),'Color','b'); hold on; 
% %     plot3(X(i,end-100:end),Y(i,end-100:end),Z(i,end-100:end),'Color','r'); hold on; 
    scatter3(X(i,:),Y(i,:).*scaleY,Z(i,:),30,cmap); hold on;
end
xlabel('X'); ylabel('Y'); zlabel('Z');
% view(-116.4279,56.0009);
%view(-162.2125,9.8783);
% view(-330.8644,17.7600);
view(33.6808,-22.8782);

figure();
scatter3(nanmean(X,1),nanmean(Y,1).*scaleY,nanmean(Z,1),30,cmap);
% plot3(nanmean(X,1),nanmean(Y,1),nanmean(Z,1),'Color','k');
% hold all
% plot3(nanmean(X(:,1:50),1),nanmean(Y(:,1:50),1),nanmean(Z(:,1:50),1),'Color','b');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Average');
% view(-116.4279,56.0009);
% view(-162.2125,9.8783);
% view(-330.8644,17.7600);
view(33.6808,-22.8782);

if flipZ==true
    Z=backupZ;
end

% load('Z:\MICROSCOPE\Kim\Final Figs\Fig1\reach trajectories\March_C\all_noOpto_reachTrajectories.mat')
% [Xout,Yout,Zout]=plotReachTrajectories(all_X,all_Y,all_Z,all_X_from_under,reachTrajTimes,100,256);
% Xnoopto=nanmean(Xout,1);
% Ynoopto=nanmean(Yout,1);
% Znoopto=nanmean(Zout,1);
% close all
% load('Z:\MICROSCOPE\Kim\Final Figs\Fig1\reach trajectories\March_C\all_opto_reachTrajectories.mat')
% [Xout,Yout,Zout]=plotReachTrajectories(all_X,all_Y,all_Z,all_X_from_under,reachTrajTimes,100,256);
% Xopto=nanmean(Xout,1);
% Yopto=nanmean(Yout,1);
% Zopto=nanmean(Zout,1);
% close all
% figure();
% colorsUpTo=size(Xnoopto,2);
% cmap=colormap(cool(colorsUpTo));
% scatter3(Xnoopto,Ynoopto.*0.25,Znoopto,30,cmap,'MarkerFaceColor','k');
% hold on;
% scatter3(Xopto,Yopto.*0.25,Zopto,30,cmap,'MarkerFaceColor','r');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Average');

end

function [X,Y,Z,X_from_under,reachBegins]=realignToReachStarts(X,Y,Z,X_from_under,Z_diff_thresh,nindsAboveForZ,realToInd)

% Z_diff_thresh=0.5
% nindsAboveForZ=10

figure(); 
plot(diff(Z,1,2)'); title('Before alignment');

% Find reach beginnings, using Z
reachBegins=nan(size(Z,1),1);
for i=1:size(Z,1)
    temp=diff(Z(i,:),1,2);
    % Find first when -diff Z > Z_diff_thresh and stays high for at least
    % nindsAboveforZ
    ishigh=-temp>Z_diff_thresh;
    reachbegin=nan;
    for j=1:length(ishigh)
        if j+nindsAboveForZ>length(ishigh)
            break
        end
        % check whether stays high
        if all(ishigh(j:j+nindsAboveForZ))
            reachbegin=j;
            break
        end
    end
    reachBegins(i)=reachbegin;
end

[X,Y,Z,X_from_under]=real(X,Y,Z,X_from_under,reachBegins,realToInd);

figure(); 
plot(diff(Z,1,2)'); title('After alignment');

end

function [X,Y,Z,X_from_under]=real(X,Y,Z,X_from_under,reachBegins,realToInd)

for i=1:size(X,1)
    if isnan(reachBegins(i))
        continue
    end
    X(i,:)=realignEachRow(X,reachBegins,i,realToInd);
    Y(i,:)=realignEachRow(Y,reachBegins,i,realToInd);
    Z(i,:)=realignEachRow(Z,reachBegins,i,realToInd);
    X_from_under(i,:)=realignEachRow(X_from_under,reachBegins,i,realToInd);
end

end

function temp=realignEachRow(currfield,reachBegins,i,realToInd)

temp=currfield(i,:);
if reachBegins(i)<realToInd
    % shift back in time
    temp=[nan(1,realToInd-reachBegins(i)) currfield(i,1:end-(realToInd-reachBegins(i)))];
elseif reachBegins(i)>realToInd
    % shift forward in time
    temp=[currfield(i,1+(reachBegins(i)-realToInd):end) nan(1,reachBegins(i)-realToInd)];
end

end

function matrix_filled=fillNans(matrix)

% Example input matrix with NaN values
% matrix = [1 2 NaN 4 NaN; 6 NaN 8 NaN 10; NaN 12 13 NaN NaN];

% Fill NaN values with linear interpolation or nearest neighboring elements
matrix_filled = matrix;

% Fill NaN values with linear interpolation or nearest neighboring elements
for i = 1:size(matrix, 1)
    row_data = matrix(i, :);
    nan_indices = find(isnan(row_data));
    non_nan_indices = find(~isnan(row_data));
    
    if numel(non_nan_indices) >= 2
        % There are non-NaN elements on both sides for interpolation
        valid_indices = [non_nan_indices(1) non_nan_indices(end)];
        
        for j = 1:numel(nan_indices)
            nan_index = nan_indices(j);
            
            if nan_index < valid_indices(1)
                % NaN value is below the range of interpolation, fill with nearest neighbor
                [~, nearest_idx] = min(abs(valid_indices(1) - nan_index));
                row_data(nan_index) = row_data(valid_indices(1) + nearest_idx);
            elseif nan_index > valid_indices(2)
                % NaN value is above the range of interpolation, fill with nearest neighbor
                [~, nearest_idx] = min(abs(valid_indices(2) - nan_index));
                row_data(nan_index) = row_data(valid_indices(2) - nearest_idx);
            else
                % Perform linear interpolation
                interp_values = row_data(valid_indices);
                interp_indices = valid_indices - valid_indices(1) + 1;
                row_data(nan_index) = interp1(interp_indices, interp_values, nan_index - valid_indices(1) + 1);
            end
        end
    elseif numel(non_nan_indices) == 1
        % Only one non-NaN element in the row, use nearest neighboring element
        nearest_index = non_nan_indices;
        row_data(nan_indices) = row_data(nearest_index);
    end
    
    matrix_filled(i, :) = row_data;
end

% Display the matrix with NaNs filled
% disp('Matrix with NaNs filled:');
% disp(matrix_filled);

end

function matrix_filtered=toss(matrix)

% Calculate the mean and standard deviation for each row
row_means = mean(matrix, 2, 'omitnan');
row_stds = std(matrix, 0, 2, 'omitnan');

% Replace outliers more than 3 standard deviations away from the mean with NaN
threshold = 3;
matrix_filtered = matrix;
for i = 1:size(matrix, 1)
    row_data = matrix(i, :);
    row_mean = row_means(i);
    row_std = row_stds(i);
    
    % Calculate the lower and upper threshold values
    lower_threshold = row_mean - threshold * row_std;
    upper_threshold = row_mean + threshold * row_std;
    
    % Replace outliers with NaN
    matrix_filtered(i, row_data < lower_threshold | row_data > upper_threshold) = NaN;
end

% Display the filtered matrix
% disp('Filtered Matrix:');
% disp(matrix_filtered);

end
