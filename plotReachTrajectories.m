function plotReachTrajectories(X,Y,Z,X_from_under,reachTrajTimes,smoobin,fps)

Xdelta_thresh=1;
atleastthismany_notnan_points=200;
scaleY=0.25;
startsinrange_X=[90 110];
startsinrange_Y=[200 450];
startsinrange_Z=[200 350];
startTime=0.25; % in seconds

% Toss outliers
X=toss(X); X=fillNans(X);
Y=toss(Y); Y=fillNans(Y);
Z=toss(Z); Z=fillNans(Z);
X_from_under=toss(X_from_under); X_from_under=fillNans(X_from_under);

% Take only reaches that start in zone defined by "startsinrange"
startInds=ceil(startTime/(1/fps));
startsatperch_X=any(X(:,1:startInds)>startsinrange_X(1),2) & any(X(:,1:startInds)<startsinrange_X(2),2);
startsatperch_Y=any(Y(:,1:startInds)>startsinrange_Y(1),2) & any(Y(:,1:startInds)<startsinrange_Y(2),2);
startsatperch_Z=any(Z(:,1:startInds)>startsinrange_Z(1),2) & any(Z(:,1:startInds)<startsinrange_Z(2),2);
startsatperch=startsatperch_X & startsatperch_Y & startsatperch_Z;
X=X(startsatperch,:);
Y=Y(startsatperch,:);
Z=Z(startsatperch,:);
X_from_under=X_from_under(startsatperch,:);

% Take only the reaches with delta X greater than Xdelta_thresh
bigenoughreach=(max(X,[],2,'omitnan')-min(X,[],2,'omitnan'))>Xdelta_thresh;
X=X(bigenoughreach,:);
Y=Y(bigenoughreach,:);
Z=Z(bigenoughreach,:);
X_from_under=X_from_under(bigenoughreach,:);

% Enough points
notnan=~isnan(X) & ~isnan(Y) & ~isnan(Z);
enoughpoints=sum(notnan,2)>atleastthismany_notnan_points;
X=X(enoughpoints,:);
Y=Y(enoughpoints,:);
Z=Z(enoughpoints,:);
X_from_under=X_from_under(enoughpoints,:);

figure(); plot(reachTrajTimes,X'); title('X');
figure(); plot(reachTrajTimes,Y'); title('Y');
figure(); plot(reachTrajTimes,Z'); title('Z');

% Smooth
if ~isempty(smoobin)
    for i=1:size(X,1)
        X(i,:)=smooth(X(i,:)',smoobin);
        Y(i,:)=smooth(Y(i,:)',smoobin);
        Z(i,:)=smooth(Z(i,:)',smoobin);
        X_from_under(i,:)=smooth(X_from_under(i,:)',smoobin);
    end
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

figure();
scatter3(nanmean(X,1),nanmean(Y,1).*scaleY,nanmean(Z,1),30,cmap);
% plot3(nanmean(X,1),nanmean(Y,1),nanmean(Z,1),'Color','k');
% hold all
% plot3(nanmean(X(:,1:50),1),nanmean(Y(:,1:50),1),nanmean(Z(:,1:50),1),'Color','b');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Average');

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
