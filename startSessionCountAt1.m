function [metadata]=startSessionCountAt1(metadata)

% for each mouse, make the first available date "session #1"
% i.e., nth_session equals 1 for this date
% shift all subsequent dates accordingly

unique_mice=unique(metadata.mouseid);

for i=1:length(unique_mice)
    currmouse=unique_mice(i);
    isthismouse=metadata.mouseid==currmouse;
    unique_sess=unique(metadata.nth_session(isthismouse));
    min_sess=nanmin(unique_sess);
    metadata.nth_session(isthismouse)=metadata.nth_session(isthismouse)-min_sess+1;
end
