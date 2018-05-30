function dts=fixChineseDVRDatetime(dts)

% Assumes format of input string is 'yyyy-MM-dd HH-mm-ss'
% Note that might be messed up by a leap year
% Each DVR is off by some fixed amount -- specify this here
datetime_adjust_years=6;
datetime_adjust_days=94;
datetime_adjust_hours=3+12;

if ~iscell(dts)
    temp=cell(1,1);
    temp{1}=dts;
    dts=temp;
end

% Convert to datetimes
for i=1:length(dts)
    dts{i}=datetime(dts{i},'InputFormat','yyyy-MM-dd HH-mm-ss');
end

% Fix date and time
for i=1:length(dts)
    dts{i}=dts{i}+days(datetime_adjust_days)+years(datetime_adjust_years)+hours(datetime_adjust_hours);
end
