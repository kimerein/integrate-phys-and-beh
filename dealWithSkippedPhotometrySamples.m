function data=dealWithSkippedPhotometrySamples(data)

notSamplingThresh=-500;

f=fieldnames(data);
for i=1:length(f)
    currf=f{i};
    temp=data.(currf);
    temp(temp<notSamplingThresh)=0;
    data.(currf)=temp;
end

% do fix for Marci's photo rig
data.distractor(data.distractor<3)=4.25;
data.distractor=4.4-data.distractor;