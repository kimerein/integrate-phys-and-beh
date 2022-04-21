function data=processPhotometry(datadir,whosRig)

% Kim's setup
% first analog channel:     photodiode green channel
% second analog channel:    photodiode red channel
% third analog channel:     copy of DAC0 out to 488 LED, DAC0 = 470nm (167Hz)
% fourth analog channel:    copy of DAC1 out to 565 LED, DAC1 = 565nm (223Hz)
% fifth analog channel:     opto
% sixth analog channel:     cue
% seventh analog channel:   distractor

% set indices to different analog inputs
switch whosRig
    case 'Kim'
        inds(1).name='green_ch'; % name of channel
        inds(1).num=1; % offset into acquired analog channel, i.e., acquisition order
        useGreenCh=true;
        inds(2).name='red_ch';
        inds(2).num=2;
        useRedCh=true;
        inds(3).name='green_mod';
        inds(3).num=3;
        inds(4).name='red_mod';
        inds(4).num=4;
        inds(5).name='opto';
        inds(5).num=5;
        inds(6).name='cue';
        inds(6).num=6;
        inds(7).name='distractor';
        inds(7).num=7;
        totalAcqChannels=7; % total number of acquired analog channels
        Fs=2000; % in Hz, sampling rate of photometry
        data.Fs=Fs;
    case 'Marci'
        inds(1).name='green_ch'; % name of channel
        inds(1).num=1; % offset into acquired analog channel, i.e., acquisition order
        useGreenCh=true;
        inds(2).name='red_ch';
        inds(2).num=2;
        useRedCh=true;
        inds(3).name='green_mod';
        inds(3).num=3;
        inds(4).name='red_mod';
        inds(4).num=4;
        inds(5).name='cue';
        inds(5).num=5;
        inds(6).name='distractor';
        inds(6).num=6;
        totalAcqChannels=6; % total number of acquired analog channels
        Fs=2000; % in Hz, sampling rate of photometry
        data.Fs=Fs;
end


% Analysis approach to extract photometry transients
settings.keepPhase='yes'; % 'yes' or 'no'
settings.fourierAmp='chronux'; % 'hilbert' or 'chronux'
settings.greenfilter_range=[1 30000]; % band pass range in Hz, green channel
settings.redfilter_range=[1 30000]; % band pass range in Hz, red channel
% necessary for chronux only
% settings.chronux.movingwin=[0.05 0.01]; % dLight1.1
% settings.chronux.movingwin=[0.4 0.01]; % dLight1.1
settings.chronux.movingwin=[0.1 0.01]; % dLight1.1, used this recently
params.tapers=[3 2]; % dLight1.1
params.Fs=Fs;
params.trialave=0;
settings.chronux.params=params;
settings.chronux.display=false;
settings.chronux.greenfilter_range=[120 200]; %[120 200]; % band pass range in Hz, green channel
settings.chronux.redfilter_range=[200 265]; %[200 265]; % band pass range in Hz, red channel
% baseline calculation or just Z-score
settings.Zscore_or_dF_F='Zscore'; % either Zscore or dF_F
settings.whichBaseline='percentile'; % can be 'percentile' or 'median'
settings.prc=10; % if using 'percentile', which prctile
settings.baselineWindow=30; % in secs, length of baseline window
settings.firstSubtractDCbaseline=true; % first subtract off DC baseline across whole trace, then get dF_F



dd=dir(datadir);
%dd=dd(1:50);
k=1;
for i=3:length(dd)
    disp(['reading ' dd(i).name]);
    a=load([datadir '\' dd(i).name]);
    flipped_DATAM=a.DATAM';
    datavec=flipped_DATAM(1:end);
    
    % set size
    if i==3
        for j=1:length(inds)
            currname=inds(j).name;
            temp=getIndsNum(inds,currname);
            data.(currname)=nan(length(3:length(dd)),length(datavec(temp-1+1:totalAcqChannels:end)));
        end
    end
    
    % read in data
    for j=1:length(inds)
        currname=inds(j).name;
        temp=getIndsNum(inds,currname);
        tempdata=data.(currname);
        tempdata(k,:)=datavec(temp-1+1:totalAcqChannels:end);
        data.(currname)=tempdata;
    end
    k=k+1;
end

% linearize data -- assumes acquisition was continuous
for i=1:length(inds)
    currname=inds(i).name;
    temp=data.(currname);
    temp=temp';
    data.(currname)=temp(1:end);
end

% deal with skipped samples
if strcmp(whosRig,'Marci')
    data=dealWithSkippedPhotometrySamples(data);
end

% process photometry
if useGreenCh==true
    % process green_ch
    switch settings.keepPhase
        case 'yes'
            data.green_ch=bandPassLFP(data.green_ch.*data.green_mod,Fs,settings.greenfilter_range(1),settings.greenfilter_range(2),0);
        case 'no'
            data.green_ch=bandPassLFP(data.green_ch,Fs,settings.greenfilter_range(1),settings.greenfilter_range(2),0);
    end
    switch settings.fourierAmp
        case 'hilbert'
            data.green_ch=real(hilbert(data.green_ch));
            data.green_time=0:1/Fs:(length(data.green_ch)-1)*(1/Fs);
        case 'chronux'
            [amp,S,t,f]=getAmpWithChronux(data.green_ch,settings.chronux.movingwin,settings.chronux.params,settings.chronux.greenfilter_range);
            if settings.chronux.display==true
                figure();
                imagesc(t,f,S');
                title('Green ch');
                figure();
                plot(t,amp);
                title('Green amp');
            end
            data.green_ch=amp;
            data.green_time=t;
    end
    data.raw_green_ch=data.green_ch;
    binsBaseWindow=floor(settings.baselineWindow/(data.green_time(2)-data.green_time(1)));
    padToLength=ceil(length(data.raw_green_ch)/binsBaseWindow)*binsBaseWindow;
    if padToLength>length(data.raw_green_ch)
        temp_data=[data.raw_green_ch data.raw_green_ch(end)*ones(1,padToLength-length(data.raw_green_ch))];
    else
        temp_data=data.raw_green_ch;
    end
    switch settings.Zscore_or_dF_F
        case 'dF_F'
            if settings.firstSubtractDCbaseline==true
                totalBase=nanmin(temp_data);
                temp_data=temp_data-totalBase;
                data.raw_green_ch=data.raw_green_ch-totalBase;
            end
            switch settings.whichBaseline
                case 'percentile'
                    temp=prctile(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow),settings.prc);
                case 'median'
                    temp=median(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow),1,'omitnan');
            end
            temp=repmat(temp,binsBaseWindow,1);
            data.green_baseline=temp(1:end);
            data.green_baseline=data.green_baseline(1:length(data.raw_green_ch));
            data.green_ch=(data.raw_green_ch-data.green_baseline)./data.green_baseline;
        case 'Zscore'
            Z=zscore(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow));
            data.green_ch=Z(1:end);
            data.green_ch=data.green_ch(1:length(data.raw_green_ch));
    end
end
if useRedCh==true
    % process red_ch
    switch settings.keepPhase
        case 'yes'
            data.red_ch=bandPassLFP(data.red_ch.*data.red_mod,Fs,settings.redfilter_range(1),settings.redfilter_range(2),0);
        case 'no'
            data.red_ch=bandPassLFP(data.red_ch,Fs,settings.redfilter_range(1),settings.redfilter_range(2),0);
    end
    switch settings.fourierAmp
        case 'hilbert'
            data.red_ch=real(hilbert(data.red_ch));
            data.red_time=0:1/Fs:(length(data.red_ch)-1)*(1/Fs);
        case 'chronux'
            [amp,S,t,f]=getAmpWithChronux(data.red_ch,settings.chronux.movingwin,settings.chronux.params,settings.chronux.redfilter_range);
            if settings.chronux.display==true
                figure();
                imagesc(t,f,S');
                title('Red ch');
                figure();
                plot(t,amp);
                title('Red amp');
            end
            data.red_ch=amp;
            data.red_time=t;
    end
    data.raw_red_ch=data.red_ch;
    binsBaseWindow=floor(settings.baselineWindow/(data.red_time(2)-data.red_time(1)));
    padToLength=ceil(length(data.raw_red_ch)/binsBaseWindow)*binsBaseWindow;
    if padToLength>length(data.raw_red_ch)
        temp_data=[data.raw_red_ch data.raw_red_ch(end)*ones(1,padToLength-length(data.raw_red_ch))];
    else
        temp_data=data.raw_red_ch;
    end
    switch settings.Zscore_or_dF_F
        case 'dF_F'
            if settings.firstSubtractDCbaseline==true
                totalBase=nanmin(temp_data);
                temp_data=temp_data-totalBase;
                data.raw_red_ch=data.raw_red_ch-totalBase;
            end
            switch settings.whichBaseline
                case 'percentile'
                    temp=prctile(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow),settings.prc);
                case 'median'
                    temp=median(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow),1,'omitnan');
            end
            temp=repmat(temp,binsBaseWindow,1);
            data.red_baseline=temp(1:end);
            data.red_baseline=data.red_baseline(1:length(data.raw_red_ch));
            data.red_ch=(data.raw_red_ch-data.red_baseline)./data.red_baseline;
        case 'Zscore'
            Z=zscore(reshape(temp_data,binsBaseWindow,length(temp_data)/binsBaseWindow));
            data.red_ch=Z(1:end);
            data.red_ch=data.red_ch(1:length(data.raw_red_ch));
    end
end




end


function temp=getIndsNum(inds,ch_name)

for i=1:length(inds)
    if strcmp(ch_name,inds(i).name)
        temp=inds(i).num;
        return 
    end
end

temp=[];       

end


function [amp,S,t,f]=getAmpWithChronux(data,movingwin,params,filter_range)

[S,t,f]=mtspecgramc(data,movingwin,params);

S=sqrt(S);

amp=nanmean(S(:,f>=filter_range(1) & f<=filter_range(2)),2)';

end