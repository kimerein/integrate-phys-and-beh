function photometry_tbt=recalcPhotometry(photometry_tbt)

useGreenCh=true;
useRedCh=true;
nanOutRedCh=true;

% baseline calculation or just Z-score
settings.Zscore_or_dF_F='Zscore'; % either Zscore or dF_F
settings.whichBaseline='percentile'; % can be 'percentile' or 'median'
settings.prc=10; % if using 'percentile', which prctile
% currently uses full trial as source of F
settings.chronux.movingwin=[0.4 0.01]; % what WAS used for moving window 

% process photometry
if useGreenCh==true
        currdata=photometry_tbt.raw_green_ch;
        for i=1:size(currdata,1)
            tempdata=currdata(i,:);
            tempdata(tempdata==0)=nan;
            switch settings.Zscore_or_dF_F
                case 'dF_F'
                    switch settings.whichBaseline
                        case 'percentile'
                            temp=prctile(tempdata(~isnan(tempdata) & tempdata~=0),settings.prc);
                        case 'median'
                            temp=median(tempdata(~isnan(tempdata) & tempdata~=0),1,'omitnan');
                    end
                    currdata(i,:)=tempdata/temp;
                case 'Zscore'
                    ftempdata=~isnan(tempdata);
                    tempdata(ftempdata)=zscore(tempdata(ftempdata));
                    currdata(i,:)=tempdata;
            end
        end
        photometry_tbt.recalc_green_ch=currdata;
end
if useRedCh==true
        currdata=photometry_tbt.raw_red_ch;
        for i=1:size(currdata,1)
            tempdata=currdata(i,:);
            tempdata(tempdata==0)=nan;
            switch settings.Zscore_or_dF_F
                case 'dF_F'
                    switch settings.whichBaseline
                        case 'percentile'
                            temp=prctile(tempdata(~isnan(tempdata) & tempdata~=0),settings.prc);
                        case 'median'
                            temp=median(tempdata(~isnan(tempdata) & tempdata~=0),1,'omitnan');
                    end
                    currdata(i,:)=tempdata/temp;
                case 'Zscore'
                    ftempdata=~isnan(tempdata);
                    tempdata(ftempdata)=zscore(tempdata(ftempdata));
                    currdata(i,:)=tempdata;
            end
        end
        photometry_tbt.recalc_red_ch=currdata;
end
if nanOutRedCh==true
        currdata=photometry_tbt.raw_red_ch;
        for i=1:size(currdata,1)
            tempdata=currdata(i,:);
            % nan out during cue
            % find cue on this trial
            cuet=photometry_tbt.cue_times(i,photometry_tbt.cue(i,:)>0.1);
            for j=1:length(cuet)
                tempdata(photometry_tbt.red_time(i,:)>cuet(j)-settings.chronux.movingwin(1) & photometry_tbt.red_time(i,:)<cuet(j)+settings.chronux.movingwin(1))=nan;
            end
            tempdata(tempdata==0)=nan;
            switch settings.Zscore_or_dF_F
                case 'dF_F'
                    switch settings.whichBaseline
                        case 'percentile'
                            temp=prctile(tempdata(~isnan(tempdata) & tempdata~=0),settings.prc);
                        case 'median'
                            temp=median(tempdata(~isnan(tempdata) & tempdata~=0),1,'omitnan');
                    end
                    currdata(i,:)=tempdata/temp;
                case 'Zscore'
                    ftempdata=~isnan(tempdata);
                    tempdata(ftempdata)=zscore(tempdata(ftempdata));
                    currdata(i,:)=tempdata;
            end
        end
        photometry_tbt.nan_out_red_ch=currdata;  
end

end