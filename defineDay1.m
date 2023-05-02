function [day1,metadata]=defineDay1(alltbt,trialTypes,metadata,isreachout_permouse,permouse_mouseid)

pelletPercThresh=61; % less than this percent pellets loaded per all wheel turns
successThresh=20; % at least this many pellets successfully
touchedPelletThresh=19; % at least this many touches of pellet
useTouchedPellet=true;
expectedNTrialsPerSess=260; % expect AT LEAST this many trials per session

day1=nan(1,length(permouse_mouseid));
metadata.sess_wrt_day1=metadata.nth_session;
for i=1:length(permouse_mouseid)
    % find day 1 for this mouse
    % using following criteria:
    %   1. Mouse successfully grabbed and consumed >30 pellets over the course of a >=45 minute session
    %   2. Pellet was present after the cue less than 50% of the time
    temp=isreachout_permouse{i};
    [~,si]=sort(temp.nth_session);
    nthsessions=temp.nth_session(si);
    hasSuccess=temp.hasSuccess(si);
    hasTouch=temp.touched_pellet(si);
    tottrials=temp.totalTrials(si);
    pelletPresentAtCue=temp.pelletPresent(si)./temp.totalTrials(si);
    pelletPerc=nan(1,length(nthsessions));
    for j=1:length(nthsessions)
        % If I only got half the trials for this session (1 vid), need to
        % multiply hasSuccess by 2
        ntrials=nansum(metadata.mouseid==permouse_mouseid(i) & metadata.nth_session==nthsessions(j));
        if tottrials(j)<expectedNTrialsPerSess/2
            % only half a session here
            hasSuccess(j)=hasSuccess(j)*2;
            hasTouch(j)=hasTouch(j)*2;
        end
        % Get pellet %
        pelletPerc(j)=mode(metadata.pelletPresentFromTable(metadata.mouseid==permouse_mouseid(i) & metadata.nth_session==nthsessions(j)));
        if isnan(pelletPerc(j))
            % some other way to estimate in what fraction of wheel turns
            % pellet was presented
            % use trial ITI to figure this out
            cueStarts=alltbt.timesFromSessionStart(:,94);
            ITIs=[diff(cueStarts); nan];
            % greater than 18 sec is two wheel turns
            thissessITIs=ITIs(metadata.mouseid==permouse_mouseid(i) & metadata.nth_session==nthsessions(j));
            pelletPerc(j)=((nansum(thissessITIs<18)./length(thissessITIs))-0.05)*100; % -5 to account for some slop in this calculation
        end
    end
    if useTouchedPellet==true
        f=find((hasTouch>=touchedPelletThresh & pelletPerc<pelletPercThresh) | (hasSuccess>=successThresh & pelletPerc<pelletPercThresh),1,'first');
    else
        f=find(hasSuccess>=successThresh & pelletPerc<pelletPercThresh,1,'first');
    end
    day1(i)=nthsessions(f);
    % put into metadata
    metadata.sess_wrt_day1(metadata.mouseid==permouse_mouseid(i))=metadata.sess_wrt_day1(metadata.mouseid==permouse_mouseid(i))-day1(i)+1;
end

end