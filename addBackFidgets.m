function beh_tbt=addBackFidgets(beh_tbt,autoReachSettings,fidget,whichMovie)

% note that movieframeinds directly references the movie frames
% fidget inds match movie frames AFTER subtraction of
% autoReachSettings.discardFirstNFrames

if ~isfield(beh_tbt,'isFidgeting')
    beh_tbt.isFidgeting=nan(size(beh_tbt.movie_distractor));
    beh_tbt.fidgetData=nan(size(beh_tbt.movie_distractor));
end

realmovieframes=[1:length(fidget.fidgetData)]+autoReachSettings.discardFirstNFrames;

for i=1:size(beh_tbt.isFidgeting,1)
    for j=1:size(beh_tbt.isFidgeting,2)
        currmovframe=beh_tbt.movieframeinds(i,j);
        if beh_tbt.from_first_video(i,j)==1
            currmovie=1;
        elseif beh_tbt.from_second_video(i,j)==1
            currmovie=2;
        elseif beh_tbt.from_third_video(i,j)==1
            currmovie=3;
        end
        if isnan(currmovframe)
            continue
        end
        if currmovie~=whichMovie
            continue
        end
        % find fidget from same movie frame
        [~,mi]=nanmin(abs(currmovframe-realmovieframes));
        beh_tbt.isFidgeting(i,j)=fidget.isFidgeting(mi);
        beh_tbt.fidgetData(i,j)=fidget.fidgetData(mi);
    end
end

