function reach_batch=reachBatchSettings()

% "Reach batch" definition
reach_batch.window=0.4; % in seconds, reaches must occur within this many seconds of each other to be in same batch
reach_batch.firstreach_type={'miss_reachStarts'}; % first reach must be one of these types
reach_batch.secondreach_type={'success_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel','miss_reachStarts_pawOnWheel'}; % second reach must be one of these types
reach_batch.take_first_or_second_type=2; % number indicates whether to convert batch to the type of the first or second reach
reach_batch.take_first_or_second_timing=1; % number indicates whether to set reach timing at the timing of the first or second reach in batch
reach_batch.removePawOnWheel=1; % if removePawOnWheel=1, will change reach from reach_pawOnWheel to just reach type
