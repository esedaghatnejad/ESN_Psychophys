function ESN_event_alignment_v2
%% load EPHYS EVENT DATA
[file_name,file_path] = uigetfile([pwd filesep '*_EVE1.mat'], 'Select EVENT DATA file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path filesep file_name]);
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = file_path;
fprintf(' --> Completed. \n')

[file_name, file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*_CMN1.mat'], 'Select _CMN1.mat file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_CMN = load([file_path filesep file_name]);
EPHYS.time = EPHYS.CH_CMN.ch_time;
EPHYS.CH_EVE.ch_time_EPHYS = EPHYS.time;
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
[file_name,file_path] = uigetfile([file_path filesep '*_REDUCED.mat'], 'Select _REDUCED file');
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path filesep file_name], 'TRIALS_DATA');
fprintf(' --> Completed. \n')

%% extract event data in ephys
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Event Data', ' ... ']);
% ch_data
% 0 : STR_TARGET_PURSUIT
% 1 : STR_TARGET_FIXATION
% 2 : DETECT_SACCADE_START
% 3 : DETECT_SACCADE_END
% 4 : END_TARGET_FIXATION
% 5 : 1Hz ossilation
% 6 : photodiode: STR_TARGET_FIXATION+DETECT_SACCADE_END
% eventId
% 1 : rising
% 0 : falling
EPHYS.CH_EVE.data = [EPHYS.CH_EVE.ch_time(:) EPHYS.CH_EVE.ch_data(:) EPHYS.CH_EVE.ch_info.eventId(:)];
fprintf(' --> Completed. \n')

%% Build state fixate from EPHYS events
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS state_fixation', ' ... ']);


time_state_str_fixation   = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 1) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
time_state_sac_detect_on  = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 2) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
time_state_sac_detect_off = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
time_state_end_fixation   = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
time_state_iti            = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1);

time_state_end_fixation
time_state_iti


eye_time            = DATA.TRIAL.eye_time;
length_time         = length(eye_time);
spike_time          = DATA.NEURAL.spike_time;
sac_str_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Saccade O')).values_times; % 'Saccade O'
sac_end_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Saccade E')).values_times; % 'Saccade E'
tgt_jmp_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Target Ti')).values_times; % 'Target Ti'
spike_time(  end+1) = max([eye_time(end), spike_time(end)])+1;
sac_str_time(end+1) = max([eye_time(end), sac_str_time(end)])+1;
sac_end_time(end+1) = max([eye_time(end), sac_end_time(end)])+1;
tgt_jmp_time(end+1) = max([eye_time(end), tgt_jmp_time(end)])+1;
spike_event         = false(length_time, 1);
sac_str_event       = false(length_time, 1);
sac_end_event       = false(length_time, 1);
tgt_jmp_event       = false(length_time, 1);
counter_spike_event   = find(spike_time    > eye_time(1), 1, 'first');
counter_sac_str_event = find(sac_str_time  > eye_time(1), 1, 'first');
counter_sac_end_event = find(sac_end_time  > eye_time(1), 1, 'first');
counter_tgt_jmp_event = find(tgt_jmp_time  > eye_time(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = eye_time(counter_time_point);
    if time_ponit_ >= spike_time(  counter_spike_event)
        spike_event(  counter_time_point) = true;
        counter_spike_event   = counter_spike_event + 1;
    end
    if time_ponit_ >= sac_str_time(counter_sac_str_event)
        sac_str_event(counter_time_point) = true;
        counter_sac_str_event = counter_sac_str_event + 1;
    end
    if time_ponit_ >= sac_end_time(counter_sac_end_event)
        sac_end_event(counter_time_point) = true;
        counter_sac_end_event = counter_sac_end_event + 1;
    end
    if time_ponit_ >= tgt_jmp_time(counter_tgt_jmp_event)
        tgt_jmp_event(counter_time_point) = true;
        counter_tgt_jmp_event = counter_tgt_jmp_event + 1;
    end
end
DATA.TRIAL.spike_event   = spike_event;
DATA.TRIAL.sac_str_event = sac_str_event;
DATA.TRIAL.sac_end_event = sac_end_event;
DATA.TRIAL.tgt_jmp_event = tgt_jmp_event;








EPHYS.state_fixation = zeros(length(EPHYS.time), 1);
state_start_time_end_target_fixation = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
state_start_time_iti                 = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1);
% make sure that the two vectors are in the right order
% the first element of fixation should be less than the first element of
% iti, this happens if recording starts at the middle of fixation state
if state_start_time_iti(1) < state_start_time_end_target_fixation(1)
    state_start_time_iti(1) = [];
end
% make sure the vectors are the same size
% this happens if recording stops at the middle of fixation state
if length(state_start_time_iti) ~= length(state_start_time_end_target_fixation)
    min_length = min([ length(state_start_time_iti),  length(state_start_time_end_target_fixation)]);
    state_start_time_iti = state_start_time_iti(1:min_length);
    state_start_time_end_target_fixation = state_start_time_end_target_fixation(1:min_length);
end

num_trials = length(state_start_time_end_target_fixation);
for counter_trial = 1 : 1 : num_trials
    inds_ = (EPHYS.time > state_start_time_end_target_fixation(counter_trial, 1)) & (EPHYS.time < state_start_time_iti(counter_trial, 1));
    EPHYS.state_fixation(inds_, 1) = 1;
end
EPHYS.trial_nums_aligned = (1 : 1 : num_trials)';
fprintf(' --> Completed. \n');

%% Build state fixate from BEHAVE
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE state_fixation', ' ... ']);
num_trials = length(BEHAVE.TRIALS_DATA.time_start);
BEHAVE.time = BEHAVE.TRIALS_DATA.time_start(1):0.001:BEHAVE.TRIALS_DATA.time_end(end);
BEHAVE.state_fixation = zeros(length(BEHAVE.time), 1);
state_start_time_end_target_fixation = nan(num_trials, 1);
state_start_time_iti                 = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    state_start_time_end_target_fixation(counter_trial, 1) = BEHAVE.TRIALS_DATA.time_state_end_fixation(counter_trial);
    state_start_time_iti(                counter_trial, 1) = BEHAVE.TRIALS_DATA.time_state_iti(counter_trial);
    inds_ = (BEHAVE.time > state_start_time_end_target_fixation(counter_trial, 1)) & (BEHAVE.time < state_start_time_iti(counter_trial, 1));
    BEHAVE.state_fixation(inds_, 1) = 1;
end
fprintf(' --> Completed. \n');

%% ALIGN resample EPHYS.state_fixation and xcorr it with BEHAVE.state_fixation
clearvars -except EPHYS BEHAVE
fprintf(['Aligning state_fixation', ' ... ']);
time_EPHYS            = EPHYS.time;
state_fixation_EPHYS  = EPHYS.state_fixation;
time_BEHAVE           =  BEHAVE.time;
state_fixation_BEHAVE = BEHAVE.state_fixation;
sample_freq_EPHYS = 30e3;
sample_freq_BEHAVE = 1e3;
[P1,Q1] = rat(sample_freq_BEHAVE/sample_freq_EPHYS);          % Rational fraction approximation
state_fixation_EPHYS_resample = resample(state_fixation_EPHYS,P1,Q1);        % Change sampling rate by rational factor
time_EPHYS_resample = resample(time_EPHYS,P1,Q1);        % Change sampling rate by rational factor
state_fixation_EPHYS_resample = state_fixation_EPHYS_resample(Q1:end-Q1); % remove the first and last Q1 points, the filter is not stable for those points
time_EPHYS_resample           = time_EPHYS_resample(          Q1:end-Q1); % remove the first and last Q1 points, the filter is not stable for those points
state_fixation_EPHYS_resample = round(state_fixation_EPHYS_resample); % rount the state_fixation to make them 0 and 1

[xcorr_value,xcorr_lag] = xcorr(state_fixation_EPHYS_resample,state_fixation_BEHAVE); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);
if sample_diff > 0
    state_fixation_EPHYS_resample = state_fixation_EPHYS_resample(abs(sample_diff):end);
    time_EPHYS_resample           = time_EPHYS_resample(          abs(sample_diff):end);
elseif sample_diff < 0
    state_fixation_BEHAVE         = state_fixation_BEHAVE(        abs(sample_diff):end);
    time_BEHAVE                   = time_BEHAVE(                  abs(sample_diff):end);
end
% make the vectors the same size
if length(state_fixation_BEHAVE) ~= length(state_fixation_EPHYS_resample)
    min_length = min([ length(state_fixation_BEHAVE),  length(state_fixation_EPHYS_resample)]);
    state_fixation_EPHYS_aligned  = state_fixation_EPHYS_resample(1:min_length);
    time_EPHY_aligned             = time_EPHYS_resample(          1:min_length);
    state_fixation_BEHAVE_aligned = state_fixation_BEHAVE(        1:min_length);
    time_BEHAVE_aligned           = time_BEHAVE(                  1:min_length);
end
EPHYS.state_fixation_aligned  = state_fixation_EPHYS_aligned;
EPHYS.time_aligned            = time_EPHY_aligned(:);
BEHAVE.state_fixation_aligned = state_fixation_BEHAVE_aligned;
BEHAVE.time_aligned           = time_BEHAVE_aligned(:);

num_trial_BEHAVE = length(BEHAVE.TRIALS_DATA.time_start);
ind_first_trial_fixation = find(diff(BEHAVE.state_fixation_aligned)==1, 1, 'first');
time_first_trial_fixation = BEHAVE.time_aligned(ind_first_trial_fixation);
trial_num_first_aligned = find(BEHAVE.TRIALS_DATA.time_start < time_first_trial_fixation, 1, 'last');
ind_last_trial_fixation = find(diff(BEHAVE.state_fixation_aligned)==1, 1, 'last');
time_last_trial_fixation = BEHAVE.time_aligned(ind_last_trial_fixation);
trial_num_last_aligned  = find(BEHAVE.TRIALS_DATA.time_start < time_last_trial_fixation, 1, 'last');
BEHAVE.trial_nums_aligned = (trial_num_first_aligned : 1 : trial_num_last_aligned)';
fprintf(' --> Completed. \n');

%% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE
ch_data = EPHYS.CH_EVE.ch_data;
ch_info = EPHYS.CH_EVE.ch_info;
ch_time = EPHYS.CH_EVE.ch_time;
data    = EPHYS.CH_EVE.data;
EPHYS_time                   = EPHYS.time;
EPHYS_state_fixation         = EPHYS.state_fixation;
EPHYS_state_fixation_aligned = EPHYS.state_fixation_aligned;
EPHYS_time_aligned           = EPHYS.time_aligned;
EPHYS_ind_fixation_begin     = find(diff(EPHYS.state_fixation_aligned)>+0.1);
EPHYS_ind_fixation_end       = find(diff(EPHYS.state_fixation_aligned)<-0.1);
BEHAVE_time                   = BEHAVE.time;
BEHAVE_state_fixation         = BEHAVE.state_fixation;
BEHAVE_state_fixation_aligned = BEHAVE.state_fixation_aligned;
BEHAVE_time_aligned           = BEHAVE.time_aligned;
BEHAVE_trial_nums_aligned     = BEHAVE.trial_nums_aligned;
BEHAVE_ind_fixation_begin     = find(diff(BEHAVE.state_fixation_aligned)>+0.1);
BEHAVE_ind_fixation_end       = find(diff(BEHAVE.state_fixation_aligned)<-0.1);
file_name = EPHYS.file_name_CH_EVE;
file_path = EPHYS.file_path_CH_EVE;
[~, file_name, ~] = fileparts(file_name);
file_name = [file_name '_aligned.mat'];
clearvars EPHYS BEHAVE
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([file_path filesep file_name], '-v7.3');
fprintf(' --> Completed. \n')
