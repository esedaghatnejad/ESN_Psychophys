function ESN_monkey_behavior_all_saccades(file_path, file_name)
%% Handle inputs
if nargin < 1
    [file_name,file_path] = uigetfile([pwd filesep '*_ANALYZED.mat'], 'Select behavior file');
end

%% Load data
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
BEHAVE = load([file_path file_name]);
BEHAVE.EXPERIMENT_PARAMS.file_name = file_name;
BEHAVE.EXPERIMENT_PARAMS.file_path = file_path;
EPHYS.CH_EVE = load([file_path file_name(1:13) '_EVE1_aligned.mat']);
fprintf(' --> Completed. \n')

%% Build BEHAVE stream & aligned
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE stream & aligned', ' ... ']);
variable_list = {'time_1K','eye_r_px_filt','eye_r_py_filt', 'eye_r_vx_filt', 'eye_r_vy_filt', 'eye_r_vm_filt', 'tgt_px', 'tgt_py'};
BEHAVE.stream = struct;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    BEHAVE.stream.(variable_name) = cell2mat(BEHAVE.TRIALS_DATA.(variable_name)(:));
end

time_1K_stream = BEHAVE.stream.time_1K(:);
time_1K_aligned     = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K(:);
length_time_ = length(time_1K_aligned);
idx_stream_to_aligned  = nan(size(time_1K_aligned));
time_1K_stream(end+1)    = max([time_1K_aligned(end), time_1K_stream(end)])+1;
counter_time_1K_stream = find(time_1K_stream >= time_1K_aligned(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = time_1K_aligned(counter_time_point);
    if time_ponit_>=time_1K_stream(counter_time_1K_stream)
        idx_stream_to_aligned(counter_time_point) = counter_time_1K_stream;
        counter_time_1K_stream = counter_time_1K_stream + 1;
    end
end

idx_bool_nan = isnan(idx_stream_to_aligned);
idx_stream_to_aligned(idx_bool_nan) = 1;
BEHAVE.aligned = struct;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    variable_data_stream = BEHAVE.stream.(variable_name);
    variable_data_aligned = variable_data_stream(idx_stream_to_aligned);
    variable_data_aligned(idx_bool_nan) = NaN;
    BEHAVE.aligned.(variable_name) = variable_data_aligned;
end

fprintf(' --> Completed. \n')

%% Build BEHAVE SACS_ALL_DATA
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE SACS_ALL_DATA', ' ... ']);
threshold = 75; % deg/s
eye_vm_ = BEHAVE.aligned.eye_r_vm_filt(:);
eye_vm_ = [eye_vm_; eye_vm_(end)];
idx_bool_rise_threshold = (eye_vm_(2:end)-threshold > 0) & (eye_vm_(1:end-1)-threshold <= 0);
idx_int_rise_threshold = find(idx_bool_rise_threshold);
num_saccades = length(idx_int_rise_threshold);
eye_velocity_trace     = BEHAVE.aligned.eye_r_vm_filt(:);
eye_velocity_trace(isnan(eye_velocity_trace)) = 0;
num_sac_datapoints = 150;
all_sac_validity   = false(1, num_saccades);
% all_sac_inds       = nan(num_sac_datapoints, num_saccades);
all_sac_ind_start  = nan(1, num_saccades);
all_sac_ind_vmax   = nan(1, num_saccades);
all_sac_ind_finish = nan(1, num_saccades);
for counter_saccade = 1 : num_saccades
    ind_ = idx_int_rise_threshold(counter_saccade);
    ind_search_begin       = ind_-50;
    ind_search_end         = ind_+100;
    if(ind_search_begin < 1); ind_search_begin=1; end
    if(ind_search_end > length(eye_velocity_trace)); ind_search_end=length(eye_velocity_trace); end
    
    params_sac.MinPeakHeight       = threshold; % deg/s
    params_sac.MinPeakProminence   = 50; % data points
    params_sac.rough_threshold     = 50.0; % deg/s
    params_sac.fine_threshold      = 20.0; % deg/s
    params_sac.sampling_freq       = 1000.0; % Hz
    params_sac.cutoff_freq         = 50.0; % Hz
    params_sac.window_half_length  = 4; % data points
    params_sac.prominence_or_first = 'first'; % which peak to select, 'prominent' or 'first'
    
    output_ = ESN_Sac_Finder(eye_velocity_trace, ...
        ind_search_begin, ind_search_end, params_sac);
    
    all_sac_validity(:,counter_saccade)   = output_.validity(:);
%     all_sac_inds(:,counter_saccade)       = output_.inds(:);
    all_sac_ind_start(:,counter_saccade)  = output_.ind_start(:);
    all_sac_ind_vmax(:,counter_saccade)   = output_.ind_vmax(:);
    all_sac_ind_finish(:,counter_saccade) = output_.ind_finish(:);
end

BEHAVE.SACS_ALL_DATA = struct;
BEHAVE.SACS_ALL_DATA.validity   = reshape(all_sac_validity(   :,all_sac_validity), 1, []);
% BEHAVE.SACS_ALL_DATA.inds       = reshape(all_sac_inds(       :,all_sac_validity), num_sac_datapoints, []);
BEHAVE.SACS_ALL_DATA.ind_start  = reshape(all_sac_ind_start(  :,all_sac_validity), 1, []);
BEHAVE.SACS_ALL_DATA.ind_vmax   = reshape(all_sac_ind_vmax(   :,all_sac_validity), 1, []);
BEHAVE.SACS_ALL_DATA.ind_finish = reshape(all_sac_ind_finish( :,all_sac_validity), 1, []);
BEHAVE.SACS_ALL_DATA = build_sac_data(BEHAVE.SACS_ALL_DATA, BEHAVE.aligned);

for counter_iteration = 1 : 5
all_sac_validity = BEHAVE.SACS_ALL_DATA.validity;
all_sac_inds = BEHAVE.SACS_ALL_DATA.inds;
all_sac_ind_start = BEHAVE.SACS_ALL_DATA.ind_start;
all_sac_ind_vmax = BEHAVE.SACS_ALL_DATA.ind_vmax;
all_sac_ind_finish = BEHAVE.SACS_ALL_DATA.ind_finish;

eye_velocity_trace     = BEHAVE.aligned.eye_r_vm_filt(:);
eye_velocity_trace(isnan(eye_velocity_trace)) = 0;
all_sac_eye_r_vm = reshape(eye_velocity_trace(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_px = BEHAVE.SACS_ALL_DATA.eye_r_px;
all_sac_eye_r_py = BEHAVE.SACS_ALL_DATA.eye_r_py;
all_sac_eye_r_amp_m = BEHAVE.SACS_ALL_DATA.eye_r_amp_m;
all_sac_eye_r_vm_max = BEHAVE.SACS_ALL_DATA.eye_r_vm_max;
all_sac_duration = BEHAVE.SACS_ALL_DATA.duration;

all_sac_validity( max(all_sac_eye_r_vm) > 1000 ) = false; % invalidate if saccade trace has velocity > 1000 deg/s
all_sac_validity( max(all_sac_eye_r_vm(1:40,:)) > 300 ) = false; % invalidate if saccade trace during 1-40ms has velocity > 300 deg/s
all_sac_validity( max(all_sac_eye_r_vm(100:end,:)) > 300 ) = false; % invalidate if saccade trace during 100-150ms has velocity > 300 deg/s
all_sac_validity( all_sac_eye_r_vm_max < threshold ) = false; % invalidate if vm_max < threshold deg/s
all_sac_validity( all_sac_duration > 100 ) = false; % invalidate if saccades duration > 100 ms
all_sac_validity( all_sac_ind_finish-all_sac_ind_vmax > 90 ) = false; % invalidate if saccades diff > 90 ms
all_sac_validity( all_sac_ind_vmax-all_sac_ind_start > 60 ) = false; % invalidate if saccades diff > 60 ms
all_sac_validity( all_sac_ind_finish-all_sac_ind_vmax < 1 ) = false; % invalidate if saccades diff < 1 ms
all_sac_validity( all_sac_ind_vmax-all_sac_ind_start < 1 ) = false; % invalidate if saccades diff < 1 ms
all_sac_validity( all_sac_eye_r_amp_m > 15 ) = false;  % invalidate if saccade amplitude > 15 deg
all_sac_validity( all_sac_eye_r_amp_m < 0.5 ) = false;  % invalidate if saccade amplitude < 0.5 deg
all_sac_validity( max(abs(all_sac_eye_r_px)) > 20 ) = false;   % invalidate if abs of eye position > 20 deg
all_sac_validity( max(abs(all_sac_eye_r_py)) > 20 ) = false;   % invalidate if abs of eye position > 20 deg

all_sac_validity( [(abs(diff(all_sac_ind_vmax)) < 5) false] ) = false; % invalidate if ind_vmax is the same
all_sac_validity( [(abs(diff(all_sac_ind_start)) < 5) false] ) = false; % invalidate if ind_start is the same
all_sac_validity( [(abs(diff(all_sac_ind_finish)) < 5) false] ) = false; % invalidate if ind_finish is the same

BEHAVE.SACS_ALL_DATA.validity = all_sac_validity;
BEHAVE.SACS_ALL_DATA = build_sac_data(BEHAVE.SACS_ALL_DATA, BEHAVE.aligned);
end
fprintf(' --> Completed. \n')

%% extract BEHAVE_ind from BEHAVE_EB_xcorr_time_1K
clearvars -except EPHYS BEHAVE
fprintf(['Extract BEHAVE_ind', ' ... ']);
num_trials = length(BEHAVE.TRIALS_DATA.time_end);
BEHAVE_time_1K     = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K;
length_time_ = length(BEHAVE_time_1K);
BEHAVE_time_cue_present    = nan(num_trials, 1);
BEHAVE_time_primSac_onset  = nan(num_trials, 1);
BEHAVE_time_primSac_vmax  = nan(num_trials, 1);
BEHAVE_time_primSac_offset = nan(num_trials, 1);
BEHAVE_time_corrSac_onset  = nan(num_trials, 1);
BEHAVE_time_corrSac_vmax  = nan(num_trials, 1);
BEHAVE_time_corrSac_offset  = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    BEHAVE_time_cue_present(counter_trial) = BEHAVE.TRIALS_DATA.time_state_cue_present{1,counter_trial}(end);
    
    ind_primSac_onset_  = (BEHAVE.SACS_PRIM_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_onset_) ~= 1
        ind_primSac_onset_ = 1;
    end
    BEHAVE_time_primSac_onset(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_onset_, counter_trial);
    
    ind_primSac_vmax_  = (BEHAVE.SACS_PRIM_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_vmax_) ~= 1
        ind_primSac_vmax_ = 60;
    end
    BEHAVE_time_primSac_vmax(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_vmax_, counter_trial);
    
    ind_primSac_offset_ = (BEHAVE.SACS_PRIM_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_offset_) ~= 1
        ind_primSac_offset_ = 150;
    end
    BEHAVE_time_primSac_offset(counter_trial) = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_offset_, counter_trial);
    
    ind_corrSac_onset_  = (BEHAVE.SACS_CORR_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_onset_) ~= 1
        ind_corrSac_onset_ = 1;
    end
    BEHAVE_time_corrSac_onset(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_onset_, counter_trial);
    
    ind_corrSac_vmax_  = (BEHAVE.SACS_CORR_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_vmax_) ~= 1
        ind_corrSac_vmax_ = 60;
    end
    BEHAVE_time_corrSac_vmax(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_vmax_, counter_trial);
    
    ind_corrSac_offset_ = (BEHAVE.SACS_CORR_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_offset_) ~= 1
        ind_corrSac_offset_ = 150;
    end
    BEHAVE_time_corrSac_offset(counter_trial) = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_offset_, counter_trial);
end

% num_trials = sum(BEHAVE_time_corrSac_offset < BEHAVE_time_1K(end));
last_trial_num = find(BEHAVE_time_corrSac_offset < BEHAVE_time_1K(end), 1, 'last');
first_trial_num = find(BEHAVE_time_cue_present > BEHAVE_time_1K(1), 1, 'first');
range_trials = (first_trial_num : last_trial_num)';
num_trials = length(range_trials);

BEHAVE_time_cue_present    = BEHAVE_time_cue_present(   range_trials);
BEHAVE_time_primSac_onset  = BEHAVE_time_primSac_onset( range_trials);
BEHAVE_time_primSac_vmax   = BEHAVE_time_primSac_vmax(  range_trials);
BEHAVE_time_primSac_offset = BEHAVE_time_primSac_offset(range_trials);
BEHAVE_time_corrSac_onset  = BEHAVE_time_corrSac_onset( range_trials);
BEHAVE_time_corrSac_vmax   = BEHAVE_time_corrSac_vmax(  range_trials);
BEHAVE_time_corrSac_offset = BEHAVE_time_corrSac_offset(range_trials);

BEHAVE_time_cue_present(end+1)    = max([BEHAVE_time_1K(end), BEHAVE_time_cue_present(end)])+1;
BEHAVE_time_primSac_onset(end+1)  = max([BEHAVE_time_1K(end), BEHAVE_time_primSac_onset(end)])+1;
BEHAVE_time_primSac_vmax(end+1)   = max([BEHAVE_time_1K(end), BEHAVE_time_primSac_vmax(end)])+1;
BEHAVE_time_primSac_offset(end+1) = max([BEHAVE_time_1K(end), BEHAVE_time_primSac_offset(end)])+1;
BEHAVE_time_corrSac_onset(end+1)  = max([BEHAVE_time_1K(end), BEHAVE_time_corrSac_onset(end)])+1;
BEHAVE_time_corrSac_vmax(end+1)   = max([BEHAVE_time_1K(end), BEHAVE_time_corrSac_vmax(end)])+1;
BEHAVE_time_corrSac_offset(end+1) = max([BEHAVE_time_1K(end), BEHAVE_time_corrSac_offset(end)])+1;
BEHAVE_ind_cue_present    = nan(num_trials, 1);
BEHAVE_ind_primSac_onset  = nan(num_trials, 1);
BEHAVE_ind_primSac_vmax   = nan(num_trials, 1);
BEHAVE_ind_primSac_offset = nan(num_trials, 1);
BEHAVE_ind_corrSac_onset  = nan(num_trials, 1);
BEHAVE_ind_corrSac_vmax   = nan(num_trials, 1);
BEHAVE_ind_corrSac_offset = nan(num_trials, 1);
counter_cue_present    = find(BEHAVE_time_cue_present    >= BEHAVE_time_1K(1), 1, 'first');
counter_primSac_onset  = find(BEHAVE_time_primSac_onset  >= BEHAVE_time_1K(1), 1, 'first');
counter_primSac_vmax   = find(BEHAVE_time_primSac_vmax   >= BEHAVE_time_1K(1), 1, 'first');
counter_primSac_offset = find(BEHAVE_time_primSac_offset >= BEHAVE_time_1K(1), 1, 'first');
counter_corrSac_onset  = find(BEHAVE_time_corrSac_onset  >= BEHAVE_time_1K(1), 1, 'first');
counter_corrSac_vmax   = find(BEHAVE_time_corrSac_vmax   >= BEHAVE_time_1K(1), 1, 'first');
counter_corrSac_offset = find(BEHAVE_time_corrSac_offset >= BEHAVE_time_1K(1), 1, 'first');
counter_trial_cue_present    = 1;
counter_trial_primSac_onset  = 1;
counter_trial_primSac_vmax   = 1;
counter_trial_primSac_offset = 1;
counter_trial_corrSac_onset  = 1;
counter_trial_corrSac_vmax   = 1;
counter_trial_corrSac_offset = 1;
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_time_1K(counter_time_point);
    if time_ponit_>=BEHAVE_time_cue_present(counter_cue_present)
        BEHAVE_ind_cue_present(counter_trial_cue_present) = counter_time_point;
        counter_cue_present = counter_cue_present + 1;
        counter_trial_cue_present = counter_trial_cue_present + 1;
    end
    if time_ponit_>=BEHAVE_time_primSac_onset(counter_primSac_onset)
        BEHAVE_ind_primSac_onset(counter_trial_primSac_onset) = counter_time_point;
        counter_primSac_onset = counter_primSac_onset + 1;
        counter_trial_primSac_onset = counter_trial_primSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_time_primSac_vmax(counter_primSac_vmax)
        BEHAVE_ind_primSac_vmax(counter_trial_primSac_vmax) = counter_time_point;
        counter_primSac_vmax = counter_primSac_vmax + 1;
        counter_trial_primSac_vmax = counter_trial_primSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_time_primSac_offset(counter_primSac_offset)
        BEHAVE_ind_primSac_offset(counter_trial_primSac_offset) = counter_time_point;
        counter_primSac_offset = counter_primSac_offset + 1;
        counter_trial_primSac_offset = counter_trial_primSac_offset + 1;
    end
    if time_ponit_>=BEHAVE_time_corrSac_onset(counter_corrSac_onset)
        BEHAVE_ind_corrSac_onset(counter_trial_corrSac_onset) = counter_time_point;
        counter_corrSac_onset = counter_corrSac_onset + 1;
        counter_trial_corrSac_onset = counter_trial_corrSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_time_corrSac_vmax(counter_corrSac_vmax)
        BEHAVE_ind_corrSac_vmax(counter_trial_corrSac_vmax) = counter_time_point;
        counter_corrSac_vmax = counter_corrSac_vmax + 1;
        counter_trial_corrSac_vmax = counter_trial_corrSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_time_corrSac_offset(counter_corrSac_offset)
        BEHAVE_ind_corrSac_offset(counter_trial_corrSac_offset) = counter_time_point;
        counter_corrSac_offset = counter_corrSac_offset + 1;
        counter_trial_corrSac_offset = counter_trial_corrSac_offset + 1;
    end
end

BEHAVE.aligned.BEHAVE_range_trials       = range_trials;
BEHAVE.aligned.BEHAVE_ind_cue_present    = BEHAVE_ind_cue_present;
BEHAVE.aligned.BEHAVE_ind_primSac_onset  = BEHAVE_ind_primSac_onset;
BEHAVE.aligned.BEHAVE_ind_primSac_vmax   = BEHAVE_ind_primSac_vmax;
BEHAVE.aligned.BEHAVE_ind_primSac_offset = BEHAVE_ind_primSac_offset;
BEHAVE.aligned.BEHAVE_ind_corrSac_onset  = BEHAVE_ind_corrSac_onset;
BEHAVE.aligned.BEHAVE_ind_corrSac_vmax   = BEHAVE_ind_corrSac_vmax;
BEHAVE.aligned.BEHAVE_ind_corrSac_offset = BEHAVE_ind_corrSac_offset;
BEHAVE.aligned.BEHAVE_validity_prim      = BEHAVE.SACS_PRIM_DATA.validity(:,range_trials)';
BEHAVE.aligned.BEHAVE_validity_corr      = BEHAVE.SACS_CORR_DATA.validity(:,range_trials)';
fprintf(' --> Completed. \n')

%% Resolve Prim Corr conflicts and misses between BEHAVE_ind and SACS_ALL_DATA
clearvars -except EPHYS BEHAVE
ind_BEHAVE_primSac = reshape(BEHAVE.aligned.BEHAVE_ind_primSac_vmax(BEHAVE.aligned.BEHAVE_validity_prim),[], 1);
idx_int_primSac = knnsearch(BEHAVE.SACS_ALL_DATA.ind_vmax(:),ind_BEHAVE_primSac);
ind_SACS_ALL_DATA_primSac = reshape(BEHAVE.SACS_ALL_DATA.ind_vmax(:,idx_int_primSac),[], 1);
ind_primSac_overlap = abs(ind_BEHAVE_primSac - ind_SACS_ALL_DATA_primSac) < 5; 

ind_BEHAVE_corrSac = reshape(BEHAVE.aligned.BEHAVE_ind_corrSac_vmax(BEHAVE.aligned.BEHAVE_validity_corr),[], 1);
idx_int_corrSac = knnsearch(BEHAVE.SACS_ALL_DATA.ind_vmax(:),ind_BEHAVE_corrSac);
ind_SACS_ALL_DATA_corrSac = reshape(BEHAVE.SACS_ALL_DATA.ind_vmax(:,idx_int_corrSac),[], 1);
ind_corrSac_overlap = abs(ind_BEHAVE_corrSac - ind_SACS_ALL_DATA_corrSac) < 5; 

BEHAVE.SACS_ALL_DATA.validity(idx_int_primSac(ind_primSac_overlap)) = false;
BEHAVE.SACS_ALL_DATA.validity(idx_int_corrSac(ind_corrSac_overlap)) = false;
BEHAVE.SACS_ALL_DATA = build_sac_data(BEHAVE.SACS_ALL_DATA, BEHAVE.aligned);

primSac_to_be_included_ind_vmax   = reshape(BEHAVE.aligned.BEHAVE_ind_primSac_vmax,  1, []);
primSac_to_be_included_ind_start  = reshape(BEHAVE.aligned.BEHAVE_ind_primSac_onset, 1, []);
primSac_to_be_included_ind_finish = reshape(BEHAVE.aligned.BEHAVE_ind_primSac_offset,1, []);

corrSac_to_be_included_ind_vmax   = reshape(BEHAVE.aligned.BEHAVE_ind_corrSac_vmax,  1, []);
corrSac_to_be_included_ind_start  = reshape(BEHAVE.aligned.BEHAVE_ind_corrSac_onset, 1, []);
corrSac_to_be_included_ind_finish = reshape(BEHAVE.aligned.BEHAVE_ind_corrSac_offset,1, []);

BEHAVE.SACS_ALL_DATA.ind_start  = [BEHAVE.SACS_ALL_DATA.ind_start  primSac_to_be_included_ind_start  corrSac_to_be_included_ind_start];
BEHAVE.SACS_ALL_DATA.ind_vmax   = [BEHAVE.SACS_ALL_DATA.ind_vmax   primSac_to_be_included_ind_vmax   corrSac_to_be_included_ind_vmax];
BEHAVE.SACS_ALL_DATA.ind_finish = [BEHAVE.SACS_ALL_DATA.ind_finish primSac_to_be_included_ind_finish corrSac_to_be_included_ind_finish];

BEHAVE.SACS_ALL_DATA.is_primSac = [false(1,length(BEHAVE.SACS_ALL_DATA.validity)) true( 1,length(ind_primSac_overlap)) false(1,length(ind_corrSac_overlap)) ];
BEHAVE.SACS_ALL_DATA.is_corrSac = [false(1,length(BEHAVE.SACS_ALL_DATA.validity)) false(1,length(ind_primSac_overlap)) true( 1,length(ind_corrSac_overlap)) ];
BEHAVE.SACS_ALL_DATA.validity   = [BEHAVE.SACS_ALL_DATA.validity true(1,length(ind_primSac_overlap)) true(1,length(ind_corrSac_overlap))];

BEHAVE.SACS_ALL_DATA = build_sac_data(BEHAVE.SACS_ALL_DATA, BEHAVE.aligned);

%% Finding saccade categories
clearvars -except EPHYS BEHAVE
fprintf(['Finding saccade categories ', ' ... ']);
all_sac_eye_r_px_finish  = BEHAVE.SACS_ALL_DATA.eye_r_px_finish;
all_sac_eye_r_py_finish  = BEHAVE.SACS_ALL_DATA.eye_r_py_finish;
all_sac_eye_r_px_start   = BEHAVE.SACS_ALL_DATA.eye_r_px_start; 
all_sac_eye_r_py_start   = BEHAVE.SACS_ALL_DATA.eye_r_py_start; 
all_sac_tgt_px_finish    = reshape(BEHAVE.aligned.tgt_px(BEHAVE.SACS_ALL_DATA.ind_finish),1, []);
all_sac_tgt_py_finish    = reshape(BEHAVE.aligned.tgt_py(BEHAVE.SACS_ALL_DATA.ind_finish),1, []);
all_sac_tgt_px_start     = reshape(BEHAVE.aligned.tgt_px(BEHAVE.SACS_ALL_DATA.ind_start),1, []);
all_sac_tgt_py_start     = reshape(BEHAVE.aligned.tgt_py(BEHAVE.SACS_ALL_DATA.ind_start),1, []);

tgtStr_px = nanmean(BEHAVE.TRIALS_DATA.start_x(:,BEHAVE.aligned.BEHAVE_range_trials));
tgtStr_py = nanmean(BEHAVE.TRIALS_DATA.start_y(:,BEHAVE.aligned.BEHAVE_range_trials));

distance_eye_to_tgt_finish    = sqrt( (all_sac_eye_r_px_finish - all_sac_tgt_px_finish).^2 + (all_sac_eye_r_py_finish - all_sac_tgt_py_finish).^2 );
distance_eye_to_tgt_start     = sqrt( (all_sac_eye_r_px_start  - all_sac_tgt_px_start).^2  + (all_sac_eye_r_py_start  - all_sac_tgt_py_start).^2 );
distance_eye_to_tgtStr_finish = sqrt( (all_sac_eye_r_px_finish - tgtStr_px).^2 + (all_sac_eye_r_py_finish - tgtStr_py).^2 );
distance_eye_to_tgtStr_start  = sqrt( (all_sac_eye_r_px_start  - tgtStr_px).^2 + (all_sac_eye_r_py_start  - tgtStr_py).^2 );

BEHAVE.SACS_ALL_DATA.is_toTgt        = (distance_eye_to_tgt_finish < 2.0);
BEHAVE.SACS_ALL_DATA.is_fromTgt      = (distance_eye_to_tgt_start < 2.0);
BEHAVE.SACS_ALL_DATA.is_lowAmp       = (BEHAVE.SACS_ALL_DATA.eye_r_amp_m < 2.0);
BEHAVE.SACS_ALL_DATA.is_toTgtStr     = (distance_eye_to_tgtStr_finish < 2.0);
BEHAVE.SACS_ALL_DATA.is_fromCenter   = (distance_eye_to_tgtStr_start < 2.0);
BEHAVE.SACS_ALL_DATA.is_aroundCenter = (distance_eye_to_tgtStr_start < 2.0) & (distance_eye_to_tgtStr_finish < 2.0);

BEHAVE.SACS_ALL_DATA.is_aroundCenter = BEHAVE.SACS_ALL_DATA.is_aroundCenter | BEHAVE.SACS_ALL_DATA.is_lowAmp;
BEHAVE.SACS_ALL_DATA.is_toTgtStr = BEHAVE.SACS_ALL_DATA.is_toTgtStr & BEHAVE.SACS_ALL_DATA.is_toTgt;
BEHAVE.SACS_ALL_DATA.is_nonTask = ~(BEHAVE.SACS_ALL_DATA.is_primSac | ...
                                    BEHAVE.SACS_ALL_DATA.is_corrSac | ...
                                    BEHAVE.SACS_ALL_DATA.is_toTgtStr | ...
                                    BEHAVE.SACS_ALL_DATA.is_fromCenter | ...
                                    BEHAVE.SACS_ALL_DATA.is_aroundCenter | ...
                                    BEHAVE.SACS_ALL_DATA.is_toTgt | ...
                                    BEHAVE.SACS_ALL_DATA.is_fromTgt);

BEHAVE.SACS_ALL_DATA.is_corrSac(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_toTgt(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_fromTgt(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_lowAmp(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_toTgtStr(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_fromCenter(BEHAVE.SACS_ALL_DATA.is_primSac) = false;
BEHAVE.SACS_ALL_DATA.is_aroundCenter(BEHAVE.SACS_ALL_DATA.is_primSac) = false;

BEHAVE.SACS_ALL_DATA.is_toTgt(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;
BEHAVE.SACS_ALL_DATA.is_fromTgt(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;
BEHAVE.SACS_ALL_DATA.is_lowAmp(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;
BEHAVE.SACS_ALL_DATA.is_toTgtStr(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;
BEHAVE.SACS_ALL_DATA.is_fromCenter(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;
BEHAVE.SACS_ALL_DATA.is_aroundCenter(BEHAVE.SACS_ALL_DATA.is_corrSac) = false;

BEHAVE.SACS_ALL_DATA.is_toTgt(BEHAVE.SACS_ALL_DATA.is_aroundCenter) = false;
BEHAVE.SACS_ALL_DATA.is_fromTgt(BEHAVE.SACS_ALL_DATA.is_aroundCenter) = false;
BEHAVE.SACS_ALL_DATA.is_toTgtStr(BEHAVE.SACS_ALL_DATA.is_aroundCenter) = false;
BEHAVE.SACS_ALL_DATA.is_fromCenter(BEHAVE.SACS_ALL_DATA.is_aroundCenter) = false;

fprintf(' --> Completed. \n')

%% Plot
clearvars -except EPHYS BEHAVE
amp_edges = -.25 : 0.5 : 15.25;
ang_edges = -191.25 : 22.5 : +191.25;
hFig = figure(1);
clf(hFig)
hold on

subplot(6,6,[1 2 7 8])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_primSac), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_primSac), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_primSac), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_primSac), 'om')
title(['primary' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_primSac)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[3 4 9 10])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), 'om')
title(['to start tgt' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_toTgtStr)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[5 6 11 12])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), 'om')
title(['low amplitude' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_aroundCenter)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[13])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_primSac), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_primSac), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[14])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_primSac), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_primSac), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

subplot(6,6,[15])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[16])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_toTgtStr), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

subplot(6,6,[17])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[18])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_aroundCenter), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

subplot(6,6,[19 20 25 26])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_corrSac), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_corrSac), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_corrSac), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_corrSac), 'om')
title(['corrective' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_corrSac)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[21 22 27 28])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), 'om')
title(['from center' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_fromCenter)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[23 24 29 30])
hold on
plot(BEHAVE.SACS_ALL_DATA.eye_r_px(:,BEHAVE.SACS_ALL_DATA.is_nonTask), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py(:,BEHAVE.SACS_ALL_DATA.is_nonTask), 'k')
plot(BEHAVE.SACS_ALL_DATA.eye_r_px_finish(:,BEHAVE.SACS_ALL_DATA.is_nonTask), ...
     BEHAVE.SACS_ALL_DATA.eye_r_py_finish(:,BEHAVE.SACS_ALL_DATA.is_nonTask), 'om')
title(['non task related' ', ' num2str(sum(BEHAVE.SACS_ALL_DATA.is_nonTask)) ' saccades'])
xlim([-15, 15])
ylim([-15, 15])

subplot(6,6,[31])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_corrSac), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_corrSac), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[32])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_corrSac), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_corrSac), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

subplot(6,6,[33])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[34])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_fromCenter), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

subplot(6,6,[35])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_nonTask), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_amp_m(:,BEHAVE.SACS_ALL_DATA.is_nonTask), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:2:14)
ylabel('Amplitude Hist')

subplot(6,6,[36])
hold on
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_nonTask), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(BEHAVE.SACS_ALL_DATA.eye_r_ang(:,BEHAVE.SACS_ALL_DATA.is_nonTask), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([-190, 190])
set(gca, 'XTick', -180:90:180)
ylabel('Angle Hist')

sgtitle(BEHAVE.EXPERIMENT_PARAMS.file_name(1:13), 'interpret', 'none');
ESN_Beautify_Plot(hFig, [18 13])

%% Save _ANALYZED.mat Data to disk
clearvars -except EPHYS BEHAVE hFig
EXPERIMENT_PARAMS = BEHAVE.EXPERIMENT_PARAMS;
TRIALS_DATA = BEHAVE.TRIALS_DATA;
SACS_PRIM_DATA = BEHAVE.SACS_PRIM_DATA;
SACS_CORR_DATA = BEHAVE.SACS_CORR_DATA;
SACS_ALL_DATA = BEHAVE.SACS_ALL_DATA;
aligned = BEHAVE.aligned;
file_path_ = BEHAVE.EXPERIMENT_PARAMS.file_path;
if ~strcmp(file_path_(end), filesep);file_path_ = [file_path_ filesep];end
file_name_ = BEHAVE.EXPERIMENT_PARAMS.file_name;
fprintf(['Saving ' file_name_ ' ...'])
save([file_path_ file_name_], ...
    'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_PRIM_DATA', 'SACS_CORR_DATA', 'SACS_ALL_DATA', 'aligned', '-v7.3');
fprintf(' --> Completed. \n')

%% Save _REDUCED.mat Data to disk
rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
    'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
    'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
    'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
    'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};
TRIALS_DATA = rmfield(TRIALS_DATA,rmfields_list);
file_name_REDUCED_ = [file_name_(1:13) '_REDUCED.mat'];
fprintf(['Saving ' file_name_REDUCED_ ' ...'])
save([file_path_ file_name_REDUCED_], ...
    'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_PRIM_DATA', 'SACS_CORR_DATA', 'SACS_ALL_DATA', '-v7.3');
fprintf(' --> Completed. \n')

%% Save Fig
file_name_plot_ = file_name_(1:end-4);
file_path_plot_ = [file_path_ '..' filesep 'analyzed_figs' filesep];
fprintf(['Saving ' file_name_plot_ ' ...'])
saveas(hFig,[file_path_plot_ file_name_plot_ '_allSaccades'], 'pdf');
saveas(hFig,[file_path_plot_ file_name_plot_ '_allSaccades'], 'png');
close(hFig)
fprintf(' --> Completed. \n')
end

function SACS_ALL_DATA = build_sac_data(SACS_ALL_DATA, aligned)
all_sac_validity   = reshape(SACS_ALL_DATA.validity,   1, []);
all_sac_ind_start  = reshape(SACS_ALL_DATA.ind_start,  1, []);
all_sac_ind_vmax   = reshape(SACS_ALL_DATA.ind_vmax,   1, []);
all_sac_ind_finish = reshape(SACS_ALL_DATA.ind_finish, 1, []);
% all_sac_inds       = reshape(SACS_ALL_DATA.inds,       num_sac_datapoints, []);

inds_span_    = ((-60+1) : 1 : (90))';
num_sac_datapoints = length(inds_span_);
length_time_      = length(aligned.time_1K);

all_sac_inds = repmat( all_sac_ind_vmax(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(all_sac_ind_vmax), 1);
all_sac_inds( all_sac_inds < 1 ) = 1;
all_sac_inds( all_sac_inds > length_time_ ) = length_time_;
all_sac_inds = all_sac_inds';

all_sac_time            = reshape(aligned.time_1K(      all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_px        = reshape(aligned.eye_r_px_filt(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_py        = reshape(aligned.eye_r_py_filt(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_vx        = reshape(aligned.eye_r_vx_filt(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_vy        = reshape(aligned.eye_r_vy_filt(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_vm        = reshape(aligned.eye_r_vm_filt(all_sac_inds),num_sac_datapoints, []);
all_sac_eye_r_vm_max    = reshape(aligned.eye_r_vm_filt(all_sac_ind_vmax),1, []);
all_sac_eye_r_px_start  = reshape(aligned.eye_r_px_filt(all_sac_ind_start),1, []);
all_sac_eye_r_px_finish = reshape(aligned.eye_r_px_filt(all_sac_ind_finish),1, []);
all_sac_eye_r_py_start  = reshape(aligned.eye_r_py_filt(all_sac_ind_start),1, []);
all_sac_eye_r_py_finish = reshape(aligned.eye_r_py_filt(all_sac_ind_finish),1, []);
all_sac_eye_r_amp_x     = reshape((all_sac_eye_r_px_finish - all_sac_eye_r_px_start),1, []);
all_sac_eye_r_amp_y     = reshape((all_sac_eye_r_py_finish - all_sac_eye_r_py_start),1, []);
all_sac_eye_r_amp_m     = reshape((sqrt(  all_sac_eye_r_amp_y.^2 + all_sac_eye_r_amp_x.^2)),1, []);
all_sac_eye_r_ang       = reshape((atan2d(all_sac_eye_r_amp_y,     all_sac_eye_r_amp_x)),   1, []);
all_sac_duration        = reshape((all_sac_ind_finish - all_sac_ind_start),1, []);

SACS_ALL_DATA.validity        = all_sac_validity(       :,all_sac_validity);
SACS_ALL_DATA.inds            = all_sac_inds(           :,all_sac_validity);
SACS_ALL_DATA.ind_start       = all_sac_ind_start(      :,all_sac_validity);
SACS_ALL_DATA.ind_vmax        = all_sac_ind_vmax(       :,all_sac_validity);
SACS_ALL_DATA.ind_finish      = all_sac_ind_finish(     :,all_sac_validity);
SACS_ALL_DATA.time            = all_sac_time(           :,all_sac_validity);
SACS_ALL_DATA.eye_r_px        = all_sac_eye_r_px(       :,all_sac_validity);
SACS_ALL_DATA.eye_r_py        = all_sac_eye_r_py(       :,all_sac_validity);
SACS_ALL_DATA.eye_r_vx        = all_sac_eye_r_vx(       :,all_sac_validity);
SACS_ALL_DATA.eye_r_vy        = all_sac_eye_r_vy(       :,all_sac_validity);
SACS_ALL_DATA.eye_r_vm        = all_sac_eye_r_vm(       :,all_sac_validity);
SACS_ALL_DATA.eye_r_vm_max    = all_sac_eye_r_vm_max(   :,all_sac_validity);
SACS_ALL_DATA.eye_r_px_start  = all_sac_eye_r_px_start( :,all_sac_validity);
SACS_ALL_DATA.eye_r_px_finish = all_sac_eye_r_px_finish(:,all_sac_validity);
SACS_ALL_DATA.eye_r_py_start  = all_sac_eye_r_py_start( :,all_sac_validity);
SACS_ALL_DATA.eye_r_py_finish = all_sac_eye_r_py_finish(:,all_sac_validity);
SACS_ALL_DATA.eye_r_amp_x     = all_sac_eye_r_amp_x(    :,all_sac_validity);
SACS_ALL_DATA.eye_r_amp_y     = all_sac_eye_r_amp_y(    :,all_sac_validity);
SACS_ALL_DATA.eye_r_amp_m     = all_sac_eye_r_amp_m(    :,all_sac_validity);
SACS_ALL_DATA.eye_r_ang       = all_sac_eye_r_ang(      :,all_sac_validity);
SACS_ALL_DATA.duration        = all_sac_duration(       :,all_sac_validity);
end


