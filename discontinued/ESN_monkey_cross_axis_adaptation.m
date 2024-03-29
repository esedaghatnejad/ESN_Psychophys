function ESN_monkey_cross_axis_adaptation(pathnames)
%% Check nargin
% if there is no inputs, then set pathnames to pwd
if nargin < 1
    pathnames = pwd;
end
%% Get list of files
% add '/' to the end of pathnames
if ~strcmp(pathnames(end), filesep)
    pathnames = [pathnames filesep];
end
% get the list of MAT files to analyze
x_axis_files = dir([pathnames 'cross_axis_adaptation_*.mat']);
x_axis_files_cell = struct2cell(x_axis_files);
filenames = x_axis_files_cell(1,:);
% if there is no cross_axis_adaptation file, then terminate the function
if isempty(filenames)
    fprintf('ERROR in sac_analyze_monkey_x_axis: no cross_axis_adaptation_*.mat has been found in path.\n');
end

SESSION_PARAMS.filenames = sort(filenames);
SESSION_PARAMS.pathnames = pathnames;

%% MAIN FILE LOOP
for counter_file = 1 : 1 : length(SESSION_PARAMS.filenames)
    fprintf('#######################################\n')
    %% Load Data
    fprintf('Loading ...\n')
    clearvars -except counter_file SESSION_PARAMS TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL;
    filename = SESSION_PARAMS.filenames{counter_file};
    pathname = SESSION_PARAMS.pathnames;
    load([pathname filename]);
    fprintf('%s: Loading complete.\n', filename);
    [~, foldername, ~] = fileparts(pathname(1:end-1)); % (1:end-1) to delete the '/' character
    SESSION_PARAMS.filename{counter_file} = filename;
    SESSION_PARAMS.pathname{counter_file} = pathname;
    SESSION_PARAMS.foldername{counter_file} = foldername;
    SESSION_PARAMS.num_trials(counter_file) = length(data.trials);
    
    %% Analyze trials
    clearvars -except counter_file SESSION_PARAMS TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL data;
    num_trials = length(data.trials);
    if ~exist('counter_file','var'); counter_file = 1; end;
    fprintf([SESSION_PARAMS.filename{counter_file} ': Analyzing TRIALS ...'])
    for counter_trial = 1 : 1 : num_trials-1
        %% Extract Trial Varibales
        clearvars -except counter_file SESSION_PARAMS ...
            TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
            TRIALS SACS_PRIM SACS_CORR data counter_trial;
        if ~exist('counter_trial','var'); counter_trial = 1; end;
        % get trial struct
        trial_struct = data.trials{1, counter_trial};
        % trial start/end time
        TRIAL.time_start = double(trial_struct.trial_start_time);
        TRIAL.time_end   = double(trial_struct.trial_end_time);
        % trial variables
        TRIAL.start_x                    = ESN_Round(trial_struct.start_x, 0.01);
        TRIAL.start_y                    = ESN_Round(trial_struct.start_y, 0.01);
        TRIAL.cue_x                      = ESN_Round(trial_struct.cue_x, 0.01);
        TRIAL.cue_y                      = ESN_Round(trial_struct.cue_y, 0.01);
        TRIAL.end_x                      = ESN_Round(trial_struct.end_x, 0.01);
        TRIAL.end_y                      = ESN_Round(trial_struct.end_y, 0.01);
        TRIAL.iss_x                      = ESN_Round(trial_struct.iss_x, 0.01);
        TRIAL.iss_y                      = ESN_Round(trial_struct.iss_y, 0.01);
        if isfield(TRIAL, 'pursuit_x')
            TRIAL.pursuit_x                  = ESN_Round(trial_struct.pursuit_x, 0.01);
            TRIAL.pursuit_y                  = ESN_Round(trial_struct.pursuit_y, 0.01);
            TRIAL.reward_area                = ESN_Round(trial_struct.reward_area, 0.01);
            TRIAL.time_pursuit               = double(trial_struct.pursuit_duration);
        else
            TRIAL.pursuit_x                  = nan;
            TRIAL.pursuit_y                  = nan;
            TRIAL.reward_area                = nan;
            TRIAL.time_pursuit               = nan;
        end
        TRIAL.time_state_str_pursuit     = double(trial_struct.state_start_time_str_target_pursuit);
        TRIAL.time_state_str_present     = double(trial_struct.state_start_time_str_target_present);
        TRIAL.time_state_str_fixation    = double(trial_struct.state_start_time_str_target_fixation);
        TRIAL.time_state_cue_present     = double(trial_struct.state_start_time_cue_target_present);
        TRIAL.time_state_sac_detect_on   = double(trial_struct.state_start_time_detect_sac_start);
        TRIAL.time_state_sac_onset       = double(trial_struct.state_start_time_saccade);
        TRIAL.time_state_sac_detect_off  = double(trial_struct.state_start_time_detect_sac_end);
        TRIAL.time_state_reward          = double(trial_struct.state_start_time_deliver_reward);
        TRIAL.time_state_end_fixation    = double(trial_struct.state_start_time_end_target_fixation);
        TRIAL.time_state_iti             = double(trial_struct.state_start_time_iti);
        TRIAL.time_state_next_trial      = double(trial_struct.state_start_time_next_trial);
        TRIAL.time_iti                   = double(trial_struct.iti);
        TRIAL.time_punishment            = double(trial_struct.punishment_time);
        TRIAL.time_fixation              = double(trial_struct.fixation_time);
        % trial start/end indices
        length_data = min([length(data.eyelink_time) length(data.t)]);
        time_array = double(data.t(1: length_data));
        TRIAL.ind_trial_str   = find(time_array>TRIAL.time_start, 1, 'first');
        TRIAL.ind_trial_end   = find(time_array<TRIAL.time_end, 1, 'last');
        TRIAL.inds_trial          = TRIAL.ind_trial_str:TRIAL.ind_trial_end;
        % trial timeseries
        TRIAL.inds_invalid   = false(1, length(TRIAL.inds_trial));
        TRIAL.time_eyelink   = double(data.eyelink_time(1, TRIAL.inds_trial));                           TRIAL.inds_invalid = isnan(TRIAL.time_eyelink) | TRIAL.inds_invalid;
        TRIAL.eye_r_px       = double(data.right_horizontal_eye(1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_px)     | TRIAL.inds_invalid;
        TRIAL.eye_r_py       = double(data.right_vertical_eye(  1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_py)     | TRIAL.inds_invalid;
        TRIAL.eye_l_px       = double(data.left_horizontal_eye( 1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_px)     | TRIAL.inds_invalid;
        TRIAL.eye_l_py       = double(data.left_vertical_eye(   1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_py)     | TRIAL.inds_invalid;
        TRIAL.eye_r_vx       = double(data.right_horizontal_eye_velocity_filtered(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vx)     | TRIAL.inds_invalid;
        TRIAL.eye_r_vy       = double(data.right_vertical_eye_velocity_filtered(  1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vy)     | TRIAL.inds_invalid;
        TRIAL.eye_l_vx       = double(data.left_horizontal_eye_velocity_filtered( 1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vx)     | TRIAL.inds_invalid;
        TRIAL.eye_l_vy       = double(data.left_vertical_eye_velocity_filtered(   1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vy)     | TRIAL.inds_invalid;
        TRIAL.eye_r_vm       = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
        TRIAL.eye_l_vm       = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
        TRIAL.time_tgt       = double(data.t(1, TRIAL.inds_trial));                                      TRIAL.inds_invalid = isnan(TRIAL.time_tgt)     | TRIAL.inds_invalid;
        TRIAL.tgt_px         = double(data.target_x(1, TRIAL.inds_trial));                               TRIAL.inds_invalid = isnan(TRIAL.tgt_px)       | TRIAL.inds_invalid;
        TRIAL.tgt_py         = double(data.target_y(1, TRIAL.inds_trial));                               TRIAL.inds_invalid = isnan(TRIAL.tgt_py)       | TRIAL.inds_invalid;
        TRIAL.reward         = double(data.reward(1, TRIAL.inds_trial));                                 TRIAL.inds_invalid = isnan(TRIAL.reward)       | TRIAL.inds_invalid;
        TRIAL.target_visible = logical(double(data.target_visible(1, TRIAL.inds_trial)));
        % correct for the bias between time_eyelink and time_tgt
        TRIAL.time_eyelink   = TRIAL.time_eyelink .* (TRIAL.time_tgt(end)-TRIAL.time_tgt(1)) ./ (TRIAL.time_eyelink(end)-TRIAL.time_eyelink(1));
        TRIAL.time_eyelink   = TRIAL.time_eyelink - TRIAL.time_eyelink(1) + TRIAL.time_tgt(1);
        TRIAL.time_1K        = TRIAL.time_eyelink(1) : 0.001 : TRIAL.time_eyelink(end);
        % make non unique points of eye traces invalid
        TRIAL.inds_invalid = ([false (diff(TRIAL.time_eyelink)==0)])       | TRIAL.inds_invalid;
        % remove invalid values
        TRIAL.time_eyelink(TRIAL.inds_invalid) = [];
        TRIAL.eye_r_px(    TRIAL.inds_invalid) = [];
        TRIAL.eye_r_py(    TRIAL.inds_invalid) = [];
        TRIAL.eye_l_px(    TRIAL.inds_invalid) = [];
        TRIAL.eye_l_py(    TRIAL.inds_invalid) = [];
        % reconstruct eye_r data
        TRIAL.eye_r_px = interp1(TRIAL.time_eyelink, TRIAL.eye_r_px, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_r_py = interp1(TRIAL.time_eyelink, TRIAL.eye_r_py, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_r_vx = diff(TRIAL.eye_r_px)./diff(TRIAL.time_1K); TRIAL.eye_r_vx=[TRIAL.eye_r_vx(1) TRIAL.eye_r_vx];
        TRIAL.eye_r_vy = diff(TRIAL.eye_r_py)./diff(TRIAL.time_1K); TRIAL.eye_r_vy=[TRIAL.eye_r_vy(1) TRIAL.eye_r_vy];
        TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
        % reconstruct eye_l data
        TRIAL.eye_l_px = interp1(TRIAL.time_eyelink, TRIAL.eye_l_px, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_l_py = interp1(TRIAL.time_eyelink, TRIAL.eye_l_py, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_l_vx = diff(TRIAL.eye_l_px)./diff(TRIAL.time_1K); TRIAL.eye_l_vx=[TRIAL.eye_l_vx(1) TRIAL.eye_l_vx];
        TRIAL.eye_l_vy = diff(TRIAL.eye_l_py)./diff(TRIAL.time_1K); TRIAL.eye_l_vy=[TRIAL.eye_l_vy(1) TRIAL.eye_l_vy];
        TRIAL.eye_l_vm = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
        TRIAL.time     = TRIAL.time_1K;
        % filter params
        sampling_freq = 1000.0;
        cutoff_freq = 100.0;
        [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
        % filter eye_r data
        TRIAL.eye_r_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_px);
        TRIAL.eye_r_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_py);
        TRIAL.eye_r_vx_filt = diff(TRIAL.eye_r_px_filt)./diff(TRIAL.time_1K); TRIAL.eye_r_vx_filt=[TRIAL.eye_r_vx_filt(1) TRIAL.eye_r_vx_filt];
        TRIAL.eye_r_vy_filt = diff(TRIAL.eye_r_py_filt)./diff(TRIAL.time_1K); TRIAL.eye_r_vy_filt=[TRIAL.eye_r_vy_filt(1) TRIAL.eye_r_vy_filt];
        TRIAL.eye_r_vm_filt = sqrt(TRIAL.eye_r_vx_filt.^2 + TRIAL.eye_r_vy_filt.^2);
        % filter eye_l data
        TRIAL.eye_l_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_px);
        TRIAL.eye_l_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_py);
        TRIAL.eye_l_vx_filt = diff(TRIAL.eye_l_px_filt)./diff(TRIAL.time_1K); TRIAL.eye_l_vx_filt=[TRIAL.eye_l_vx_filt(1) TRIAL.eye_l_vx_filt];
        TRIAL.eye_l_vy_filt = diff(TRIAL.eye_l_py_filt)./diff(TRIAL.time_1K); TRIAL.eye_l_vy_filt=[TRIAL.eye_l_vy_filt(1) TRIAL.eye_l_vy_filt];
        TRIAL.eye_l_vm_filt = sqrt(TRIAL.eye_l_vx_filt.^2 + TRIAL.eye_l_vy_filt.^2);
        % trial state indices
        TRIAL.ind_state_str_pursuit    = find(TRIAL.time>TRIAL.time_state_str_pursuit(end), 1, 'first');
        TRIAL.ind_state_str_present    = find(TRIAL.time>TRIAL.time_state_str_present(end), 1, 'first');
        TRIAL.ind_state_str_fixation   = find(TRIAL.time>TRIAL.time_state_str_fixation(end), 1, 'first');
        TRIAL.ind_state_cue_present    = find(TRIAL.time>TRIAL.time_state_cue_present(end), 1, 'first');
        TRIAL.ind_state_sac_detect_on  = find(TRIAL.time>TRIAL.time_state_sac_detect_on(end), 1, 'first');
        TRIAL.ind_state_sac_onset      = find(TRIAL.time>TRIAL.time_state_sac_onset(end), 1, 'first');
        TRIAL.ind_state_sac_detect_off = find(TRIAL.time>TRIAL.time_state_sac_detect_off(end), 1, 'first');
        TRIAL.ind_state_reward         = find(TRIAL.time>TRIAL.time_state_reward(end), 1, 'first');
        TRIAL.ind_state_end_fixation   = find(TRIAL.time>TRIAL.time_state_end_fixation(end), 1, 'first');
        TRIAL.ind_state_iti            = find(TRIAL.time>TRIAL.time_state_iti(end), 1, 'first');
        TRIAL.ind_state_next_trial     = find(TRIAL.time>TRIAL.time_state_next_trial(end), 1, 'first');
        
        %% Extract Primary Sac
        clearvars -except counter_file SESSION_PARAMS ...
            TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
            TRIALS SACS_PRIM SACS_CORR data counter_trial ...
            TRIAL SAC_PRIM SAC_CORR
        
        trial_eye_velocity_trace        = TRIAL.eye_r_vm_filt;
        ind_search_begin_sac_prim       = TRIAL.ind_state_cue_present;
        ind_search_end_sac_prim         = TRIAL.ind_state_end_fixation;
        
        params_prim.MinPeakHeight       = 150.0; % deg/s
        params_prim.MinPeakProminence   = 100; % data points
        params_prim.rough_threshold     = 50.0; % deg/s
        params_prim.fine_threshold      = 20.0; % deg/s
        params_prim.sampling_freq       = 1000.0; % Hz
        params_prim.cutoff_freq         = 50.0; % Hz
        params_prim.window_half_length  = 4; % data points
        params_prim.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
        
        output_ = ESN_Sac_Finder(trial_eye_velocity_trace, ...
                                 ind_search_begin_sac_prim, ind_search_end_sac_prim, params_prim);
        
        validity_sac_prim   = output_.validity;
        inds_sac_prim       = output_.inds;
        ind_sac_prim_start  = output_.ind_start;
        ind_sac_prim_vmax   = output_.ind_vmax;
        ind_sac_prim_finish = output_.ind_finish;
        
        %% Save Primary Sac data to SAC_PRIM
        SAC_PRIM.validity                  = validity_sac_prim;
        SAC_PRIM.inds                      = inds_sac_prim;
        SAC_PRIM.ind_start                 = ind_sac_prim_start;
        SAC_PRIM.ind_vmax                  = ind_sac_prim_vmax;
        SAC_PRIM.ind_finish                = ind_sac_prim_finish;
        
        SAC_PRIM.time                      = TRIAL.time_1K( SAC_PRIM.inds);
        SAC_PRIM.eye_r_px                  = TRIAL.eye_r_px(SAC_PRIM.inds);
        SAC_PRIM.eye_r_py                  = TRIAL.eye_r_py(SAC_PRIM.inds);
        SAC_PRIM.eye_r_vx                  = TRIAL.eye_r_vx(SAC_PRIM.inds);
        SAC_PRIM.eye_r_vy                  = TRIAL.eye_r_vy(SAC_PRIM.inds);
        SAC_PRIM.eye_r_vm                  = TRIAL.eye_r_vm(SAC_PRIM.inds);
        SAC_PRIM.eye_r_vm_max              = TRIAL.eye_r_vm(SAC_PRIM.ind_vmax);
        SAC_PRIM.eye_r_px_centered         = SAC_PRIM.eye_r_px - TRIAL.start_x;
        SAC_PRIM.eye_r_py_centered         = SAC_PRIM.eye_r_py - TRIAL.start_y;
        
        SAC_PRIM.eye_r_px_start            = TRIAL.eye_r_px(SAC_PRIM.ind_start);
        SAC_PRIM.eye_r_px_finish           = TRIAL.eye_r_px(SAC_PRIM.ind_finish);
        SAC_PRIM.eye_r_py_start            = TRIAL.eye_r_py(SAC_PRIM.ind_start);
        SAC_PRIM.eye_r_py_finish           = TRIAL.eye_r_py(SAC_PRIM.ind_finish);
        SAC_PRIM.eye_r_px_start_centered   = SAC_PRIM.eye_r_px_start  - TRIAL.start_x;
        SAC_PRIM.eye_r_px_finish_centered  = SAC_PRIM.eye_r_px_finish - TRIAL.start_x;
        SAC_PRIM.eye_r_py_start_centered   = SAC_PRIM.eye_r_py_start  - TRIAL.start_y;
        SAC_PRIM.eye_r_py_finish_centered  = SAC_PRIM.eye_r_py_finish - TRIAL.start_y;
        SAC_PRIM.eye_r_amp_x               = (SAC_PRIM.eye_r_px_finish - SAC_PRIM.eye_r_px_start);
        SAC_PRIM.eye_r_amp_y               = (SAC_PRIM.eye_r_py_finish - SAC_PRIM.eye_r_py_start);
        SAC_PRIM.eye_r_amp_m               = (sqrt(SAC_PRIM.eye_r_amp_x^2+SAC_PRIM.eye_r_amp_y^2));
        SAC_PRIM.reaction                  = SAC_PRIM.ind_start - TRIAL.ind_state_cue_present;
        
        %% Extract Corrective Sac
        clearvars -except counter_file SESSION_PARAMS ...
            TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
            TRIALS SACS_PRIM SACS_CORR data counter_trial ...
            TRIAL SAC_PRIM SAC_CORR
        
        trial_eye_velocity_trace       = TRIAL.eye_r_vm_filt;
        ind_search_begin_sac_corr      = SAC_PRIM.ind_finish + 20;
        ind_search_end_sac_corr        = TRIAL.ind_state_iti;
        
        params_corr.MinPeakHeight      = 110.0; % deg/s
        params_corr.MinPeakProminence  = 80; % data points
        params_corr.rough_threshold    = 50.0; % deg/s
        params_corr.fine_threshold     = 20.0; % deg/s
        params_corr.sampling_freq      = 1000.0; % Hz
        params_corr.cutoff_freq        = 75.0; % Hz
        params_corr.window_half_length = 4; % data points
        params_corr.prominence_or_first = 'first'; % which peak to select, 'prominent' or 'first'
        
        output_ = ESN_Sac_Finder(trial_eye_velocity_trace, ...
                                 ind_search_begin_sac_corr, ind_search_end_sac_corr, params_corr);
        
        validity_sac_corr   = output_.validity;
        inds_sac_corr       = output_.inds;
        ind_sac_corr_start  = output_.ind_start;
        ind_sac_corr_vmax   = output_.ind_vmax;
        ind_sac_corr_finish = output_.ind_finish;
        
        %% Save Corrective Sac data to SAC_CORR
        SAC_CORR.validity                  = validity_sac_corr;
        SAC_CORR.inds                      = inds_sac_corr;
        SAC_CORR.ind_start                 = ind_sac_corr_start;
        SAC_CORR.ind_vmax                  = ind_sac_corr_vmax;
        SAC_CORR.ind_finish                = ind_sac_corr_finish;
        
        SAC_CORR.time                      = TRIAL.time_1K( SAC_CORR.inds);
        SAC_CORR.eye_r_px                  = TRIAL.eye_r_px(SAC_CORR.inds);
        SAC_CORR.eye_r_py                  = TRIAL.eye_r_py(SAC_CORR.inds);
        SAC_CORR.eye_r_vx                  = TRIAL.eye_r_vx(SAC_CORR.inds);
        SAC_CORR.eye_r_vy                  = TRIAL.eye_r_vy(SAC_CORR.inds);
        SAC_CORR.eye_r_vm                  = TRIAL.eye_r_vm(SAC_CORR.inds);
        SAC_CORR.eye_r_vm_max              = TRIAL.eye_r_vm(SAC_CORR.ind_vmax);
        SAC_CORR.eye_r_px_centered         = SAC_CORR.eye_r_px - TRIAL.start_x;
        SAC_CORR.eye_r_py_centered         = SAC_CORR.eye_r_py - TRIAL.start_y;
        
        SAC_CORR.eye_r_px_start            = TRIAL.eye_r_px(SAC_CORR.ind_start);
        SAC_CORR.eye_r_px_finish           = TRIAL.eye_r_px(SAC_CORR.ind_finish);
        SAC_CORR.eye_r_py_start            = TRIAL.eye_r_py(SAC_CORR.ind_start);
        SAC_CORR.eye_r_py_finish           = TRIAL.eye_r_py(SAC_CORR.ind_finish);
        SAC_CORR.eye_r_px_start_centered   = SAC_CORR.eye_r_px_start  - TRIAL.start_x;
        SAC_CORR.eye_r_px_finish_centered  = SAC_CORR.eye_r_px_finish - TRIAL.start_x;
        SAC_CORR.eye_r_py_start_centered   = SAC_CORR.eye_r_py_start  - TRIAL.start_y;
        SAC_CORR.eye_r_py_finish_centered  = SAC_CORR.eye_r_py_finish - TRIAL.start_y;
        SAC_CORR.eye_r_amp_x               = (SAC_CORR.eye_r_px_finish - SAC_CORR.eye_r_px_start);
        SAC_CORR.eye_r_amp_y               = (SAC_CORR.eye_r_py_finish - SAC_CORR.eye_r_py_start);
        SAC_CORR.eye_r_amp_m               = (sqrt(SAC_CORR.eye_r_amp_x^2+SAC_CORR.eye_r_amp_y^2));
        SAC_CORR.reaction                  = SAC_CORR.ind_start - SAC_PRIM.ind_finish;
        
        %% Compute the start position Bias
        inds_start_fixation = TRIAL.ind_state_cue_present - round(TRIAL.time_fixation * 1000) : 1 : TRIAL.ind_state_cue_present;
        state_start_eye_r_px_start_fixation = TRIAL.eye_l_px_filt(inds_start_fixation);
        state_start_eye_r_py_start_fixation = TRIAL.eye_l_py_filt(inds_start_fixation);
        tgt_start_x = TRIAL.start_x;
        tgt_start_y = TRIAL.start_y;
        start_x_bias = mean(state_start_eye_r_px_start_fixation) - tgt_start_x;
        start_y_bias = mean(state_start_eye_r_py_start_fixation) - tgt_start_y;
        TRIAL.start_x_bias = start_x_bias;
        TRIAL.start_y_bias = start_y_bias;
        
        %% Build TRIALS, SACS_PRIM, SACS_CORR
        TRIALS(counter_trial) = TRIAL;
        SACS_PRIM(counter_trial) = SAC_PRIM;
        SACS_CORR(counter_trial) = SAC_CORR;
        
        % print a dot every 20 trials
        if rem(counter_trial, 20) == 0
            fprintf('.');
        end
        
    end
    fprintf(' --> Completed. \n')
    
    %% Arrange 'TRIALS_DATA'
    clearvars -except counter_file SESSION_PARAMS ...
        TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
        TRIALS SACS_PRIM SACS_CORR ...
        TRIALS_DATA SACS_PRIM_DATA SACS_CORR_DATA
    fprintf([SESSION_PARAMS.filename{counter_file} ': Arranging TRIALS_DATA ...'])
    clearvars('TRIALS_DATA'); TRIALS_DATA = struct;
    field_names_TRIALS = fieldnames(TRIALS);
    for counter_fields = 1 : 1 : length(field_names_TRIALS)
        for counter_trials = 1 : 1 : length(TRIALS)
            variable_TRIALS_ = TRIALS(counter_trials).(field_names_TRIALS{counter_fields});
            % handling an error which the variable_TRIALS_ was []
            if isempty(variable_TRIALS_)
                variable_TRIALS_ = nan;
            end
            variable_TRIALS_ = variable_TRIALS_(:);
            if max(size(variable_TRIALS_)) > 1
                variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
            end
            % the field does not exist in TRIALS_DATA
            if ~isfield(TRIALS_DATA, field_names_TRIALS{counter_fields})
                TRIALS_DATA.(field_names_TRIALS{counter_fields}) = [];
            end
            variable_TRIALS_DATA_ = TRIALS_DATA.(field_names_TRIALS{counter_fields});
            % variable_TRIALS_ is cell array
            if iscell(variable_TRIALS_)
                % variable_TRIALS_DATA_ is cell array
                % variables are compatible (both are cell), add new data
                if iscell(variable_TRIALS_DATA_)
                    variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
                end
                % variable_TRIALS_DATA_ is matrix array
                % convert variable_TRIALS_DATA_ to cell, and add new data
                if isnumeric(variable_TRIALS_DATA_)
                    variable_TRIALS_DATA_ = num2cell(variable_TRIALS_DATA_);
                    variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
                end
            end
            % variable_TRIALS_ is matrix array
            if isnumeric(variable_TRIALS_)
                % variable_TRIALS_DATA_ is cell array
                % convert variable_TRIALS_ to cell, and add new data
                if iscell(variable_TRIALS_DATA_)
                    variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
                    variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
                end
                % variable_TRIALS_DATA_ is matrix array
                % variables are compatible (both are matrix), add new data
                if isnumeric(variable_TRIALS_DATA_)
                    variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
                end
            end
            TRIALS_DATA.(field_names_TRIALS{counter_fields}) = variable_TRIALS_DATA_;
        end
        fprintf('.');
    end
    fprintf(' --> Completed. \n')
    
    %% Arrange 'SACS_PRIM_DATA'
    clearvars -except counter_file SESSION_PARAMS ...
        TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
        TRIALS SACS_PRIM SACS_CORR ...
        TRIALS_DATA SACS_PRIM_DATA SACS_CORR_DATA
    fprintf([SESSION_PARAMS.filename{counter_file} ': Arranging SACS_PRIM_DATA ...'])
    clearvars('SACS_PRIM_DATA'); SACS_PRIM_DATA = struct;
    field_names_SACS_PRIM_DATA = fieldnames(SACS_PRIM);
    SACS_PRIM_DATA_cell = struct2cell(SACS_PRIM);
    for counter_fields = 1 : 1 : length(field_names_SACS_PRIM_DATA)
        SACS_PRIM_DATA_field_cell = SACS_PRIM_DATA_cell(counter_fields,:,:);
        SACS_PRIM_DATA_field_cell = reshape(SACS_PRIM_DATA_field_cell, 1, []);
        nz = max(cellfun(@numel,SACS_PRIM_DATA_field_cell));
        SACS_PRIM_DATA_field_mat = cell2mat(cellfun(@(x) vertcat(double(x(:)),NaN(nz-numel(x), 1)),SACS_PRIM_DATA_field_cell,'uni',false));
        SACS_PRIM_DATA.(field_names_SACS_PRIM_DATA{counter_fields}) = SACS_PRIM_DATA_field_mat;
    end
    SACS_PRIM_DATA.validity  = logical(SACS_PRIM_DATA.validity);
    fprintf(' --> Completed. \n')
    
    %% Arrange 'SACS_CORR_DATA'
    clearvars -except counter_file SESSION_PARAMS ...
        TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
        TRIALS SACS_PRIM SACS_CORR ...
        TRIALS_DATA SACS_PRIM_DATA SACS_CORR_DATA
    fprintf([SESSION_PARAMS.filename{counter_file} ': Arranging SACS_CORR_DATA ...'])
    clearvars('SACS_CORR_DATA'); SACS_CORR_DATA = struct;
    field_names_SACS_CORR_DATA = fieldnames(SACS_CORR);
    SACS_CORR_DATA_cell = struct2cell(SACS_CORR);
    for counter_fields = 1 : 1 : length(field_names_SACS_CORR_DATA)
        SACS_CORR_DATA_field_cell = SACS_CORR_DATA_cell(counter_fields,:,:);
        SACS_CORR_DATA_field_cell = reshape(SACS_CORR_DATA_field_cell, 1, []);
        nz = max(cellfun(@numel,SACS_CORR_DATA_field_cell));
        SACS_CORR_DATA_field_mat = cell2mat(cellfun(@(x) vertcat(double(x(:)),NaN(nz-numel(x), 1)),SACS_CORR_DATA_field_cell,'uni',false));
        SACS_CORR_DATA.(field_names_SACS_CORR_DATA{counter_fields}) = SACS_CORR_DATA_field_mat;
    end
    SACS_CORR_DATA.validity  = logical(SACS_CORR_DATA.validity);
    fprintf(' --> Completed. \n')
    
    %% Build TRIALS_DATA_ALL, SACS_PRIM_DATA_ALL, SACS_CORR_DATA_ALL
    TRIALS_DATA_ALL(counter_file) = TRIALS_DATA;
    SACS_PRIM_DATA_ALL(counter_file) = SACS_PRIM_DATA;
    SACS_CORR_DATA_ALL(counter_file) = SACS_CORR_DATA;
    fprintf('#######################################\n')
    
end

%% Arrange 'TRIALS_SESSION'
clearvars -except SESSION_PARAMS ...
    TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
    TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION
fprintf('Arranging TRIALS_SESSION ...')
clearvars('TRIALS_SESSION'); TRIALS_SESSION = struct;
field_names_TRIALS_DATA_ALL = fieldnames(TRIALS_DATA_ALL);
for counter_fields = 1 : 1 : length(field_names_TRIALS_DATA_ALL)
    for counter_files = 1 : 1 : length(TRIALS_DATA_ALL)
        variable_TRIALS_DATA_ALL_ = TRIALS_DATA_ALL(counter_files).(field_names_TRIALS_DATA_ALL{counter_fields});
        % the field does not exist in TRIALS_SESSION
        if ~isfield(TRIALS_SESSION, field_names_TRIALS_DATA_ALL{counter_fields})
            TRIALS_SESSION.(field_names_TRIALS_DATA_ALL{counter_fields}) = [];
        end
        variable_TRIALS_SESSION_ = TRIALS_SESSION.(field_names_TRIALS_DATA_ALL{counter_fields});
        variable_TRIALS_SESSION_ = horzcat(variable_TRIALS_SESSION_, variable_TRIALS_DATA_ALL_);
        TRIALS_SESSION.(field_names_TRIALS_DATA_ALL{counter_fields}) = variable_TRIALS_SESSION_;
    end
    fprintf('.');
end
fprintf(' --> Completed. \n')

%% Arrange 'SACS_PRIM_SESSION'
clearvars -except SESSION_PARAMS ...
    TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
    TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION
fprintf('Arranging SACS_PRIM_SESSION ...')
clearvars('SACS_PRIM_SESSION'); SACS_PRIM_SESSION = struct;
field_names_SACS_PRIM_DATA_ALL = fieldnames(SACS_PRIM_DATA_ALL);
for counter_fields = 1 : 1 : length(field_names_SACS_PRIM_DATA_ALL)
    for counter_files = 1 : 1 : length(SACS_PRIM_DATA_ALL)
        variable_SACS_PRIM_DATA_ALL_ = SACS_PRIM_DATA_ALL(counter_files).(field_names_SACS_PRIM_DATA_ALL{counter_fields});
        % the field does not exist in SACS_PRIM_SESSION
        if ~isfield(SACS_PRIM_SESSION, field_names_SACS_PRIM_DATA_ALL{counter_fields})
            SACS_PRIM_SESSION.(field_names_SACS_PRIM_DATA_ALL{counter_fields}) = [];
        end
        variable_SACS_PRIM_SESSION_ = SACS_PRIM_SESSION.(field_names_SACS_PRIM_DATA_ALL{counter_fields});
        variable_SACS_PRIM_SESSION_ = horzcat(variable_SACS_PRIM_SESSION_, variable_SACS_PRIM_DATA_ALL_);
        SACS_PRIM_SESSION.(field_names_SACS_PRIM_DATA_ALL{counter_fields}) = variable_SACS_PRIM_SESSION_;
    end
    fprintf('.');
end
SACS_PRIM_SESSION.validity  = logical(SACS_PRIM_SESSION.validity);
fprintf(' --> Completed. \n')

%% Arrange 'SACS_CORR_SESSION'
clearvars -except SESSION_PARAMS ...
    TRIALS_DATA_ALL SACS_PRIM_DATA_ALL SACS_CORR_DATA_ALL ...
    TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION
fprintf('Arranging SACS_CORR_SESSION ...')
clearvars('SACS_CORR_SESSION'); SACS_CORR_SESSION = struct;
field_names_SACS_CORR_DATA_ALL = fieldnames(SACS_CORR_DATA_ALL);
for counter_fields = 1 : 1 : length(field_names_SACS_CORR_DATA_ALL)
    for counter_files = 1 : 1 : length(SACS_CORR_DATA_ALL)
        variable_SACS_CORR_DATA_ALL_ = SACS_CORR_DATA_ALL(counter_files).(field_names_SACS_CORR_DATA_ALL{counter_fields});
        % the field does not exist in SACS_CORR_SESSION
        if ~isfield(SACS_CORR_SESSION, field_names_SACS_CORR_DATA_ALL{counter_fields})
            SACS_CORR_SESSION.(field_names_SACS_CORR_DATA_ALL{counter_fields}) = [];
        end
        variable_SACS_CORR_SESSION_ = SACS_CORR_SESSION.(field_names_SACS_CORR_DATA_ALL{counter_fields});
        variable_SACS_CORR_SESSION_ = horzcat(variable_SACS_CORR_SESSION_, variable_SACS_CORR_DATA_ALL_);
        SACS_CORR_SESSION.(field_names_SACS_CORR_DATA_ALL{counter_fields}) = variable_SACS_CORR_SESSION_;
    end
    fprintf('.');
end
SACS_CORR_SESSION.validity  = logical(SACS_CORR_SESSION.validity);
fprintf(' --> Completed. \n')

%% Finding different sections of the experiment
clearvars -except SESSION_PARAMS TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION
iss_m = sqrt(TRIALS_SESSION.iss_x.^2 + TRIALS_SESSION.iss_y.^2);
inds_no_perturb = iss_m == 0;
inds_perturb = iss_m ~= 0;

trials_num_baseline  = [];
trials_num_adapt_1   = [];
trials_num_washout_1 = [];
trials_num_adapt_2   = [];
trials_num_washout_2 = [];
trials_num_adapt_3   = [];
trials_num_washout_3 = [];

num_trials = length(iss_m);
inds_switch_perturbation = find(diff(inds_no_perturb));
num_switch_perturbation = length(inds_switch_perturbation);
if     (num_switch_perturbation == 0)
    % no swith happened, either all the trials are perturb or no_perturb
    if sum(inds_no_perturb) > sum(inds_perturb)
        % all the trials are baseline
        trials_num_baseline  = 1 : num_trials;
    else
        % all the trials are adapt_1
        trials_num_adapt_1   = 1 : num_trials;
    end
    
elseif (num_switch_perturbation == 1)
    % one switch happened, either from base to adapt_1 or
    % from adapt_1 to washout_1
    if sum( inds_no_perturb(1:inds_switch_perturbation(1)) ) > sum( inds_perturb(1:inds_switch_perturbation(1)) )
        % first no_perturb (baseline), then perturb (adapt_1)
        % first part is baseline
        trials_num_baseline  = ( 1 : inds_switch_perturbation(1) );
        % second (last) part is adapt_1
        trials_num_adapt_1   = ( inds_switch_perturbation(1)+1 : num_trials );
    else
        % no_perturb happened after perturb
        % first part is adapt_1
        trials_num_adapt_1   = ( 1 : inds_switch_perturbation(1) );
        % second (last) part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(1)+1 : num_trials );
    end
    
elseif (num_switch_perturbation == 2)
    % two switches happened
    if sum( inds_no_perturb(1:inds_switch_perturbation(1)) ) > sum( inds_perturb(1:inds_switch_perturbation(1)) )
        % first part is baseline
        trials_num_baseline  = ( 1 : inds_switch_perturbation(1) );
        % second part is adapt_1
        trials_num_adapt_1   = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third (last) part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(2)+1 : num_trials );
    else
        % first part is adapt_1
        trials_num_adapt_1   = ( 1 : inds_switch_perturbation(1) );
        % second part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third (last) part is adapt_2
        trials_num_adapt_2   = ( inds_switch_perturbation(2)+1 : num_trials );
    end
    
elseif (num_switch_perturbation == 3)
    % three switches happened
    if sum( inds_no_perturb(1:inds_switch_perturbation(1)) ) > sum( inds_perturb(1:inds_switch_perturbation(1)) )
        % first part is baseline
        trials_num_baseline  = ( 1 : inds_switch_perturbation(1) );
        % second part is adapt_1
        trials_num_adapt_1   = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(2)+1 : inds_switch_perturbation(3) );
        % forth (last) part is adapt_2
        trials_num_adapt_2   = ( inds_switch_perturbation(3)+1 : num_trials );
    else
        % first part is adapt_1
        trials_num_adapt_1   = ( 1 : inds_switch_perturbation(1) );
        % second part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third part is adapt_2
        trials_num_adapt_2   = ( inds_switch_perturbation(2)+1 : inds_switch_perturbation(3) );
        % forth (last) part is washout_2
        trials_num_washout_2 = ( inds_switch_perturbation(3)+1 : num_trials );
    end
elseif (num_switch_perturbation == 4)
    % four switches happened
    if sum( inds_no_perturb(1:inds_switch_perturbation(1)) ) > sum( inds_perturb(1:inds_switch_perturbation(1)) )
        % first part is baseline
        trials_num_baseline  = ( 1 : inds_switch_perturbation(1) );
        % second part is adapt_1
        trials_num_adapt_1   = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(2)+1 : inds_switch_perturbation(3) );
        % forth part is adapt_2
        trials_num_adapt_2   = ( inds_switch_perturbation(3)+1 : inds_switch_perturbation(4) );
        % fifth (last) part is washout_2
        trials_num_washout_2 = ( inds_switch_perturbation(4)+1 : num_trials );
    else
        % first part is adapt_1
        trials_num_adapt_1   = ( 1 : inds_switch_perturbation(1) );
        % second part is washout_1
        trials_num_washout_1 = ( inds_switch_perturbation(1)+1 : inds_switch_perturbation(2) );
        % third part is adapt_2
        trials_num_adapt_2   = ( inds_switch_perturbation(2)+1 : inds_switch_perturbation(3) );
        % forth part is washout_2
        trials_num_washout_2 = ( inds_switch_perturbation(3)+1 : inds_switch_perturbation(4) );
        % fifth (last) part is adapt_3
        trials_num_adapt_3   = ( inds_switch_perturbation(4)+1 : num_trials );
        fprintf('SECTION DETECTION: ADAPT_3 DETECTED, CONSIDER MODIFYING INPUT FILES.\n');
    end
    
else
    fprintf('SECTION DETECTION: MORE THAN 4 SWITCHES HAPPENED, CONSIDER MODIFYING INPUT FILES.\n');
end

TRIALS_SESSION.trials_num           = 1 : num_trials;
TRIALS_SESSION.trials_num_baseline  = trials_num_baseline;
TRIALS_SESSION.trials_num_adapt_1   = trials_num_adapt_1;
TRIALS_SESSION.trials_num_washout_1 = trials_num_washout_1;
TRIALS_SESSION.trials_num_adapt_2   = trials_num_adapt_2;
TRIALS_SESSION.trials_num_washout_2 = trials_num_washout_2;
TRIALS_SESSION.trials_num_adapt_3   = trials_num_adapt_3;
TRIALS_SESSION.trials_num_washout_3 = trials_num_washout_3;

%% LEARNING
clearvars -except SESSION_PARAMS TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION
bin_size_ = 10;
trial_num_edges = 1 : bin_size_ : 1101;
inds_validity = SACS_PRIM_SESSION.validity;
inds_r_dir = (TRIALS_SESSION.cue_x - TRIALS_SESSION.start_x) > 0;
inds_l_dir = (TRIALS_SESSION.cue_x - TRIALS_SESSION.start_x) < 0;
iss_x = TRIALS_SESSION.iss_x;
inds_iss_x_negative = iss_x < 0;
if sum(inds_iss_x_negative(inds_r_dir)) > sum(inds_iss_x_negative(inds_l_dir))
    inds_u_dir = inds_l_dir;
    inds_d_dir = inds_r_dir;
else
    inds_u_dir = inds_r_dir;
    inds_d_dir = inds_l_dir;
end
sac_prim_px   = SACS_PRIM_SESSION.eye_r_px_finish_centered;
sac_prim_px(~inds_validity) = nan;
sac_prim_px_ = align_trials(TRIALS_SESSION, sac_prim_px);
inds_r_dir_  = align_trials(TRIALS_SESSION, inds_r_dir) == 1;
inds_l_dir_  = align_trials(TRIALS_SESSION, inds_l_dir) == 1;
inds_u_dir_  = align_trials(TRIALS_SESSION, inds_u_dir) == 1;
inds_d_dir_  = align_trials(TRIALS_SESSION, inds_d_dir) == 1;
[sac_prim_px_r_dir_, trials_num_r_dir_] = directional_data(sac_prim_px_, inds_r_dir_);
[sac_prim_px_l_dir_, trials_num_l_dir_] = directional_data(sac_prim_px_, inds_l_dir_);
[sac_prim_px_u_dir_, trials_num_u_dir_] = directional_data(sac_prim_px_, inds_u_dir_);
[sac_prim_px_d_dir_, trials_num_d_dir_] = directional_data(sac_prim_px_, inds_d_dir_);
sac_prim_px_ = nanmean([sac_prim_px_r_dir_; sac_prim_px_l_dir_; sac_prim_px_u_dir_; sac_prim_px_d_dir_]);
trials_num_  = nanmean([trials_num_r_dir_;  trials_num_l_dir_;  trials_num_u_dir_;  trials_num_d_dir_ ]);
[trials_num_r_dir_binned_, sac_prim_px_r_dir_binned_, ~, ~] = ESN_Bin(trials_num_r_dir_, sac_prim_px_r_dir_, trial_num_edges);
[trials_num_l_dir_binned_, sac_prim_px_l_dir_binned_, ~, ~] = ESN_Bin(trials_num_l_dir_, sac_prim_px_l_dir_, trial_num_edges);
[trials_num_u_dir_binned_, sac_prim_px_u_dir_binned_, ~, ~] = ESN_Bin(trials_num_u_dir_, sac_prim_px_u_dir_, trial_num_edges);
[trials_num_d_dir_binned_, sac_prim_px_d_dir_binned_, ~, ~] = ESN_Bin(trials_num_d_dir_, sac_prim_px_d_dir_, trial_num_edges);

LEARNING.sac_prim_px_ = sac_prim_px_;
LEARNING.trials_num_  = trials_num_;
LEARNING.inds_r_dir_  = inds_r_dir_;
LEARNING.inds_l_dir_  = inds_l_dir_;
LEARNING.inds_u_dir_  = inds_u_dir_;
LEARNING.inds_d_dir_  = inds_d_dir_;
LEARNING.sac_prim_px_r_dir_ = sac_prim_px_r_dir_;
LEARNING.sac_prim_px_l_dir_ = sac_prim_px_l_dir_;
LEARNING.sac_prim_px_u_dir_ = sac_prim_px_u_dir_;
LEARNING.sac_prim_px_d_dir_ = sac_prim_px_d_dir_;
LEARNING.trials_num_r_dir_  = trials_num_r_dir_;
LEARNING.trials_num_l_dir_  = trials_num_l_dir_;
LEARNING.trials_num_u_dir_  = trials_num_u_dir_;
LEARNING.trials_num_d_dir_  = trials_num_d_dir_;
LEARNING.sac_prim_px_r_dir_binned_  = sac_prim_px_r_dir_binned_';
LEARNING.sac_prim_px_l_dir_binned_  = sac_prim_px_l_dir_binned_';
LEARNING.sac_prim_px_u_dir_binned_  = sac_prim_px_u_dir_binned_';
LEARNING.sac_prim_px_d_dir_binned_  = sac_prim_px_d_dir_binned_';
LEARNING.trials_num_r_dir_binned_   = trials_num_r_dir_binned_';
LEARNING.trials_num_l_dir_binned_   = trials_num_l_dir_binned_';
LEARNING.trials_num_u_dir_binned_   = trials_num_u_dir_binned_';
LEARNING.trials_num_d_dir_binned_   = trials_num_d_dir_binned_';

%% Build Reduced_sized_struct
rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
        'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
        'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
        'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
        'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};

Reduced_size_struct = struct;
Reduced_size_struct.TRIALS_SESSION    = TRIALS_SESSION;
Reduced_size_struct.SACS_PRIM_SESSION = SACS_PRIM_SESSION;
Reduced_size_struct.SACS_CORR_SESSION = SACS_CORR_SESSION;
Reduced_size_struct.LEARNING          = LEARNING;
Reduced_size_struct.SESSION_PARAMS    = SESSION_PARAMS;

Reduced_size_struct.TRIALS_SESSION = rmfield(Reduced_size_struct.TRIALS_SESSION,rmfields_list);

%% save: TRIALS_SESSION, SACS_PRIM_SESSION, SACS_CORR_SESSION, LEARNING to disk
clearvars -except SESSION_PARAMS TRIALS_SESSION SACS_PRIM_SESSION SACS_CORR_SESSION LEARNING Reduced_size_struct
filename = 'SESSION_DATA';
pathname = SESSION_PARAMS.pathname{end};
foldername = SESSION_PARAMS.foldername{end};
fprintf(['Saving ' filename '_' foldername ' file ...'])
save([pathname filename '_' foldername '.mat'], 'SESSION_PARAMS', 'TRIALS_SESSION', ...
    'SACS_PRIM_SESSION', 'SACS_CORR_SESSION', 'LEARNING', ...
    'Reduced_size_struct', '-v7.3');
fprintf(' --> Completed. \n')

%% Terminate the program if running in shell mode
if ~usejava('desktop')
    return;
end

%% plot-02: Reaction time histogram
h_fig_ = figure(2);
clf(h_fig_)
hold on
cmap_lines_ = lines(7);
hist_edges_reaction_time = 0:5:300;
ylim([0 0.15])
xlim([min(hist_edges_reaction_time) max(hist_edges_reaction_time)])

hist_X_PRIM = SACS_PRIM_SESSION.reaction(SACS_PRIM_SESSION.validity);
hist_X_CORR = SACS_CORR_SESSION.reaction(SACS_CORR_SESSION.validity);
[hist_N_PRIM,hist_edges_PRIM] = histcounts(hist_X_PRIM,hist_edges_reaction_time, 'Normalization', 'probability');
[hist_N_CORR,hist_edges_CORR] = histcounts(hist_X_CORR,hist_edges_reaction_time, 'Normalization', 'probability');
[~, ind_mode_PRIM] = max(hist_N_PRIM);
[~, ind_mode_CORR] = max(hist_N_CORR);
react_time_mode_PRIM = mean([hist_edges_PRIM(ind_mode_PRIM) hist_edges_PRIM(ind_mode_PRIM+1)]);
react_time_mode_CORR = mean([hist_edges_CORR(ind_mode_CORR) hist_edges_CORR(ind_mode_CORR+1)]);

plot(repmat(react_time_mode_CORR,1,2), ylim, '-', 'linewidth', 2, 'color', cmap_lines_(2,:));
plot(repmat(react_time_mode_PRIM,1,2), ylim, '-', 'linewidth', 2, 'color', cmap_lines_(1,:));
h_hist_CORR_ = histogram(hist_X_CORR, hist_edges_reaction_time, ...
    'Normalization', 'probability', 'FaceColor', cmap_lines_(2,:));
h_hist_PRIM_ = histogram(hist_X_PRIM, hist_edges_reaction_time, ...
    'Normalization', 'probability', 'FaceColor', cmap_lines_(1,:));

xlabel('Reaction Time (ms)')
ylabel('Probability')
legend([h_hist_PRIM_ h_hist_CORR_], {'Primary', 'Corrective'})

ESN_Beautify_Plot

%% save: Reaction time histogram
filename = 'fig2_reactiontime';
pathname = SESSION_PARAMS.pathname{end};
foldername = SESSION_PARAMS.foldername{end};
fprintf(['Saving ' filename '_' foldername ' file ...'])
saveas(h_fig_, [pathname filename '_' foldername], 'png');
saveas(h_fig_, [pathname filename '_' foldername], 'pdf');
fprintf(' --> Completed. \n')

%% plot-01: Learning curve
h_fig_ = figure(1);
num_trials = length(SACS_PRIM_SESSION.eye_r_py_finish_centered);
trials_num = TRIALS_SESSION.trials_num;
inds_prim_r_dir = (TRIALS_SESSION.cue_x - TRIALS_SESSION.start_x) > 0;
inds_prim_l_dir = (TRIALS_SESSION.cue_x - TRIALS_SESSION.start_x) < 0;
sac_prim_px_finish = SACS_PRIM_SESSION.eye_r_px_finish_centered;
sac_prim_px_finish_r_dir = sac_prim_px_finish(inds_prim_r_dir);
sac_prim_px_finish_l_dir = sac_prim_px_finish(inds_prim_l_dir);
inds_validity = SACS_PRIM_SESSION.validity;
inds_validity_r_dir = inds_validity(inds_prim_r_dir);
inds_validity_l_dir = inds_validity(inds_prim_l_dir);
trials_num_r_dir = trials_num(inds_prim_r_dir);
trials_num_l_dir = trials_num(inds_prim_l_dir);
sac_prim_px_finish_r_dir(~inds_validity_r_dir) = nan;
sac_prim_px_finish_l_dir(~inds_validity_l_dir) = nan;
sac_prim_px_finish_r_dir = ESN_Outlier(sac_prim_px_finish_r_dir, 25);
sac_prim_px_finish_l_dir = ESN_Outlier(sac_prim_px_finish_l_dir, 25);

iss_x = TRIALS_SESSION.iss_x;
inds_iss_x_negative = iss_x < 0;
if sum(inds_iss_x_negative(inds_prim_r_dir)) > sum(inds_iss_x_negative(inds_prim_l_dir))
    %sac_prim_y_r_dir = -1 .* sac_prim_y_r_dir;
    r_sign_ = '-';
    l_sign_ = '+';
else
    %sac_prim_y_l_dir = -1 .* sac_prim_y_l_dir;
    r_sign_ = '+';
    l_sign_ = '-';
end

inds_nan_r_dir = isnan(sac_prim_px_finish_r_dir);
inds_nan_l_dir = isnan(sac_prim_px_finish_l_dir);
sac_prim_px_finish_r_dir_ = sac_prim_px_finish_r_dir;
sac_prim_px_finish_l_dir_ = sac_prim_px_finish_l_dir;
sac_prim_px_finish_r_dir_(inds_nan_r_dir) = [];
sac_prim_px_finish_l_dir_(inds_nan_l_dir) = [];
trials_num_r_dir_ = trials_num_r_dir;
trials_num_l_dir_ = trials_num_l_dir;
trials_num_r_dir_(inds_nan_r_dir) = [];
trials_num_l_dir_(inds_nan_l_dir) = [];

trial_num_edges = 1 : 40 : num_trials+1;
[trials_num_r_dir_binned, sac_prim_px_finish_r_dir_binned, ~, ~] = ESN_Bin(trials_num_r_dir_, sac_prim_px_finish_r_dir_, trial_num_edges);
[trials_num_l_dir_binned, sac_prim_px_finish_l_dir_binned, ~, ~] = ESN_Bin(trials_num_l_dir_, sac_prim_px_finish_l_dir_, trial_num_edges);

mean_sac_prim_px_finish = nanmean([nanmean(sac_prim_px_finish_r_dir), nanmean(sac_prim_px_finish_l_dir)]);

clf(h_fig_)
hold on
plot(trials_num(TRIALS_SESSION.trials_num_baseline), mean_sac_prim_px_finish, '.k')
plot(trials_num(TRIALS_SESSION.trials_num_adapt_1),  mean_sac_prim_px_finish, '.r')
plot(trials_num(TRIALS_SESSION.trials_num_washout_1),mean_sac_prim_px_finish, '.b')
h_plot_r = plot(trials_num_r_dir_binned, sac_prim_px_finish_r_dir_binned, '.-');
h_plot_l = plot(trials_num_l_dir_binned, sac_prim_px_finish_l_dir_binned, '.-');
legend([h_plot_r h_plot_l], {['right (' r_sign_ ')'], ['left (' l_sign_ ')']}, 'Location', 'northoutside', 'Orientation', 'horizontal')
xlabel('Trial Number (#)')
ylabel('Prim Sac Vertical Position (deg)')
% ylim([-.05 .3])
ESN_Beautify_Plot

%% save: Learning Plot
filename = 'fig1_learning';
pathname = SESSION_PARAMS.pathname{end};
foldername = SESSION_PARAMS.foldername{end};
fprintf(['Saving ' filename '_' foldername ' file ...'])
saveas(h_fig_, [pathname filename '_' foldername], 'png');
saveas(h_fig_, [pathname filename '_' foldername], 'pdf');
fprintf(' --> Completed. \n')

%% plot-03: vertical velocity trajectory
% variables from the "learning_curve" plot are necessary for this plot
h_fig_ = figure(3);
clf(h_fig_)
inds_validity_r_dir = ~inds_nan_r_dir;
inds_validity_l_dir = ~inds_nan_l_dir;
sac_prim_vx = SACS_PRIM_SESSION.eye_r_vx;
sac_prim_vx_r_dir = sac_prim_vx(:,inds_prim_r_dir);
sac_prim_vx_l_dir = sac_prim_vx(:,inds_prim_l_dir);
sac_prim_vy = SACS_PRIM_SESSION.eye_r_vy;
sac_prim_vy_r_dir = sac_prim_vy(:,inds_prim_r_dir);
sac_prim_vy_l_dir = sac_prim_vy(:,inds_prim_l_dir);

trials_of_interest = TRIALS_SESSION.trials_num_baseline;
inds_trials_of_interest_r_dir = (trials_num_r_dir >= trials_of_interest(1)) & (trials_num_r_dir <= trials_of_interest(end));
inds_trials_of_interest_l_dir = (trials_num_l_dir >= trials_of_interest(1)) & (trials_num_l_dir <= trials_of_interest(end));

subplot(2,4,1)
hold on
plot(     sac_prim_vx_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir));
plot(mean(sac_prim_vx_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,2)
hold on
plot(     sac_prim_vx_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir));
plot(mean(sac_prim_vx_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,5)
hold on
plot(     sac_prim_vy_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir));
plot(mean(sac_prim_vy_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,6)
hold on
plot(     sac_prim_vy_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir));
plot(mean(sac_prim_vy_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir), 2), 'k', 'linewidth', 2);

trials_of_interest = TRIALS_SESSION.trials_num_washout_1;
inds_trials_of_interest_r_dir = (trials_num_r_dir >= trials_of_interest(1)) & (trials_num_r_dir <= trials_of_interest(end));
inds_trials_of_interest_l_dir = (trials_num_l_dir >= trials_of_interest(1)) & (trials_num_l_dir <= trials_of_interest(end));

subplot(2,4,3)
hold on
plot(     sac_prim_vx_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir));
plot(mean(sac_prim_vx_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,4)
hold on
plot(     sac_prim_vx_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir));
plot(mean(sac_prim_vx_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,7)
hold on
plot(     sac_prim_vy_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir));
plot(mean(sac_prim_vy_r_dir(:, inds_validity_r_dir&inds_trials_of_interest_r_dir), 2), 'k', 'linewidth', 2);
subplot(2,4,8)
hold on
plot(     sac_prim_vy_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir));
plot(mean(sac_prim_vy_l_dir(:, inds_validity_l_dir&inds_trials_of_interest_l_dir), 2), 'k', 'linewidth', 2);

end

function output_ = align_trials(TRIALS_SESSION, input_)
% this function will align trials so that
% trials 001-100:  baseline  (100 trials)
% trials 101-800:  adapt_1   (700 trials)
% trials 801-1100: washout_1 (300 trials)
trials_num           = TRIALS_SESSION.trials_num;
trials_num_baseline  = trials_num(TRIALS_SESSION.trials_num_baseline);
trials_num_adapt_1   = trials_num(TRIALS_SESSION.trials_num_adapt_1);
trials_num_washout_1 = trials_num(TRIALS_SESSION.trials_num_washout_1);
input_baseline       = input_(trials_num_baseline);
input_adapt_1        = input_(trials_num_adapt_1);
input_washout_1      = input_(trials_num_washout_1);
output_              = nan(1, 1100);
output_(1+000:min([length(input_baseline),  100])+000) = input_baseline( 1:min([length(input_baseline),  100]));
output_(1+100:min([length(input_adapt_1),   700])+100) = input_adapt_1(  1:min([length(input_adapt_1),   700]));
output_(1+800:min([length(input_washout_1), 300])+800) = input_washout_1(1:min([length(input_washout_1), 300]));
end

function [output_, trials_] = directional_data(input_, inds_direction_)
% input data should be in aligned format so trials 101-200 be the first 100 trials of the adaptation
% phase
output_ = input_;
output_(~inds_direction_) = nan;
output_ = ESN_Outlier(output_, 25);
bias_ = nanmean(output_(101:200));
output_ = output_ - bias_;
trials_ = 1:1:1100; 
trials_(isnan(output_)) = nan;
end


