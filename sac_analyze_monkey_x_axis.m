%% Commits
% 20180620: change the variable names to make the code compatible with the
% new version of the experiment
% 20180801: start changing the code to analyze saccade_adaptation, cross_axis_adaptation experiment
% 20180805: for arranging the TRIALS_DATA started using cell arrays instead of putting nans and
%           making the sizes the same. This approach resulted in saving space during saving to disc.
% 20180808: adding new cell "Finding different sections of the experiment"
% 20180809: fix an issue for sac_prim detection, instead of
%           length(TRIAL.inds_trial) use length(TRIAL.time_1K)
% 20180810: make the code compatible with the previous versions to analyze
%           random corrective task

%% CLEAR
clc; clear;close all;

%% MAIN FILE LOOP
% get the list of MAT files to analyze
[filenames, pathnames] = uigetfile('*.mat', 'Pick list of MAT files', 'MultiSelect', 'on');
SESSION_PARAMS.filenames = sort(filenames);
SESSION_PARAMS.pathnames = pathnames;
for counter_file = 1 : 1 : length(SESSION_PARAMS.filenames)
    fprintf('#######################################\n')
    %% Load Data
    fprintf('Loading ...\n')
    clearvars -except counter_file SESSION_PARAMS TRIALS_DATA_ALL SACS_PRIM_DATA_ALL;
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
    clearvars -except counter_file SESSION_PARAMS TRIALS_DATA_ALL SACS_PRIM_DATA_ALL data;
    num_trials = length(data.trials);
    fprintf([SESSION_PARAMS.filename{counter_file} ': Analyzing TRIALS ...'])
    for counter_trial = 1 : 1 : num_trials-1
        %% Extract Trial Varibales
        clearvars -except counter_file SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
            data TRIALS SACS_PRIM counter_trial;
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
        TRIAL.inds_invalid = false(1, length(TRIAL.inds_trial));
        TRIAL.time     = double(data.t(1, TRIAL.inds_trial));                                      TRIAL.inds_invalid = isnan(TRIAL.time)     | TRIAL.inds_invalid;
        TRIAL.eye_r_px = double(data.right_horizontal_eye(1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_px) | TRIAL.inds_invalid;
        TRIAL.eye_r_py = double(data.right_vertical_eye(  1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_py) | TRIAL.inds_invalid;
        TRIAL.eye_l_px = double(data.left_horizontal_eye( 1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_px) | TRIAL.inds_invalid;
        TRIAL.eye_l_py = double(data.left_vertical_eye(   1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_py) | TRIAL.inds_invalid;
        TRIAL.eye_r_vx = double(data.right_horizontal_eye_velocity_filtered(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vx) | TRIAL.inds_invalid;
        TRIAL.eye_r_vy = double(data.right_vertical_eye_velocity_filtered(  1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vy) | TRIAL.inds_invalid;
        TRIAL.eye_l_vx = double(data.left_horizontal_eye_velocity_filtered( 1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vx) | TRIAL.inds_invalid;
        TRIAL.eye_l_vy = double(data.left_vertical_eye_velocity_filtered(   1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vy) | TRIAL.inds_invalid;
        TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
        TRIAL.eye_l_vm = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
        TRIAL.tgt_px   = double(data.target_x(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.tgt_px) | TRIAL.inds_invalid;
        TRIAL.tgt_py   = double(data.target_y(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.tgt_py) | TRIAL.inds_invalid;
        TRIAL.target_visible = logical(double(data.target_visible(1, TRIAL.inds_trial)));
        TRIAL.reward         = double(data.reward(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.reward) | TRIAL.inds_invalid;
        % remove invalid values
        TRIAL.time(    TRIAL.inds_invalid) = [];
        TRIAL.eye_r_px(TRIAL.inds_invalid) = [];
        TRIAL.eye_r_py(TRIAL.inds_invalid) = [];
        TRIAL.eye_l_px(TRIAL.inds_invalid) = [];
        TRIAL.eye_l_py(TRIAL.inds_invalid) = [];
        % reconstruct eye_r data
        TRIAL.time_1K  = TRIAL.time(1) : 0.001 : TRIAL.time(end);
        TRIAL.eye_r_px = interp1(TRIAL.time, TRIAL.eye_r_px, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_r_py = interp1(TRIAL.time, TRIAL.eye_r_py, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_r_vx = diff(TRIAL.eye_r_px)./diff(TRIAL.time_1K); TRIAL.eye_r_vx=[TRIAL.eye_r_vx(1) TRIAL.eye_r_vx];
        TRIAL.eye_r_vy = diff(TRIAL.eye_r_py)./diff(TRIAL.time_1K); TRIAL.eye_r_vy=[TRIAL.eye_r_vy(1) TRIAL.eye_r_vy];
        TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
        % reconstruct eye_l data
        TRIAL.eye_l_px = interp1(TRIAL.time, TRIAL.eye_l_px, TRIAL.time_1K, 'linear', 'extrap');
        TRIAL.eye_l_py = interp1(TRIAL.time, TRIAL.eye_l_py, TRIAL.time_1K, 'linear', 'extrap');
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
        clearvars -except counter_file SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
            data TRIALS SACS_PRIM counter_trial TRIAL
        % filter params
        sampling_freq = 1000.0;
        cutoff_freq = 100.0;
        [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
        % search slot for primary saccade
        sac_inds_search_slot = TRIAL.ind_state_cue_present : 1 : TRIAL.ind_state_iti;
        % extract primary sac
        sac_analyze_flag  = true;
        sac_validity      = true;
        sac_vm            = TRIAL.eye_r_vm_filt(1, sac_inds_search_slot);
        sac_vm_filt_heavy = filtfilt(b_butter,a_butter,sac_vm);
        [~, ind_sac_vmax] = findpeaks(sac_vm_filt_heavy, 'MinPeakProminence',75, 'MinPeakHeight', 100);
        % peaks happen very close to each other
        if(sum(diff(ind_sac_vmax)<80))
            sac_validity         = false;
        end
        [sac_vmax, ind_sac_vmax] = findpeaks(sac_vm_filt_heavy, 'MinPeakProminence',100,'SortStr','descend', 'NPeaks', 1, 'MinPeakHeight', 150);
        
        if(isempty(sac_vmax))
            sac_validity         = false;
        end
        
        if(~sac_validity)
            ind_sac_vmax         = round((TRIAL.ind_state_sac_onset+TRIAL.ind_state_reward)/2);
            inds_sac             = (ind_sac_vmax -61) : 1 : min([(ind_sac_vmax -61+149), length(TRIAL.time_1K)]);
            ind_sac_start        = ind_sac_vmax - 50;
            ind_sac_finish       = ind_sac_vmax + 50;
            sac_vmax             = TRIAL.eye_r_vm_filt(ind_sac_vmax);
            sac_analyze_flag     = false;
        end
        
        % find the ind_sac_vmax in original time series (non heavily filtered)
        % re-define the search slot
        if(sac_analyze_flag)
            if ((ind_sac_vmax+10) > (length(sac_vm))) || ((ind_sac_vmax-10) < 1)
                ind_sac_vmax_ = 9; % handling an error in which the max happened in the end
            else
                [~, ind_sac_vmax_] = max(sac_vm(ind_sac_vmax-10:ind_sac_vmax+10));
            end
            ind_sac_vmax_   = ind_sac_vmax_ + ((ind_sac_vmax-10) - 1);
            ind_sac_vmax    = ind_sac_vmax_ + (sac_inds_search_slot(1) - 1);
            % put the v_max at index 100, narrow the search to 300ms
            sac_inds_search_slot = (ind_sac_vmax - 100) : 1 : min([(ind_sac_vmax - 100 + 299), length(TRIAL.time_1K)]);
            sac_vm               = TRIAL.eye_r_vm_filt(1, sac_inds_search_slot);
            sac_vm_filt_heavy    = filtfilt(b_butter,a_butter,sac_vm);
        end
        
        % find sac begining based on 50deg/s threshold
        sac_begining_flag = true;
        if(sac_analyze_flag)
            % find the index that sac_vm_filt_heavy is more than 50deg/s
            ind_sac_start_50        = find(sac_vm_filt_heavy(1:100)<50, 1, 'last');
            % if the data is too noisy set the ind_sac_start to 50 and make the saccade invalid
            if(isempty(ind_sac_start_50))
                ind_sac_start     = 50;
                sac_validity      = false;
                sac_begining_flag = false;
            end
        end
        
        % find sac begining based on 20deg/s threshold
        ind_sac_start_20 = [];
        if (sac_begining_flag && sac_analyze_flag)
            ind_sac_start_20    = find(sac_vm_filt_heavy(1:ind_sac_start_50)<20, 1, 'last');
        end
        if(isempty(ind_sac_start_20) && sac_begining_flag && sac_analyze_flag)
            % if the sac start is between 20deg/s and 50deg/s then
            % find the ind_sac_start from original data based on 50deg/s threshold
            window_half_length = 10;
            ind_sac_start_       = find(sac_vm( max([(ind_sac_start_50-window_half_length), 1]) : min([(ind_sac_start_50+window_half_length), length(sac_vm)]) )<50, 1, 'last') - 1 + (max([(ind_sac_start_50-window_half_length), 1]) );
            if(isempty(ind_sac_start_))
                % if for whatever reason the original data is noisy but heavily filtered data is OK
                % use the index from heavily filtered data
                ind_sac_start    = ind_sac_start_50;
            else
                ind_sac_start    = ind_sac_start_;
            end
        end
        if((~isempty(ind_sac_start_20)) && sac_begining_flag && sac_analyze_flag)
            % find the ind_sac_start from original data based on 20deg/s threshold
            window_half_length = 10;
            ind_sac_start_       = find(sac_vm(max([(ind_sac_start_20-window_half_length), 1]) : min([(ind_sac_start_20+window_half_length), length(sac_vm)]) )<20, 1, 'last') - 1 + (max([(ind_sac_start_20-window_half_length), 1]) );
            if(isempty(ind_sac_start_))
                % if for whatever reason the original data is noisy but heavily filtered data is OK
                % use the index from heavily filtered data
                ind_sac_start    = ind_sac_start_20;
            else
                ind_sac_start    = ind_sac_start_;
            end
        end
        
        % find sac ending based on 50deg/s threshold
        sac_ending_flag = true;
        if(sac_analyze_flag)
            % find the index that sac_vm_filt_heavy is less than 50deg/s
            ind_sac_finish_50       = find(sac_vm_filt_heavy(100:end)<50, 1, 'first') - 1 + 100;
            if(isempty(ind_sac_finish_50))
                % if the data is too noisy set the ind_sac_finish to 150 and make the saccade invalid
                ind_sac_finish      = 150;
                sac_validity        = false;
                sac_ending_flag     = false;
            end
        end
        
        % find sac ending based on 20deg/s threshold
        ind_sac_finish_20 = [];
        if ( sac_ending_flag && sac_analyze_flag)
            ind_sac_finish_20   = find(sac_vm_filt_heavy( 100 : min([(ind_sac_finish_50+50), length(sac_vm_filt_heavy)]) )<20, 1, 'first') - 1 + 100;
        end
        if(isempty(ind_sac_finish_20) && sac_ending_flag && sac_analyze_flag)
            % if the sac end is between 20deg/s and 50deg/s then
            % find the ind_sac_finish from original data based on 50deg/s threshold
            window_half_length = 10;
            ind_sac_finish_       = find(sac_vm(max([(ind_sac_finish_50-window_half_length), 1]) : min([(ind_sac_finish_50+window_half_length), length(sac_vm)]) )<50, 1, 'first') - 1 + (max([(ind_sac_finish_50-window_half_length), 1]) );
            if(isempty(ind_sac_finish_))
                % if for whatever reason the original data is noisy but heavily filtered data is OK
                % use the index from heavily filtered data
                ind_sac_finish    = ind_sac_finish_50;
            else
                ind_sac_finish    = ind_sac_finish_;
            end
        end
        if((~isempty(ind_sac_finish_20)) && sac_ending_flag && sac_analyze_flag)
            % find the ind_sac_finish from original data based on 20deg/s threshold
            window_half_length = 10;
            ind_sac_finish_       = find(sac_vm(max([(ind_sac_finish_20-window_half_length), 1]) : min([(ind_sac_finish_20+window_half_length), length(sac_vm)]) )<20, 1, 'first') - 1 + (max([(ind_sac_finish_20-window_half_length), 1]));
            if(isempty(ind_sac_finish_))
                % if for whatever reason the original data is noisy but heavily filtered data is OK
                % use the index from heavily filtered data
                ind_sac_finish    = ind_sac_finish_20;
            else
                ind_sac_finish    = ind_sac_finish_;
            end
        end
        
        if(sac_analyze_flag)
            ind_sac_start    = ind_sac_vmax - 1 - 100 + ind_sac_start;
            ind_sac_finish   = ind_sac_vmax - 1 - 100 + ind_sac_finish;
            inds_sac         = (ind_sac_vmax -61) : 1 : min([(ind_sac_vmax -61+149), length(TRIAL.time_1K)]);
        end
        
        SAC_PRIM.validity    = sac_validity;
        SAC_PRIM.inds        = inds_sac;
        SAC_PRIM.ind_start   = ind_sac_start;
        SAC_PRIM.ind_vmax    = ind_sac_vmax;
        SAC_PRIM.ind_finish  = ind_sac_finish;
        
        SAC_PRIM.time                  = TRIAL.time_1K( inds_sac);
        SAC_PRIM.eye_r_px              = TRIAL.eye_r_px(inds_sac);
        SAC_PRIM.eye_r_py              = TRIAL.eye_r_py(inds_sac);
        SAC_PRIM.eye_r_vx              = TRIAL.eye_r_vx(inds_sac);
        SAC_PRIM.eye_r_vy              = TRIAL.eye_r_vy(inds_sac);
        SAC_PRIM.eye_r_vm              = TRIAL.eye_r_vm(inds_sac);
        SAC_PRIM.eye_r_vm_max          = TRIAL.eye_r_vm(ind_sac_vmax);
        SAC_PRIM.eye_r_px_centered     = SAC_PRIM.eye_r_px - TRIAL.start_x;
        SAC_PRIM.eye_r_py_centered     = SAC_PRIM.eye_r_py - TRIAL.start_y;
        
        SAC_PRIM.eye_r_px_start            = TRIAL.eye_r_px(ind_sac_start);
        SAC_PRIM.eye_r_px_finish           = TRIAL.eye_r_px(ind_sac_finish);
        SAC_PRIM.eye_r_py_start            = TRIAL.eye_r_py(ind_sac_start);
        SAC_PRIM.eye_r_py_finish           = TRIAL.eye_r_py(ind_sac_finish);
        SAC_PRIM.eye_r_px_start_centered   = SAC_PRIM.eye_r_px_start  - TRIAL.start_x;
        SAC_PRIM.eye_r_px_finish_centered  = SAC_PRIM.eye_r_px_finish - TRIAL.start_x;
        SAC_PRIM.eye_r_py_start_centered   = SAC_PRIM.eye_r_py_start  - TRIAL.start_y;
        SAC_PRIM.eye_r_py_finish_centered  = SAC_PRIM.eye_r_py_finish - TRIAL.start_y;
        SAC_PRIM.eye_r_amp_x               = (SAC_PRIM.eye_r_px_finish - SAC_PRIM.eye_r_px_start);
        SAC_PRIM.eye_r_amp_y               = (SAC_PRIM.eye_r_py_finish - SAC_PRIM.eye_r_py_start);
        SAC_PRIM.eye_r_amp_m               = (sqrt(SAC_PRIM.eye_r_amp_x^2+SAC_PRIM.eye_r_amp_y^2));
        SAC_PRIM.reaction                  = SAC_PRIM.ind_start - TRIAL.ind_state_cue_present;
        
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
        
        %% Build TRIAL and SAC_PRIM
        TRIALS(counter_trial) = TRIAL;
        SACS_PRIM(counter_trial) = SAC_PRIM;
        
    end
    fprintf(' --> Completed. \n')
    
    %% Arrange 'TRIALS_DATA'
    clearvars -except counter_file SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
        TRIALS SACS_PRIM TRIALS_DATA SACS_PRIM_DATA
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
    clearvars -except counter_file SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
        TRIALS SACS_PRIM TRIALS_DATA SACS_PRIM_DATA
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
    
    %% Build TRIALS_DATA_ALL and SACS_PRIM_DATA_ALL
    TRIALS_DATA_ALL(counter_file) = TRIALS_DATA;
    SACS_PRIM_DATA_ALL(counter_file) = SACS_PRIM_DATA;
    fprintf('#######################################\n')

    
end

%% Arrange 'TRIALS_SESSION'
clearvars -except SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
    TRIALS_SESSION SACS_PRIM_SESSION
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
clearvars -except SESSION_PARAMS  TRIALS_DATA_ALL SACS_PRIM_DATA_ALL ...
    TRIALS_SESSION SACS_PRIM_SESSION
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

%% Finding different sections of the experiment
clearvars -except SESSION_PARAMS TRIALS_SESSION SACS_PRIM_SESSION
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

%% save: TRIALS_SESSION and SACS_PRIM_SESSION to disk
clearvars -except SESSION_PARAMS TRIALS_SESSION SACS_PRIM_SESSION
filename = 'SESSION_DATA';
pathname = SESSION_PARAMS.pathname{end};
foldername = SESSION_PARAMS.foldername{end};
fprintf(['Saving ' filename '_' foldername ' file ...'])
save([pathname filename '_' foldername '.mat'], 'SESSION_PARAMS', 'TRIALS_SESSION', 'SACS_PRIM_SESSION', '-v7.3');
fprintf(' --> Completed. \n')

%% plot-02: Reaction time histogram
h_fig_ = figure(2);
clf(h_fig_)
trials_reaction_time_edges = 0:10:400;
histogram(SACS_PRIM_SESSION.reaction, trials_reaction_time_edges)
xlabel('Reaction Time (ms)')
ylabel('Frequency')
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
sac_prim_py_finish = SACS_PRIM_SESSION.eye_r_py_finish_centered;
sac_prim_py_finish_r_dir = sac_prim_py_finish(inds_prim_r_dir);
sac_prim_py_finish_l_dir = sac_prim_py_finish(inds_prim_l_dir);
inds_validity = SACS_PRIM_SESSION.validity;
inds_validity_r_dir = inds_validity(inds_prim_r_dir);
inds_validity_l_dir = inds_validity(inds_prim_l_dir);
trials_num_r_dir = trials_num(inds_prim_r_dir);
trials_num_l_dir = trials_num(inds_prim_l_dir);
sac_prim_py_finish_r_dir(~inds_validity_r_dir) = nan;
sac_prim_py_finish_l_dir(~inds_validity_l_dir) = nan;
sac_prim_py_finish_r_dir = ESN_Outlier(sac_prim_py_finish_r_dir, 25);
sac_prim_py_finish_l_dir = ESN_Outlier(sac_prim_py_finish_l_dir, 25);

iss_y = TRIALS_SESSION.iss_y;
inds_iss_y_negative = iss_y < 0;
if sum(inds_iss_y_negative(inds_prim_r_dir)) > sum(inds_iss_y_negative(inds_prim_l_dir))
    %sac_prim_y_r_dir = -1 .* sac_prim_y_r_dir;
    r_sign_ = '-';
    l_sign_ = '+';
else
    %sac_prim_y_l_dir = -1 .* sac_prim_y_l_dir;
    r_sign_ = '+';
    l_sign_ = '-';
end

inds_nan_r_dir = isnan(sac_prim_py_finish_r_dir);
inds_nan_l_dir = isnan(sac_prim_py_finish_l_dir);
sac_prim_py_finish_r_dir_ = sac_prim_py_finish_r_dir;
sac_prim_py_finish_l_dir_ = sac_prim_py_finish_l_dir;
sac_prim_py_finish_r_dir_(inds_nan_r_dir) = [];
sac_prim_py_finish_l_dir_(inds_nan_l_dir) = [];
trials_num_r_dir_ = trials_num_r_dir;
trials_num_l_dir_ = trials_num_l_dir;
trials_num_r_dir_(inds_nan_r_dir) = [];
trials_num_l_dir_(inds_nan_l_dir) = [];

trial_num_edges = 1 : 40 : num_trials+1;
[trials_num_r_dir_binned, sac_prim_py_finish_r_dir_binned, ~, ~] = ESN_BINNING(trials_num_r_dir_, sac_prim_py_finish_r_dir_, trial_num_edges);
[trials_num_l_dir_binned, sac_prim_py_finish_l_dir_binned, ~, ~] = ESN_BINNING(trials_num_l_dir_, sac_prim_py_finish_l_dir_, trial_num_edges);

mean_sac_prim_py_finish = nanmean([nanmean(sac_prim_py_finish_r_dir), nanmean(sac_prim_py_finish_l_dir)]);

clf(h_fig_)
hold on
plot(trials_num(TRIALS_SESSION.trials_num_baseline), mean_sac_prim_py_finish, '.k')
plot(trials_num(TRIALS_SESSION.trials_num_adapt_1),  mean_sac_prim_py_finish, '.r')
plot(trials_num(TRIALS_SESSION.trials_num_washout_1),mean_sac_prim_py_finish, '.b')
h_plot_r = plot(trials_num_r_dir_binned, sac_prim_py_finish_r_dir_binned, '.-');
h_plot_l = plot(trials_num_l_dir_binned, sac_prim_py_finish_l_dir_binned, '.-');
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


