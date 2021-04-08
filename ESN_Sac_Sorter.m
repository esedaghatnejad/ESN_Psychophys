function [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = ESN_Sac_Sorter(TRIALS_DATA, EXPERIMENT_PARAMS, flag_session_figure)
%% Handle inputs
if nargin < 3
    flag_session_figure = true;
end

if nargin < 2
    EXPERIMENT_PARAMS.file_name = '20202020_2020';
end

if nargin < 1
    SACS_ALL_DATA = [];
    TRIALS_DATA = [];
    EXPERIMENT_PARAMS = [];
    return;
end

EXPERIMENT_PARAMS.flag_session_figure = flag_session_figure;

%% Define saccade tags
EXPERIMENT_PARAMS.sac_tag_list = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    ... % Add more tags here, do not reorder or change the tags defined above.
    };

%% Analyze trials
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA;
num_trials = length(TRIALS_DATA.time_1K);
fprintf([EXPERIMENT_PARAMS.file_name ': Analyzing TRIALS ...'])
for counter_trial = 1 : 1 : num_trials
    %% Clear variables
    clearvars -except EXPERIMENT_PARAMS TRIALS_DATA counter_trial num_trials SACS_ALL
    
    %% Extract Trial Varibales
    TRIAL.start_x       = TRIALS_DATA.start_x(1,counter_trial);
    TRIAL.start_y       = TRIALS_DATA.start_y(1,counter_trial);
    TRIAL.cue_x         = TRIALS_DATA.cue_x(  1,counter_trial);
    TRIAL.cue_y         = TRIALS_DATA.cue_y(  1,counter_trial);
    TRIAL.end_x         = TRIALS_DATA.end_x(  1,counter_trial);
    TRIAL.end_y         = TRIALS_DATA.end_y(  1,counter_trial);
    TRIAL.iss_x         = TRIALS_DATA.iss_x(  1,counter_trial);
    TRIAL.iss_y         = TRIALS_DATA.iss_y(  1,counter_trial);
    
    TRIAL.time_start                 = TRIALS_DATA.time_start(               1,counter_trial);
    TRIAL.time_end                   = TRIALS_DATA.time_end(                 1,counter_trial);
    variable_list = {'time_state_str_pursuit','time_state_str_present','time_state_str_fixation','time_state_cue_present',...
        'time_state_sac_detect_on','time_state_sac_onset','time_state_sac_detect_off',...
        'time_state_reward','time_state_end_fixation','time_state_iti','time_state_next_trial'};
    for counter_variable = 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        if isa(TRIALS_DATA.(variable_name),'double')
            TRIAL.(variable_name)     = TRIALS_DATA.(variable_name)(   1,counter_trial);
        elseif isa(TRIALS_DATA.(variable_name),'cell')
            TRIAL.(variable_name)     = TRIALS_DATA.(variable_name){   1,counter_trial};
        end
    end
%     TRIAL.time_state_str_pursuit     = TRIALS_DATA.time_state_str_pursuit{   1,counter_trial};
%     TRIAL.time_state_str_present     = TRIALS_DATA.time_state_str_present{   1,counter_trial};
%     TRIAL.time_state_str_fixation    = TRIALS_DATA.time_state_str_fixation{  1,counter_trial};
%     TRIAL.time_state_cue_present     = TRIALS_DATA.time_state_cue_present{   1,counter_trial};
%     TRIAL.time_state_sac_detect_on   = TRIALS_DATA.time_state_sac_detect_on{ 1,counter_trial};
%     TRIAL.time_state_sac_onset       = TRIALS_DATA.time_state_sac_onset{     1,counter_trial};
%     TRIAL.time_state_sac_detect_off  = TRIALS_DATA.time_state_sac_detect_off{1,counter_trial};
%     TRIAL.time_state_reward          = TRIALS_DATA.time_state_reward{        1,counter_trial};
%     TRIAL.time_state_end_fixation    = TRIALS_DATA.time_state_end_fixation(  1,counter_trial);
%     TRIAL.time_state_iti             = TRIALS_DATA.time_state_iti(           1,counter_trial);
%     TRIAL.time_state_next_trial      = TRIALS_DATA.time_state_next_trial(    1,counter_trial);
    TRIAL.time_iti                   = TRIALS_DATA.time_iti(                 1,counter_trial);
    TRIAL.time_punishment            = TRIALS_DATA.time_punishment(          1,counter_trial);
    TRIAL.time_fixation              = TRIALS_DATA.time_fixation(            1,counter_trial);
    TRIAL.time_pursuit               = TRIALS_DATA.time_pursuit(             1,counter_trial);
    TRIAL.time_pursuit(isnan(TRIAL.time_pursuit)) = 0.200;
    
    TRIAL.time_1K       = TRIALS_DATA.time_1K{      1,counter_trial};
    TRIAL.tgt_px        = TRIALS_DATA.tgt_px{       1,counter_trial};
    TRIAL.tgt_py        = TRIALS_DATA.tgt_py{       1,counter_trial};
    TRIAL.reward        = TRIALS_DATA.reward{       1,counter_trial};
    TRIAL.target_visible = TRIALS_DATA.target_visible{1,counter_trial};
    TRIAL.eye_r_px_filt = TRIALS_DATA.eye_r_px_filt{1,counter_trial};
    TRIAL.eye_r_py_filt = TRIALS_DATA.eye_r_py_filt{1,counter_trial};
    TRIAL.eye_r_vx_filt = TRIALS_DATA.eye_r_vx_filt{1,counter_trial};
    TRIAL.eye_r_vy_filt = TRIALS_DATA.eye_r_vy_filt{1,counter_trial};
    TRIAL.eye_r_vm_filt = TRIALS_DATA.eye_r_vm_filt{1,counter_trial};
    
    %% Extract Saccades
    clearvars -except EXPERIMENT_PARAMS TRIALS_DATA counter_trial num_trials TRIAL SACS_ALL
    
    threshold = 75; % deg/s
    eye_vm_ = TRIAL.eye_r_vm_filt;
    eye_vm_ = [eye_vm_; eye_vm_(end)];
    idx_bool_rise_threshold = (eye_vm_(2:end)-threshold > 0) & (eye_vm_(1:end-1)-threshold <= 0);
    idx_int_rise_threshold = find(idx_bool_rise_threshold);
    num_saccades = length(idx_int_rise_threshold);
    eye_velocity_trace     = TRIAL.eye_r_vm_filt;
    eye_velocity_trace(isnan(eye_velocity_trace)) = 0;
    all_sac_validity   = false(1, num_saccades);
    all_sac_ind_onset  = nan(1, num_saccades);
    all_sac_ind_vmax   = nan(1, num_saccades);
    all_sac_ind_offset = nan(1, num_saccades);
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
        all_sac_ind_onset(:,counter_saccade)  = output_.ind_start(:);
        all_sac_ind_vmax(:,counter_saccade)   = output_.ind_vmax(:);
        all_sac_ind_offset(:,counter_saccade) = output_.ind_finish(:);
    end
    
    inds_span_    = ((-60+1) : 1 : (90))';
    num_sac_datapoints = length(inds_span_);
    length_time_      = length(TRIAL.eye_r_vm_filt);
    
    all_sac_inds = repmat( all_sac_ind_vmax(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(all_sac_ind_vmax), 1);
    all_sac_inds( all_sac_inds < 1 ) = 1;
    all_sac_inds( all_sac_inds > length_time_ ) = length_time_;
    all_sac_inds = all_sac_inds';
    %
    for counter_iteration = 1 : 5
        all_sac_ind_onset_  = all_sac_ind_onset(:,all_sac_validity);
        all_sac_ind_vmax_   = all_sac_ind_vmax(:,all_sac_validity);
        all_sac_ind_offset_ = all_sac_ind_offset(:,all_sac_validity);
        all_sac_inds_       = all_sac_inds(:,all_sac_validity);
        all_sac_validity_   = true(size(all_sac_ind_onset_));
        
        eye_velocity_trace_     = eye_velocity_trace;
        all_sac_eye_r_vm_ = reshape(eye_velocity_trace_(all_sac_inds_),num_sac_datapoints, []);
        all_sac_eye_r_px_ = reshape(TRIAL.eye_r_px_filt(all_sac_inds_),num_sac_datapoints, []);
        all_sac_eye_r_py_ = reshape(TRIAL.eye_r_py_filt(all_sac_inds_),num_sac_datapoints, []);
        all_sac_eye_r_vm_max_ = reshape(eye_velocity_trace_(all_sac_ind_vmax_), 1, []);
        all_sac_duration_ = all_sac_ind_offset_ - all_sac_ind_onset_;
        all_sac_eye_r_amp_m_ = sqrt(...
            ((TRIAL.eye_r_px_filt(all_sac_ind_offset_) - TRIAL.eye_r_px_filt(all_sac_ind_onset_)).^2) + ...
            ((TRIAL.eye_r_py_filt(all_sac_ind_offset_) - TRIAL.eye_r_py_filt(all_sac_ind_onset_)).^2) );
        all_sac_eye_r_amp_m_ = reshape(all_sac_eye_r_amp_m_, 1, []);
        
        all_sac_validity_( max(all_sac_eye_r_vm_) > 1200 ) = false; % invalidate if saccade trace has velocity > 1000 deg/s
        all_sac_validity_( max(all_sac_eye_r_vm_(1:40,:)) > 300 ) = false; % invalidate if saccade trace during 1-40ms has velocity > 300 deg/s
        all_sac_validity_( max(all_sac_eye_r_vm_(100:end,:)) > 300 ) = false; % invalidate if saccade trace during 100-150ms has velocity > 300 deg/s
        all_sac_validity_( all_sac_eye_r_vm_max_ < threshold ) = false; % invalidate if vm_max < threshold deg/s
        all_sac_validity_( all_sac_duration_ > 100 ) = false; % invalidate if saccades duration > 100 ms
        all_sac_validity_( all_sac_ind_offset_-all_sac_ind_vmax_ > 90 ) = false; % invalidate if deceleration > 90 ms
        all_sac_validity_( all_sac_ind_vmax_-all_sac_ind_onset_ > 60 ) = false; % invalidate if acceleration > 60 ms
        all_sac_validity_( all_sac_ind_offset_-all_sac_ind_vmax_ < 1 ) = false; % invalidate if deceleration < 1 ms
        all_sac_validity_( all_sac_ind_vmax_-all_sac_ind_onset_ < 1 ) = false; % invalidate if acceleration < 1 ms
        all_sac_validity_( all_sac_eye_r_amp_m_ > 20 ) = false;  % invalidate if saccade amplitude > 20 deg
%         all_sac_validity_( all_sac_eye_r_amp_m_ < 0.5 ) = false;  % invalidate if saccade amplitude < 0.5 deg
        all_sac_validity_( max(abs(all_sac_eye_r_px_)) > 20 ) = false;   % invalidate if abs of eye position > 20 deg
        all_sac_validity_( max(abs(all_sac_eye_r_py_)) > 20 ) = false;   % invalidate if abs of eye position > 20 deg
        
        all_sac_validity_( [(abs(diff(all_sac_ind_vmax_)) < 5) false] ) = false; % invalidate if ind_vmax is the same
        all_sac_validity_( [(abs(diff(all_sac_ind_onset_)) < 5) false] ) = false; % invalidate if ind_onset is the same
        all_sac_validity_( [(abs(diff(all_sac_ind_offset_)) < 5) false] ) = false; % invalidate if ind_offset is the same
        
        all_sac_ind_onset  = all_sac_ind_onset_(:,all_sac_validity_);
        all_sac_ind_vmax   = all_sac_ind_vmax_(:,all_sac_validity_);
        all_sac_ind_offset = all_sac_ind_offset_(:,all_sac_validity_);
        all_sac_inds       = all_sac_inds_(:,all_sac_validity_);
        all_sac_validity   = true(size(all_sac_ind_onset));
    end
    
    %% Init SACS_ALL and add common params
    SACS_ALL_TRIAL = struct;
    SACS_ALL_TRIAL.validity    = reshape(all_sac_validity,1, []);
    SACS_ALL_TRIAL.trial_num   = reshape(ones(size(all_sac_validity))*counter_trial, 1, []);
    SACS_ALL_TRIAL.tag         = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.count       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.time_onset  = reshape(TRIAL.time_1K(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.time_vmax   = reshape(TRIAL.time_1K(all_sac_ind_vmax),  1, []);
    SACS_ALL_TRIAL.time_offset = reshape(TRIAL.time_1K(all_sac_ind_offset),1, []);
    
    SACS_ALL_TRIAL.time_visual      = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.time_visual_2    = nan(size(all_sac_validity)); % this is a special scenario and applies to back_center_prim tag 7
    SACS_ALL_TRIAL.visual_px_onset  = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_px_offset = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_py_onset  = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_py_offset = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.reaction         = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_x     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_y     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_m     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_ang       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_start       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_finish      = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_ang         = nan(size(all_sac_validity));
    
    SACS_ALL_TRIAL.duration         = reshape((SACS_ALL_TRIAL.time_offset - SACS_ALL_TRIAL.time_onset) * 1000.0, 1, []);
    SACS_ALL_TRIAL.eye_r_vm_max     = reshape(TRIAL.eye_r_vm_filt(all_sac_ind_vmax),   1, []);
    SACS_ALL_TRIAL.eye_r_px_onset   = reshape(TRIAL.eye_r_px_filt(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.eye_r_px_offset  = reshape(TRIAL.eye_r_px_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_r_py_onset   = reshape(TRIAL.eye_r_py_filt(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.eye_r_py_offset  = reshape(TRIAL.eye_r_py_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_r_amp_x      = reshape((SACS_ALL_TRIAL.eye_r_px_offset - SACS_ALL_TRIAL.eye_r_px_onset), 1, []);
    SACS_ALL_TRIAL.eye_r_amp_y      = reshape((SACS_ALL_TRIAL.eye_r_py_offset - SACS_ALL_TRIAL.eye_r_py_onset), 1, []);
    SACS_ALL_TRIAL.eye_r_amp_m      = reshape(sqrt((SACS_ALL_TRIAL.eye_r_amp_x.^2) + (SACS_ALL_TRIAL.eye_r_amp_y.^2)), 1, []);
    SACS_ALL_TRIAL.eye_r_ang        = reshape((atan2d(SACS_ALL_TRIAL.eye_r_amp_y, SACS_ALL_TRIAL.eye_r_amp_x)),        1, []);
    SACS_ALL_TRIAL.tgt_px_onset     = reshape(TRIAL.tgt_px(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.tgt_px_offset    = reshape(TRIAL.tgt_px(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.tgt_py_onset     = reshape(TRIAL.tgt_py(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.tgt_py_offset    = reshape(TRIAL.tgt_py(all_sac_ind_offset), 1, []);
    
    %% Tag saccades
    threshold_pos = 1.5; % deg
    threshold_ang = 45.0; % deg
    
    % Tag the last prim_success
    % search the period after the last cue_present till the next_trial
    time_start_search  = TRIAL.time_state_cue_present(end);
    time_finish_search = TRIAL.time_state_next_trial(end);
    idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');
    if ~isempty(idx_sac)
        SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
        SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.start_x;
        SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.start_y;
        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.cue_x;
        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.cue_y;
        SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
        SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
        SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
        SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
        SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
        SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
            ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
            ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
        SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
            ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
            ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
        diff_finish_tgt_end = sqrt( ...
            ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - TRIAL.end_x).^2) + ...
            ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - TRIAL.end_y).^2) );
        SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
            ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) * SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
              (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) * SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
            / (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) / (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) );
        validity_sac = true;
        validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos);
        validity_sac = validity_sac && ( (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos) || (diff_finish_tgt_end < threshold_pos) );
        validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_ang(   idx_sac) < threshold_ang);
        
        if validity_sac
            SACS_ALL_TRIAL.tag(idx_sac) = 1; % 'prim_success' tag 1
        elseif ( SACS_ALL_TRIAL.diff_ang(idx_sac) < threshold_ang ) && (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos)
            SACS_ALL_TRIAL.tag(idx_sac) = 2; % 'prim_attempt' tag 2
        else
            SACS_ALL_TRIAL.tag(idx_sac) = 3; % 'prim_fail' tag 3
        end
    end
    
    % Tag the prim_fail
    % search the period after cue_present till the next str_present
    % exclude the last cue_present since we already took care of it
    length_cue_presentation = length(TRIAL.time_state_cue_present) - 1;
    if length_cue_presentation > 0
    for counter_cue_pres = 1 : length_cue_presentation
        time_start_search  = TRIAL.time_state_cue_present(counter_cue_pres);
        time_finish_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation > time_start_search, 1, 'first') );
        idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');
        if ~isempty(idx_sac)
            SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
            SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.start_x;
            SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.start_y;
            SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.cue_x;
            SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.cue_y;
            SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
            SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
            SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
            SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
            diff_finish_tgt_end = sqrt( ...
                ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - TRIAL.end_x).^2) + ...
                ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - TRIAL.end_y).^2) );
            SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
                ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                  (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) );
            validity_sac = true;
            validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos);
            validity_sac = validity_sac && ( (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos) || (diff_finish_tgt_end < threshold_pos) );
            validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_ang(   idx_sac) < threshold_ang);

            if validity_sac
                SACS_ALL_TRIAL.tag(idx_sac) = 1; % 'prim_success' tag 1
            elseif ( SACS_ALL_TRIAL.diff_ang(idx_sac) < threshold_ang ) && (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos)
                SACS_ALL_TRIAL.tag(idx_sac) = 2; % 'prim_attempt' tag 2
            else
                SACS_ALL_TRIAL.tag(idx_sac) = 3; % 'prim_fail' tag 3
            end
        end
    end
    end
    
    % Tag the corr_success and corr_fail
    % search the period after a success sac_prim offset till the next str_present
    % or till the next_trial
    idx_prim_success = find(SACS_ALL_TRIAL.tag == 1); % find prim_success
    if ~isempty(idx_prim_success)
        for counter_prim = 1 : length(idx_prim_success)
            idx_prim = idx_prim_success(counter_prim);
            time_last_cue_pres = TRIAL.time_state_cue_present(end);
            time_start_search = SACS_ALL_TRIAL.time_offset(idx_prim);
            if time_start_search < time_last_cue_pres
                time_finish_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation > time_start_search, 1, 'first') );
            else
                time_finish_search = TRIAL.time_state_next_trial(end);
            end
            idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');
            if ~isempty(idx_sac)
                SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
                % we will use tgt_cue and tgt_end for visual params, but we will use
                % the prim offset for diff values
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.cue_x;
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.cue_y;
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.end_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.end_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                % we will use the prim offset for diff values
                visual_amp_x_prim = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.eye_r_px_offset(idx_prim));
                visual_amp_y_prim = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.eye_r_py_offset(idx_prim));
                visual_amp_m_prim = sqrt((visual_amp_x_prim.^2) + (visual_amp_y_prim.^2));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.eye_r_px_offset( idx_prim)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.eye_r_py_offset( idx_prim)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
                    ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) .* visual_amp_x_prim) + ...
                      (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) .* visual_amp_y_prim) ) ...
                    ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) ./ (visual_amp_m_prim) );
                validity_sac = true;
                validity_sac = validity_sac && (SACS_ALL_TRIAL.visual_amp_m(idx_sac) > eps);
                validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos);
                validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos);
                validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_ang(   idx_sac) < threshold_ang);
                if validity_sac
                    SACS_ALL_TRIAL.tag(idx_sac) = 4; % 'corr_success' tag 4
                else
                    SACS_ALL_TRIAL.tag(idx_sac) = 5; % 'corr_fail' tag 5
                end
            end
        end
    end
    
    % Tag back_center
    % search the period before cue_present till the previous str_present (str_pursuit)
    % tag them back_center_irrelev, later we will revisit these tags.
    for counter_cue_pres = 1 : length(TRIAL.time_state_cue_present)
        time_finish_search  = TRIAL.time_state_cue_present(counter_cue_pres);
        if counter_cue_pres==1
            time_start_search = TRIAL.time_start;
        else
            time_start_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation < time_finish_search, 1, 'last') ) - TRIAL.time_pursuit; % str_pursuit is the shiftet version of str_present
        end
        idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'last');
        if ~isempty(idx_sac)
            if (SACS_ALL_TRIAL.tag(idx_sac) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac) <= 3)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3
                % if the saccade is tagged as a prim sac, do not re-tag it
                % and skip to next saccade.
                continue;
            end
            SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
            SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_px_onset(idx_sac); % this entry is for the sake of filling up the nan values
            SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_py_onset(idx_sac); % this entry is for the sake of filling up the nan values
            SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
            SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
            SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
            validity_sac = true;
            validity_sac = validity_sac && (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos);

            if validity_sac
                SACS_ALL_TRIAL.tag(idx_sac) = 8; % 'back_center_irrelev' tag 8
            else
                SACS_ALL_TRIAL.tag(idx_sac) = 10; % 'other_irrelev' tag 10
            end
        end
    end
    
    % Re-tag the back_center_irrelev
    % review tags and based on the events before these tags, re-tag them as
    % back_center_success or back_center_prim
    idx_back_center = find(SACS_ALL_TRIAL.tag == 8); % find back_center_irrelev
    if ~isempty(idx_back_center)
        for counter_back_center = 1 : length(idx_back_center)
            idx_sac = idx_back_center(counter_back_center);
            if (idx_sac == 1)
                % this is the 1st sac of the trial, i.e., back_center_success
                SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_start;
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_px_onset(idx_sac);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_py_onset(idx_sac);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
                    ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                      (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                    ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) );
                SACS_ALL_TRIAL.tag(idx_sac) = 6; % 'back_center_success' tag 6
                continue;
            end
            if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 5)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'corr_success' tag 4 % 'cord_fail' tag 5
                % If the previous sac is a prim or corr saccade, then this is a back_center_prim.
                % The time_visual in case of back_center_prim, can be considered as when the animal observes the target
                % after the str_pursuit which is the time_visual. However, another definition can be when the boop got played 
                % at the end of prim saccade. To cover this scenario we added time_visual_2.
                SACS_ALL_TRIAL.time_visual_2(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                % In case of corrective saccades the time_visual should be set to the offset of corr sac and not str_pursuit
                if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 4) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 5)
                    % 'corr_success' tag 4 % 'cord_fail' tag 5
                    SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                end
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_px_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_py_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
                    ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                      (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                    ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) );
                SACS_ALL_TRIAL.tag(idx_sac) = 7; % 'back_center_prim' tag 7
                continue;
            end
        end
    end
    
    % Tag the potential back_center_prim which happen after prim or corr
    % but did not trigger the "time_state_cue_present"
    idx_prim_corr = find((SACS_ALL_TRIAL.tag >= 1) & (SACS_ALL_TRIAL.tag <= 5));
    % find 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'corr_success' tag 4 % 'cord_fail' tag 5
    if ~isempty(idx_prim_corr)
        for counter_prim_corr = 1 : length(idx_prim_corr)
            idx_sac = idx_prim_corr(counter_prim_corr)+1;
            if idx_sac > length(SACS_ALL_TRIAL.tag)
                % the idx is out of range. break the for loop.
                break;
            end
            if ~isnan(SACS_ALL_TRIAL.tag(idx_sac))
                % the saccade has been tagged before, go to the next
                % saccade
                continue;
            end
            diff_finish = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - TRIAL.start_x).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - TRIAL.start_y).^2) );
            if diff_finish < threshold_pos
                % this is a missed back_center_prim and should be tagged.
                
                
                % The time_visual in case of back_center_prim, can be considered as when the animal observes the target
                % after the str_pursuit which is the time_visual. However, another definition can be when the boop got played 
                % at the end of prim saccade. To cover this scenario we added time_visual_2.
                SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                SACS_ALL_TRIAL.time_visual_2(idx_sac)    = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 3)
                    % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3
                    SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1) + TRIAL.time_punishment;
                end
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_px_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_r_py_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_r_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_r_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = acosd( ...
                    ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                      (SACS_ALL_TRIAL.eye_r_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                    ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) );
                SACS_ALL_TRIAL.tag(idx_sac) = 7; % 'back_center_prim' tag 7
            end
        end
    end
    
    % Tag target_irrelev. These are the remaining saccades that the saccade offset has landed on the target.
    diff_finish = sqrt( ...
        ((SACS_ALL_TRIAL.eye_r_px_offset - SACS_ALL_TRIAL.tgt_px_offset).^2) + ...
        ((SACS_ALL_TRIAL.eye_r_py_offset - SACS_ALL_TRIAL.tgt_py_offset).^2) );
    idx_target_irrelev = isnan(SACS_ALL_TRIAL.tag) & (diff_finish < threshold_pos);
    SACS_ALL_TRIAL.visual_px_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_r_px_onset(idx_target_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_py_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_r_py_onset(idx_target_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_px_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_px_offset( idx_target_irrelev);
    SACS_ALL_TRIAL.visual_py_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_py_offset( idx_target_irrelev);
    SACS_ALL_TRIAL.diff_finish(     idx_target_irrelev) = diff_finish(idx_target_irrelev);
    SACS_ALL_TRIAL.tag(             idx_target_irrelev) = 9; % 'target_irrelev' tag 9
    
    % Tag the remaining saccades as other_irrelev
    idx_other_irrelev = isnan(SACS_ALL_TRIAL.tag);
    SACS_ALL_TRIAL.visual_px_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_r_px_onset( idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_py_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_r_py_onset( idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_px_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_r_px_offset(idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_py_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_r_py_offset(idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.diff_finish(     idx_other_irrelev) = 0.0;
    SACS_ALL_TRIAL.tag(             idx_other_irrelev) = 10; % 'other_irrelev' tag 10
    
    % Fill up the nan values
    idx_nan_visual_values = isnan(SACS_ALL_TRIAL.visual_amp_x);
    SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values)     = (SACS_ALL_TRIAL.visual_px_offset(idx_nan_visual_values) - SACS_ALL_TRIAL.visual_px_onset(idx_nan_visual_values));
    SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values)     = (SACS_ALL_TRIAL.visual_py_offset(idx_nan_visual_values) - SACS_ALL_TRIAL.visual_py_onset(idx_nan_visual_values));
    SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values).^2));
    SACS_ALL_TRIAL.visual_ang(idx_nan_visual_values)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values), SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values));
    SACS_ALL_TRIAL.diff_start(idx_nan_visual_values)  = sqrt( ...
        ((SACS_ALL_TRIAL.eye_r_px_onset( idx_nan_visual_values) - SACS_ALL_TRIAL.visual_px_onset( idx_nan_visual_values)).^2) + ...
        ((SACS_ALL_TRIAL.eye_r_py_onset( idx_nan_visual_values) - SACS_ALL_TRIAL.visual_py_onset( idx_nan_visual_values)).^2) );
    SACS_ALL_TRIAL.diff_ang(idx_nan_visual_values) = acosd( ...
        ( (SACS_ALL_TRIAL.eye_r_amp_x(idx_nan_visual_values) .* SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values)) + ...
          (SACS_ALL_TRIAL.eye_r_amp_y(idx_nan_visual_values) .* SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values)) ) ...
        ./ (SACS_ALL_TRIAL.eye_r_amp_m(idx_nan_visual_values)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values)) );
    
    
    % Set the reaction time of target_irrelev and other_irrelev as the time
    % from the offset of the previous saccade.
    idx_irrelev = find((SACS_ALL_TRIAL.tag==9) | (SACS_ALL_TRIAL.tag==10));
     % 'target_irrelev' tag 9 % 'other_irrelev' tag 10
    if ~isempty(idx_irrelev)
        for counter_target_irrelev = 1 : length(idx_irrelev)
            idx_sac = idx_irrelev(counter_target_irrelev);
            if idx_sac == 1
                SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_start;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            else
                SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_onset(idx_sac-1);
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            end
        end
    end
    
    % Fill up the count values
    prim_tags = (SACS_ALL_TRIAL.tag==1) | (SACS_ALL_TRIAL.tag==2); % 'prim_success' tag 1 % 'prim_attempt' tag 2
    if sum(prim_tags) > 0
        SACS_ALL_TRIAL.count(prim_tags) = 1 : sum(prim_tags);
    end
    corr_tags = (SACS_ALL_TRIAL.tag==4) | (SACS_ALL_TRIAL.tag==5); % 'corr_success' tag 4 % 'cord_fail' tag 5
    if sum(corr_tags) > 0
        SACS_ALL_TRIAL.count(corr_tags) = 1 : sum(corr_tags);
    end
    back_center_tags = (SACS_ALL_TRIAL.tag==6) | (SACS_ALL_TRIAL.tag==7); % 'back_center_success' tag 6 % 'back_center_prim' tag 7
    if sum(back_center_tags) > 0
        SACS_ALL_TRIAL.count(back_center_tags) = 1 : sum(back_center_tags);
    end
    irrelev_tags = (SACS_ALL_TRIAL.tag==3) | (SACS_ALL_TRIAL.tag==8) | (SACS_ALL_TRIAL.tag==9) | (SACS_ALL_TRIAL.tag==10);
    % 'prim_fail' tag 3 % 'back_center_irrelev' tag 8 % 'target_irrelev' tag 9 % 'other_irrelev' tag 10
    if sum(irrelev_tags) > 0
        SACS_ALL_TRIAL.count(irrelev_tags) = 1 : sum(irrelev_tags);
    end
    
    %% Plot trial
    %{
    trial_num = counter_trial;
    
    hFig_ = figure(1);
    clf(hFig_);
    
    hAx_(1) = subplot(3,1,1);
    hold on;
    plot(TRIAL.time_1K, TRIAL.eye_r_px_filt)
    plot(TRIAL.time_1K, TRIAL.tgt_px)
    plot(TRIAL.time_1K(all_sac_ind_onset(all_sac_validity)), TRIAL.eye_r_px_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
    plot(TRIAL.time_1K(all_sac_ind_offset(all_sac_validity)), TRIAL.eye_r_px_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
    plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
    plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
    plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
    plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
    plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
    plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
    ylim([-10 10])
%     xlabel('Time (s)')
%     set(gca,'xtick',[])
    ylabel('Horz. Eye (deg)')
    title(['Trial no.: ' num2str(trial_num) ', (191111-142741)'])
    
    hAx_(2) = subplot(3,1,2);
    hold on;
    plot(TRIAL.time_1K, TRIAL.eye_r_py_filt)
    plot(TRIAL.time_1K, TRIAL.tgt_py)
    plot(TRIAL.time_1K(all_sac_ind_onset(all_sac_validity)), TRIAL.eye_r_py_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
    plot(TRIAL.time_1K(all_sac_ind_offset(all_sac_validity)), TRIAL.eye_r_py_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
    plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
    plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
    plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
    plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
    plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
    plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
    ylim([-10 10])
%     xlabel('Time (s)')
%     set(gca,'xtick',[])
    ylabel('Vert. Eye (deg)')
    
    hAx_(3) = subplot(3,1,3);
    hold on;
    plot(TRIAL.time_1K, TRIAL.eye_r_vm_filt)
    plot(TRIAL.time_1K(all_sac_ind_onset(all_sac_validity)), TRIAL.eye_r_vm_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
    plot(TRIAL.time_1K(all_sac_ind_vmax(all_sac_validity)), TRIAL.eye_r_vm_filt(all_sac_ind_vmax(all_sac_validity)), 'or')
    plot(TRIAL.time_1K(all_sac_ind_offset(all_sac_validity)), TRIAL.eye_r_vm_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
    plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
    plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
    plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
    plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
    plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
    plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
    ylim([0 1000])
    xlabel('Time (s)')
    ylabel('Eye Vel. (deg/s)')
    
    for counter_sac = 1 : length(SACS_ALL_TRIAL.validity)
        tag_ = SACS_ALL_TRIAL.tag(counter_sac);
        if     tag_==1
            color_ = 'r';
        elseif tag_==2
            color_ = [0.6350    0.0780    0.1840]; % 'r';
        elseif tag_==3
            color_ = [0.4 0.4 0.4];
        elseif tag_==4
            color_ = 'b';
        elseif tag_==5
            color_ = [0    0.4470    0.7410]; % 'b';
        elseif tag_==6
            color_ = [0.4940    0.1840    0.5560]; % 'm';
        elseif tag_==7
            color_ = [0.4940    0.1840    0.5560]; % 'm';
        elseif tag_==8
            color_ = [0.6 0.6 0.6];
        elseif tag_==9
            color_ = [0.6 0.6 0.6];
        elseif tag_==10
            color_ = [0.6 0.6 0.6];
        else
            color_ = [0.6 0.6 0.6];
        end
        label_ = [EXPERIMENT_PARAMS.sac_tag_list{tag_} '_' num2str(SACS_ALL_TRIAL.count(counter_sac))];
        xline(hAx_(1), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
        xline(hAx_(2), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
        xline(hAx_(3), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
%         linkaxes(hAx_,'x')
    end
    w = waitforbuttonpress;
%     ESN_Beautify_Plot(hFig_, [8,4], 12)
    %}
    
    %% Save Saccades
    SACS_ALL(counter_trial) = SACS_ALL_TRIAL;
    
    %% print a dot every 20 trials
    if rem(counter_trial, 20) == 0
        fprintf('.');
    end
    
end
fprintf(' --> Completed. \n')

%% Arrange 'SACS_ALL_DATA'
clearvars -except EXPERIMENT_PARAMS ...
    TRIALS_DATA SACS_ALL;
fprintf([EXPERIMENT_PARAMS.file_name ': Arranging SACS_ALL_DATA ...']);
clearvars('SACS_ALL_DATA'); SACS_ALL_DATA = struct;
field_names_SACS_ALL_DATA = fieldnames(SACS_ALL);
SACS_ALL_DATA_cell = struct2cell(SACS_ALL);
for counter_fields = 1 : 1 : length(field_names_SACS_ALL_DATA)
    SACS_ALL_DATA_field_cell = SACS_ALL_DATA_cell(counter_fields,:,:);
    SACS_ALL_DATA_field_cell = reshape(SACS_ALL_DATA_field_cell, 1, []);
    SACS_ALL_DATA_field_mat = cell2mat(SACS_ALL_DATA_field_cell);
    SACS_ALL_DATA.(field_names_SACS_ALL_DATA{counter_fields}) = SACS_ALL_DATA_field_mat;
end
fprintf(' --> Completed. \n');

%% Add saccade trajectory to the SACS_ALL_DATA
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA SACS_ALL_DATA;
fprintf([EXPERIMENT_PARAMS.file_name ': Adding saccade trajectory ...']);
num_sacs = length(SACS_ALL_DATA.validity);
SACS_ALL_DATA.eye_r_px = nan(150, num_sacs);
SACS_ALL_DATA.eye_r_py = nan(150, num_sacs);
for counter_sac = 1 : num_sacs
    sac_time_vmax = SACS_ALL_DATA.time_vmax(counter_sac);
    sac_trial_num = SACS_ALL_DATA.trial_num(counter_sac);
    trial_time_1K       = TRIALS_DATA.time_1K{      1,sac_trial_num};
    trial_eye_r_px_filt = TRIALS_DATA.eye_r_px_filt{1,sac_trial_num};
    trial_eye_r_py_filt = TRIALS_DATA.eye_r_py_filt{1,sac_trial_num};
    ind_vmax = find(trial_time_1K >= sac_time_vmax, 1, 'first');
    if isempty(ind_vmax)
        continue;
    end
    inds_span_    = ((-60+1) : 1 : (90))';
    num_sac_datapoints = length(inds_span_);
    length_time_      = length(trial_time_1K);
    
    sac_inds = repmat( ind_vmax, 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_vmax), 1);
    sac_inds( sac_inds < 1 ) = 1;
    sac_inds( sac_inds > length_time_ ) = length_time_;
    sac_inds = reshape(sac_inds,num_sac_datapoints, 1);
    SACS_ALL_DATA.eye_r_px(:,counter_sac) = reshape(trial_eye_r_px_filt(sac_inds), num_sac_datapoints, 1);
    SACS_ALL_DATA.eye_r_py(:,counter_sac) = reshape(trial_eye_r_py_filt(sac_inds), num_sac_datapoints, 1);
end
fprintf(' --> Completed. \n');

%% Plot SACS_ALL_DATA
if ~EXPERIMENT_PARAMS.flag_session_figure
    return;
end

% Set parameters
amp_edges = -.25 : 0.5 : 15.25;
ang_edges = (-pi-(pi/16)) : (pi/8) : (pi-(pi/16));
react_edges = -12.5: 25 : 512.5;
num_row = 9;
num_col = 9;

% Init plot
hFig = figure(1);
clf(hFig)
hold on

% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for counter_tag = 1 : 9
if counter_tag == 8
    idx_tag = (SACS_ALL_DATA.tag==8) | (SACS_ALL_DATA.tag==9);
    title_ = [EXPERIMENT_PARAMS.sac_tag_list{8} ' & ' EXPERIMENT_PARAMS.sac_tag_list{9}];
elseif counter_tag == 9
    idx_tag = (SACS_ALL_DATA.tag==10);
    title_ = EXPERIMENT_PARAMS.sac_tag_list{10};
else
    idx_tag = (SACS_ALL_DATA.tag==counter_tag);
    title_ = EXPERIMENT_PARAMS.sac_tag_list{counter_tag};
end

axes_minor_nums = reshape(1:num_row*num_col, num_row, num_col)';
axes_main_row = floor((counter_tag - 1) / 3)+1;
axes_main_col = mod(counter_tag, 3); if (axes_main_col==0); axes_main_col=3; end;
row1_ = ((axes_main_row-1)*3)+1;
row2_ = ((axes_main_row-1)*3)+2;
row3_ = ((axes_main_row-1)*3)+3;
col1_ = ((axes_main_col-1)*3)+1;
col2_ = ((axes_main_col-1)*3)+2;
col3_ = ((axes_main_col-1)*3)+3;

axes_trace = [axes_minor_nums(row1_,col1_), axes_minor_nums(row1_,col2_), axes_minor_nums(row1_,col3_)...
              axes_minor_nums(row2_,col1_), axes_minor_nums(row2_,col2_), axes_minor_nums(row2_,col3_) ];
axes_amp_dis = axes_minor_nums(row3_,col1_);
axes_ang_dis = axes_minor_nums(row3_,col2_);
axes_react_dis = axes_minor_nums(row3_,col3_);

subplot(num_row,num_col,axes_trace)
hold on
plot(SACS_ALL_DATA.eye_r_px(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py(:,idx_tag), 'k')
plot(SACS_ALL_DATA.eye_r_px_offset(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py_offset(:,idx_tag), 'om')
title([title_ ': ' num2str(sum(idx_tag)) ' sac'], 'interpret', 'none');
xlim([-17, 17])
ylim([-15, 15])
% axis equal;

subplot(num_row,num_col,axes_amp_dis)
hold on
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:3:15)
ylabel('Amplitude')

subplot(num_row,num_col,axes_ang_dis)
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
hold on
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
set(gca, 'ThetaTick', [])
set(gca, 'RTick', [])
set(gca, 'Title', [])

subplot(num_row,num_col,axes_react_dis)
hold on
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 500]/1000)
set(gca, 'XTick', (0:200:500)/1000)
ylabel('Reaction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgtitle([EXPERIMENT_PARAMS.file_name(1:13) ', ' ...
    'trial: ' num2str(length(TRIALS_DATA.time_1K)) ', ' ...
    'sac: ' num2str(length(SACS_ALL_DATA.validity)) ', ' ...
    'dur: ' num2str((TRIALS_DATA.time_end(end) - TRIALS_DATA.time_start(1))/60,3) 'min' ...
    ], ...
    'interpret', 'none');
ESN_Beautify_Plot(hFig, [13 13], 8)

end