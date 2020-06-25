%% Load CH_DATA_ALL and CH_DATA_PSORT
clear;
fprintf(['Loading CH_DATA_ALL and CH_DATA_PSORT', ' ... ']);
file_path = pwd;
[file_name,file_path] = uigetfile([file_path filesep '*.mat'], 'Select converted smr file');
[~, file_name, ~]  = fileparts(file_name);
load([file_path filesep file_name '.mat'])
CH_DATA_PSORT = Psort_read_psort([file_path filesep file_name '.psort']);

CH_DATA_PSORT.file_name = file_name;
CH_DATA_PSORT.file_path = file_path;
fprintf(' --> Completed. \n')

%% Build DATA
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
fprintf(['Building DATA', ' ... ']);
clearvars DATA;
DATA.ch_data = CH_DATA_PSORT.topLevel_data.ch_data;
DATA.time_50K = CH_DATA_PSORT.topLevel_data.ch_time;
DATA.ss_index = logical(CH_DATA_PSORT.topLevel_data.ss_index);

CH_DATA_ALL_cell = struct2cell(CH_DATA_ALL);
% field_names = CH_DATA_ALL_cell(4,:);
ch_number = cell2mat(CH_DATA_ALL_cell(1,:));
DATA.eye_px = CH_DATA_ALL(ch_number==2).values; % CH_DATA_ALL(ismember(field_names,'H Eye')).values; % 'HE'
DATA.eye_py = CH_DATA_ALL(ch_number==3).values; % CH_DATA_ALL(ismember(field_names,'V Eye')).values; % 'VE'
DATA.tgt_px = CH_DATA_ALL(ch_number==4).values; % CH_DATA_ALL(ismember(field_names,'H Targ')).values; % 'HT'
DATA.tgt_py = CH_DATA_ALL(ch_number==5).values; % CH_DATA_ALL(ismember(field_names,'V Targ')).values; % 'VT'
if sum(ch_number==10) > 0
    DATA.tgt_visible = (CH_DATA_ALL(ch_number==10).values > 10);
else
    DATA.tgt_visible = true(size(DATA.tgt_px));
end
length_time_1K = length(DATA.eye_px);
DATA.time_1K = linspace(CH_DATA_ALL(ch_number==2).offset, ...
                   CH_DATA_ALL(ch_number==2).max_time, length_time_1K)';

% filter params
sampling_freq = 1 / (CH_DATA_ALL(ch_number==2).interval); % 'HE'
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye data
DATA.eye_px_filt = filtfilt(b_butter,a_butter,DATA.eye_px);
DATA.eye_py_filt = filtfilt(b_butter,a_butter,DATA.eye_py);
DATA.eye_vx_filt = diff(DATA.eye_px_filt)./diff(DATA.time_1K); DATA.eye_vx_filt=[DATA.eye_vx_filt(1); DATA.eye_vx_filt];
DATA.eye_vy_filt = diff(DATA.eye_py_filt)./diff(DATA.time_1K); DATA.eye_vy_filt=[DATA.eye_vy_filt(1); DATA.eye_vy_filt];
DATA.eye_vm_filt = sqrt(DATA.eye_vx_filt.^2 + DATA.eye_vy_filt.^2);
DATA.eye_ax_filt = diff(DATA.eye_vx_filt)./diff(DATA.time_1K); DATA.eye_ax_filt=[DATA.eye_ax_filt(1); DATA.eye_ax_filt];
DATA.eye_ay_filt = diff(DATA.eye_vy_filt)./diff(DATA.time_1K); DATA.eye_ay_filt=[DATA.eye_ay_filt(1); DATA.eye_ay_filt];
DATA.eye_am_filt = sqrt(DATA.eye_ax_filt.^2 + DATA.eye_ay_filt.^2);

DATA.tgt_px_filt = ESN_Round(DATA.tgt_px, 1);
DATA.tgt_py_filt = ESN_Round(DATA.tgt_py, 1);

% Firing Rate
sample_freq = 50e3;
norm_dist_sigma = 2.5e-3; % seconds
norm_dist_x = -(4*norm_dist_sigma*sample_freq):1:(4*norm_dist_sigma*sample_freq);
norm_dist_y = normpdf(norm_dist_x,0, norm_dist_sigma*sample_freq);
time_50K = DATA.time_50K;
time_1K = DATA.time_1K;
event_spike_50K = DATA.ss_index;
time_spike_50K = time_50K(event_spike_50K);
instant_firing_rate = 1./diff(time_spike_50K);
instant_firing_rate_50K = zeros(size(time_50K));
index_spike_50K = find(event_spike_50K);
index_spike_50K = [1; index_spike_50K(:); length(event_spike_50K)];
instant_firing_rate = [0; instant_firing_rate(:); 0];
for counter_spike = 1 : 1 : length(index_spike_50K)-1
    ind_start = index_spike_50K(counter_spike);
    ind_finish = index_spike_50K(counter_spike+1);
    instant_firing_rate_50K(ind_start:ind_finish) = instant_firing_rate(counter_spike);
end
instant_firing_rate_50K_conv = conv(instant_firing_rate_50K, norm_dist_y, 'same');
instant_firing_rate_1K_conv = interp1(time_50K, instant_firing_rate_50K_conv, time_1K);
DATA.firing_rate_50K = instant_firing_rate_50K_conv;
DATA.firing_rate_1K = instant_firing_rate_1K_conv;

% Detect 40ms target flash
edge_rise = (diff(DATA.tgt_visible) == 1);
edge_rise = [false; edge_rise(:);];
edge_fall = (diff(DATA.tgt_visible) == -1);
edge_fall = [edge_fall(:); false;];
edge_rise_ind = find(edge_rise);
edge_fall_ind = find(edge_fall);
if ~(isempty(edge_fall_ind) || isempty(edge_rise_ind))
if edge_fall_ind(1) < edge_rise_ind(1)
    edge_fall_ind(1) = [];
end
if length(edge_fall_ind) > length(edge_rise_ind)
    edge_fall_ind = edge_fall_ind(1:length(edge_rise_ind));
end
if length(edge_rise_ind) > length(edge_fall_ind)
    edge_rise_ind = edge_rise_ind(1:length(edge_fall_ind));
end
time_diff_switch = DATA.time_1K(edge_fall_ind) - DATA.time_1K(edge_rise_ind);
target_flash_45ms = edge_rise_ind(time_diff_switch < 0.045);
DATA.target_flash_45ms = target_flash_45ms;
else
target_flash_45ms = false(size(DATA.tgt_visible));
DATA.target_flash_45ms = target_flash_45ms;
end
% Build event_tgt_jump_1K
% if sum(ch_number==32) > 0
%     Events_tgt_jump = (CH_DATA_ALL(ch_number==32).values_codes(:,1)==1);
%     DATA.time_tgt_jump_1K = CH_DATA_ALL(ch_number==32).values_times(Events_tgt_jump,1);
% else
    ind_visible = find(DATA.tgt_visible);
    tgt_px_visible = DATA.tgt_px(ind_visible);
    tgt_py_visible = DATA.tgt_py(ind_visible);
    
    tgt_jump_x_1 = [( abs(diff(tgt_px_visible)) > 1.0 ); false];
    tgt_jump_x_2 = [false; ( abs(diff(tgt_px_visible)) > 1.0 )];
    tgt_jump_x = tgt_jump_x_2;
    tgt_jump_x(tgt_jump_x_1&tgt_jump_x_2) = false;
    tgt_jump_y_1 = [( abs(diff(tgt_py_visible)) > 1.0 ); false];
    tgt_jump_y_2 = [false; ( abs(diff(tgt_py_visible)) > 1.0 )];
    tgt_jump_y = tgt_jump_y_2;
    tgt_jump_y(tgt_jump_y_1&tgt_jump_y_2) = false;
    tgt_jump_xy = tgt_jump_x | tgt_jump_y;
    ind_tgt_jump = ind_visible(tgt_jump_xy);
    event_tgt_jump_1K = false(size(DATA.tgt_visible));
    event_tgt_jump_1K(ind_tgt_jump) = true;
    event_tgt_jump_1K(target_flash_45ms) = false;
    DATA.time_tgt_jump_1K = DATA.time_1K(event_tgt_jump_1K);
    
    delta_tgt_jump_x = tgt_px_visible(tgt_jump_xy) - tgt_px_visible([tgt_jump_xy(2:end); false;]);
    delta_tgt_jump_y = tgt_py_visible(tgt_jump_xy) - tgt_py_visible([tgt_jump_xy(2:end); false;]);
    delta_tgt_jump_x_1K = nan(size(DATA.tgt_visible));
    delta_tgt_jump_x_1K(ind_tgt_jump) = delta_tgt_jump_x;
    delta_tgt_jump_x_1K(target_flash_45ms) = nan;
    delta_tgt_jump_y_1K = nan(size(DATA.tgt_visible));
    delta_tgt_jump_y_1K(ind_tgt_jump) = delta_tgt_jump_y;
    delta_tgt_jump_y_1K(target_flash_45ms) = nan;
    DATA.delta_tgt_jump_x_1K = delta_tgt_jump_x_1K;
    DATA.delta_tgt_jump_y_1K = delta_tgt_jump_y_1K;
% end

fprintf(' --> Completed. \n')

%% Finding events in time_1K
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
fprintf(['Finding events in time_1K', ' ... ']);
% Finding events in time_1K
time_reference = DATA.time_1K;
time_spike_1K = DATA.time_50K(DATA.ss_index);
time_tgt_jump_1K = DATA.time_tgt_jump_1K;

length_time = length(time_reference);
variable_list = {'spike_1K', 'tgt_jump_1K'};
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time_' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time_'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event_'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter_' variable_name ' = ' 'counter_temp_' ';']);
end

for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    %{
    % Do not use 'eval', it is very slow, instead use the actual variables
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ ...
            'if ( time_ponit_ >= time_' variable_name '(  counter_' variable_name ') );', ...
                'event_' variable_name '(    counter_time_point) = true;', ...
                'counter_' variable_name '   = counter_' variable_name '   + 1;', ...
            'end;', ...
        ]);
    end
    % Here is the template for actual variables
    if time_ponit_ >= time_VARIABLE(  counter_VARIABLE)
        event_VARIABLE(    counter_time_point) = true;
        counter_VARIABLE   = counter_VARIABLE   + 1;
    end
    %}
    if time_ponit_ >= time_spike_1K(  counter_spike_1K)
        event_spike_1K(    counter_time_point) = true;
        counter_spike_1K   = counter_spike_1K   + 1;
    end
    if time_ponit_ >= time_tgt_jump_1K(  counter_tgt_jump_1K)
        event_tgt_jump_1K(    counter_time_point) = true;
        counter_tgt_jump_1K   = counter_tgt_jump_1K   + 1;
    end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'DATA.event_' variable_name ' = ' 'event_' variable_name ';']);
end

fprintf(' --> Completed. \n')

%% Extract Directions
fprintf(['Extract Directions', ' ... ']);
ind_tgt_jump_1K = find(DATA.event_tgt_jump_1K);
% ind_tgt_jump_1K(ind_tgt_jump_1K<6) = 6;
delta_tgt_jump_x = DATA.delta_tgt_jump_x_1K(~isnan(DATA.delta_tgt_jump_x_1K)); % DATA.tgt_px(ind_tgt_jump_1K) - DATA.tgt_px(ind_tgt_jump_1K-1);
delta_tgt_jump_y = DATA.delta_tgt_jump_y_1K(~isnan(DATA.delta_tgt_jump_y_1K)); % DATA.tgt_py(ind_tgt_jump_1K) - DATA.tgt_py(ind_tgt_jump_1K-1);
angle_tgt_jump = atan2d(delta_tgt_jump_y, delta_tgt_jump_x);
events_000 = (angle_tgt_jump > -45) & (angle_tgt_jump <= +45);
events_090 = (angle_tgt_jump > +45) & (angle_tgt_jump <= +135);
events_270 = (angle_tgt_jump > -135) & (angle_tgt_jump <= -45);
events_180 = ((angle_tgt_jump > +135) & (angle_tgt_jump <= 225)) | ((angle_tgt_jump > -225) & (angle_tgt_jump <= -135));
DATA.tgt_jump_amplitude = sqrt(delta_tgt_jump_x.^2 + delta_tgt_jump_y.^2);
DATA.tgt_jump_angle = angle_tgt_jump;
DATA.trial_num_000 = find(events_000);
DATA.trial_num_090 = find(events_090);
DATA.trial_num_180 = find(events_180);
DATA.trial_num_270 = find(events_270);
fprintf(' --> Completed. \n')

%% Detect saccades
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
fprintf(['Building DATA_SAC', ' ... ']);
DATA_SAC = struct;
ind_tgt_jump_1K = find(DATA.event_tgt_jump_1K);
ind_tgt_jump_1K = [ind_tgt_jump_1K(:); length(DATA.event_tgt_jump_1K)];
num_trials = sum(DATA.event_tgt_jump_1K);

% sac_start_coarse = (DATA.eye_vm_filt(1:end-1) <= 100.0) & (DATA.eye_vm_filt(2:end) > 100.0);
% sac_start_coarse = [sac_start_coarse; false];
% sac_start_coarse(end) = true;
% num_sac = sum(sac_start_coarse) - 1;
% ind_sac_start_coarse = find(sac_start_coarse);
for counter_trial = 1 : num_trials
    eye_velocity_trace     = DATA.eye_vm_filt;
    ind_search_begin       = ind_tgt_jump_1K(counter_trial);
    ind_search_end         = ind_tgt_jump_1K(counter_trial+1);
    
    params_sac.MinPeakHeight       = 150.0; % deg/s
    params_sac.MinPeakProminence   = 100; % data points
    params_sac.rough_threshold     = 50.0; % deg/s
    params_sac.fine_threshold      = 20.0; % deg/s
    params_sac.sampling_freq       = 1000.0; % Hz
    params_sac.cutoff_freq         = 50.0; % Hz
    params_sac.window_half_length  = 4; % data points
    params_sac.prominence_or_first = 'first'; % which peak to select, 'prominent' or 'first'
    
    output_ = ESN_Sac_Finder(eye_velocity_trace, ...
        ind_search_begin, ind_search_end, params_sac);
    
    validity_sac   = output_.validity(:);
    inds_sac       = output_.inds(:);
    ind_sac_start  = output_.ind_start(:);
    ind_sac_vmax   = output_.ind_vmax(:);
    ind_sac_finish = output_.ind_finish(:);
    
    %% Save Primary Sac data to SAC_PRIM
    DATA_SAC.validity(counter_trial,:)                  = validity_sac;
    DATA_SAC.inds(counter_trial,:)                      = inds_sac;
    DATA_SAC.ind_start(counter_trial,:)                 = ind_sac_start;
    DATA_SAC.ind_vmax(counter_trial,:)                  = ind_sac_vmax;
    DATA_SAC.ind_finish(counter_trial,:)                = ind_sac_finish;
    
    DATA_SAC.time(counter_trial,:)                    = DATA.time_1K(    DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_px(counter_trial,:)                  = DATA.eye_px_filt(DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_py(counter_trial,:)                  = DATA.eye_py_filt(DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_vx(counter_trial,:)                  = DATA.eye_vx_filt(DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_vy(counter_trial,:)                  = DATA.eye_vy_filt(DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_vm(counter_trial,:)                  = DATA.eye_vm_filt(DATA_SAC.inds(counter_trial,:));
    DATA_SAC.eye_vm_max(counter_trial,:)              = DATA.eye_vm_filt(DATA_SAC.ind_vmax(counter_trial,:));
%     DATA_SAC.eye_px_centered(counter_sac,:)         = DATA_SAC.eye_px(counter_sac,:) - DATA.start_x;
%     DATA_SAC.eye_py_centered(counter_sac,:)         = DATA_SAC.eye_py(counter_sac,:) - DATA.start_y;
    
    DATA_SAC.eye_px_start(counter_trial,:)            = DATA.eye_px_filt(DATA_SAC.ind_start(counter_trial,:));
    DATA_SAC.eye_px_finish(counter_trial,:)           = DATA.eye_px_filt(DATA_SAC.ind_finish(counter_trial,:));
    DATA_SAC.eye_py_start(counter_trial,:)            = DATA.eye_py_filt(DATA_SAC.ind_start(counter_trial,:));
    DATA_SAC.eye_py_finish(counter_trial,:)           = DATA.eye_py_filt(DATA_SAC.ind_finish(counter_trial,:));
%     DATA_SAC.eye_px_start_centered(counter_sac,:)   = DATA_SAC.eye_px_start(counter_sac,:)  - DATA.start_x;
%     DATA_SAC.eye_px_finish_centered(counter_sac,:)  = DATA_SAC.eye_px_finish(counter_sac,:) - DATA.start_x;
%     DATA_SAC.eye_py_start_centered(counter_sac,:)   = DATA_SAC.eye_py_start(counter_sac,:)  - DATA.start_y;
%     DATA_SAC.eye_py_finish_centered(counter_sac,:)  = DATA_SAC.eye_py_finish(counter_sac,:) - DATA.start_y;
    DATA_SAC.eye_amp_x(counter_trial,:)               = (DATA_SAC.eye_px_finish(counter_trial,:) - DATA_SAC.eye_px_start(counter_trial,:));
    DATA_SAC.eye_amp_y(counter_trial,:)               = (DATA_SAC.eye_py_finish(counter_trial,:) - DATA_SAC.eye_py_start(counter_trial,:));
    DATA_SAC.eye_amp_m(counter_trial,:)               = (sqrt(DATA_SAC.eye_amp_x(counter_trial,:)^2+DATA_SAC.eye_amp_y(counter_trial,:)^2));
%     DATA_SAC.reaction(counter_sac,:)                  = DATA_SAC.ind_start - TRIAL.ind_state_cue_present;
    DATA_SAC.duration(counter_trial,:)                  = DATA_SAC.ind_finish(counter_trial,:) - DATA_SAC.ind_start(counter_trial,:);
end
% DATA_SAC.validity(isoutlier(DATA_SAC.duration)) = false;
% DATA_SAC.validity(isoutlier(DATA_SAC.eye_vm_max)) = false;

DATA.event_sac_start_1K = false(size(DATA.time_1K));
DATA.event_sac_start_1K(DATA_SAC.ind_start) = true;
DATA.event_sac_finish_1K = false(size(DATA.time_1K));
DATA.event_sac_finish_1K(DATA_SAC.ind_finish) = true;
DATA.event_sac_vmax_1K = false(size(DATA.time_1K));
DATA.event_sac_vmax_1K(DATA_SAC.ind_vmax) = true;
fprintf(' --> Completed. \n')

%% Finding events in time_50K
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
fprintf(['Finding events in time_50K', ' ... ']);

% Finding events in time_50K
time_reference = DATA.time_50K;
time_tgt_jump_50K   = DATA.time_1K(DATA.event_tgt_jump_1K);
time_sac_start_50K  = DATA.time_1K(DATA.event_sac_start_1K);
time_sac_finish_50K = DATA.time_1K(DATA.event_sac_finish_1K);
time_sac_vmax_50K   = DATA.time_1K(DATA.event_sac_vmax_1K);

length_time = length(time_reference);
variable_list = {'tgt_jump_50K', 'sac_start_50K', 'sac_finish_50K', 'sac_vmax_50K'};
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time_' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time_'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event_'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter_' variable_name ' = ' 'counter_temp_' ';']);
end

for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    %{
    % Do not use 'eval', it is very slow, instead use the actual variables
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ ...
            'if ( time_ponit_ >= time_' variable_name '(  counter_' variable_name ') );', ...
                'event_' variable_name '(    counter_time_point) = true;', ...
                'counter_' variable_name '   = counter_' variable_name '   + 1;', ...
            'end;', ...
        ]);
    end
    % Here is the template for actual variables
    if time_ponit_ >= time_VARIABLE(  counter_VARIABLE)
        event_VARIABLE(    counter_time_point) = true;
        counter_VARIABLE   = counter_VARIABLE   + 1;
    end
    %}
    if time_ponit_ >= time_tgt_jump_50K(  counter_tgt_jump_50K )
        event_tgt_jump_50K( counter_time_point ) = true;
        counter_tgt_jump_50K = counter_tgt_jump_50K   + 1;
    end
    if time_ponit_ >= time_sac_start_50K(  counter_sac_start_50K)
        event_sac_start_50K(    counter_time_point) = true;
        counter_sac_start_50K   = counter_sac_start_50K   + 1;
    end
    if time_ponit_ >= time_sac_finish_50K(  counter_sac_finish_50K)
        event_sac_finish_50K(    counter_time_point) = true;
        counter_sac_finish_50K   = counter_sac_finish_50K   + 1;
    end
    if time_ponit_ >= time_sac_vmax_50K(  counter_sac_vmax_50K)
        event_sac_vmax_50K(    counter_time_point) = true;
        counter_sac_vmax_50K   = counter_sac_vmax_50K   + 1;
    end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'DATA.event_' variable_name ' = ' 'event_' variable_name ';']);
end

fprintf(' --> Completed. \n')

%% Build Raster Data
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
fprintf(['Building Raster Data', ' ... ']);
length_time_      = length(DATA.time_1K);
ind_cue_onset_1K = find(DATA.event_tgt_jump_1K);
ind_sac_onset_1K = (DATA_SAC.ind_start);
inds_span_cue_onset    = ((-300+1) : 1 : (300))';
inds_span_sac_onset  = ((-300+1) : 1 : (300))';
inds_cue_onset_1K = repmat( ind_cue_onset_1K(:), 1, length(inds_span_cue_onset)) + repmat(inds_span_cue_onset(:)', length(ind_cue_onset_1K), 1);
inds_cue_onset_1K( inds_cue_onset_1K < 1 ) = 1;
inds_cue_onset_1K( inds_cue_onset_1K > length_time_ ) = length_time_;

inds_sac_onset_1K = repmat( ind_sac_onset_1K(:), 1, length(inds_span_sac_onset)) + repmat(inds_span_sac_onset(:)', length(ind_sac_onset_1K), 1);
inds_sac_onset_1K( inds_sac_onset_1K < 1 ) = 1;
inds_sac_onset_1K( inds_sac_onset_1K > length_time_ ) = length_time_;

DATA.inds_span_cue_onset = inds_span_cue_onset;
DATA.inds_span_sac_onset = inds_span_sac_onset;
DATA.inds_cue_onset_1K = inds_cue_onset_1K;
DATA.inds_sac_onset_1K = inds_sac_onset_1K;
fprintf(' --> Completed. \n')

%% Plot Raster data sac_onset
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
clearvars hAx firing_range
% plot cue_onset
hFig = figure(1);
clf(hFig);

hAx(1) = subplot(3,3,5);
hAx_ = hAx(1);
inds_sac_onset_1K = DATA.inds_sac_onset_1K;
tgt_jump_amplitude = DATA.tgt_jump_amplitude;
[~, ind_] = sort(tgt_jump_amplitude);
inds_sac_onset_1K = inds_sac_onset_1K(ind_,:);
inds_span = DATA.inds_span_sac_onset;
ESN_plot_raster(hAx_, DATA, inds_sac_onset_1K, inds_span);
axes(hAx_);
title('All Directions')

hAx(2) = subplot(3,3,6);
hAx_ = hAx(2);
if length(DATA.trial_num_000) > 1
inds_sac_onset_1K = DATA.inds_sac_onset_1K(DATA.trial_num_000, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_000, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_sac_onset_1K = inds_sac_onset_1K(ind_,:);
inds_span = DATA.inds_span_sac_onset;
ESN_plot_raster(hAx_, DATA, inds_sac_onset_1K, inds_span);
end
axes(hAx_);
title('Right')
xlabel('Saccade Onset (ms)')
yyaxis right;
ylabel('Firing Rate (Hz)')

hAx(3) = subplot(3,3,2);
hAx_ = hAx(3);
if length(DATA.trial_num_090) > 1
inds_sac_onset_1K = DATA.inds_sac_onset_1K(DATA.trial_num_090, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_090, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_sac_onset_1K = inds_sac_onset_1K(ind_,:);
inds_span = DATA.inds_span_sac_onset;
ESN_plot_raster(hAx_, DATA, inds_sac_onset_1K, inds_span);
end
axes(hAx_);
ylabel('Trial')
yyaxis right;
ylabel('Firing Rate (Hz)')

hAx(4) = subplot(3,3,4);
hAx_ = hAx(4);
if length(DATA.trial_num_180) > 1
inds_sac_onset_1K = DATA.inds_sac_onset_1K(DATA.trial_num_180, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_180, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_sac_onset_1K = inds_sac_onset_1K(ind_,:);
inds_span = DATA.inds_span_sac_onset;
ESN_plot_raster(hAx_, DATA, inds_sac_onset_1K, inds_span);
end
axes(hAx_);
title('Left')
ylabel('Trial')
xlabel('Saccade Onset (ms)')

hAx(5) = subplot(3,3,8);
hAx_ = hAx(5);
if length(DATA.trial_num_270) > 1
inds_sac_onset_1K = DATA.inds_sac_onset_1K(DATA.trial_num_270, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_270, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_sac_onset_1K = inds_sac_onset_1K(ind_,:);
inds_span = DATA.inds_span_sac_onset;
ESN_plot_raster(hAx_, DATA, inds_sac_onset_1K, inds_span);
end
axes(hAx_);
title('Down')
ylabel('Trial')
xlabel('Saccade Onset (ms)')
yyaxis right;
ylabel('Firing Rate (Hz)')

sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');

ESN_Beautify_Plot(hFig, [11 8.5])
%%
fprintf(['Saving plots', ' ...'])
file_name = CH_DATA_PSORT.file_name;
file_path = CH_DATA_PSORT.file_path;
saveas(hFig,[file_path filesep file_name '_modulation_sac_onset.png'], 'png');
saveas(hFig,[file_path filesep file_name '_modulation_sac_onset.pdf'], 'pdf');
close(hFig)
fprintf(' --> Completed. \n')
%}

%% Plot Raster data cue_onset
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
% plot cue_onset
hFig = figure(2);
clf(hFig);

hAx(1) = subplot(3,3,5);
hAx_ = hAx(1);
inds_cue_onset_1K = DATA.inds_cue_onset_1K;
tgt_jump_amplitude = DATA.tgt_jump_amplitude;
[~, ind_] = sort(tgt_jump_amplitude);
inds_cue_onset_1K = inds_cue_onset_1K(ind_,:);
inds_span = DATA.inds_span_cue_onset;
ESN_plot_raster(hAx_, DATA, inds_cue_onset_1K, inds_span);
axes(hAx_);
title('All Directions')

hAx(2) = subplot(3,3,6);
hAx_ = hAx(2);
if length(DATA.trial_num_000) > 1
inds_cue_onset_1K = DATA.inds_cue_onset_1K(DATA.trial_num_000, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_000, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_cue_onset_1K = inds_cue_onset_1K(ind_,:);
inds_span = DATA.inds_span_cue_onset;
ESN_plot_raster(hAx_, DATA, inds_cue_onset_1K, inds_span);
end
axes(hAx_);
title('Right')
xlabel('Cue Onset (ms)')
yyaxis right;
ylabel('Firing Rate (Hz)')

hAx(3) = subplot(3,3,2);
hAx_ = hAx(3);
if length(DATA.trial_num_090) > 1
inds_cue_onset_1K = DATA.inds_cue_onset_1K(DATA.trial_num_090, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_090, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_cue_onset_1K = inds_cue_onset_1K(ind_,:);
inds_span = DATA.inds_span_cue_onset;
ESN_plot_raster(hAx_, DATA, inds_cue_onset_1K, inds_span);
end
axes(hAx_);
ylabel('Trial')
yyaxis right;
ylabel('Firing Rate (Hz)')

hAx(4) = subplot(3,3,4);
hAx_ = hAx(4);
if length(DATA.trial_num_180) > 1
inds_cue_onset_1K = DATA.inds_cue_onset_1K(DATA.trial_num_180, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_180, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_cue_onset_1K = inds_cue_onset_1K(ind_,:);
inds_span = DATA.inds_span_cue_onset;
ESN_plot_raster(hAx_, DATA, inds_cue_onset_1K, inds_span);
end
axes(hAx_);
title('Left')
ylabel('Trial')
xlabel('Cue Onset (ms)')

hAx(5) = subplot(3,3,8);
hAx_ = hAx(5);
if length(DATA.trial_num_270) > 1
inds_cue_onset_1K = DATA.inds_cue_onset_1K(DATA.trial_num_270, :);
tgt_jump_amplitude = DATA.tgt_jump_amplitude(DATA.trial_num_270, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_cue_onset_1K = inds_cue_onset_1K(ind_,:);
inds_span = DATA.inds_span_cue_onset;
ESN_plot_raster(hAx_, DATA, inds_cue_onset_1K, inds_span);
end
axes(hAx_);
title('Down')
ylabel('Trial')
xlabel('Cue Onset (ms)')
yyaxis right;
ylabel('Firing Rate (Hz)')

sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');

ESN_Beautify_Plot(hFig, [11 8.5])
%%
fprintf(['Saving plots', ' ...'])
file_name = CH_DATA_PSORT.file_name;
file_path = CH_DATA_PSORT.file_path;
saveas(hFig,[file_path filesep file_name '_modulation_cue_onset.png'], 'png');
saveas(hFig,[file_path filesep file_name '_modulation_cue_onset.pdf'], 'pdf');
close(hFig)
fprintf(' --> Completed. \n')
%}

%% Plot Eye and Tgt raw data
clearvars -except CH_DATA_ALL CH_DATA_PSORT DATA DATA_SAC;
Line_Color = lines(7);
hFig = figure(3);
clf(hFig);

hold on
[tgt_px_1,index_tgt_px_1,~] = unique(DATA.tgt_px);
tgt_py_1 = DATA.tgt_py(index_tgt_px_1);
[tgt_py_2,index_tgt_py_2,~] = unique(DATA.tgt_py);
tgt_px_2 = DATA.tgt_px(index_tgt_py_2);

tgt_jump_px = DATA.tgt_px(DATA.event_tgt_jump_1K);
tgt_jump_py = DATA.tgt_py(DATA.event_tgt_jump_1K);
[tgt_jump_px_1,index_tgt_jump_px_1,~] = unique(tgt_jump_px);
tgt_jump_py_1 = tgt_jump_py(index_tgt_jump_px_1);
[tgt_jump_py_2,index_tgt_jump_py_2,~] = unique(tgt_jump_py);
tgt_jump_px_2 = tgt_jump_px(index_tgt_jump_py_2);

plot(tgt_px_1, tgt_py_1,'.', 'Color', Line_Color(7,:))
plot(tgt_px_2, tgt_py_2,'.', 'Color', Line_Color(7,:))
plot(tgt_jump_px_1, tgt_jump_py_1,'o', 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(tgt_jump_px_2, tgt_jump_py_2,'o', 'LineWidth', 1, 'Color', Line_Color(7,:))
range_tgt_x = max(DATA.tgt_px) - min(DATA.tgt_px);
range_tgt_y = max(DATA.tgt_py) - min(DATA.tgt_py);
xlim([(min([DATA.tgt_px(:); -5])-(0.05*range_tgt_x)), (max([DATA.tgt_px(:); +5])+(0.05*range_tgt_x))])
ylim([(min([DATA.tgt_py(:); -5])-(0.05*range_tgt_y)), (max([DATA.tgt_py(:); +5])+(0.05*range_tgt_y))])
xlabel('H Tgt (deg)')
ylabel('V Tgt (deg)')

title(CH_DATA_PSORT.file_name, 'Interpreter', 'none');

ESN_Beautify_Plot(hFig, [4 4])
%%
fprintf(['Saving plots', ' ...'])
file_name = CH_DATA_PSORT.file_name;
file_path = CH_DATA_PSORT.file_path;
saveas(hFig,[file_path filesep file_name '_modulation.png'], 'png');
saveas(hFig,[file_path filesep file_name '_modulation.pdf'], 'pdf');
close(hFig)
fprintf(' --> Completed. \n')

%% Plot Raster data
%{
Line_Color = lines(7);
% plot cue_onset
hFig = figure(1);
clf(hFig);
hAx(1) = subplot(1,2,1);
yyaxis left;
hold on
train_data_logic = DATA.event_spike_1K(DATA.inds_cue_onset_1K);
inds_span = DATA.inds_span_cue_onset;
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', Line_Color(1,:))
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trial Number')
% hAx(2) = subplot(3,2,5);
yyaxis right;
hold on
plot(inds_span, mean(train_data_logic)*1e3,'-', 'Color', Line_Color(2,:))
plot(inds_span, mean(DATA.firing_rate_1K(DATA.inds_cue_onset_1K)),'-', 'LineWidth', 2, 'Color', Line_Color(7,:))
% xlim([min(inds_span)-1 max(inds_span)+1])
ylabel('Firing rate (Hz)')
set(gca, 'YColor', Line_Color(7,:))
yyaxis left;
xlabel('Cue Onset')

hAx(2) = subplot(1,2,2);
yyaxis left;
hold on
train_data_logic = DATA.event_spike_1K(DATA.inds_sac_onset_1K);
inds_span = DATA.inds_span_sac_onset;
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', Line_Color(1,:))
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trial Number')
% hAx(4) = subplot(3,2,6);
yyaxis right;
hold on
plot(inds_span, mean(train_data_logic)*1e3,'-', 'Color', Line_Color(2,:))
plot(inds_span, mean(DATA.firing_rate_1K(DATA.inds_sac_onset_1K)),'-', 'LineWidth', 2, 'Color', Line_Color(7,:))
% xlim([min(inds_span)-1 max(inds_span)+1])
ylabel('Firing rate (Hz)')
set(gca, 'YColor', Line_Color(7,:))
yyaxis left;
xlabel('Sac Onset')

sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');

ESN_Beautify_Plot(hFig, [8.5 11])
set(hAx(1), 'XGrid', 'on', 'YGrid', 'off')
set(hAx(2), 'XGrid', 'on', 'YGrid', 'off')
% set(hAx(2), 'XGrid', 'on', 'YGrid', 'on')
% set(hAx(4), 'XGrid', 'on', 'YGrid', 'on')
fprintf(['Saving plots', ' ...'])
% file_name = CH_DATA_PSORT.file_name;
% file_path = CH_DATA_PSORT.file_path;
% saveas(hFig,[file_path filesep file_name '_modulation.png'], 'png');
% close(hFig)
fprintf(' --> Completed. \n')
%}

%% Save mat file
fprintf(['Saving mat file', ' ...'])
% file_name = CH_DATA_PSORT.file_name;
% file_path = CH_DATA_PSORT.file_path;
% save([file_path filesep file_name '_modulation.mat'], 'DATA', 'DATA_SAC', '-v7.3')
fprintf(' --> Completed. \n')
fprintf(['### ' file_name ' ###\n'])

%% Plot Raw Data
%{
clearvars h_axes
hFig = figure(2);
clf(hFig)
h_axes(1) = subplot(3,1,1);
hold on
plot(DATA.time_50K, DATA.ch_data, 'k');
plot(DATA.time_50K(DATA.ss_index), DATA.ch_data(DATA.ss_index), '*b');
ylabel('Ephys')

sample_freq = 50e3;
norm_dist_sigma = 2.5e-3; % seconds
norm_dist_x = -(4*norm_dist_sigma*sample_freq):1:(4*norm_dist_sigma*sample_freq);
norm_dist_y = normpdf(norm_dist_x,0, norm_dist_sigma*sample_freq);
time_50K = DATA.time_50K;
time_1K = DATA.time_1K;
event_spike_50K = DATA.ss_index;
time_spike_50K = time_50K(event_spike_50K);
instant_firing_rate = 1./diff(time_spike_50K);
instant_firing_rate_50K = zeros(size(time_50K));
index_spike_50K = find(event_spike_50K);
index_spike_50K = [1; index_spike_50K(:); length(event_spike_50K)];
instant_firing_rate = [0; instant_firing_rate(:); 0];
for counter_spike = 1 : 1 : length(index_spike_50K)-1
    ind_start = index_spike_50K(counter_spike);
    ind_finish = index_spike_50K(counter_spike+1);
    instant_firing_rate_50K(ind_start:ind_finish) = instant_firing_rate(counter_spike);
end
instant_firing_rate_50K_conv = conv(instant_firing_rate_50K, norm_dist_y, 'same');
instant_firing_rate_1K_conv = interp1(time_50K, instant_firing_rate_50K_conv, time_1K);

h_axes(2) = subplot(3,1,2);
hold on
plot(DATA.time_50K, instant_firing_rate_50K_conv, 'k');
ylabel('Firing Rate')

h_axes(3) = subplot(3,1,3);
hold on
plot(DATA.time_1K, DATA.eye_vm_filt, 'k');
plot(DATA.time_1K(DATA.event_sac_1K), DATA.eye_vm_filt(DATA.event_sac_1K), '*r');
plot(DATA.time_1K(DATA.event_tgt_jump_1K), DATA.eye_vm_filt(DATA.event_tgt_jump_1K), 'og');
ylabel('Eye Velocity')

linkaxes(h_axes,'x');
%}

%% Plot Eye and Tgt data 3 row version
%{
clearvars h_axes
clf(figure(11))
h_axes(1) = subplot(3,1,1);
hold on
plot(DATA.time_1K( DATA.tgt_visible), DATA.tgt_px_filt( DATA.tgt_visible),'.')
plot(DATA.time_1K, DATA.eye_px_filt, 'LineWidth', 1)
plot(DATA.time_1K(DATA.event_tgt_jump_1K), DATA.tgt_px_filt(DATA.event_tgt_jump_1K), '*')
% plot(DATA.time_1K(DATA.target_flash_45ms), DATA.tgt_px_filt(DATA.target_flash_45ms), 'o')
plot(DATA.time_1K(DATA_SAC.ind_start(DATA_SAC.validity)), DATA.eye_px_filt(DATA_SAC.ind_start(DATA_SAC.validity)), 'o')
ylabel('H Target (deg)')
range_tgt_x = max(DATA.tgt_px) - min(DATA.tgt_px);
ylim([(min([DATA.tgt_px(:); -5])-(0.05*range_tgt_x)), (max([DATA.tgt_px(:); +5])+(0.05*range_tgt_x))])


h_axes(2) = subplot(3,1,2);
hold on
plot(DATA.time_1K( DATA.tgt_visible), DATA.tgt_py_filt( DATA.tgt_visible),'.')
% plot(DATA.time_1K(~DATA.tgt_visible), DATA.tgt_py_filt(~DATA.tgt_visible),'.r')
plot(DATA.time_1K, DATA.eye_py_filt, 'LineWidth', 1)
plot(DATA.time_1K(DATA.event_tgt_jump_1K), DATA.tgt_py_filt(DATA.event_tgt_jump_1K), '*')
% plot(DATA.time_1K(DATA_SAC.ind_start(DATA_SAC.validity)), DATA.eye_py_filt(DATA_SAC.ind_start(DATA_SAC.validity)), 'o')
ylabel('V Target (deg)')
range_tgt_y = max(DATA.tgt_py) - min(DATA.tgt_py);
ylim([(min([DATA.tgt_py(:); -5])-(0.05*range_tgt_y)), (max([DATA.tgt_py(:); +5])+(0.05*range_tgt_y))])

h_axes(3) = subplot(3,1,3);
hold on
% plot(DATA.time_1K, DATA.tgt_py_filt,'.')
% plot(DATA.time_1K, DATA.eye_py_filt)
% plot(DATA.time_1K(DATA.tgt_jump_event), DATA.tgt_py_filt(DATA.tgt_jump_event), 'o')
plot(DATA.time_1K, DATA.eye_vm_filt, '-', 'LineWidth', 1)
% plot(DATA.time_1K(DATA_SAC.ind_start(DATA_SAC.validity)), DATA.eye_vm_filt(DATA_SAC.ind_start(DATA_SAC.validity)), 'o')
ylabel('Eye Velocity (deg/s)')
linkaxes(h_axes,'x');
ESN_Beautify_Plot

% plot-12
hFig = figure(12);
clf(hFig);
z_range = [0 250];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
trial_num_dir = DATA.trial_num_180;
inds_alignment_type = DATA.inds_cue_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_cue_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.firing_rate_1K(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Cue Onset')
title('Left')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
trial_num_dir = DATA.trial_num_180;
inds_alignment_type = DATA.inds_sac_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_sac_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.firing_rate_1K(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Sac Onset')
title('Left')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
trial_num_dir = DATA.trial_num_000;
inds_alignment_type = DATA.inds_cue_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_cue_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.firing_rate_1K(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Cue Onset')
title('Right')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
trial_num_dir = DATA.trial_num_000;
inds_alignment_type = DATA.inds_sac_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_sac_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.firing_rate_1K(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Sac Onset')
title('Right')
colormap hot

sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');


% plot-13
hFig = figure(13);
clf(hFig);
z_range = [0 600];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
trial_num_dir = DATA.trial_num_180;
inds_alignment_type = DATA.inds_cue_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_cue_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.eye_vm_filt(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Cue Onset')
title('Left')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
trial_num_dir = DATA.trial_num_180;
inds_alignment_type = DATA.inds_sac_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_sac_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.eye_vm_filt(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Sac Onset')
title('Left')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
trial_num_dir = DATA.trial_num_000;
inds_alignment_type = DATA.inds_cue_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_cue_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.eye_vm_filt(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Cue Onset')
title('Right')
colormap hot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
trial_num_dir = DATA.trial_num_000;
inds_alignment_type = DATA.inds_sac_onset_1K(trial_num_dir, :);
inds_span = DATA.inds_span_sac_onset;

tgt_jump_amplitude = DATA.tgt_jump_amplitude(trial_num_dir, :);
[~, ind_] = sort(tgt_jump_amplitude);
inds_alignment_type = inds_alignment_type(ind_,:);
% train_data_logic = DATA.event_spike_1K(inds_alignment_type);
firing_rate = DATA.eye_vm_filt(inds_alignment_type);
% firing_rate = firing_rate - repmat(mean(firing_rate(1:100), 2), 1, size(firing_rate,2));
h_plot = surf(inds_span, 1:size(firing_rate,1), firing_rate);
view(2)
zlim(z_range)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(firing_rate,1)+3)])
h_plot.EdgeColor = 'none';
xlabel('Sac Onset')
title('Right')
colormap hot

sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');


%}

%% function ESN_arrange_raster_plot_data
function [x_axis, y_axis] = ESN_arrange_raster_plot_data(train_data_logic, x_axis_values, line_half_len)
if nargin < 2
    x_axis_values = 1 : 1 : size(train_data_logic, 2);
    line_half_len = 0.5;
end
if nargin < 3
    line_half_len = 0.5;
end
train_data_logic = train_data_logic > 0.1;
train_data_row_number = nan(size(train_data_logic));
for counter_row = 1 : size(train_data_logic, 1)
    train_data_row_number(counter_row, train_data_logic(counter_row,:)) = counter_row;
end
train_data_col_number = repmat(x_axis_values(:)', size(train_data_logic,1), 1);
x_axis = [train_data_col_number(:)'; train_data_col_number(:)'; nan(length(train_data_col_number(:)), 1)'];
y_axis = [(train_data_row_number(:)-line_half_len)'; (train_data_row_number(:)+line_half_len)'; nan(length(train_data_row_number(:)), 1)'];
x_axis = x_axis(:);
y_axis = y_axis(:);
end

%% function ESN_arrange_raster_plot_data
function ESN_plot_raster(hAx, DATA, inds_event, inds_span)
raster_data_spike      = DATA.event_spike_1K(     inds_event);
raster_data_sac_start  = DATA.event_sac_start_1K( inds_event);
raster_data_sac_finish = DATA.event_sac_finish_1K(inds_event);
raster_data_sac_vmax   = DATA.event_sac_vmax_1K(  inds_event);
firing_rate            = DATA.firing_rate_1K(     inds_event);

Line_Color = lines(7);
% hAx(1) = subplot(3,3,5);
axes(hAx)
yyaxis left;
hold on
% train_data_logic = DATA.event_spike_1K(DATA.inds_sac_onset_1K);
% inds_span = DATA.inds_span_sac_onset;

[x_axis_data, y_axis_data] = ESN_arrange_raster_plot_data(raster_data_spike,      inds_span, 0.5);
plot(x_axis_data(:), y_axis_data(:), 'LineWidth', 2, 'Color', Line_Color(1,:))

[x_axis_data, y_axis_data] = ESN_arrange_raster_plot_data(raster_data_sac_start,  inds_span, 1);
plot(x_axis_data(:), y_axis_data(:), 'LineWidth', 3, 'Color', Line_Color(3,:))

[x_axis_data, y_axis_data] = ESN_arrange_raster_plot_data(raster_data_sac_finish, inds_span, 1);
plot(x_axis_data(:), y_axis_data(:), 'LineWidth', 3, 'Color', Line_Color(4,:))

[x_axis_data, y_axis_data] = ESN_arrange_raster_plot_data(raster_data_sac_vmax,   inds_span, 1);
plot(x_axis_data(:), y_axis_data(:), 'LineWidth', 3, 'Color', Line_Color(2,:))

xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(raster_data_spike,1)+3)])
set(gca, 'YColor', [0 0 0.5])
% ylabel('Trial Number')
% hAx(4) = subplot(3,2,6);
yyaxis right;
hold on
% plot(inds_span, mean(train_data_logic)*1e3,'-', 'Color', Line_Color(2,:))
plot(inds_span, mean(firing_rate),'-', 'LineWidth', 2, 'Color', Line_Color(7,:))
% xlim([min(inds_span)-1 max(inds_span)+1])
% ylabel('Firing rate (Hz)')
set(gca, 'YColor', Line_Color(7,:))
ylim([0 200])
yyaxis left;
% xlabel('Sac Onset')
end
