function ESN_event_alignment_v4
%% clear
clear;

%% load EPHYS EVENT DATA
file_path = pwd;
[file_name,file_path] = uigetfile([file_path filesep 'all_channels.events'], 'Select all_channels.events file');
fprintf(['Loading ', file_name, ' ... ']);
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path filesep file_name]);
EPHYS.CH_EVE.ch_data = ch_data;
EPHYS.CH_EVE.ch_time = ch_time;
EPHYS.CH_EVE.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = file_path;
fprintf(' --> Completed. \n')

%% load ADC1 and ADC2
[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*_ADC1.continuous'], 'Select ADC1.continuous file');
fprintf(['Loading ', file_name, ' ... ']);
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path filesep file_name]);
EPHYS.ADC1.ch_data = ch_data;
EPHYS.ADC1.ch_time = ch_time;
EPHYS.ADC1.ch_info = ch_info;
fprintf(' --> Completed. \n')

[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*_ADC2.continuous'], 'Select ADC2.continuous file');
fprintf(['Loading ', file_name, ' ... ']);
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path filesep file_name]);
EPHYS.ADC2.ch_data = ch_data;
EPHYS.ADC2.ch_time = ch_time;
EPHYS.ADC2.ch_info = ch_info;
fprintf(' --> Completed. \n')

EPHYS.time_30K   = EPHYS.ADC1.ch_time;

%% load BEHAVE DATA
[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep 'random_corrective_saccades_*.mat'], 'Select random_corrective_saccades file');
fprintf(['Loading ', file_name, ' ... ']);
data = load([file_path filesep file_name],'data');
BEHAVE = data.data;
fprintf(' --> Completed. \n')

%% Filter BEHAVE Eye Data
clearvars -except EPHYS BEHAVE
fprintf(['Filter BEHAVE Eye Data', ' ... ']);

min_length = min([ ...
    length(BEHAVE.eyelink_time),...
    length(BEHAVE.right_horizontal_eye),...
    length(BEHAVE.right_vertical_eye),...
    length(BEHAVE.t),...
    ]);
% timeseries
inds_invalid   = false(min_length, 1);
time_eyelink   = double(BEHAVE.eyelink_time(1:min_length)');          inds_invalid = isnan(time_eyelink) | inds_invalid;
eye_r_px       = double(BEHAVE.right_horizontal_eye(1:min_length)');  inds_invalid = isnan(eye_r_px)     | inds_invalid;
eye_r_py       = double(BEHAVE.right_vertical_eye(1:min_length)');    inds_invalid = isnan(eye_r_py)     | inds_invalid;
time_tgt       = double(BEHAVE.t(1:min_length)');                     inds_invalid = isnan(time_tgt)     | inds_invalid;
% correct for the bias between time_eyelink and time_tgt
time_eyelink   = time_eyelink .* (time_tgt(end)-time_tgt(1)) ./ (time_eyelink(end)-time_eyelink(1));
time_eyelink   = time_eyelink - time_eyelink(1) + time_tgt(1);
time_1K        = (time_eyelink(1) : 0.001 : time_eyelink(end))';
% make non unique points of eye traces invalid
inds_invalid = ([false; (diff(time_eyelink)==0)]) | inds_invalid;
inds_invalid = (abs(eye_r_px) > 15.0) | inds_invalid;
inds_invalid = (abs(eye_r_py) > 15.0) | inds_invalid;
% remove invalid values
time_eyelink(inds_invalid) = [];
eye_r_px(    inds_invalid) = [];
eye_r_py(    inds_invalid) = [];
% reconstruct eye_r data
eye_r_px = interp1(time_eyelink, eye_r_px, time_1K, 'linear', 'extrap');
eye_r_py = interp1(time_eyelink, eye_r_py, time_1K, 'linear', 'extrap');
% filter params
sampling_freq = 1000.0;
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye_r data
eye_r_px_filt = filtfilt(b_butter,a_butter,eye_r_px);
eye_r_py_filt = filtfilt(b_butter,a_butter,eye_r_py);
eye_r_vx_filt = diff(eye_r_px_filt)./diff(time_1K); eye_r_vx_filt=[eye_r_vx_filt(1); eye_r_vx_filt];
eye_r_vy_filt = diff(eye_r_py_filt)./diff(time_1K); eye_r_vy_filt=[eye_r_vy_filt(1); eye_r_vy_filt];
eye_r_vm_filt = sqrt(eye_r_vx_filt.^2 + eye_r_vy_filt.^2);
eye_r_vm_filt = abs(filtfilt(b_butter,a_butter,eye_r_vm_filt));

BEHAVE.time_1K  = time_1K;
BEHAVE.eye_r_px_filt = eye_r_px_filt;
BEHAVE.eye_r_py_filt = eye_r_py_filt;
BEHAVE.eye_r_vx_filt = eye_r_vx_filt;
BEHAVE.eye_r_vy_filt = eye_r_vy_filt;
BEHAVE.eye_r_vm_filt = eye_r_vm_filt;

fprintf(' --> Completed. \n')

%% Filter EPHYS Eye Data
clearvars -except EPHYS BEHAVE
fprintf(['Filter EPHYS Eye Data', ' ... ']);
% timeseries
inds_invalid   = false(length(EPHYS.time_30K), 1);
time_30K       = double(EPHYS.time_30K);
time_ephys     = double(EPHYS.time_30K);             inds_invalid = isnan(time_ephys)   | inds_invalid;
eye_r_px       = double(EPHYS.ADC2.ch_data) * 5.0;   inds_invalid = isnan(eye_r_px)     | inds_invalid;
eye_r_py       = double(EPHYS.ADC1.ch_data) * 5.0;   inds_invalid = isnan(eye_r_py)     | inds_invalid;
% make non unique points of eye traces invalid
inds_invalid = ([false; (diff(time_ephys)==0)]) | inds_invalid;
inds_invalid = (abs(eye_r_px) > 15.0) | inds_invalid;
inds_invalid = (abs(eye_r_py) > 15.0) | inds_invalid;
% remove invalid values
time_ephys(  inds_invalid) = [];
eye_r_px(    inds_invalid) = [];
eye_r_py(    inds_invalid) = [];
% reconstruct eye_r data
eye_r_px = interp1(time_ephys, eye_r_px, time_30K, 'linear', 'extrap');
eye_r_py = interp1(time_ephys, eye_r_py, time_30K, 'linear', 'extrap');
% filter params
sampling_freq = 30000.0;
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye_r data
eye_r_px_filt = filtfilt(b_butter,a_butter,eye_r_px);
eye_r_py_filt = filtfilt(b_butter,a_butter,eye_r_py);
eye_r_vx_filt = diff(eye_r_px_filt)./diff(time_30K); eye_r_vx_filt=[eye_r_vx_filt(1); eye_r_vx_filt];
eye_r_vy_filt = diff(eye_r_py_filt)./diff(time_30K); eye_r_vy_filt=[eye_r_vy_filt(1); eye_r_vy_filt];
eye_r_vm_filt = sqrt(eye_r_vx_filt.^2 + eye_r_vy_filt.^2);
eye_r_vm_filt = abs(filtfilt(b_butter,a_butter,eye_r_vm_filt));

EPHYS.time_30K  = time_30K;
EPHYS.eye_r_px_filt = eye_r_px_filt;
EPHYS.eye_r_py_filt = eye_r_py_filt;
EPHYS.eye_r_vx_filt = eye_r_vx_filt;
EPHYS.eye_r_vy_filt = eye_r_vy_filt;
EPHYS.eye_r_vm_filt = eye_r_vm_filt;

fprintf(' --> Completed. \n')

%% extract event data in ephys
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Event Data', ' ... ']);
% ch_data
% 2 : STR_TARGET_PURSUIT 0
% 3 : STR_TARGET_FIXATION 1
% 4 : DETECT_SACCADE_START 2
% - : DETECT_SACCADE_END 3
% 4 : END_TARGET_FIXATION 4
% 5 : 1Hz ossilation 5
% 6 : photodiode: STR_TARGET_FIXATION+DETECT_SACCADE_END
% eventId
% 1 : rising
% 0 : falling
EPHYS.CH_EVE.data = [EPHYS.CH_EVE.ch_time(:) EPHYS.CH_EVE.ch_data(:) EPHYS.CH_EVE.ch_info.eventId(:)];
fprintf(' --> Completed. \n')

%% Build BEHAVE Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE Alignment events', ' ... ']);

time_reference      = BEHAVE.time_1K;
length_time         = length(time_reference);
time_state_str_fixation = [];
time_state_sac_detect_on = [];
time_state_sac_detect_off = [];
time_state_end_fixation = [];
time_state_iti = [];
for counter_trial = 1 : length(BEHAVE.trials)-1
    time_state_str_fixation = [time_state_str_fixation; BEHAVE.trials{counter_trial}.state_start_time_str_target_fixation(:)];
    time_state_sac_detect_on = [time_state_sac_detect_on; BEHAVE.trials{counter_trial}.state_start_time_detect_sac_start(:)];
    time_state_sac_detect_off = [time_state_sac_detect_off; BEHAVE.trials{counter_trial}.state_start_time_detect_sac_end(:)];
    time_state_end_fixation = [time_state_end_fixation; BEHAVE.trials{counter_trial}.state_start_time_end_target_fixation(:)];
    time_state_iti = [time_state_iti; BEHAVE.trials{counter_trial}.state_start_time_iti(:)];
end

variable_list = {'_state_str_fixation','_state_sac_detect_on','_state_sac_detect_off','_state_end_fixation','_state_iti'};

BEHAVE.Alignment.time_1K = time_reference;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'BEHAVE.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
end
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_state_str_fixation(  counter_state_str_fixation)
        event_state_str_fixation(    counter_time_point) = true;
        counter_state_str_fixation   = counter_state_str_fixation   + 1;
    end
    if time_ponit_ >= time_state_sac_detect_on( counter_state_sac_detect_on)
        event_state_sac_detect_on(   counter_time_point) = true;
        counter_state_sac_detect_on  = counter_state_sac_detect_on  + 1;
    end
    if time_ponit_ >= time_state_sac_detect_off(counter_state_sac_detect_off)
        event_state_sac_detect_off(  counter_time_point) = true;
        counter_state_sac_detect_off = counter_state_sac_detect_off + 1;
    end
    if time_ponit_ >= time_state_end_fixation(  counter_state_end_fixation)
        event_state_end_fixation(    counter_time_point) = true;
        counter_state_end_fixation   = counter_state_end_fixation   + 1;
    end
    if time_ponit_ >= time_state_iti(           counter_state_iti)
        event_state_iti(             counter_time_point) = true;
        counter_state_iti            = counter_state_iti            + 1;
    end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'BEHAVE.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
    length_time = length(event_temp_);
    inds_span_ = (0) : 1 : (100);
    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name '_100' ' = ' 'event_temp_'   ';']);
end

event_state_combined = ...
    double(event_state_str_fixation_100)   .* 1 + ...
    double(event_state_sac_detect_on_100)  .* 2 + ...
    double(event_state_sac_detect_off_100) .* 3;% + ...
%     double(event_state_end_fixation_100)   .* 4 + ...
%     double(event_state_iti_100)            .* 5 ;
event_state_combined(event_state_combined>3) = 3;
event_photodiode_combined = ...
    double(event_state_str_fixation_100)   .* 1 + ...
    double(event_state_sac_detect_off_100) .* 1 ;%+ ...
    %double(event_state_sac_detect_on_100)  .* 2 + ...
    %double(event_state_end_fixation_100)   .* 2 ;

BEHAVE.Alignment.event_state_combined      = event_state_combined;
BEHAVE.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'state_str_fixation: 1 , ', ...
'state_sac_detect_on: 2 , ', ...
'state_sac_detect_off: 3 , '];%, ...
% 'state_end_fixation: 4 , ', ...
% 'state_iti: 5', ...
% ];
photodiode_description = [
'state_str_fixation: 1 , ', ...
'state_sac_detect_off: 1 , ', ...
];
BEHAVE.Alignment.state_description      = state_description;
BEHAVE.Alignment.photodiode_description = photodiode_description;

BEHAVE.Alignment.eye_r_vm_filt = interp1(BEHAVE.time_1K, BEHAVE.eye_r_vm_filt, time_reference, 'linear', 'extrap');

fprintf(' --> Completed. \n');

%% Build EPHYS Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Alignment events', ' ... ']);

time_reference      = ( ESN_Round(EPHYS.time_30K(1),0.001) : 0.001 : ESN_Round(EPHYS.time_30K(end),0.001) )';
length_time         = length(time_reference);
time_state_str_fixation   = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_state_sac_detect_on  = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_state_sac_detect_off = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
% time_state_end_fixation   = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
% time_state_iti            = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_photodiode_rise      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_photodiode_fall      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
variable_list = {'_state_str_fixation','_state_sac_detect_on','_state_sac_detect_off','_photodiode_rise','_photodiode_fall'};

if isempty(time_photodiode_rise)
    time_photodiode_rise = [time_state_str_fixation; time_state_sac_detect_off];
    time_photodiode_rise = sort(time_photodiode_rise);
end

if isempty(time_photodiode_fall)
    time_photodiode_fall = [time_state_str_fixation; time_state_sac_detect_off];
    time_photodiode_fall = sort(time_photodiode_fall);
end

EPHYS.Alignment.time_1K = time_reference;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
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
    if time_ponit_ >= time_state_str_fixation(  counter_state_str_fixation)
        event_state_str_fixation(    counter_time_point) = true;
        counter_state_str_fixation   = counter_state_str_fixation   + 1;
    end
    if time_ponit_ >= time_state_sac_detect_on( counter_state_sac_detect_on)
        event_state_sac_detect_on(   counter_time_point) = true;
        counter_state_sac_detect_on  = counter_state_sac_detect_on  + 1;
    end
    if time_ponit_ >= time_state_sac_detect_off(counter_state_sac_detect_off)
        event_state_sac_detect_off(  counter_time_point) = true;
        counter_state_sac_detect_off = counter_state_sac_detect_off + 1;
    end
%     if time_ponit_ >= time_state_end_fixation(  counter_state_end_fixation)
%         event_state_end_fixation(    counter_time_point) = true;
%         counter_state_end_fixation   = counter_state_end_fixation   + 1;
%     end
%     if time_ponit_ >= time_state_iti(           counter_state_iti)
%         event_state_iti(             counter_time_point) = true;
%         counter_state_iti            = counter_state_iti            + 1;
%     end
    if time_ponit_ >= time_photodiode_rise(     counter_photodiode_rise)
        event_photodiode_rise(       counter_time_point) = true;
        counter_photodiode_rise      = counter_photodiode_rise      + 1;
    end
    if time_ponit_ >= time_photodiode_fall(     counter_photodiode_fall)
        event_photodiode_fall(       counter_time_point) = true;
        counter_photodiode_fall      = counter_photodiode_fall      + 1;
    end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
    length_time = length(event_temp_);
    inds_span_ = (0) : 1 : (100);
    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name '_100' ' = ' 'event_temp_'   ';']);
end

event_state_combined = ...
    double(event_state_str_fixation_100)   .* 1 + ...
    double(event_state_sac_detect_on_100)  .* 2 + ...
    double(event_state_sac_detect_off_100) .* 3 ;%+ ...
%     double(event_state_end_fixation_100)   .* 4 + ...
%     double(event_state_iti_100)            .* 5 ;
event_state_combined(event_state_combined>3) = 3;
event_photodiode_combined = ...
    double(event_photodiode_rise_100)  .* 1;% + ...
    %double(event_photodiode_fall_100)  .* 2 ;

EPHYS.Alignment.event_state_combined      = event_state_combined;
EPHYS.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'state_str_fixation: 1 , ', ...
'state_sac_detect_on: 2 , ', ...
'state_sac_detect_off: 3 , '];%, ...
% 'state_end_fixation: 4 , ', ...
% 'state_iti: 5', ...
% ];
photodiode_description = [
'state_str_fixation: 1 , ', ...
'state_sac_detect_off: 1 , ', ...
];
EPHYS.Alignment.state_description      = state_description;
EPHYS.Alignment.photodiode_description = photodiode_description;

EPHYS.Alignment.eye_r_vm_filt = interp1(EPHYS.time_30K, EPHYS.eye_r_vm_filt, time_reference, 'linear', 'extrap');

fprintf(' --> Completed. \n');


%% ALIGN EPHYS and BEHAVE state_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE state_combined', ' ... ']);
EPHYS_time_1K              = EPHYS.Alignment.time_1K;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_1K             = BEHAVE.Alignment.time_1K;
EPHYS_eye_r_vm_filt  = EPHYS.Alignment.eye_r_vm_filt;
BEHAVE_eye_r_vm_filt = BEHAVE.Alignment.eye_r_vm_filt;
EPHYS_eye_r_vm_filt(EPHYS_eye_r_vm_filt < 100) = 0;
BEHAVE_eye_r_vm_filt(BEHAVE_eye_r_vm_filt < 100) = 0;
EPHYS_eye_r_vm_filt(EPHYS_eye_r_vm_filt > 1500) = 1500;
BEHAVE_eye_r_vm_filt(BEHAVE_eye_r_vm_filt > 1500) = 1500;
% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_eye_r_vm_filt, BEHAVE_eye_r_vm_filt); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

EPHYS_state_combined       = EPHYS.Alignment.event_state_combined;
BEHAVE_state_combined      = BEHAVE.Alignment.event_state_combined;

if  sample_diff > 0
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_state_combined(  abs(sample_diff):end);
    EPHYS_EB_xcorr_time_1K            = EPHYS_time_1K(            abs(sample_diff):end);
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_state_combined;
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_time_1K;
elseif sample_diff < 0
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_state_combined;
    EPHYS_EB_xcorr_time_1K            = EPHYS_time_1K;
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_state_combined( abs(sample_diff):end);
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_time_1K(           abs(sample_diff):end);
end
% state_combined: make the vectors the same size
if length(BEHAVE_EB_xcorr_state_combined_1K) ~= length(EPHYS_EB_xcorr_state_combined_1K)
    min_length = min([ length(BEHAVE_EB_xcorr_state_combined_1K),  length(EPHYS_EB_xcorr_state_combined_1K)]);
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_EB_xcorr_state_combined_1K(  1:min_length);
    EPHYS_EB_xcorr_time_1K            = EPHYS_EB_xcorr_time_1K(            1:min_length);
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_EB_xcorr_state_combined_1K( 1:min_length);
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_EB_xcorr_time_1K(           1:min_length);
end

% accounting for INCORRECT_SACCADE which is recorded in EPHYS but is
% missing in BEHAVE
EPHYS_EB_xcorr_state_combined_1K((EPHYS_EB_xcorr_state_combined_1K == 3) & (BEHAVE_EB_xcorr_state_combined_1K == 0)) = 0;

% %% state_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_EB_xcorr_state_combined_1K_filt  = filtfilt(b_butter,a_butter,EPHYS_EB_xcorr_state_combined_1K);
BEHAVE_EB_xcorr_state_combined_1K_filt = filtfilt(b_butter,a_butter,BEHAVE_EB_xcorr_state_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_EB_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(EPHYS_EB_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_EB_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_EB_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_EB_state_combined_chunk  = EPHYS_EB_xcorr_state_combined_1K_filt(inds_chunk);
    BEHAVE_EB_state_combined_chunk = BEHAVE_EB_xcorr_state_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_EB_state_combined_chunk,BEHAVE_EB_state_combined_chunk, 15, 'absolute');  % allow upto 15ms warp
    EPHYS_EB_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_EB_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_EB_inds_DTW  = cell2mat(EPHYS_EB_inds_DTW);
BEHAVE_EB_inds_DTW = cell2mat(BEHAVE_EB_inds_DTW);
EPHYS_EB_inds      = ( 1 : 1 : length(EPHYS_EB_xcorr_time_1K ) )';
BEHAVE_EB_inds     = ( 1 : 1 : length(BEHAVE_EB_xcorr_time_1K) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
EB_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_EB_inds));
EB_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_EB_inds));
for counter_ind = 1 : 1 : length(EPHYS_EB_inds_DTW)
    ind_EPHYS_EB_DTW  = EPHYS_EB_inds_DTW(counter_ind);
    ind_BEHAVE_EB_DTW = BEHAVE_EB_inds_DTW(counter_ind);
    EB_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_EB_DTW) = ind_EPHYS_EB_DTW;
    EB_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_EB_DTW)  = ind_BEHAVE_EB_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_1K = EPHYS_EB_xcorr_time_1K;
time_EPHYS_EB_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_1K(end)])+1;
event_EPHYS_EB_xcorr_30K       = nan(length(EPHYS_EB_xcorr_time_1K), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_1K(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_30K(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_1K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_1K = EPHYS_EB_xcorr_time_1K;
time_EPHYS_EB_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_1K(end)])+1;
event_EPHYS_EB_xcorr_1K       = nan(length(EPHYS_EB_xcorr_time_1K), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_1K(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_1K(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_1K(:);
length_time         = length(time_reference);
time_BEHAVE_EB_xcorr_1K = BEHAVE_EB_xcorr_time_1K;
time_BEHAVE_EB_xcorr_1K(end+1) = max([time_reference(end), time_BEHAVE_EB_xcorr_1K(end)])+1;
event_BEHAVE_EB_xcorr_1K       = nan(length(BEHAVE_EB_xcorr_time_1K), 1);
counter_BEHAVE_EB_xcorr     = find(time_BEHAVE_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_EB_xcorr_1K(  counter_BEHAVE_EB_xcorr)
        event_BEHAVE_EB_xcorr_1K(    counter_BEHAVE_EB_xcorr) = counter_time_point;
        counter_BEHAVE_EB_xcorr   = counter_BEHAVE_EB_xcorr   + 1;
    end
end

EPHYS_EB_xcorr_ind_30K   = event_EPHYS_EB_xcorr_30K;
EPHYS_EB_xcorr_ind_1K    = event_EPHYS_EB_xcorr_1K;
BEHAVE_EB_xcorr_ind_1K   = event_BEHAVE_EB_xcorr_1K;
EPHYS_EB_aligned_ind_30K = EPHYS_EB_xcorr_ind_30K(EB_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_EB_aligned_ind_1K  = EPHYS_EB_xcorr_ind_1K( EB_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_EB_aligned_ind_1K = BEHAVE_EB_xcorr_ind_1K(EB_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_30K          = EPHYS_EB_aligned_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K           = EPHYS_EB_aligned_ind_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_aligned_ind_1K          = BEHAVE_EB_aligned_ind_1K;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_BEHAVE_to_EPHYS = EB_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_EPHYS_to_BEHAVE = EB_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K            = EPHYS_EB_xcorr_time_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K           = BEHAVE_EB_xcorr_time_1K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_30K            = EPHYS_EB_xcorr_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_1K             = EPHYS_EB_xcorr_ind_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K            = BEHAVE_EB_xcorr_ind_1K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_state_combined_1K  = EPHYS_EB_xcorr_state_combined_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_EB_xcorr_state_combined_1K;
% EPHYS.CH_EVE.align_states.state_description            = BEHAVE.Alignment.state_description;
fprintf(' --> Completed. \n');

%% ALIGN EPHYS and BEHAVE photodiode_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE photodiode_combined', ' ... ']);
EPHYS_time_1K              = EPHYS.Alignment.time_1K;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_1K             = BEHAVE.Alignment.time_1K;
EPHYS_photodiode_combined  = EPHYS.Alignment.event_photodiode_combined;
BEHAVE_photodiode_combined = BEHAVE.Alignment.event_photodiode_combined;

% photodiode_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_photodiode_combined+1, BEHAVE_photodiode_combined+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);
if sample_diff > 0
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_photodiode_combined(  abs(sample_diff):end);
    EPHYS_PD_xcorr_time_1K                 = EPHYS_time_1K(                 abs(sample_diff):end);
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_photodiode_combined;
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K;
elseif sample_diff < 0
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_photodiode_combined;
    EPHYS_PD_xcorr_time_1K                 = EPHYS_time_1K;
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_photodiode_combined( abs(sample_diff):end);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K(                abs(sample_diff):end);
end
% photodiode_combined: make the vectors the same size
if length(BEHAVE_PD_xcorr_photodiode_combined_1K) ~= length(EPHYS_PD_xcorr_photodiode_combined_1K)
    min_length = min([ length(BEHAVE_PD_xcorr_photodiode_combined_1K),  length(EPHYS_PD_xcorr_photodiode_combined_1K)]);
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_PD_xcorr_photodiode_combined_1K(  1:min_length);
    EPHYS_PD_xcorr_time_1K                 = EPHYS_PD_xcorr_time_1K(                 1:min_length);
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_PD_xcorr_photodiode_combined_1K( 1:min_length);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_1K(                1:min_length);
end

% photodiode_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_PD_xcorr_photodiode_combined_1K_filt  = filtfilt(b_butter,a_butter,EPHYS_PD_xcorr_photodiode_combined_1K);
BEHAVE_PD_xcorr_photodiode_combined_1K_filt = filtfilt(b_butter,a_butter,BEHAVE_PD_xcorr_photodiode_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_PD_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(EPHYS_PD_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_PD_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_PD_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_PD_photodiode_combined_chunk  = EPHYS_PD_xcorr_photodiode_combined_1K_filt(inds_chunk);
    BEHAVE_PD_photodiode_combined_chunk = BEHAVE_PD_xcorr_photodiode_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_PD_photodiode_combined_chunk,BEHAVE_PD_photodiode_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    EPHYS_PD_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_PD_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_PD_inds_DTW  = cell2mat(EPHYS_PD_inds_DTW);
BEHAVE_PD_inds_DTW = cell2mat(BEHAVE_PD_inds_DTW);
EPHYS_PD_inds      = ( 1 : 1 : length(EPHYS_PD_xcorr_time_1K ) )';
BEHAVE_PD_inds     = ( 1 : 1 : length(BEHAVE_PD_xcorr_time_1K) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
PD_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_PD_inds));
PD_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_PD_inds));
for counter_ind = 1 : 1 : length(EPHYS_PD_inds_DTW)
    ind_EPHYS_PD_DTW  = EPHYS_PD_inds_DTW(counter_ind);
    ind_BEHAVE_PD_DTW = BEHAVE_PD_inds_DTW(counter_ind);
    PD_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_PD_DTW) = ind_EPHYS_PD_DTW;
    PD_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_PD_DTW)  = ind_BEHAVE_PD_DTW;
end


time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_1K = EPHYS_PD_xcorr_time_1K;
time_EPHYS_PD_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_1K(end)])+1;
event_EPHYS_PD_xcorr_30K       = nan(length(EPHYS_PD_xcorr_time_1K), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_1K(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_30K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_1K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_1K = EPHYS_PD_xcorr_time_1K;
time_EPHYS_PD_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_1K(end)])+1;
event_EPHYS_PD_xcorr_1K       = nan(length(EPHYS_PD_xcorr_time_1K), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_1K(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_1K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_1K(:);
length_time         = length(time_reference);
time_BEHAVE_PD_xcorr_1K = BEHAVE_PD_xcorr_time_1K;
time_BEHAVE_PD_xcorr_1K(end+1) = max([time_reference(end), time_BEHAVE_PD_xcorr_1K(end)])+1;
event_BEHAVE_PD_xcorr_1K       = nan(length(BEHAVE_PD_xcorr_time_1K), 1);
counter_BEHAVE_PD_xcorr     = find(time_BEHAVE_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_PD_xcorr_1K(  counter_BEHAVE_PD_xcorr)
        event_BEHAVE_PD_xcorr_1K(    counter_BEHAVE_PD_xcorr) = counter_time_point;
        counter_BEHAVE_PD_xcorr   = counter_BEHAVE_PD_xcorr   + 1;
    end
end

EPHYS_PD_xcorr_ind_30K   = event_EPHYS_PD_xcorr_30K;
EPHYS_PD_xcorr_ind_1K    = event_EPHYS_PD_xcorr_1K;
BEHAVE_PD_xcorr_ind_1K   = event_BEHAVE_PD_xcorr_1K;
EPHYS_PD_aligned_ind_30K = EPHYS_PD_xcorr_ind_30K(PD_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_PD_aligned_ind_1K  = EPHYS_PD_xcorr_ind_1K( PD_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_PD_aligned_ind_1K = BEHAVE_PD_xcorr_ind_1K(PD_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_30K               = EPHYS_PD_aligned_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K                = EPHYS_PD_aligned_ind_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_aligned_ind_1K               = BEHAVE_PD_aligned_ind_1K;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_BEHAVE_to_EPHYS = PD_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_EPHYS_to_BEHAVE = PD_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_time_1K                 = EPHYS_PD_xcorr_time_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_1K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_30K                 = EPHYS_PD_xcorr_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_1K                  = EPHYS_PD_xcorr_ind_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_ind_1K                 = BEHAVE_PD_xcorr_ind_1K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_PD_xcorr_photodiode_combined_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_PD_xcorr_photodiode_combined_1K;
EPHYS.CH_EVE.align_photodiode.photodiode_description              = BEHAVE.Alignment.photodiode_description;
fprintf(' --> Completed. \n');

%% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE
EPHYS_time_1K    = EPHYS.Alignment.time_1K;
EPHYS_time_30K   = EPHYS.time_30K;
BEHAVE_time_1K   = BEHAVE.Alignment.time_1K;
align_photodiode = EPHYS.CH_EVE.align_photodiode;
align_states     = EPHYS.CH_EVE.align_states;

file_name = EPHYS.file_name_CH_EVE;
file_path = EPHYS.file_path_CH_EVE;
[~, file_name, ~] = fileparts(file_name);
file_name = [file_name '_EVE1_aligned.mat'];
clearvars EPHYS BEHAVE
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([file_path filesep file_name], '-v7.3');
fprintf(' --> Completed. \n')

end
