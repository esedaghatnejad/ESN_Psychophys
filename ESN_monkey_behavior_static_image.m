function ESN_monkey_behavior_static_image
%% load BEHAVE DATA
file_path = pwd;
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[file_name,file_path] = uigetfile([file_path '*.mat'], 'Select _corrective_saccades file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
data = load([file_path file_name], 'data');
BEHAVE.data = data.data;
fprintf(' --> Completed. \n')

BEHAVE.EXPERIMENT_PARAMS.mat_FileName = file_name;
BEHAVE.EXPERIMENT_PARAMS.mat_PathName = file_path;
[~, BEHAVE.EXPERIMENT_PARAMS.file_name, ~] = fileparts(file_name);
[~, foldername, ~] = fileparts(file_path(1:end-1)); % (1:end-1) to delete the filesep character
BEHAVE.EXPERIMENT_PARAMS.folder_name = foldername;

%% load EVE1_aligned
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[file_name,file_path] = uigetfile([file_path '*' '_EVE1_aligned.mat'], 'Select _EVE1_aligned file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
EPHYS.CH_EVE = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% Build BEHAVE stream
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE stream ', ' ... ']);
BEHAVE.stream = struct;
% start/end indices
BEHAVE.stream.inds_trial          = 1 : 1 : length(BEHAVE.data.eyelink_time);
% timeseries
BEHAVE.stream.inds_invalid   = false(1, length(BEHAVE.stream.inds_trial));
BEHAVE.stream.time_eyelink   = double(BEHAVE.data.eyelink_time(1, BEHAVE.stream.inds_trial));                           BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.time_eyelink) | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_r_px       = double(BEHAVE.data.right_horizontal_eye(1, BEHAVE.stream.inds_trial));                   BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_r_px)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_r_py       = double(BEHAVE.data.right_vertical_eye(  1, BEHAVE.stream.inds_trial));                   BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_r_py)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_l_px       = double(BEHAVE.data.left_horizontal_eye( 1, BEHAVE.stream.inds_trial));                   BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_l_px)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_l_py       = double(BEHAVE.data.left_vertical_eye(   1, BEHAVE.stream.inds_trial));                   BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_l_py)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_r_vx       = double(BEHAVE.data.right_horizontal_eye_velocity_filtered(1, BEHAVE.stream.inds_trial)); BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_r_vx)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_r_vy       = double(BEHAVE.data.right_vertical_eye_velocity_filtered(  1, BEHAVE.stream.inds_trial)); BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_r_vy)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_l_vx       = double(BEHAVE.data.left_horizontal_eye_velocity_filtered( 1, BEHAVE.stream.inds_trial)); BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_l_vx)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_l_vy       = double(BEHAVE.data.left_vertical_eye_velocity_filtered(   1, BEHAVE.stream.inds_trial)); BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.eye_l_vy)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.eye_r_vm       = sqrt(BEHAVE.stream.eye_r_vx.^2 + BEHAVE.stream.eye_r_vy.^2);
BEHAVE.stream.eye_l_vm       = sqrt(BEHAVE.stream.eye_l_vx.^2 + BEHAVE.stream.eye_l_vy.^2);
BEHAVE.stream.time_tgt       = double(BEHAVE.data.t(1, BEHAVE.stream.inds_trial));                                      BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.time_tgt)     | BEHAVE.stream.inds_invalid;
BEHAVE.stream.tgt_px         = double(BEHAVE.data.target_x(1, BEHAVE.stream.inds_trial));                               BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.tgt_px)       | BEHAVE.stream.inds_invalid;
BEHAVE.stream.tgt_py         = double(BEHAVE.data.target_y(1, BEHAVE.stream.inds_trial));                               BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.tgt_py)       | BEHAVE.stream.inds_invalid;
BEHAVE.stream.reward         = double(BEHAVE.data.reward(1, BEHAVE.stream.inds_trial));                                 BEHAVE.stream.inds_invalid = isnan(BEHAVE.stream.reward)       | BEHAVE.stream.inds_invalid;
BEHAVE.stream.target_visible = logical(double(BEHAVE.data.target_visible(1, BEHAVE.stream.inds_trial)));
% correct for the bias between time_eyelink and time_tgt
BEHAVE.stream.time_eyelink   = BEHAVE.stream.time_eyelink .* (BEHAVE.stream.time_tgt(end)-BEHAVE.stream.time_tgt(1)) ./ (BEHAVE.stream.time_eyelink(end)-BEHAVE.stream.time_eyelink(1));
BEHAVE.stream.time_eyelink   = BEHAVE.stream.time_eyelink - BEHAVE.stream.time_eyelink(1) + BEHAVE.stream.time_tgt(1);
BEHAVE.stream.time_1K        = BEHAVE.stream.time_eyelink(1) : 0.001 : BEHAVE.stream.time_eyelink(end);
% make non unique points of eye traces invalid
BEHAVE.stream.inds_invalid = ([false (diff(BEHAVE.stream.time_eyelink)==0)])       | BEHAVE.stream.inds_invalid;
% remove invalid values
BEHAVE.stream.time_eyelink(BEHAVE.stream.inds_invalid) = [];
BEHAVE.stream.eye_r_px(    BEHAVE.stream.inds_invalid) = [];
BEHAVE.stream.eye_r_py(    BEHAVE.stream.inds_invalid) = [];
BEHAVE.stream.eye_l_px(    BEHAVE.stream.inds_invalid) = [];
BEHAVE.stream.eye_l_py(    BEHAVE.stream.inds_invalid) = [];
% reconstruct eye_r data
BEHAVE.stream.eye_r_px = interp1(BEHAVE.stream.time_eyelink, BEHAVE.stream.eye_r_px, BEHAVE.stream.time_1K, 'linear', 'extrap');
BEHAVE.stream.eye_r_py = interp1(BEHAVE.stream.time_eyelink, BEHAVE.stream.eye_r_py, BEHAVE.stream.time_1K, 'linear', 'extrap');
BEHAVE.stream.eye_r_vx = diff(BEHAVE.stream.eye_r_px)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_r_vx=[BEHAVE.stream.eye_r_vx(1) BEHAVE.stream.eye_r_vx];
BEHAVE.stream.eye_r_vy = diff(BEHAVE.stream.eye_r_py)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_r_vy=[BEHAVE.stream.eye_r_vy(1) BEHAVE.stream.eye_r_vy];
BEHAVE.stream.eye_r_vm = sqrt(BEHAVE.stream.eye_r_vx.^2 + BEHAVE.stream.eye_r_vy.^2);
% reconstruct eye_l data
BEHAVE.stream.eye_l_px = interp1(BEHAVE.stream.time_eyelink, BEHAVE.stream.eye_l_px, BEHAVE.stream.time_1K, 'linear', 'extrap');
BEHAVE.stream.eye_l_py = interp1(BEHAVE.stream.time_eyelink, BEHAVE.stream.eye_l_py, BEHAVE.stream.time_1K, 'linear', 'extrap');
BEHAVE.stream.eye_l_vx = diff(BEHAVE.stream.eye_l_px)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_l_vx=[BEHAVE.stream.eye_l_vx(1) BEHAVE.stream.eye_l_vx];
BEHAVE.stream.eye_l_vy = diff(BEHAVE.stream.eye_l_py)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_l_vy=[BEHAVE.stream.eye_l_vy(1) BEHAVE.stream.eye_l_vy];
BEHAVE.stream.eye_l_vm = sqrt(BEHAVE.stream.eye_l_vx.^2 + BEHAVE.stream.eye_l_vy.^2);
BEHAVE.stream.time     = BEHAVE.stream.time_1K;
% filter params
sampling_freq = 1000.0;
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye_r data
BEHAVE.stream.eye_r_px_filt = filtfilt(b_butter,a_butter,BEHAVE.stream.eye_r_px);
BEHAVE.stream.eye_r_py_filt = filtfilt(b_butter,a_butter,BEHAVE.stream.eye_r_py);
BEHAVE.stream.eye_r_vx_filt = diff(BEHAVE.stream.eye_r_px_filt)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_r_vx_filt=[BEHAVE.stream.eye_r_vx_filt(1) BEHAVE.stream.eye_r_vx_filt];
BEHAVE.stream.eye_r_vy_filt = diff(BEHAVE.stream.eye_r_py_filt)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_r_vy_filt=[BEHAVE.stream.eye_r_vy_filt(1) BEHAVE.stream.eye_r_vy_filt];
BEHAVE.stream.eye_r_vm_filt = sqrt(BEHAVE.stream.eye_r_vx_filt.^2 + BEHAVE.stream.eye_r_vy_filt.^2);
% filter eye_l data
BEHAVE.stream.eye_l_px_filt = filtfilt(b_butter,a_butter,BEHAVE.stream.eye_l_px);
BEHAVE.stream.eye_l_py_filt = filtfilt(b_butter,a_butter,BEHAVE.stream.eye_l_py);
BEHAVE.stream.eye_l_vx_filt = diff(BEHAVE.stream.eye_l_px_filt)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_l_vx_filt=[BEHAVE.stream.eye_l_vx_filt(1) BEHAVE.stream.eye_l_vx_filt];
BEHAVE.stream.eye_l_vy_filt = diff(BEHAVE.stream.eye_l_py_filt)./diff(BEHAVE.stream.time_1K); BEHAVE.stream.eye_l_vy_filt=[BEHAVE.stream.eye_l_vy_filt(1) BEHAVE.stream.eye_l_vy_filt];
BEHAVE.stream.eye_l_vm_filt = sqrt(BEHAVE.stream.eye_l_vx_filt.^2 + BEHAVE.stream.eye_l_vy_filt.^2);

fprintf(' --> Completed. \n')

%% Build BEHAVE aligned
fprintf(['Building BEHAVE aligned', ' ... ']);
variable_list = {'time_1K','eye_r_px_filt','eye_r_py_filt', 'eye_r_vx_filt', 'eye_r_vy_filt', 'eye_r_vm_filt', 'tgt_px', 'tgt_py'};
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

tgtStr_px = 0.0;
tgtStr_py = 0.0;

distance_eye_to_tgt_finish    = sqrt( (all_sac_eye_r_px_finish - all_sac_tgt_px_finish).^2 + (all_sac_eye_r_py_finish - all_sac_tgt_py_finish).^2 );
distance_eye_to_tgt_start     = sqrt( (all_sac_eye_r_px_start  - all_sac_tgt_px_start).^2  + (all_sac_eye_r_py_start  - all_sac_tgt_py_start).^2 );
distance_eye_to_tgtStr_finish = sqrt( (all_sac_eye_r_px_finish - tgtStr_px).^2 + (all_sac_eye_r_py_finish - tgtStr_py).^2 );
distance_eye_to_tgtStr_start  = sqrt( (all_sac_eye_r_px_start  - tgtStr_px).^2 + (all_sac_eye_r_py_start  - tgtStr_py).^2 );

BEHAVE.SACS_ALL_DATA.is_primSac      = false(size(BEHAVE.SACS_ALL_DATA.validity));
BEHAVE.SACS_ALL_DATA.is_corrSac      = false(size(BEHAVE.SACS_ALL_DATA.validity));
BEHAVE.SACS_ALL_DATA.is_toTgt        = (distance_eye_to_tgt_finish < 2.0);
BEHAVE.SACS_ALL_DATA.is_fromTgt      = (distance_eye_to_tgt_start < 2.0);
BEHAVE.SACS_ALL_DATA.is_lowAmp       = (BEHAVE.SACS_ALL_DATA.eye_r_amp_m < 2.0);
BEHAVE.SACS_ALL_DATA.is_toTgtStr     = (distance_eye_to_tgtStr_finish < 2.0);
BEHAVE.SACS_ALL_DATA.is_fromCenter   = (distance_eye_to_tgtStr_start < 2.0);
BEHAVE.SACS_ALL_DATA.is_aroundCenter = (distance_eye_to_tgtStr_start < 2.0) & (distance_eye_to_tgtStr_finish < 2.0);
BEHAVE.SACS_ALL_DATA.is_all          = BEHAVE.SACS_ALL_DATA.validity;
BEHAVE.SACS_ALL_DATA.is_notLowAmp    = BEHAVE.SACS_ALL_DATA.validity;

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

BEHAVE.SACS_ALL_DATA.is_notLowAmp(BEHAVE.SACS_ALL_DATA.is_lowAmp) = false;

fprintf(' --> Completed. \n')

%% Save _ANALYZED.mat Data to disk
clearvars -except EPHYS BEHAVE hFig
EXPERIMENT_PARAMS = BEHAVE.EXPERIMENT_PARAMS;
SACS_ALL_DATA = BEHAVE.SACS_ALL_DATA;
aligned = BEHAVE.aligned;
file_path_ = BEHAVE.EXPERIMENT_PARAMS.mat_PathName;
if ~strcmp(file_path_(end), filesep);file_path_ = [file_path_ filesep];end
file_name_ = [EXPERIMENT_PARAMS.file_name '_ANALYZED.mat'];
fprintf(['Saving ' file_name_ ' ...'])
save([file_path_ file_name_], ...
    'EXPERIMENT_PARAMS', 'SACS_ALL_DATA', 'aligned', '-v7.3');
fprintf(' --> Completed. \n')

end

%% function build_sac_data(SACS_ALL_DATA, aligned)
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


