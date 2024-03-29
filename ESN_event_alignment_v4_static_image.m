function ESN_event_alignment_v4_static_image
%% clear
clear;
EPHYS.debug_figures = false;

%% load ADC1 and ADC2
file_path = pwd;
[file_name,file_path] = uigetfile([file_path filesep '*_ADC1.continuous'], 'Select ADC1.continuous file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path file_name]);
EPHYS.ADC1.ch_data = ch_data;
EPHYS.ADC1.ch_time = ch_time;
EPHYS.ADC1.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = file_path;
fprintf(' --> Completed. \n')

[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*_ADC2.continuous'], 'Select ADC2.continuous file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path file_name]);
EPHYS.ADC2.ch_data = ch_data;
EPHYS.ADC2.ch_time = ch_time;
EPHYS.ADC2.ch_info = ch_info;
fprintf(' --> Completed. \n')

EPHYS.time_30K   = EPHYS.ADC1.ch_time;

%% load BEHAVE DATA
[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.mat'], 'Select _corrective_saccades file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
data = load([file_path file_name],'data');
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

%% Build EPHYS Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Alignment events', ' ... ']);

time_reference      = ( ESN_Round(EPHYS.time_30K(1),0.001) : 0.001 : ESN_Round(EPHYS.time_30K(end),0.001) )';
EPHYS.Alignment.time_1K = time_reference;
EPHYS.Alignment.eye_r_vm_filt = interp1(EPHYS.time_30K, EPHYS.eye_r_vm_filt, time_reference, 'linear', 'extrap');

fprintf(' --> Completed. \n');

%% Build BEHAVE Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE Alignment events', ' ... ']);

time_reference      = BEHAVE.time_1K;
BEHAVE.Alignment.time_1K = time_reference;
BEHAVE.Alignment.eye_r_vm_filt = interp1(BEHAVE.time_1K, BEHAVE.eye_r_vm_filt, time_reference, 'linear', 'extrap');

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

EPHYS_state_combined       = EPHYS_eye_r_vm_filt;
BEHAVE_state_combined      = BEHAVE_eye_r_vm_filt;

% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_eye_r_vm_filt, BEHAVE_eye_r_vm_filt); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

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

% %% state_combined: find the Dynamic Time Warp (DTW) between 2 time series
EPHYS_EB_xcorr_state_combined_1K_filt  = EPHYS_EB_xcorr_state_combined_1K;
BEHAVE_EB_xcorr_state_combined_1K_filt = BEHAVE_EB_xcorr_state_combined_1K;
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
%     [~,ix,iy] = dtw(EPHYS_EB_state_combined_chunk,BEHAVE_EB_state_combined_chunk, 15, 'absolute');  % allow upto 15ms warp
    ix = 1 : 1 : length(EPHYS_EB_state_combined_chunk);
    iy = 1 : 1 : length(BEHAVE_EB_state_combined_chunk);
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
fprintf(' --> Completed. \n');

if EPHYS.debug_figures
clf(figure(1));subplot(2,1,1);plot(EPHYS_eye_r_vm_filt);subplot(2,1,2);plot(BEHAVE_eye_r_vm_filt);
clf(figure(2));subplot(2,1,1);plot(EPHYS_state_combined, '.-');subplot(2,1,2);plot(BEHAVE_state_combined, '.-');
clf(figure(3));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_1K, '.-');
clf(figure(4));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_1K), '.-';
clf(figure(5));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_1K(EB_ind_convert_from_BEHAVE_to_EPHYS), '.-');
end
disp(['length EPHYS: ' num2str( EPHYS.Alignment.time_1K(end)-EPHYS.Alignment.time_1K(1) )])
disp(['length BEHAVE: ' num2str( BEHAVE.Alignment.time_1K(end)-BEHAVE.Alignment.time_1K(1) )])
disp(['xcorr diff: ' num2str( sample_diff )])

%% ALIGN EPHYS and BEHAVE photodiode_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE photodiode_combined', ' ... ']);
EPHYS_time_1K              = EPHYS.Alignment.time_1K;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_1K             = BEHAVE.Alignment.time_1K;
EPHYS_eye_r_vm_filt  = EPHYS.Alignment.eye_r_vm_filt;
BEHAVE_eye_r_vm_filt = BEHAVE.Alignment.eye_r_vm_filt;
EPHYS_eye_r_vm_filt(EPHYS_eye_r_vm_filt < 100) = 0;
BEHAVE_eye_r_vm_filt(BEHAVE_eye_r_vm_filt < 100) = 0;
EPHYS_eye_r_vm_filt(EPHYS_eye_r_vm_filt > 1500) = 1500;
BEHAVE_eye_r_vm_filt(BEHAVE_eye_r_vm_filt > 1500) = 1500;

EPHYS_photodiode_combined       = EPHYS_eye_r_vm_filt;
BEHAVE_photodiode_combined      = BEHAVE_eye_r_vm_filt;

% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_eye_r_vm_filt, BEHAVE_eye_r_vm_filt); % cross-correlate signals with each other
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
EPHYS_PD_xcorr_photodiode_combined_1K_filt  = EPHYS_PD_xcorr_photodiode_combined_1K;
BEHAVE_PD_xcorr_photodiode_combined_1K_filt = BEHAVE_PD_xcorr_photodiode_combined_1K;
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
%     [~,ix,iy] = dtw(EPHYS_PD_photodiode_combined_chunk,BEHAVE_PD_photodiode_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    ix = 1 : 1 : length(EPHYS_PD_photodiode_combined_chunk);
    iy = 1 : 1 : length(BEHAVE_PD_photodiode_combined_chunk);
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

if EPHYS.debug_figures
clf(figure(2));subplot(2,1,1);plot(EPHYS_photodiode_combined, '.-');subplot(2,1,2);plot(BEHAVE_photodiode_combined, '.-');
clf(figure(3));subplot(2,1,1);plot(EPHYS_PD_xcorr_photodiode_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_PD_xcorr_photodiode_combined_1K, '.-');
clf(figure(4));subplot(2,1,1);plot(EPHYS_PD_xcorr_photodiode_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_PD_xcorr_photodiode_combined_1K), '.-';
clf(figure(5));subplot(2,1,1);plot(EPHYS_PD_xcorr_photodiode_combined_1K, '.-');subplot(2,1,2);plot(BEHAVE_PD_xcorr_photodiode_combined_1K(PD_ind_convert_from_BEHAVE_to_EPHYS), '.-');
end

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
