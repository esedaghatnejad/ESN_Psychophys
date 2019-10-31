%% function ESN_plot_neural_modulation_v2(num_data_set)
function ESN_plot_neural_modulation_v2(num_data_set)
if nargin < 1
    num_data_set = 1;
end
%% Build EPHYS_ and BEHAVE_ for each single dataset
clearvars EPHYS BEHAVE
for counter_dataset = 1 : 1 : num_data_set
    [EPHYS_(counter_dataset), BEHAVE_(counter_dataset)] = build_EPHYS_BEHAVE_single_dataset;
end

%% Plot-1 SS & CS train cue_present 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_cue_present_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_cue_present_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'prim');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 1;
plot_data.xlabel_text_raster_    = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_cue_present;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

[~, file_name, ~]  = fileparts(EPHYS_(1).CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Plot-2 SS & CS train primSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_primSac_onset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_primSac_onset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'prim');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 2;
plot_data.xlabel_text_raster_    = {'Time relative to prim sac onset (ms)', 'Directions based on primary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_onset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

[~, file_name, ~]  = fileparts(EPHYS_(1).CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Plot-3 SS & CS train primSac_offset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_primSac_offset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_primSac_offset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'corr');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 3;
plot_data.xlabel_text_raster_    = {'Time relative to prim sac offset (ms)', 'Directions based on secondary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_offset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

[~, file_name, ~]  = fileparts(EPHYS_(1).CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Plot-4 SS & CS train corrSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_corrSac_onset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_corrSac_onset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'corr');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 4;
plot_data.xlabel_text_raster_    = {'Time relative to corr sac onset (ms)', 'Directions based on secondary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_corrSac_onset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

[~, file_name, ~]  = fileparts(EPHYS_(1).CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Plot-5 Neural Properties
CH_sorted_ = concatenate_dataset(EPHYS_, 'CH_sorted', @horzcat);
EPHYS.CH_sorted.SS_data   = concatenate_dataset(CH_sorted_.SS_data, [], @vertcat);
EPHYS.CH_sorted.CS_data   = concatenate_dataset(CH_sorted_.CS_data, [], @vertcat);
EPHYS.CH_sorted.Corr_data = concatenate_dataset(CH_sorted_.Corr_data, [], @vertcat);

if isfield(EPHYS.CH_sorted.SS_data, 'SS_waveform_hipass')
    EPHYS.CH_sorted.SS_data.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform_hipass;
    EPHYS.CH_sorted.CS_data.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform_hipass;
end

% figure
plot_data.fig_num_ = 5;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform)+std(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform)-std(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])

plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform)+std(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform)-std(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])

xlabel('Time (ms)')
ylabel('Voltage (uv)')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

SSxSS_AUTO = EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_sorted.Corr_data.SS_inds_span);
else
    inds_span =      EPHYS.CH_sorted.Corr_data.SS_inds_span;
end
bin_size_time = mean(EPHYS.CH_sorted.Corr_data.SS_bin_size_time);
if (~isempty(SSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(SSxSS_AUTO, 1) > 1
        prob_value_ = mean(SSxSS_AUTO);
    else
        prob_value_ =      SSxSS_AUTO;
    end
    prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(SSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(SSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(SSxSS_AUTO));
    y_axis_mean_ = nan(size(SSxSS_AUTO));
    y_axis_stdv_ = nan(size(SSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])

CSxSS_AUTO    = EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_sorted.Corr_data.CS_inds_span);
else
    inds_span      = EPHYS.CH_sorted.Corr_data.CS_inds_span;
end
bin_size_time = mean(EPHYS.CH_sorted.Corr_data.CS_bin_size_time);
if (~isempty(CSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxSS_AUTO, 1) > 1
        prob_value_ = mean(CSxSS_AUTO);
    else
        prob_value_ =      CSxSS_AUTO;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxSS_AUTO));
    y_axis_mean_ = nan(size(CSxSS_AUTO));
    y_axis_stdv_ = nan(size(CSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

ylim([0 inf])
xlabel('Time (ms)')
ylabel('Probability')

% subplot(2,2,3) ISI SS
plot_handle_(3) = subplot(2,2,3);
hold on

edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(EPHYS.CH_sorted.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.1 0.1 0.9]);
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(EPHYS.CH_sorted.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.9 0.1 0.1]);
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')

ESN_Beautify_Plot
hFig = fig_handle_(plot_data.fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

[~, file_name, ~]  = fileparts(EPHYS_(1).CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Save Fig
EPHYS.CH_sorted_file_name = EPHYS_(1).CH_sorted_file_name;
EPHYS.CH_sorted_file_path = EPHYS_(1).CH_sorted_file_path;
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
if num_data_set > 1
    file_name = [file_name '_combine_' num2str(num_data_set)];
end

fprintf(['Saving plots', ' ...'])
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path file_name '_modulation_cue_present'], 'pdf');
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path file_name '_modulation_cue_present'], 'png');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_onset'], 'pdf');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_onset'], 'png');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_offset'], 'pdf');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_offset'], 'png');
saveas(fig_handle_(4),[EPHYS.CH_sorted_file_path file_name '_modulation_corrSac_onset'], 'pdf');
saveas(fig_handle_(4),[EPHYS.CH_sorted_file_path file_name '_modulation_corrSac_onset'], 'png');
saveas(fig_handle_(5),[EPHYS.CH_sorted_file_path file_name '_properties'], 'pdf');
saveas(fig_handle_(5),[EPHYS.CH_sorted_file_path file_name '_properties'], 'png');
close(fig_handle_(1))
close(fig_handle_(2))
close(fig_handle_(3))
close(fig_handle_(4))
close(fig_handle_(5))
fprintf(' --> Completed. \n')

%% Report Properties
for counter_dataset = 1 : 1 : num_data_set
    [~, file_name, ~]  = fileparts(EPHYS_(counter_dataset).CH_sorted_file_name);
    duration = (EPHYS_(counter_dataset).CH_EVE.EPHYS_time_15K(end)-EPHYS_(counter_dataset).CH_EVE.EPHYS_time_15K(1));
    numCS = length(EPHYS_(counter_dataset).CH_sorted.CS_data.CS_ind);
    freqCS = numCS/duration;
    numSS = length(EPHYS_(counter_dataset).CH_sorted.SS_data.SS_ind);
    freqSS = numSS/duration;
    numTrial = length(BEHAVE_(counter_dataset).TRIALS_DATA.time_end);
    fprintf(['*******************************************' '\n'])
    fprintf([file_name '\n'])
    fprintf([       'dur'   '\t'        'numCS'   '\t'        'freqCS'   '\t'        'numSS'   '\t'        'freqSS'   '\t'        'numTrial'   '\n'])
    fprintf([num2str(duration/60,'%.1f') '\t'  num2str(numCS,'%.0f') '\t' num2str(freqCS,'%.2f') '\t' num2str(numSS,'%.0f') '\t' num2str(freqSS,'%.2f') '\t' num2str(numTrial,'%.0f') '\n'])
end

end

%% function ESN_raster_plot_axes
function [x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, x_axis_values, line_half_len)
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

%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
smooth_data_ = smooth(data_, 21, 'sgolay', 2);
end

%% function build_EPHYS_BEHAVE_single_dataset
function [EPHYS, BEHAVE] = build_EPHYS_BEHAVE_single_dataset
%% load EPHYS EVENT DATA
[file_name,file_path] = uigetfile([pwd filesep '*_aligned.mat'], 'Select EVENT DATA file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path filesep file_name]);
EPHYS.CH_EVE.EPHYS_time_15K = EPHYS.CH_EVE.EPHYS_time_15K(:);
EPHYS.CH_EVE.EPHYS_time_1K  = EPHYS.CH_EVE.EPHYS_time_1K(:);
EPHYS.CH_EVE.BEHAVE_time_1K = EPHYS.CH_EVE.BEHAVE_time_1K(:);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
% [file_name,file_path] = uigetfile([file_path filesep '*_REDUCED.mat'], 'Select _REDUCED file');
[file_name,file_path] = uigetfile([file_path filesep '*_ANALYZED.mat'], 'Select _ANALYZED file');
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% load EPHYS sorted DATA
[file_name,file_path] = uigetfile([file_path filesep '*_sorted*.mat'], 'Select _sorted* file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE
fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data.CS_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_sorted.Corr_data.CS_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_sorted.Corr_data.SS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO    = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO    = Corr_data.CS_CSxSS_AUTO;

fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

fprintf(' --> Completed. \n')

%% Build BEHAVE_eye_r_vm_filt
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE_eye_r_vm_filt', ' ... ']);
eye_r_vm_filt = cell2mat(BEHAVE.TRIALS_DATA.eye_r_vm_filt(:));
time_1K_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.time_1K(:));
BEHAVE_time_1K     = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_ = length(BEHAVE_time_1K);
BEHAVE_eye_r_vm_filt_1K = nan(size(BEHAVE_time_1K));
time_1K_cell2mat(end+1)    = max([BEHAVE_time_1K(end), time_1K_cell2mat(end)])+1;
counter_time_1K_cell2mat = find(time_1K_cell2mat >= BEHAVE_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_time_1K(counter_time_point);
    if time_ponit_>=time_1K_cell2mat(counter_time_1K_cell2mat)
        BEHAVE_eye_r_vm_filt_1K(counter_time_point) = eye_r_vm_filt(counter_time_1K_cell2mat);
        counter_time_1K_cell2mat = counter_time_1K_cell2mat + 1;
    end
end
EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_1K = BEHAVE_eye_r_vm_filt_1K;
fprintf(' --> Completed. \n')

%% Build Raster Data (inds_span)
clearvars -except EPHYS BEHAVE
fprintf(['Building Raster Plot Data', ' ... ']);
num_trials = length(BEHAVE.TRIALS_DATA.time_end);
BEHAVE_EB_xcorr_time_1K     = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K;
length_time_ = length(BEHAVE_EB_xcorr_time_1K);
BEHAVE_EB_xcorr_time_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_onset  = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    BEHAVE_EB_xcorr_time_cue_present(counter_trial) = BEHAVE.TRIALS_DATA.time_state_cue_present{1,counter_trial}(end);
    ind_primSac_onset_  = (BEHAVE.SACS_PRIM_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_onset_) ~= 1
        ind_primSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_primSac_onset(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_onset_, counter_trial);
    ind_primSac_offset_ = (BEHAVE.SACS_PRIM_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_offset_) ~= 1
        ind_primSac_offset_ = 150;
    end
    BEHAVE_EB_xcorr_time_primSac_offset(counter_trial) = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_offset_, counter_trial);
    ind_corrSac_onset_  = (BEHAVE.SACS_CORR_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_onset_) ~= 1
        ind_corrSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_corrSac_onset(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_onset_, counter_trial);
end
BEHAVE_EB_xcorr_time_cue_present(end+1)    = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_cue_present(end)])+1;
BEHAVE_EB_xcorr_time_primSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_primSac_onset(end)])+1;
BEHAVE_EB_xcorr_time_primSac_offset(end+1) = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_primSac_offset(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_corrSac_onset(end)])+1;
BEHAVE_EB_xcorr_ind_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_onset  = nan(num_trials, 1);
counter_cue_present    = find(BEHAVE_EB_xcorr_time_cue_present    >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_primSac_onset  = find(BEHAVE_EB_xcorr_time_primSac_onset  >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_primSac_offset = find(BEHAVE_EB_xcorr_time_primSac_offset >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_corrSac_onset  = find(BEHAVE_EB_xcorr_time_corrSac_onset  >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_trial_cue_present    = 1;
counter_trial_primSac_onset  = 1;
counter_trial_primSac_offset = 1;
counter_trial_corrSac_onset  = 1;
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_EB_xcorr_time_1K(counter_time_point);
    if time_ponit_>=BEHAVE_EB_xcorr_time_cue_present(counter_cue_present)
        BEHAVE_EB_xcorr_ind_cue_present(counter_trial_cue_present) = counter_time_point;
        counter_cue_present = counter_cue_present + 1;
        counter_trial_cue_present = counter_trial_cue_present + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_onset(counter_primSac_onset)
        BEHAVE_EB_xcorr_ind_primSac_onset(counter_trial_primSac_onset) = counter_time_point;
        counter_primSac_onset = counter_primSac_onset + 1;
        counter_trial_primSac_onset = counter_trial_primSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_offset(counter_primSac_offset)
        BEHAVE_EB_xcorr_ind_primSac_offset(counter_trial_primSac_offset) = counter_time_point;
        counter_primSac_offset = counter_primSac_offset + 1;
        counter_trial_primSac_offset = counter_trial_primSac_offset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_onset(counter_corrSac_onset)
        BEHAVE_EB_xcorr_ind_corrSac_onset(counter_trial_corrSac_onset) = counter_time_point;
        counter_corrSac_onset = counter_corrSac_onset + 1;
        counter_trial_corrSac_onset = counter_trial_corrSac_onset + 1;
    end
end
% convert xcorr to aligned for EPHYS. We find the events on BEHAVE and then
% should convert it to EPHYS for SS & CS (EPHYS related events). We should not convert BEHAVE for
% BEHAVE related events.
EPHYS_EB_aligned_ind_cue_present_1K    = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);
EPHYS_EB_aligned_ind_primSac_onset_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_primSac_onset);
EPHYS_EB_aligned_ind_primSac_offset_1K = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_primSac_offset);
EPHYS_EB_aligned_ind_corrSac_onset_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_onset);
inds_span_cue_present    = ((-100+1) : 1 : (300))';
inds_span_primSac_onset  = ((-300+1) : 1 : (100))';
inds_span_primSac_offset = ((-100+1) : 1 : (300))';
inds_span_corrSac_onset  = ((-300+1) : 1 : (100))';
% Build EPHYS_EB_aligned_inds
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_      = length(EPHYS_time_1K);

EPHYS_EB_aligned_inds_cue_present_1K = repmat( EPHYS_EB_aligned_ind_cue_present_1K(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
EPHYS_EB_aligned_inds_cue_present_1K( EPHYS_EB_aligned_inds_cue_present_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_cue_present_1K( EPHYS_EB_aligned_inds_cue_present_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_onset_1K = repmat( EPHYS_EB_aligned_ind_primSac_onset_1K(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
EPHYS_EB_aligned_inds_primSac_onset_1K( EPHYS_EB_aligned_inds_primSac_onset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_onset_1K( EPHYS_EB_aligned_inds_primSac_onset_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_offset_1K = repmat( EPHYS_EB_aligned_ind_primSac_offset_1K(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
EPHYS_EB_aligned_inds_primSac_offset_1K( EPHYS_EB_aligned_inds_primSac_offset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_offset_1K( EPHYS_EB_aligned_inds_primSac_offset_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_onset_1K = repmat( EPHYS_EB_aligned_ind_corrSac_onset_1K(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
EPHYS_EB_aligned_inds_corrSac_onset_1K( EPHYS_EB_aligned_inds_corrSac_onset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_onset_1K( EPHYS_EB_aligned_inds_corrSac_onset_1K > length_time_ ) = length_time_;

% We should not convert BEHAVE for BEHAVE related events.
BEHAVE_EB_aligned_ind_cue_present_1K    = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);
BEHAVE_EB_aligned_ind_primSac_onset_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_primSac_onset);
BEHAVE_EB_aligned_ind_primSac_offset_1K = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_primSac_offset);
BEHAVE_EB_aligned_ind_corrSac_onset_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_onset);
% Build BEHAVE_EB_aligned_inds
BEHAVE_time_1K     = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_      = length(BEHAVE_time_1K);

BEHAVE_EB_aligned_inds_cue_present_1K = repmat( BEHAVE_EB_aligned_ind_cue_present_1K(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
BEHAVE_EB_aligned_inds_cue_present_1K( BEHAVE_EB_aligned_inds_cue_present_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_cue_present_1K( BEHAVE_EB_aligned_inds_cue_present_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_onset_1K = repmat( BEHAVE_EB_aligned_ind_primSac_onset_1K(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
BEHAVE_EB_aligned_inds_primSac_onset_1K( BEHAVE_EB_aligned_inds_primSac_onset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_onset_1K( BEHAVE_EB_aligned_inds_primSac_onset_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_offset_1K = repmat( BEHAVE_EB_aligned_ind_primSac_offset_1K(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
BEHAVE_EB_aligned_inds_primSac_offset_1K( BEHAVE_EB_aligned_inds_primSac_offset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_offset_1K( BEHAVE_EB_aligned_inds_primSac_offset_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_onset_1K = repmat( BEHAVE_EB_aligned_ind_corrSac_onset_1K(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
BEHAVE_EB_aligned_inds_corrSac_onset_1K( BEHAVE_EB_aligned_inds_corrSac_onset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_onset_1K( BEHAVE_EB_aligned_inds_corrSac_onset_1K > length_time_ ) = length_time_;

EPHYS.CH_EVE.EPHYS_EB_aligned_ind_cue_present_1K     = EPHYS_EB_aligned_ind_cue_present_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_onset_1K   = EPHYS_EB_aligned_ind_primSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_offset_1K  = EPHYS_EB_aligned_ind_primSac_offset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_onset_1K   = EPHYS_EB_aligned_ind_corrSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_cue_present_1K    = EPHYS_EB_aligned_inds_cue_present_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_onset_1K  = EPHYS_EB_aligned_inds_primSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_offset_1K = EPHYS_EB_aligned_inds_primSac_offset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_onset_1K  = EPHYS_EB_aligned_inds_corrSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_cue_present_1K     = BEHAVE_EB_aligned_ind_cue_present_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_onset_1K   = BEHAVE_EB_aligned_ind_primSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_offset_1K  = BEHAVE_EB_aligned_ind_primSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_onset_1K   = BEHAVE_EB_aligned_ind_corrSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_cue_present_1K    = BEHAVE_EB_aligned_inds_cue_present_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_onset_1K  = BEHAVE_EB_aligned_inds_primSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_offset_1K = BEHAVE_EB_aligned_inds_primSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_onset_1K  = BEHAVE_EB_aligned_inds_corrSac_onset_1K;
EPHYS.CH_EVE.inds_span_cue_present      = inds_span_cue_present(:)';
EPHYS.CH_EVE.inds_span_primSac_onset    = inds_span_primSac_onset(:)';
EPHYS.CH_EVE.inds_span_primSac_offset   = inds_span_primSac_offset(:)';
EPHYS.CH_EVE.inds_span_corrSac_onset    = inds_span_corrSac_onset(:)';

fprintf(' --> Completed. \n')

end

%% function concatenate_dataset
function upper_field_struct = concatenate_dataset(dataset_, upper_field_name, horz_OR_vert)
if ~isempty(upper_field_name)
dataset = struct(upper_field_name,struct());
field_names_ = fieldnames(dataset_(1).(upper_field_name));
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(upper_field_name).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset.(upper_field_name), field_names_{counter_fields})
            dataset.(upper_field_name).(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(upper_field_name).(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(upper_field_name).(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset.(upper_field_name);
elseif isempty(upper_field_name)
dataset = struct();
field_names_ = fieldnames(dataset_);
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset, field_names_{counter_fields})
            dataset.(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset;
end
end

%% function plot_rasters_data
function fig_handle_ = plot_rasters_data(raster_data, plot_data)
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_CS_probab_ = plot_data.xlabel_text_CS_probab_;
inds_span              = plot_data.inds_span;

train_data_logic_SS_top = raster_data.train_data_logic_SS_top;
train_data_logic_CS_top = raster_data.train_data_logic_CS_top;
velocity_data_top = raster_data.velocity_data_top;
% velocity_data_mean_top = raster_data.velocity_data_mean_top;
% velocity_data_stdv_top = raster_data.velocity_data_stdv_top;

train_data_logic_SS_left = raster_data.train_data_logic_SS_left;
train_data_logic_CS_left = raster_data.train_data_logic_CS_left;
velocity_data_left = raster_data.velocity_data_left;
% velocity_data_mean_left = raster_data.velocity_data_mean_left;
% velocity_data_stdv_left = raster_data.velocity_data_stdv_left;

train_data_logic_SS_right = raster_data.train_data_logic_SS_right;
train_data_logic_CS_right = raster_data.train_data_logic_CS_right;
velocity_data_right = raster_data.velocity_data_right;
% velocity_data_mean_right = raster_data.velocity_data_mean_right;
% velocity_data_stdv_right = raster_data.velocity_data_stdv_right;

train_data_logic_SS_down = raster_data.train_data_logic_SS_down;
train_data_logic_CS_down = raster_data.train_data_logic_CS_down;
velocity_data_down = raster_data.velocity_data_down;
% velocity_data_mean_down = raster_data.velocity_data_mean_down;
% velocity_data_stdv_down = raster_data.velocity_data_stdv_down;


% % CS Probab
range_inds_probability = 101:300;
prob_top   = sum( sum(train_data_logic_CS_top(  :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_top, 1);
prob_left  = sum( sum(train_data_logic_CS_left( :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_left, 1);
prob_right = sum( sum(train_data_logic_CS_right(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_right, 1);
prob_down  = sum( sum(train_data_logic_CS_down( :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_down, 1);
prob_amplitude = [prob_right prob_top prob_left prob_down prob_right];
% % plot xlim and ylim
range_SS_Firing = [0 200];
range_vm        = [0 600];
% % xlabel text
% xlabel_text_raster_ = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
% xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% % figure
% fig_num_ = 1;
fig_handle_ = figure(fig_num_);
clf(fig_handle_)

% % Top Plot Raster
subplot(3,3,2)
hold on
% train_data_logic_SS_top = SS_train_aligned(inds_event(inds_top,:));
[x_axis_SS_top, y_axis_SS_top] = ESN_raster_plot_axes(train_data_logic_SS_top, inds_span, 0.5);
plot(x_axis_SS_top(:), y_axis_SS_top(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_top = CS_train_aligned(inds_event(inds_top,:));
[x_axis_CS_top, y_axis_CS_top] = ESN_raster_plot_axes(train_data_logic_CS_top, inds_span, 3);
plot(x_axis_CS_top(:), y_axis_CS_top(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_top,1)+3)])
ylabel('Trials')
title('Top')

% % Top Plot Velocity
subplot(3,3,1)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean_top = nanmean(velocity_data_top);
velocity_data_stdv_top = nanstd( velocity_data_top);
plot(inds_span, velocity_data_mean_top+velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_top-velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_top,                        '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_top,:));
firing_SS_top = mean(train_data_logic_SS_top) * 1000;
plot(inds_span, ESN_smooth(firing_SS_top), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', [0.1 0.1 0.9])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')
set(gca, 'YColor', [0.1 0.9 0.1])

% % Left Plot Raster
subplot(3,3,4)
hold on
% train_data_logic_SS_left = SS_train_aligned(inds_event(inds_left,:));
[x_axis_SS_left, y_axis_SS_left] = ESN_raster_plot_axes(train_data_logic_SS_left, inds_span, 0.5);
plot(x_axis_SS_left(:), y_axis_SS_left(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_left = CS_train_aligned(inds_event(inds_left,:));
[x_axis_CS_left, y_axis_CS_left] = ESN_raster_plot_axes(train_data_logic_CS_left, inds_span, 3);
plot(x_axis_CS_left(:), y_axis_CS_left(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_left,1)+3)])
ylabel('Trials')
title('Left')

% % Left Plot Velocity
subplot(3,3,7)
hold on
yyaxis left;
% velocity_data_left = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean_left = nanmean(velocity_data_left);
velocity_data_stdv_left = nanstd( velocity_data_left);
plot(inds_span, velocity_data_mean_left+velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_left-velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_left,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_left,:));
firing_SS_left = mean(train_data_logic_SS_left) * 1000;
plot(inds_span, ESN_smooth(firing_SS_left), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', [0.1 0.1 0.9])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')
set(gca, 'YColor', [0.1 0.9 0.1])

% % Right Plot Raster
subplot(3,3,6)
hold on
% train_data_logic_SS_right = SS_train_aligned(inds_event(inds_right,:));
[x_axis_SS_right, y_axis_SS_right] = ESN_raster_plot_axes(train_data_logic_SS_right, inds_span, 0.5);
plot(x_axis_SS_right(:), y_axis_SS_right(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_right = CS_train_aligned(inds_event(inds_right,:));
[x_axis_CS_right, y_axis_CS_right] = ESN_raster_plot_axes(train_data_logic_CS_right, inds_span, 3);
plot(x_axis_CS_right(:), y_axis_CS_right(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_right,1)+3)])
ylabel('Trials')
title('Right')

% % Right Plot Velocity
subplot(3,3,3)
hold on
yyaxis left;
% velocity_data_right = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean_right = nanmean(velocity_data_right);
velocity_data_stdv_right = nanstd( velocity_data_right);
plot(inds_span, velocity_data_mean_right+velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_right-velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_right,                          '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_right,:));
firing_SS_right = mean(train_data_logic_SS_right) * 1000;
plot(inds_span, ESN_smooth(firing_SS_right), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', [0.1 0.1 0.9])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')
set(gca, 'YColor', [0.1 0.9 0.1])

% % Down Plot Raster
subplot(3,3,8)
hold on
% train_data_logic_SS_down = SS_train_aligned(inds_event(inds_down,:));
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic_SS_down, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_down = CS_train_aligned(inds_event(inds_down,:));
[x_axis_CS_down, y_axis_CS_down] = ESN_raster_plot_axes(train_data_logic_CS_down, inds_span, 3);
plot(x_axis_CS_down(:), y_axis_CS_down(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_down,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% % Down Plot Velocity
subplot(3,3,9)
hold on
yyaxis left;
% velocity_data_down = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean_down = nanmean(velocity_data_down);
velocity_data_stdv_down = nanstd( velocity_data_down);
plot(inds_span, velocity_data_mean_down+velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_down-velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_down,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_down,:));
firing_SS_down = mean(train_data_logic_SS_down) * 1000;
plot(inds_span, ESN_smooth(firing_SS_down), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', [0.1 0.1 0.9])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')
set(gca, 'YColor', [0.1 0.9 0.1])

% % Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* prob_amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* prob_amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);
set(gca, 'YColor', [0.9 0.1 0.1])
set(gca, 'XColor', [0.9 0.1 0.1])

ESN_Beautify_Plot
hFig = fig_handle_;
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');
end

%% function single_dataset_raster
function raster_data = single_dataset_raster(EPHYS, BEHAVE, EPHYS_inds_event, BEHAVE_inds_event, prim_OR_corr)
% inds of interest
% % CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_1K;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_1K;
% % eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_1K;
% % check validity of trials
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
% % inds directions
if contains(prim_OR_corr, 'prim')
inds_top   = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) <  0) & ...
             inds_valid;
elseif contains(prim_OR_corr, 'corr')
inds_top   = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) <  0) & ...
             inds_valid;
end

% Top Plot Raster
train_data_logic_SS_top = SS_train_aligned(EPHYS_inds_event(inds_top,:));
train_data_logic_CS_top = CS_train_aligned(EPHYS_inds_event(inds_top,:));
% % Top Plot Velocity
velocity_data_top = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_top,:));

% % Left Plot Raster
train_data_logic_SS_left = SS_train_aligned(EPHYS_inds_event(inds_left,:));
train_data_logic_CS_left = CS_train_aligned(EPHYS_inds_event(inds_left,:));
% % Left Plot Velocity
velocity_data_left = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_left,:));

% % Right Plot Raster
train_data_logic_SS_right = SS_train_aligned(EPHYS_inds_event(inds_right,:));
train_data_logic_CS_right = CS_train_aligned(EPHYS_inds_event(inds_right,:));
% % Right Plot Velocity
velocity_data_right = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_right,:));

% % Down Plot Raster
train_data_logic_SS_down = SS_train_aligned(EPHYS_inds_event(inds_down,:));
train_data_logic_CS_down = CS_train_aligned(EPHYS_inds_event(inds_down,:));
% % Down Plot Velocity
velocity_data_down = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_down,:));

if size(EPHYS_inds_event, 2) == size(train_data_logic_SS_top, 2)
raster_data.train_data_logic_SS_top = train_data_logic_SS_top;
raster_data.train_data_logic_CS_top = train_data_logic_CS_top;
raster_data.velocity_data_top       = velocity_data_top;
else
raster_data.train_data_logic_SS_top = train_data_logic_SS_top';
raster_data.train_data_logic_CS_top = train_data_logic_CS_top';
raster_data.velocity_data_top       = velocity_data_top';
end

if size(EPHYS_inds_event, 2) == size(train_data_logic_SS_left, 2)
raster_data.train_data_logic_SS_left = train_data_logic_SS_left;
raster_data.train_data_logic_CS_left = train_data_logic_CS_left;
raster_data.velocity_data_left       = velocity_data_left;
else
raster_data.train_data_logic_SS_left = train_data_logic_SS_left';
raster_data.train_data_logic_CS_left = train_data_logic_CS_left';
raster_data.velocity_data_left       = velocity_data_left';
end

if size(EPHYS_inds_event, 2) == size(train_data_logic_SS_right, 2)
raster_data.train_data_logic_SS_right = train_data_logic_SS_right;
raster_data.train_data_logic_CS_right = train_data_logic_CS_right;
raster_data.velocity_data_right       = velocity_data_right;
else
raster_data.train_data_logic_SS_right = train_data_logic_SS_right';
raster_data.train_data_logic_CS_right = train_data_logic_CS_right';
raster_data.velocity_data_right       = velocity_data_right';
end

if size(EPHYS_inds_event, 2) == size(train_data_logic_SS_down, 2)
raster_data.train_data_logic_SS_down = train_data_logic_SS_down;
raster_data.train_data_logic_CS_down = train_data_logic_CS_down;
raster_data.velocity_data_down       = velocity_data_down;
else
raster_data.train_data_logic_SS_down = train_data_logic_SS_down';
raster_data.train_data_logic_CS_down = train_data_logic_CS_down';
raster_data.velocity_data_down       = velocity_data_down';
end


end

%% function ESN_correlogram
function Corr_data = ESN_correlogram(SS_time, CS_time)
bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

if (~isempty(CS_time)) && (~isempty(SS_time))
    ch_time_min = min([SS_time(1) CS_time(1)]);
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(CS_time))
    ch_time_min = min(  CS_time(1) );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max(  CS_time(end) ) + 2.0;
    
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(SS_time))
    ch_time_min = min( SS_time(1)  );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max( SS_time(end)  ) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
end

% SSxSS_AUTO
if (~isempty(SS_time))
    CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
    % SSxSS correlogram
    SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
    ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
else
    SSxSS_AUTO       = false(0, length(inds_span(:)'));
    ss_inds_span     = nan(0, length(inds_span(:)'));
    ss_bin_size_time = nan(0, 1);
end

% CSxSS_WITHIN
if (~isempty(CS_time)) && (~isempty(SS_time))
    CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
    % CSxSS correlogram
    CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
    cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
else
    CSxSS_AUTO       = false(0, length(inds_span(:)'));
    cs_inds_span     = nan(0, length(inds_span(:)'));
    cs_bin_size_time = nan(0, 1);
end

Corr_data = struct;
Corr_data.CS_inds_span     = cs_inds_span;
Corr_data.CS_bin_size_time = cs_bin_size_time;
Corr_data.SS_inds_span     = ss_inds_span;
Corr_data.SS_bin_size_time = ss_bin_size_time;
Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
end
