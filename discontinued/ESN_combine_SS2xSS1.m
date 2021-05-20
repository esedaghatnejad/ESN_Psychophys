function ESN_combine_SS2xSS1
%% Load CH1 dataset
[file_name,file_path] = uigetfile([pwd filesep '*_sorted*.mat'], 'Select SS-1 file');
fprintf(['Loading ', file_name, ' as CH1 ... ']);
EPHYS.CH_1_sorted = load([file_path filesep file_name], 'CH_data', 'SS_data', 'CS_data', 'Corr_data');
EPHYS.CH_1_sorted_file_name = file_name;
EPHYS.CH_1_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% Load CH2 dataset
[file_name,file_path] = uigetfile([file_path filesep '*_sorted*.mat'], 'Select SS-2 file');
fprintf(['Loading ', file_name, ' as CH2 ... ']);
EPHYS.CH_2_sorted = load([file_path filesep file_name], 'CH_data', 'SS_data', 'CS_data', 'Corr_data');
EPHYS.CH_2_sorted_file_name = file_name;
EPHYS.CH_2_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% File Name
file_name_date = EPHYS.CH_1_sorted_file_name(1:13);
file_name_CH1  = EPHYS.CH_1_sorted_file_name(15:16);
file_name_CH2  = EPHYS.CH_2_sorted_file_name(15:16);
file_name_base = [file_name_date '_CH' file_name_CH1 'x' 'CH' file_name_CH2];
EPHYS.file_name_base = file_name_base;

%% SS2xSS1_waveform & SS1xSS2_waveform
EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform = EPHYS.CH_1_sorted.CH_data.CH_hipass(EPHYS.CH_2_sorted.SS_data.SS_inds);
EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform = EPHYS.CH_1_sorted.CH_data.CH_hipass(EPHYS.CH_2_sorted.CS_data.CS_inds);
EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform = EPHYS.CH_2_sorted.CH_data.CH_hipass(EPHYS.CH_1_sorted.SS_data.SS_inds);
EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform = EPHYS.CH_2_sorted.CH_data.CH_hipass(EPHYS.CH_1_sorted.CS_data.CS_inds);

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE
fprintf(['Building SSxSS_CROSS PROBABILITY ' ' ...'])
SS_1_time   = EPHYS.CH_1_sorted.SS_data.SS_time;
SS_2_time   = EPHYS.CH_2_sorted.SS_data.SS_time;
CS_1_time   = EPHYS.CH_1_sorted.CS_data.CS_time;
CS_2_time   = EPHYS.CH_2_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_1_time, SS_2_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_1_sorted.Corr_data.SS2xSS1_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_1_sorted.Corr_data.SS2xSS1_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_1_sorted.Corr_data.SS2xSS1_CROSS         = Corr_data.CS_CSxSS_AUTO;

Corr_data = ESN_correlogram(SS_2_time, SS_1_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_2_sorted.Corr_data.SS1xSS2_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_2_sorted.Corr_data.SS1xSS2_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_2_sorted.Corr_data.SS1xSS2_CROSS         = Corr_data.CS_CSxSS_AUTO;

Corr_data = ESN_correlogram(SS_1_time, CS_2_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_1_sorted.Corr_data.CS2xSS1_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_1_sorted.Corr_data.CS2xSS1_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_1_sorted.Corr_data.CS2xSS1_CROSS         = Corr_data.CS_CSxSS_AUTO;

Corr_data = ESN_correlogram(SS_2_time, CS_1_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_2_sorted.Corr_data.CS1xSS2_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_2_sorted.Corr_data.CS1xSS2_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_2_sorted.Corr_data.CS1xSS2_CROSS         = Corr_data.CS_CSxSS_AUTO;

Corr_data = ESN_correlogram(CS_1_time, CS_2_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_AUTO          = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_1_sorted.Corr_data.CS2xCS1_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_1_sorted.Corr_data.CS2xCS1_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_1_sorted.Corr_data.CS2xCS1_CROSS         = Corr_data.CS_CSxSS_AUTO;

Corr_data = ESN_correlogram(CS_2_time, CS_1_time); % ESN_correlogram(SS_time, CS_time)
EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_AUTO          = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_2_sorted.Corr_data.CS1xCS2_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_2_sorted.Corr_data.CS1xCS2_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_2_sorted.Corr_data.CS1xCS2_CROSS         = Corr_data.CS_CSxSS_AUTO;


fprintf(' --> Completed. \n')

%% Plot-1 SSxSS AUTO and CROSS
% figure
plot_data.fig_num_ = 1;
sub_plot_num_row = 2;
sub_plot_num_col = 5;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(sub_plot_num_row,sub_plot_num_col,1);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS_waveform)+std(EPHYS.CH_1_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS_waveform)-std(EPHYS.CH_1_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH1 Voltage (uv)')
title('SS1 Triggered')

plot_handle_(2) = subplot(sub_plot_num_row,sub_plot_num_col,2);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform)+std(EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform)-std(EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.SS_data.SS2xCH1_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH1 Voltage (uv)')
title('SS2 Triggered')

max_YLim = max([plot_handle_(1).YLim(2) plot_handle_(2).YLim(2)]);
min_YLim = min([plot_handle_(1).YLim(1) plot_handle_(2).YLim(1)]);
set(plot_handle_(1), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(2), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(3) = subplot(sub_plot_num_row,sub_plot_num_col,3);
hold on
SSxSS_AUTO = EPHYS.CH_1_sorted.Corr_data.Corr_SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.Corr_SS_inds_span);
else
    inds_span =      EPHYS.CH_1_sorted.Corr_data.Corr_SS_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.Corr_SS_bin_size_time);
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
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('SS1 Probability')
title('SS1 Triggered')

plot_handle_(4) = subplot(sub_plot_num_row,sub_plot_num_col,[4 5]);
hold on
SS2xSS1_CROSS    = EPHYS.CH_1_sorted.Corr_data.SS2xSS1_CROSS;
if size(SS2xSS1_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.SS2xSS1_inds_span);
else
    inds_span      = EPHYS.CH_1_sorted.Corr_data.SS2xSS1_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.SS2xSS1_bin_size_time);
if (~isempty(SS2xSS1_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(SS2xSS1_CROSS, 1) > 1
        prob_value_ = mean(SS2xSS1_CROSS);
    else
        prob_value_ =      SS2xSS1_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(SS2xSS1_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(SS2xSS1_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(SS2xSS1_CROSS));
    y_axis_mean_ = nan(size(SS2xSS1_CROSS));
    y_axis_stdv_ = nan(size(SS2xSS1_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS1 Probability')
title('SS2 Triggered')
max_YLim = max([plot_handle_(3).YLim(2) plot_handle_(4).YLim(2)]);
set(plot_handle_(3), 'YLim', [0 max_YLim]);
set(plot_handle_(4), 'YLim', [0 max_YLim]);

% subplot(2,2,1) Waveform
plot_handle_(5) = subplot(sub_plot_num_row,sub_plot_num_col,6);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS_waveform)+std(EPHYS.CH_2_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS_waveform)-std(EPHYS.CH_2_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH2 Voltage (uv)')
title('SS2 Triggered')

plot_handle_(6) = subplot(sub_plot_num_row,sub_plot_num_col,7);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform)+std(EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform)-std(EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.SS_data.SS1xCH2_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH2 Voltage (uv)')
title('SS1 Triggered')

max_YLim = max([plot_handle_(5).YLim(2) plot_handle_(6).YLim(2)]);
min_YLim = min([plot_handle_(5).YLim(1) plot_handle_(6).YLim(1)]);
set(plot_handle_(5), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(6), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(7) = subplot(sub_plot_num_row,sub_plot_num_col,8);
hold on
SSxSS_AUTO = EPHYS.CH_2_sorted.Corr_data.Corr_SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.Corr_SS_inds_span);
else
    inds_span =      EPHYS.CH_2_sorted.Corr_data.Corr_SS_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.Corr_SS_bin_size_time);
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
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('SS2 Probability')
title('SS2 Triggered')

plot_handle_(8) = subplot(sub_plot_num_row,sub_plot_num_col,[9 10]);
hold on
SS1xSS2_CROSS    = EPHYS.CH_2_sorted.Corr_data.SS1xSS2_CROSS;
if size(SS1xSS2_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.SS1xSS2_inds_span);
else
    inds_span      = EPHYS.CH_2_sorted.Corr_data.SS1xSS2_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.SS1xSS2_bin_size_time);
if (~isempty(SS1xSS2_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(SS1xSS2_CROSS, 1) > 1
        prob_value_ = mean(SS1xSS2_CROSS);
    else
        prob_value_ =      SS1xSS2_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(SS1xSS2_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(SS1xSS2_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(SS1xSS2_CROSS));
    y_axis_mean_ = nan(size(SS1xSS2_CROSS));
    y_axis_stdv_ = nan(size(SS1xSS2_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS2 Probability')
title('SS1 Triggered')

max_YLim = max([plot_handle_(7).YLim(2) plot_handle_(8).YLim(2)]);
set(plot_handle_(7), 'YLim', [0 max_YLim]);
set(plot_handle_(8), 'YLim', [0 max_YLim]);

ESN_Beautify_Plot
hFig = fig_handle_(plot_data.fig_num_);
figure_size  = [13.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(fig_handle_(plot_data.fig_num_), [EPHYS.file_name_base '_SSxSS'], 'Interpreter', 'none');

%% Plot-2 CSxSS AUTO and CROSS
% figure
plot_data.fig_num_ = 2;
sub_plot_num_row = 2;
sub_plot_num_col = 5;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(sub_plot_num_row,sub_plot_num_col,1);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform)+std(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform)-std(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH1 Voltage (uv)')
title('CS1 Triggered')

plot_handle_(2) = subplot(sub_plot_num_row,sub_plot_num_col,2);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform)+std(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform)-std(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH1 Voltage (uv)')
title('CS2 Triggered')

max_YLim = max([plot_handle_(1).YLim(2) plot_handle_(2).YLim(2)]);
min_YLim = min([plot_handle_(1).YLim(1) plot_handle_(2).YLim(1)]);
set(plot_handle_(1), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(2), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(3) = subplot(sub_plot_num_row,sub_plot_num_col,3);
hold on
CSxSS_AUTO = EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.Corr_CS_inds_span);
else
    inds_span =      EPHYS.CH_1_sorted.Corr_data.Corr_CS_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.Corr_CS_bin_size_time);
if (~isempty(CSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxSS_AUTO, 1) > 1
        prob_value_ = mean(CSxSS_AUTO);
    else
        prob_value_ =      CSxSS_AUTO;
    end
    %prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxSS_AUTO));
    y_axis_mean_ = nan(size(CSxSS_AUTO));
    y_axis_stdv_ = nan(size(CSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('SS1 Probability')
title('CS1 Triggered')

plot_handle_(4) = subplot(sub_plot_num_row,sub_plot_num_col,[4 5]);
hold on
CS2xSS1_CROSS    = EPHYS.CH_1_sorted.Corr_data.CS2xSS1_CROSS;
if size(CS2xSS1_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.CS2xSS1_inds_span);
else
    inds_span      = EPHYS.CH_1_sorted.Corr_data.CS2xSS1_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.CS2xSS1_bin_size_time);
if (~isempty(CS2xSS1_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CS2xSS1_CROSS, 1) > 1
        prob_value_ = mean(CS2xSS1_CROSS);
    else
        prob_value_ =      CS2xSS1_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CS2xSS1_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CS2xSS1_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CS2xSS1_CROSS));
    y_axis_mean_ = nan(size(CS2xSS1_CROSS));
    y_axis_stdv_ = nan(size(CS2xSS1_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS1 Probability')
title('CS2 Triggered')
max_YLim = max([plot_handle_(3).YLim(2) plot_handle_(4).YLim(2)]);
set(plot_handle_(3), 'YLim', [0 max_YLim]);
set(plot_handle_(4), 'YLim', [0 max_YLim]);

% subplot(2,2,1) Waveform
plot_handle_(5) = subplot(sub_plot_num_row,sub_plot_num_col,6);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform)+std(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform)-std(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH2 Voltage (uv)')
title('CS2 Triggered')

plot_handle_(6) = subplot(sub_plot_num_row,sub_plot_num_col,7);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform)+std(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform)-std(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH2 Voltage (uv)')
title('CS1 Triggered')

max_YLim = max([plot_handle_(5).YLim(2) plot_handle_(6).YLim(2)]);
min_YLim = min([plot_handle_(5).YLim(1) plot_handle_(6).YLim(1)]);
set(plot_handle_(5), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(6), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(7) = subplot(sub_plot_num_row,sub_plot_num_col,8);
hold on
CSxSS_AUTO = EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.Corr_CS_inds_span);
else
    inds_span =      EPHYS.CH_2_sorted.Corr_data.Corr_CS_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.Corr_CS_bin_size_time);
if (~isempty(CSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxSS_AUTO, 1) > 1
        prob_value_ = mean(CSxSS_AUTO);
    else
        prob_value_ =      CSxSS_AUTO;
    end
    %prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxSS_AUTO));
    y_axis_mean_ = nan(size(CSxSS_AUTO));
    y_axis_stdv_ = nan(size(CSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('SS2 Probability')
title('CS2 Triggered')

plot_handle_(8) = subplot(sub_plot_num_row,sub_plot_num_col,[9 10]);
hold on
CS1xSS2_CROSS    = EPHYS.CH_2_sorted.Corr_data.CS1xSS2_CROSS;
if size(CS1xSS2_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.CS1xSS2_inds_span);
else
    inds_span      = EPHYS.CH_2_sorted.Corr_data.CS1xSS2_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.CS1xSS2_bin_size_time);
if (~isempty(CS1xSS2_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CS1xSS2_CROSS, 1) > 1
        prob_value_ = mean(CS1xSS2_CROSS);
    else
        prob_value_ =      CS1xSS2_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CS1xSS2_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CS1xSS2_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CS1xSS2_CROSS));
    y_axis_mean_ = nan(size(CS1xSS2_CROSS));
    y_axis_stdv_ = nan(size(CS1xSS2_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS2 Probability')
title('CS1 Triggered')

max_YLim = max([plot_handle_(7).YLim(2) plot_handle_(8).YLim(2)]);
set(plot_handle_(7), 'YLim', [0 max_YLim]);
set(plot_handle_(8), 'YLim', [0 max_YLim]);

ESN_Beautify_Plot
hFig = fig_handle_(plot_data.fig_num_);
figure_size  = [13.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(fig_handle_(plot_data.fig_num_), [EPHYS.file_name_base '_CSxSS'], 'Interpreter', 'none');

%% Plot-3 CSxCS AUTO and CROSS
% figure
plot_data.fig_num_ = 3;
sub_plot_num_row = 2;
sub_plot_num_col = 5;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(sub_plot_num_row,sub_plot_num_col,1);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform)+std(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform)-std(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH1 Voltage (uv)')
title('CS1 Triggered')

plot_handle_(2) = subplot(sub_plot_num_row,sub_plot_num_col,2);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform)+std(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform)-std(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_1_sorted.CS_data.CS2xCH1_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH1 Voltage (uv)')
title('CS2 Triggered')

max_YLim = max([plot_handle_(1).YLim(2) plot_handle_(2).YLim(2)]);
min_YLim = min([plot_handle_(1).YLim(1) plot_handle_(2).YLim(1)]);
set(plot_handle_(1), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(2), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(3) = subplot(sub_plot_num_row,sub_plot_num_col,3);
hold on
CSxCS_AUTO = EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_AUTO;
if size(CSxCS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_inds_span);
else
    inds_span =      EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.Corr_CS_CSxCS_bin_size_time);
if (~isempty(CSxCS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxCS_AUTO, 1) > 1
        prob_value_ = mean(CSxCS_AUTO);
    else
        prob_value_ =      CSxCS_AUTO;
    end
    prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxCS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxCS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxCS_AUTO));
    y_axis_mean_ = nan(size(CSxCS_AUTO));
    y_axis_stdv_ = nan(size(CSxCS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('CS1 Probability')
title('CS1 Triggered')
if plot_handle_(3).YLim(2) == 1
    plot_handle_(3).YLim(2) = 0.001;
end

plot_handle_(4) = subplot(sub_plot_num_row,sub_plot_num_col,[4 5]);
hold on
CS2xCS1_CROSS    = EPHYS.CH_1_sorted.Corr_data.CS2xCS1_CROSS;
if size(CS2xCS1_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_1_sorted.Corr_data.CS2xCS1_inds_span);
else
    inds_span      = EPHYS.CH_1_sorted.Corr_data.CS2xCS1_inds_span;
end
bin_size_time = mean(EPHYS.CH_1_sorted.Corr_data.CS2xCS1_bin_size_time);
if (~isempty(CS2xCS1_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CS2xCS1_CROSS, 1) > 1
        prob_value_ = mean(CS2xCS1_CROSS);
    else
        prob_value_ =      CS2xCS1_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CS2xCS1_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CS2xCS1_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CS2xCS1_CROSS));
    y_axis_mean_ = nan(size(CS2xCS1_CROSS));
    y_axis_stdv_ = nan(size(CS2xCS1_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS1 Probability')
title('CS2 Triggered')
max_YLim = max([plot_handle_(3).YLim(2) plot_handle_(4).YLim(2)]);
set(plot_handle_(3), 'YLim', [0 max_YLim]);
set(plot_handle_(4), 'YLim', [0 max_YLim]);

% subplot(2,2,1) Waveform
plot_handle_(5) = subplot(sub_plot_num_row,sub_plot_num_col,6);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform)+std(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform)-std(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time (ms)')
ylabel('CH2 Voltage (uv)')
title('CS2 Triggered')

plot_handle_(6) = subplot(sub_plot_num_row,sub_plot_num_col,7);
hold on
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform)+std(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform)-std(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_2_sorted.CS_data.CS1xCH2_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
xlabel('Time (ms)')
% ylabel('CH2 Voltage (uv)')
title('CS1 Triggered')

max_YLim = max([plot_handle_(5).YLim(2) plot_handle_(6).YLim(2)]);
min_YLim = min([plot_handle_(5).YLim(1) plot_handle_(6).YLim(1)]);
set(plot_handle_(5), 'YLim', [min_YLim max_YLim]);
set(plot_handle_(6), 'YLim', [min_YLim max_YLim]);

% subplot(2,2,2) Probablities
plot_handle_(7) = subplot(sub_plot_num_row,sub_plot_num_col,8);
hold on
CSxCS_AUTO = EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_AUTO;
if size(CSxCS_AUTO, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_inds_span);
else
    inds_span =      EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.Corr_CS_CSxCS_bin_size_time);
if (~isempty(CSxCS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxCS_AUTO, 1) > 1
        prob_value_ = mean(CSxCS_AUTO);
    else
        prob_value_ =      CSxCS_AUTO;
    end
    prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxCS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxCS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxCS_AUTO));
    y_axis_mean_ = nan(size(CSxCS_AUTO));
    y_axis_stdv_ = nan(size(CSxCS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% ylim([0 inf])
xlabel('Time (ms)')
ylabel('CS2 Probability')
title('CS2 Triggered')
if plot_handle_(7).YLim(2) == 1
    plot_handle_(7).YLim(2) = 0.001;
end

plot_handle_(8) = subplot(sub_plot_num_row,sub_plot_num_col,[9 10]);
hold on
CS1xCS2_CROSS    = EPHYS.CH_2_sorted.Corr_data.CS1xCS2_CROSS;
if size(CS1xCS2_CROSS, 1) > 1
    inds_span = mean(EPHYS.CH_2_sorted.Corr_data.CS1xCS2_inds_span);
else
    inds_span      = EPHYS.CH_2_sorted.Corr_data.CS1xCS2_inds_span;
end
bin_size_time = mean(EPHYS.CH_2_sorted.Corr_data.CS1xCS2_bin_size_time);
if (~isempty(CS1xCS2_CROSS))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CS1xCS2_CROSS, 1) > 1
        prob_value_ = mean(CS1xCS2_CROSS);
    else
        prob_value_ =      CS1xCS2_CROSS;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CS1xCS2_CROSS, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CS1xCS2_CROSS, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CS1xCS2_CROSS));
    y_axis_mean_ = nan(size(CS1xCS2_CROSS));
    y_axis_stdv_ = nan(size(CS1xCS2_CROSS));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

%ylim([0 inf])
xlabel('Time (ms)')
% ylabel('SS2 Probability')
title('CS1 Triggered')

max_YLim = max([plot_handle_(7).YLim(2) plot_handle_(8).YLim(2)]);
set(plot_handle_(7), 'YLim', [0 max_YLim]);
set(plot_handle_(8), 'YLim', [0 max_YLim]);

ESN_Beautify_Plot
hFig = fig_handle_(plot_data.fig_num_);
figure_size  = [13.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(fig_handle_(plot_data.fig_num_), [EPHYS.file_name_base '_CSxCS'], 'Interpreter', 'none');

%% Save combined file
file_name = EPHYS.file_name_base;
[save_file_name,save_file_path] = uiputfile([EPHYS.CH_1_sorted_file_path filesep file_name], 'Select where to save the _sorted mat-file.');
fprintf(['Saving ' save_file_name ' ... ']);
if isequal(save_file_name,0)
    fprintf(' --> Cancelled. \n');
    return;
end
saveas(fig_handle_(1),[save_file_path filesep save_file_name '_SSxSS'], 'pdf');
saveas(fig_handle_(1),[save_file_path filesep save_file_name '_SSxSS'], 'png');
saveas(fig_handle_(2),[save_file_path filesep save_file_name '_CSxSS'], 'pdf');
saveas(fig_handle_(2),[save_file_path filesep save_file_name '_CSxSS'], 'png');
saveas(fig_handle_(3),[save_file_path filesep save_file_name '_CSxCS'], 'pdf');
saveas(fig_handle_(3),[save_file_path filesep save_file_name '_CSxCS'], 'png');

CH_1.CS_data   = EPHYS.CH_1_sorted.CS_data;
CH_1.SS_data   = EPHYS.CH_1_sorted.SS_data;
CH_1.Corr_data = EPHYS.CH_1_sorted.Corr_data;
CH_1.file_name = EPHYS.CH_1_sorted_file_name;
CH_1.file_path = EPHYS.CH_1_sorted_file_path;
CH_2.CS_data   = EPHYS.CH_2_sorted.CS_data;
CH_2.SS_data   = EPHYS.CH_2_sorted.SS_data;
CH_2.Corr_data = EPHYS.CH_2_sorted.Corr_data;
CH_2.file_name = EPHYS.CH_2_sorted_file_name;
CH_2.file_path = EPHYS.CH_2_sorted_file_path;

save([save_file_path filesep save_file_name], 'CH_1','CH_2','-v7.3');

fprintf(' --> Completed. \n');

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
