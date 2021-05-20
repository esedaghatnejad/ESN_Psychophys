function temp_ESN_cFN_combine(num_mat_files)
%% Load DATA
fprintf(['Loading DATA', ' ... ']);
file_path = pwd;
clearvars DATA;
for counter_file = 1 : num_mat_files
    [file_name,file_path] = uigetfile([file_path filesep '*_modulation.mat'], 'Select _modulation file');
    DATA_ = load([file_path filesep file_name], 'DATA');
    DATA(counter_file) = DATA_.DATA;
end
fprintf(' --> Completed. \n')

%% Combine DATA
sample_freq = 50e3;
norm_dist_sigma = 2.5e-3; % seconds
norm_dist_x = -(4*norm_dist_sigma*sample_freq):1:(4*norm_dist_sigma*sample_freq);
norm_dist_y = normpdf(norm_dist_x,0, norm_dist_sigma*sample_freq);
event_spike_1K_cue_onset = [];
event_spike_1K_sac_onset = [];
ifr_spike_1K_cue_onset = [];
ifr_spike_1K_sac_onset = [];
for counter_file = 1 : num_mat_files
    train_data_logic = DATA(counter_file).event_spike_1K(DATA(counter_file).inds_cue_onset_1K);
    event_spike_1K_cue_onset = [event_spike_1K_cue_onset; train_data_logic];
    
    train_data_logic = DATA(counter_file).event_spike_1K(DATA(counter_file).inds_sac_onset_1K);
    event_spike_1K_sac_onset = [event_spike_1K_sac_onset; train_data_logic];
    
    time_50K = DATA(counter_file).time_50K;
    time_1K = DATA(counter_file).time_1K;
    event_spike_50K = DATA(counter_file).ss_index;
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
%     instant_firing_rate = [instant_firing_rate; instant_firing_rate(end);];
%     instant_firing_rate_50K = zeros(size(time_50K));
%     instant_firing_rate_50K(event_spike_50K) = instant_firing_rate;
    instant_firing_rate_50K_conv = conv(instant_firing_rate_50K, norm_dist_y, 'same');
    instant_firing_rate_1K_conv = interp1(time_50K, instant_firing_rate_50K_conv, time_1K);
    
    ifr_data = instant_firing_rate_1K_conv(DATA(counter_file).inds_cue_onset_1K);
    ifr_spike_1K_cue_onset = [ifr_spike_1K_cue_onset; ifr_data];
    
    ifr_data = instant_firing_rate_1K_conv(DATA(counter_file).inds_sac_onset_1K);
    ifr_spike_1K_sac_onset = [ifr_spike_1K_sac_onset; ifr_data];
end

%% Plot Raster data
% plot cue_onset
hFig = figure(1);
clf(hFig);
subplot(3,2,[1 3])
hold on
train_data_logic = event_spike_1K_cue_onset;
inds_span = DATA(1).inds_span_cue_onset;
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trial Number')
subplot(3,2,5)
hold on
% plot(inds_span, mean(train_data_logic)*1e3)
plot(inds_span, mean(ifr_spike_1K_cue_onset), 'LineWidth', 2)
grid('on')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([0 250])
ylabel('Firing Rate')
xlabel('Cue Onset')

subplot(3,2,[2 4])
hold on
train_data_logic = event_spike_1K_sac_onset;
inds_span = DATA(1).inds_span_sac_onset;
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Saccade Number')
subplot(3,2,6)
hold on
% plot(inds_span, mean(train_data_logic)*1e3)
plot(inds_span, mean(ifr_spike_1K_sac_onset), 'LineWidth', 2)
grid('on')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([0 250])
ylabel('Firing Rate')
xlabel('Sac Onset')

% sgtitle(hFig, CH_DATA_PSORT.file_name, 'Interpreter', 'none');

% ESN_Beautify_Plot(hFig, [8.5 11])
% fprintf(['Saving plots', ' ...'])
% file_name = CH_DATA_PSORT.file_name;
% file_path = CH_DATA_PSORT.file_path;
% saveas(hFig,[file_path filesep file_name '_modulation.png'], 'png');
% close(hFig)
% fprintf(' --> Completed. \n')

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
