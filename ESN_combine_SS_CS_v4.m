function ESN_combine_SS_CS_v4
%% Load SS dataset
file_path = [pwd filesep];
[file_name, file_path] = uigetfile([file_path '*.psort'], 'Select SS file');
fprintf(['Loading ', file_name, ' ... ']);
DATA_PSORT_SS = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

% [file_name,file_path] = uigetfile([pwd filesep '*_sorted*.mat'], 'Select SS file');
% fprintf(['Loading ', file_name, ' as SS ... ']);
% EPHYS.CH_SS_sorted = load([file_path filesep file_name], 'CH_data', 'SS_data');
% EPHYS.CH_sorted_file_name = file_name;
% EPHYS.CH_sorted_file_path = file_path;
% fprintf(' --> Completed. \n')

%% Load CS dataset
[file_name, file_path] = uigetfile([file_path '*.psort'], 'Select SS file');
fprintf(['Loading ', file_name, ' ... ']);
DATA_PSORT_CS = Psort_read_psort([file_path file_name]);
fprintf(' --> Completed. \n')

% [file_name,file_path] = uigetfile([file_path filesep '*_sorted*.mat'], 'Select CS file');
% fprintf(['Loading ', file_name, ' as CS ... ']);
% EPHYS.CH_CS_sorted = load([file_path filesep file_name], 'CS_data');
% fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted from DATA_PSORT
waveform_inds_span = ((-60+1) : 1 : (120));

ch_data = double(DATA_PSORT_SS.topLevel_data.ch_data);
ch_time = double(DATA_PSORT_SS.topLevel_data.ch_time);
EPHYS.CH_sorted.CH_data.ch_data = ch_data;
EPHYS.CH_sorted.CH_data.ch_time = ch_time;

ch_data = double(DATA_PSORT_SS.topLevel_data.ch_data);
ch_time = double(DATA_PSORT_SS.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT_SS.topLevel_data.ss_index)));
SS_time = ch_time(SS_index);
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
EPHYS.CH_sorted.SS_data.SS_ind  = SS_index;
EPHYS.CH_sorted.SS_data.SS_inds = SS_inds;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;

ch_data = double(DATA_PSORT_CS.topLevel_data.ch_data);
ch_time = double(DATA_PSORT_CS.topLevel_data.ch_time);
CS_index = find(logical(double(DATA_PSORT_CS.topLevel_data.cs_index)));
CS_time = ch_time(CS_index);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
CS_waveform = ch_data(CS_inds);
EPHYS.CH_sorted.CS_data.CS_ind  = CS_index;
EPHYS.CH_sorted.CS_data.CS_inds = CS_inds;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;

%% Override CS_waveform
EPHYS.CH_sorted.CS_data.CS_waveform = EPHYS.CH_sorted.CH_data.ch_data(EPHYS.CH_sorted.CS_data.CS_inds);

%% Build CS_pca
[~, pca_score, ~] = pca(EPHYS.CH_sorted.CS_data.CS_waveform);
EPHYS.CH_sorted.CS_data.CS_pca_1 = pca_score(:,1);
EPHYS.CH_sorted.CS_data.CS_pca_2 = pca_score(:,2);

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

%% Plot-5 Neural Properties
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

[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
sgtitle(fig_handle_(plot_data.fig_num_), file_name, 'Interpreter', 'none');

%% Save combined file
[save_file_name,save_file_path] = uiputfile([EPHYS.CH_sorted_file_path filesep EPHYS.CH_sorted_file_name], 'Select where to save the _sorted mat-file.');
fprintf(['Saving ' save_file_name ' ... ']);
if isequal(save_file_name,0)
    fprintf(' --> Cancelled. \n');
    return;
end
CH_data = EPHYS.CH_sorted.CH_data;
CS_data = EPHYS.CH_sorted.CS_data;
SS_data = EPHYS.CH_sorted.SS_data;
Corr_data = EPHYS.CH_sorted.Corr_data;
save([save_file_path filesep save_file_name], 'CH_data','CS_data','SS_data','Corr_data','-v7.3');
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
