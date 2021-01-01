%% function ESN_plot_neural_modulation_v3(num_data_set)
function ESN_plot_neural_modulation_v4_static_image

%% load EPHYS sorted DATA
file_path = pwd;
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[file_name, file_path] = uigetfile([file_path '*.psort'], 'Select psort file');

fprintf(['Loading ', file_name, ' ... ']);
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
file_name = EPHYS.CH_sorted_file_name(1:13);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[file_name,file_path] = uigetfile([file_path file_name '_EVE1_aligned.mat'], 'Select EVENT DATA file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);
EPHYS.CH_EVE.EPHYS_time_30K = reshape(EPHYS.CH_EVE.EPHYS_time_30K,[], 1);
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
file_name = EPHYS.CH_sorted_file_name(1:13);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[file_name,file_path] = uigetfile([file_path file_name '_ANALYZED.mat'], 'Select _ANALYZED file');
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

waveform_inds_span = ((-60+1) : 1 : (120));
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

SS_waveform  = reshape(SS_waveform,[], length(waveform_inds_span));
CS_waveform  = reshape(CS_waveform,[], length(waveform_inds_span));

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE TRAIN_DATA
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
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K; % EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K; % 
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

%% Building TRAIN_DATA
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building TRAIN_DATA', ' ... ']);
TRAIN_DATA = struct;
inds_span    = ((-300+1) : 1 : (300))';
amp_edges = [-0.5 2 4 6 8 10 50];
ang_edges = -202.5 : 45 : +202.5;
TRAIN_DATA.inds_span = inds_span;
TRAIN_DATA.amp_edges = amp_edges;
TRAIN_DATA.ang_edges = ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;

ang_edges_last_bin_id = length(ang_edges) - 1;

[~, ~, amp_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_amp_m, amp_edges);
[~, ~, ang_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_ang, ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

variable_list = {'notLowAmp', 'lowAmp', 'all'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};


for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = BEHAVE.SACS_ALL_DATA.(['is_' variable_name]);
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
    ind_indType_name = BEHAVE.SACS_ALL_DATA.(['ind_' indType_name]);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(0, length(inds_span));
            velocity_data = zeros(0, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = spike_train_data(inds_train);
        velocity_data = velocity_trace_data(inds_velocity);
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

fprintf(' --> Completed. \n')

%% Add Neural_Properties
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Add Neural_Properties', ' ... ']);
TRAIN_DATA.Neural_Properties = struct();
if ~isempty(EPHYS.CH_sorted.SS_data.SS_time)
    TRAIN_DATA.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.SS_data.SS_time(end)) - (EPHYS.CH_sorted.SS_data.SS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = TRAIN_DATA.Neural_Properties.SS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
    TRAIN_DATA.Neural_Properties.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform;
else
    TRAIN_DATA.Neural_Properties.SS_num = 0;
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.CS_data.CS_time(end)) - (EPHYS.CH_sorted.CS_data.CS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.SS_time = [];
    TRAIN_DATA.Neural_Properties.SS_waveform = [];
end

if ~isempty(EPHYS.CH_sorted.CS_data.CS_time)
    TRAIN_DATA.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    TRAIN_DATA.Neural_Properties.CS_firing_rate = TRAIN_DATA.Neural_Properties.CS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
    TRAIN_DATA.Neural_Properties.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform;
else
    TRAIN_DATA.Neural_Properties.CS_num = 0;
    TRAIN_DATA.Neural_Properties.CS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.CS_time = [];
    TRAIN_DATA.Neural_Properties.CS_waveform = [];
end
TRAIN_DATA.Neural_Properties.Corr_data = EPHYS.CH_sorted.Corr_data;
fprintf(' --> Completed. \n')

%% avg_over_amplitude
TRAIN_DATA_ang = avg_over_amplitude(TRAIN_DATA);

%% Plot-5 Neural Properties
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
Neural_Properties_data.file_name = file_name;
Neural_Properties_data.file_path = EPHYS.CH_sorted_file_path;

Neural_Properties_data.CH_sorted.SS_data   = EPHYS.CH_sorted.SS_data;
Neural_Properties_data.CH_sorted.CS_data   = EPHYS.CH_sorted.CS_data;
Neural_Properties_data.CH_sorted.Corr_data = EPHYS.CH_sorted.Corr_data;

% figure
Line_Color = lines(7);
fig_num_ = 5;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)+std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)-std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(1,:))

plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)+std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)-std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(7,:))

xlabel('Time (ms)')
ylabel('Voltage (uv)')
title('CS & SS Waveform')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

SSxSS_AUTO = Neural_Properties_data.CH_sorted.Corr_data.SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span);
else
    inds_span =      Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_bin_size_time);
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
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(1,:))

CSxSS_AUTO    = Neural_Properties_data.CH_sorted.Corr_data.CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span);
else
    inds_span      = Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_bin_size_time);
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
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(7,:))

ylim([0 inf])
xlabel('Time (ms)')
ylabel('Probability')
title('X Probability')

% subplot(2,2,3) ISI SS
plot_handle_(3) = subplot(2,2,3);
hold on

edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(Neural_Properties_data.CH_sorted.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(1,:));
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')
title('SS Inter-Spike-Interval')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(Neural_Properties_data.CH_sorted.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(7,:));
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')
title('CS Inter-Spike-Interval')

ESN_Beautify_Plot(fig_handle_(fig_num_), [8.0 8.0])

sgtitle(fig_handle_(fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-2 SS & CS train primSac_onset 
fig_num_ = 2;
variable_name  = 'notLowAmp';
% spikeType_name = 'SS';
indType_name   = 'start';
num_ang_bin = 8;
inds_span    = ((-300+1) : 1 : (300));
subplot_fig_num = [4, 7, 8, 9, 6, 3, 2, 1];

range_SS_Firing = [0 200];
% range_vm        = [0 600];
Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_SS_firing = [0    0.3    0.5];
color_CS = Line_Color(7,:);
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

for counter_ang = 1 : 1 : num_ang_bin
subplot(3,3,subplot_fig_num(counter_ang))
hold on
train_data_logic_SS_ = TRAIN_DATA_ang.(variable_name).(['SS' '_train_' indType_name]){1, counter_ang};
train_data_logic_CS_ = TRAIN_DATA_ang.(variable_name).(['CS' '_train_' indType_name]){1, counter_ang};

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)

xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])
end

% % CS Probab
range_inds_probability = 101:300;
prob_amplitude = zeros(1, 9);
ang_order = [5, 6, 7, 8, 1, 2, 3, 4];
for counter_ang = 1 : 1 : num_ang_bin
train_data_logic_CS_ = TRAIN_DATA_ang.(variable_name).(['CS' '_train_' indType_name]){1, counter_ang};
prob_  = sum( sum(train_data_logic_CS_(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_, 1);
prob_amplitude(1, ang_order(counter_ang)) = prob_;
end
prob_amplitude(1, 9) = prob_amplitude(1, 1);

subplot(3,3,5)
hold on
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(45) cosd(90) cosd(135) cosd(180) cosd(225) cosd(270) cosd(315) cosd(0)] .* prob_amplitude;
y_axis = [sind(0) sind(45) sind(90) sind(135) sind(180) sind(225) sind(270) sind(315) sind(0)] .* prob_amplitude;
x_axis = x_axis(~isnan(x_axis)); y_axis = y_axis(~isnan(y_axis));
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', color_CS)
axis equal
set(gca, 'YColor', color_CS)
set(gca, 'XColor', color_CS)

% %% xlabel
xlabel_text_raster_    = {'Time relative to sac onset (ms)', 'Directions based on sac trajectory'};
xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};

subplot(3,3,5)
xlabel(xlabel_text_CS_probab_);

subplot(3,3,8)
xlabel(xlabel_text_raster_);

ESN_Beautify_Plot(fig_handle_(fig_num_), [8.0 8.0])
sgtitle(fig_handle_(fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Save Fig
response_save_fig = questdlg('Do you want to save the figures?',...
    'Question Dialog','Yes','No','Yes');
if contains(response_save_fig, 'Yes')
fprintf(['Saving plots', ' ...'])
save_file_path = uigetdir(EPHYS.CH_sorted_file_path, 'Select where to save the figures.');
if ~isequal(save_file_path,0)
saveas(fig_handle_(2),[save_file_path filesep file_name '_modulation_sac_onset'], 'pdf');
saveas(fig_handle_(2),[save_file_path filesep file_name '_modulation_sac_onset'], 'png');
saveas(fig_handle_(5),[save_file_path filesep file_name '_properties'], 'pdf');
saveas(fig_handle_(5),[save_file_path filesep file_name '_properties'], 'png');
close(fig_handle_(2))
close(fig_handle_(5))
end
fprintf(' --> Completed. \n')
end % if contains(response_save_fig, 'Yes')

%% Report Properties
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
duration = (EPHYS.CH_EVE.EPHYS_time_30K(end)-EPHYS.CH_EVE.EPHYS_time_30K(1));
numCS = length(EPHYS.CH_sorted.CS_data.CS_ind);
freqCS = numCS/duration;
numSS = length(EPHYS.CH_sorted.SS_data.SS_ind);
freqSS = numSS/duration;
numSac = sum(BEHAVE.SACS_ALL_DATA.is_notLowAmp);
fprintf(['*******************************************' '\n'])
fprintf([file_name '\n'])
fprintf([       'dur'   '\t'        'numCS'   '\t'        'freqCS'   '\t'        'numSS'   '\t'        'freqSS'   '\t'        'numSac'   '\n'])
fprintf([num2str(duration/60,'%.1f') '\t'  num2str(numCS,'%.0f') '\t' num2str(freqCS,'%.2f') '\t' num2str(numSS,'%.0f') '\t' num2str(freqSS,'%.2f') '\t' num2str(numSac,'%.0f') '\n'])

end % function ESN_plot_neural_modulation_v4_static_image

%% function avg_over_amplitude
function TRAIN_DATA_ang = avg_over_amplitude(TRAIN_DATA)
%% Loop over pCells
fprintf(['Building TRAIN_DATA_ang', ' ... ']);
clearvars TRAIN_DATA_ang
variable_list = {'notLowAmp', 'lowAmp', 'all'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(TRAIN_DATA.all.SS_train_start{1,1}, 1);
span_width = size(TRAIN_DATA.all.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_ang = 1 : num_ang_bin
        train_data_all_ang = zeros(0, span_width);
        velocity_data_all_ang = zeros(0, span_width);
        num_saccades_all_ang = 0;
        for counter_amp = 1 : num_amp_bin
            train_data = TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
            train_data_all_ang = vertcat(train_data_all_ang, train_data);
            velocity_data_all_ang = vertcat(velocity_data_all_ang, velocity_data);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        TRAIN_DATA_ang.(variable_name).([spikeType_name '_train_' indType_name]){1, counter_ang} = train_data_all_ang;
        TRAIN_DATA_ang.(variable_name).([spikeType_name '_velocity_' indType_name]){1, counter_ang} = velocity_data_all_ang;
        TRAIN_DATA_ang.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
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
