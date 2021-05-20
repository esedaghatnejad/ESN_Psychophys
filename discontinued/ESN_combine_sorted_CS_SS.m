function ESN_combine_sorted_CS_SS
%% load EPHYS sorted DATA
[file_name, file_path] = uigetfile([pwd filesep '*_CH*.mat'], 'Select _CH#.mat file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_0 = load([file_path filesep file_name]);
EPHYS.CH_file_name = file_name;
EPHYS.CH_file_path = file_path;
[~, file_name, ~]  = fileparts(EPHYS.CH_file_name);
EPHYS.file_name = file_name;
fprintf(' --> Completed. \n')

[file_name,file_path] = uigetfile([EPHYS.CH_file_path filesep EPHYS.file_name '_CS.csv'], 'Select _CS.csv file');
fprintf(['Loading ', file_name, ' ... ']);
sample_freq = 30e3;
EPHYS.CS_data.CS_ind   = readmatrix([file_path filesep file_name]);
EPHYS.CS_data.CS_time  = EPHYS.CS_data.CS_ind/sample_freq + EPHYS.CH_0.ch_time(1) - (1/sample_freq);
fprintf(' --> Completed. \n')

[file_name,file_path] = uigetfile([EPHYS.CH_file_path filesep EPHYS.file_name '_SS.csv'], 'Select _SS.csv file');
fprintf(['Loading ', file_name, ' ... ']);
sample_freq = 30e3;
EPHYS.SS_data.SS_ind   = readmatrix([file_path filesep file_name]);
EPHYS.SS_data.SS_time  = EPHYS.SS_data.SS_ind/sample_freq + EPHYS.CH_0.ch_time(1) - (1/sample_freq);
fprintf(' --> Completed. \n')

[file_name, file_path] = uigetfile([EPHYS.CH_file_path filesep '*_CMN1.mat'], 'Select _CMN1.mat file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_CMN = load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% Filter Neural Data.
fprintf([EPHYS.file_name ': Filtering', ' ... ']);
clearvars -except EPHYS
sample_freq = 30e3; %Hz
Wn_hipass = [300 3000]/(sample_freq/2);
[b_hipass,a_hipass] = butter(4, Wn_hipass);
[b_lopass,a_lopass] = butter(4, Wn_hipass(1), 'low');

EPHYS.CH_data.ch_time     = EPHYS.CH_0.ch_time(:);
EPHYS.CH_data.ch_info     = EPHYS.CH_0.ch_info;
EPHYS.CH_data.ch_data_cmn = EPHYS.CH_0.ch_data(:) - EPHYS.CH_CMN.ch_data(:);

% band pass filter the data
EPHYS.CH_data.ch_data_hipass = filtfilt( b_hipass, a_hipass, EPHYS.CH_data.ch_data_cmn(:));
% low pass filter and then detrend the data
EPHYS.CH_data.ch_data_lopass = filtfilt( b_lopass, a_lopass, EPHYS.CH_data.ch_data_cmn(:));
EPHYS.CH_data.ch_data_lopass = detrend(  EPHYS.CH_data.ch_data_lopass);

EPHYS = rmfield(EPHYS,{'CH_0', 'CH_CMN'});
fprintf(' --> Completed. \n')

%% Init
clearvars -except EPHYS
EPHYS.CS_data.CS_num_total         = length(EPHYS.CS_data.CS_ind);
EPHYS.CS_data.CS_inds              = nan(EPHYS.CS_data.CS_num_total, 180);
EPHYS.CS_data.CS_waveform   = nan(EPHYS.CS_data.CS_num_total, 180);

EPHYS.SS_data.SS_num_total         = length(EPHYS.SS_data.SS_ind);
EPHYS.SS_data.SS_inds              = nan(EPHYS.SS_data.SS_num_total, 180);
EPHYS.SS_data.SS_waveform   = nan(EPHYS.SS_data.SS_num_total, 180);

%% CS Align
fprintf([EPHYS.file_name ': Aligning CS', ' ... ']);
clearvars -except EPHYS
data_length          = length(EPHYS.CH_data.ch_time);
CS_num_total         = EPHYS.CS_data.CS_num_total;
CS_inds              = nan(CS_num_total, 180);
CS_waveform   = nan(CS_num_total, 180);

for counter_CS = 1 : 1 : CS_num_total
    CS_ind_  = EPHYS.CS_data.CS_ind(counter_CS);
    CS_inds_ = (CS_ind_ - 59) : 1 : (CS_ind_ + 120);
    CS_inds_(CS_inds_<1) = 1;
    CS_inds_(CS_inds_>data_length) = data_length;
    CS_waveform_ = EPHYS.CH_data.ch_data_hipass(CS_inds_, 1)';
    CS_inds(             counter_CS, :) = CS_inds_;
    CS_waveform(  counter_CS, :) = CS_waveform_;
    
    if mod(counter_CS, 1000) == 0
        disp([num2str(counter_CS) ' / ' num2str(CS_num_total)])
    end
end

EPHYS.CS_data.CS_inds              = CS_inds;
EPHYS.CS_data.CS_waveform   = CS_waveform;

fprintf(' --> Completed. \n')

%% SS Align
fprintf([EPHYS.file_name ': Aligning SS', ' ... ']);
clearvars -except EPHYS
data_length          = length(EPHYS.CH_data.ch_time);
SS_num_total         = EPHYS.SS_data.SS_num_total;
SS_inds              = nan(SS_num_total, 180);
SS_waveform   = nan(SS_num_total, 180);

for counter_SS = 1 : 1 : SS_num_total
    SS_ind_  = EPHYS.SS_data.SS_ind(counter_SS);
    SS_inds_ = (SS_ind_ - 59) : 1 : (SS_ind_ + 120);
    SS_inds_(SS_inds_<1) = 1;
    SS_inds_(SS_inds_>data_length) = data_length;
    SS_waveform_ = EPHYS.CH_data.ch_data_hipass(SS_inds_, 1)';
    SS_inds(             counter_SS, :) = SS_inds_;
    SS_waveform(  counter_SS, :) = SS_waveform_;
    
    if mod(counter_SS, 10000) == 0
        disp([num2str(counter_SS) ' / ' num2str(SS_num_total)])
    end
end

EPHYS.SS_data.SS_inds              = SS_inds;
EPHYS.SS_data.SS_waveform   = SS_waveform;
fprintf(' --> Completed. \n')

%% build SSxSS AUTO PROBABILITY
fprintf([EPHYS.file_name ': Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])

bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = (-span_window_size_half+1) : 1 : (span_window_size_half);

CH__.CH_data.ch_time =  EPHYS.CH_data.ch_time - EPHYS.CH_data.ch_time(1);
CH__.SS_data.SS_time =  EPHYS.SS_data.SS_time - EPHYS.CH_data.ch_time(1);
CH__.CS_data.CS_time =  EPHYS.CS_data.CS_time - EPHYS.CH_data.ch_time(1);
CH__.CH_data.ch_time = ESN_Round(CH__.CH_data.ch_time, bin_size_time);
CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
CH__.CH_data.ch_time_reconstruct = CH__.CH_data.ch_time(1) : bin_size_time : CH__.CH_data.ch_time(end);

% SSxSS_AUTO
CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind, 1, length(inds_span)) + repmat(inds_span, length(CH__.SS_data.SS_ind), 1);
CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );

CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
CH__.SS_data.SS_event_trace( 1   ) = false;
CH__.SS_data.SS_event_trace( end ) = false;

CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
EPHYS.SS_data.SSxSS_AUTO = CH__.SS_data.SS_event_reconstruct;

% CSxSS_WITHIN
CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind, 1, length(inds_span)) + repmat(inds_span, length(CH__.CS_data.CS_ind), 1);
CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );

CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
CH__.SS_data.SS_event_trace( 1   ) = false;
CH__.SS_data.SS_event_trace( end ) = false;

CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
EPHYS.SS_data.CSxSS_AUTO = CH__.SS_data.SS_event_reconstruct;

EPHYS.inds_span = inds_span;
EPHYS.bin_size_time = bin_size_time;
fprintf(' --> Completed. \n')

%% Save CH_data CS_data SS_data
fprintf([EPHYS.file_name ': ' 'Saving ' EPHYS.file_name '_sorted.mat' ' ...'])
clearvars -except EPHYS
CH_data = EPHYS.CH_data;
CS_data = EPHYS.CS_data;
SS_data = EPHYS.SS_data;
save([EPHYS.CH_file_path EPHYS.file_name '_sorted.mat'], 'CH_data', 'CS_data', 'SS_data', '-v7.3');
fprintf(' --> Completed. \n')

%% Plot Neural Properties
fig_handle_(1) = figure(1);
clf(fig_handle_(1))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(EPHYS.SS_data.SS_waveform)+std(EPHYS.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.SS_data.SS_waveform)-std(EPHYS.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])

plot((1:180)/30-2, mean(EPHYS.CS_data.CS_waveform)+std(EPHYS.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CS_data.CS_waveform)-std(EPHYS.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])

xlabel('Time (ms)')
ylabel('Voltage (uv)')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

x_axis_ = linspace((-50+(EPHYS.bin_size_time*1000)), 50, length(EPHYS.inds_span))';

prob_value_ = mean(EPHYS.SS_data.SSxSS_AUTO);
prob_value_(round(length(EPHYS.inds_span)/2)) = 0;
y_axis_mean_ = prob_value_;
y_axis_stdv_ = sqrt(size(EPHYS.SS_data.SSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(EPHYS.SS_data.SSxSS_AUTO, 1);
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])

prob_value_ = mean(EPHYS.SS_data.CSxSS_AUTO);
y_axis_mean_ = prob_value_;
y_axis_stdv_ = sqrt(size(EPHYS.SS_data.CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(EPHYS.SS_data.CSxSS_AUTO, 1);
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
ISI_SS = diff(EPHYS.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.1 0.1 0.9]);
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(EPHYS.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.9 0.1 0.1]);
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')

%% Save Fig
fprintf([EPHYS.file_name ': Saving plot CS SS', ' ...'])

ESN_Beautify_Plot
hFig = fig_handle_(1);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

saveas(fig_handle_(1),[EPHYS.CH_file_path EPHYS.file_name '_plot'], 'pdf');
saveas(fig_handle_(1),[EPHYS.CH_file_path EPHYS.file_name '_plot'], 'png');
close(fig_handle_(1))

fprintf(' --> Completed. \n')
%{
%% Scratch, 
%% load data
clearvars -except EPHYS
[file_name, file_path] = uigetfile([pwd filesep '*_CH*.mat'], 'Select _CH#.mat file');
EPHYS.file_name = file_name;
EPHYS.file_path = file_path;
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_0 = load([file_path filesep file_name]);
EPHYS.CH_data.ch_time        = EPHYS.CH_0.CH_data.ch_data_time(:);
EPHYS.CH_data.ch_data_cmn    = EPHYS.CH_0.CH_data.ch_data_cmn(:);
EPHYS.CH_data.ch_data_hipass = EPHYS.CH_0.CH_data.ch_data_hipass(:);
EPHYS.CH_data.ch_data_lopass = EPHYS.CH_0.CH_data.ch_data_lopass(:);
EPHYS.CS_data = EPHYS.CH_0.CS_data;
EPHYS.CS_data.CS_waveform = EPHYS.CH_0.CS_data.CS_waveform;
EPHYS.SS_data = EPHYS.CH_0.SS_data;
EPHYS.SS_data.SS_waveform = EPHYS.CH_0.SS_data.SS_waveform;
fprintf(' --> Completed. \n')

%% Interactive Plot for CS 
clearvars -except EPHYS
h_fig = figure(1);
clf(h_fig);
h_ax_pca  = subplot(1, 2, 1);
h_ax_wave = subplot(1, 2, 2);
hold(h_ax_pca, 'on')
hold(h_ax_wave, 'on')
[EPHYS.CS_data.CS_pca_coeff, EPHYS.CS_data.CS_pca_score, EPHYS.CS_data.CS_pca_latent] = pca(EPHYS.CS_data.CS_waveform);
pca_1 = EPHYS.CS_data.CS_pca_score(:,1);
pca_2 = EPHYS.CS_data.CS_pca_score(:,2);
cs_wave_y_axis = EPHYS.CS_data.CS_waveform';
cs_wave_x_axis = linspace(-2, 4, size(cs_wave_y_axis, 1))';
num_cs = length(pca_1);
is_inside_ = true(num_cs, 1);

if sum(is_inside_) == length(is_inside_)
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'b');
elseif sum(is_inside_) == 0
    h_plot_pca_in   = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
else
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
end

vertices_pca = [(nanmean(pca_1)),               (nanmean(pca_2)+nanstd(pca_2));...
                (nanmean(pca_1)+nanstd(pca_1)), (nanmean(pca_2));...
                (nanmean(pca_1)+nanstd(pca_1)), (nanmean(pca_2)-nanstd(pca_2));...
                (nanmean(pca_1)-nanstd(pca_1)), (nanmean(pca_2)-nanstd(pca_2));...
                (nanmean(pca_1)-nanstd(pca_1)), (nanmean(pca_2));...
                ];
vertices_wave= [(nanmean(cs_wave_x_axis)),                                 (nanmean(nanmean(cs_wave_y_axis))+nanmean(nanstd(cs_wave_y_axis)));...
                (nanmean(cs_wave_x_axis)+nanmean(nanstd(cs_wave_x_axis))), (nanmean(nanmean(cs_wave_y_axis)));...
                (nanmean(cs_wave_x_axis)+nanmean(nanstd(cs_wave_x_axis))), (nanmean(nanmean(cs_wave_y_axis))-nanmean(nanstd(cs_wave_y_axis)));...
                (nanmean(cs_wave_x_axis)-nanmean(nanstd(cs_wave_x_axis))), (nanmean(nanmean(cs_wave_y_axis))-nanmean(nanstd(cs_wave_y_axis)));...
                (nanmean(cs_wave_x_axis)-nanmean(nanstd(cs_wave_x_axis))), (nanmean(nanmean(cs_wave_y_axis)));...
                ];

h_polygon_pca  = drawpolygon(h_ax_pca,  'Position', vertices_pca);
h_polygon_wave = drawpolygon(h_ax_wave, 'Position', vertices_wave);

%% Update based on PCA_Ellipse
hold(h_ax_pca, 'on')
hold(h_ax_wave, 'on')
is_inside_ = inpolygon(pca_1,pca_2,h_polygon_pca.Position(:,1),h_polygon_pca.Position(:,2));
delete([h_plot_pca_in;  h_plot_pca_out ])
delete([h_plot_wave_in; h_plot_wave_out])
if sum(is_inside_) == length(is_inside_)
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'b');
elseif sum(is_inside_) == 0
    h_plot_pca_in   = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
else
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
end

hold(h_ax_pca, 'off')
hold(h_ax_wave, 'off')

%% Update based on Wave_Ellipse
hold(h_ax_pca, 'on')
hold(h_ax_wave, 'on')
for counter_cs_wave = 1 : 1 : size(cs_wave_y_axis, 2)
    is_inside_wave_ = inpolygon(cs_wave_x_axis,cs_wave_y_axis(:,counter_cs_wave),...
        h_polygon_wave.Position(:,1), h_polygon_wave.Position(:,2));
    is_inside_(counter_cs_wave) = sum(is_inside_wave_) > 0;
end
delete([h_plot_pca_in;  h_plot_pca_out ])
delete([h_plot_wave_in; h_plot_wave_out])
if sum(is_inside_) == length(is_inside_)
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'b');
elseif sum(is_inside_) == 0
    h_plot_pca_in   = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
else
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
end

hold(h_ax_pca, 'off')
hold(h_ax_wave, 'off')

%% Delete selected points, Delete Selected is_inside_
hold(h_ax_pca, 'on')
hold(h_ax_wave, 'on')
pca_1 = pca_1(~is_inside_);
pca_2 = pca_2(~is_inside_);
cs_wave_y_axis = cs_wave_y_axis(:,~is_inside_);
num_cs = length(pca_1);
is_inside_ = false(num_cs, 1);

delete([h_plot_pca_in;  h_plot_pca_out ])
delete([h_plot_wave_in; h_plot_wave_out])
if sum(is_inside_) == length(is_inside_)
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'b');
elseif sum(is_inside_) == 0
    h_plot_pca_in   = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
else
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
end

hold(h_ax_pca, 'off')
hold(h_ax_wave, 'off')

%% Keep selected points, Delete Selected ~is_inside_
hold(h_ax_pca, 'on')
hold(h_ax_wave, 'on')
pca_1 = pca_1(is_inside_);
pca_2 = pca_2(is_inside_);
cs_wave_y_axis = cs_wave_y_axis(:,is_inside_);
num_cs = length(pca_1);
is_inside_ = true(num_cs, 1);

delete([h_plot_pca_in;  h_plot_pca_out ])
delete([h_plot_wave_in; h_plot_wave_out])
if sum(is_inside_) == length(is_inside_)
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'b');
elseif sum(is_inside_) == 0
    h_plot_pca_in   = plot(h_ax_pca,  nan(size(pca_1)),          nan(size(pca_2)),         '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, nan(size(cs_wave_x_axis)), nan(size(cs_wave_y_axis)), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
else
    h_plot_pca_in   = plot(h_ax_pca,  pca_1( is_inside_, :),     pca_2( is_inside_, :), '.r');
    h_plot_pca_out  = plot(h_ax_pca,  pca_1(~is_inside_, :),     pca_2(~is_inside_, :), '.b');
    h_plot_wave_in  = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:, is_inside_), 'r');
    h_plot_wave_out = plot(h_ax_wave, cs_wave_x_axis,            cs_wave_y_axis(:,~is_inside_), 'b');
end

hold(h_ax_pca, 'off')
hold(h_ax_wave, 'off')

%% K-means clustering CS
clearvars -except EPHYS
h_fig = figure(1);
clf(h_fig);
h_ax_wave = subplot(1, 1, 1);
hold(h_ax_wave, 'on')
cs_wave_y_axis = EPHYS.CS_data.CS_waveform;
cs_wave_x_axis = linspace(-2, 4, size(cs_wave_y_axis, 2))';
num_clusters = 4;
idx_cluster = kmeans(cs_wave_y_axis,num_clusters);
h_plot = [];
for counter_cluster = 1 : num_clusters
    h_plot_ = plot(h_ax_wave, cs_wave_x_axis, mean(cs_wave_y_axis(idx_cluster==counter_cluster, :))');
    h_plot = [h_plot; h_plot_];
end
hold(h_ax_wave, 'off')

%% K-means clustering SS
clearvars -except EPHYS
h_fig = figure(1);
clf(h_fig);
h_ax_wave = subplot(1, 1, 1);
hold(h_ax_wave, 'on')
ss_wave_y_axis = EPHYS.SS_data.SS_waveform;
ss_wave_x_axis = linspace(-2, 4, size(ss_wave_y_axis, 2))';
num_clusters = 10;
idx_cluster = kmeans(ss_wave_y_axis,num_clusters);
h_plot = [];
for counter_cluster = 1 : num_clusters
    h_plot_ = plot(h_ax_wave, ss_wave_x_axis, mean(ss_wave_y_axis(idx_cluster==counter_cluster, :))');
    h_plot = [h_plot; h_plot_];
end
hold(h_ax_wave, 'off')

%%
%}



