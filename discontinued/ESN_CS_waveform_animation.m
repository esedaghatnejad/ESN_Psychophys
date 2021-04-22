%% Load data
%{
clc;
clear;
fprintf(['Loading the mat file ' ' ...'])
% file_path_ = '/home/kkarbasi/cssorter/125d_data_sorted/2019-04/2019-04-24/2019-04-24_15-33-57/analyzed_data';
% file_name_ = '190424_153357_05_sorted.mat';
% file_path_ = '/home/kkarbasi/cssorter/125d_data_sorted/2019-04/2019-04-10/2019-04-10_14-45-41/analyzed_data';
% file_name_ =  '190410_144541_04_sorted.mat';
file_path_ = '/home/kkarbasi/cssorter/125d_data_sorted/2019-04/2019-04-03/2019-04-03_15-33-06/analyzed_data';
file_name_ =  '190403_153306_03_sorted';
CH_DATA = load([file_path_ filesep file_name_]);
fprintf(' --> Completed. \n')
%}

%% Clear vars
clc; close all;
clearvars -except CH_DATA BEHAVE

%% Extract variables
fprintf(['Extract variables to CH__ ' ' ...'])
CH__.CH_data.ch_time     = CH_DATA.CH_data.CH_time;
CH__.CH_data.CH_data_cmn = CH_DATA.CH_data.CH_data_cmn;
CH__.CS_data.CS_ind      = CH_DATA.CS_data.CS_ind;
fprintf(' --> Completed. \n')

%% Correlogram
fprintf(['Correlogram ' ' ...'])
SS_time   = CH_DATA.SS_data.SS_time;
CS_time   = CH_DATA.CS_data.CS_time;
CH__.Corr_data = ESN_correlogram(SS_time, CS_time); % ESN_correlogram(SS_time, CS_time) 
% Corr_data.SS_inds_span; Corr_data.SS_bin_size_time; Corr_data.SS_SSxSS_AUTO; Corr_data.CS_inds_span; Corr_data.CS_bin_size_time; Corr_data.CS_CSxSS_AUTO;

CS_inds_span  = CH__.Corr_data.CS_inds_span(1, :);
CS_CSxSS_AUTO = CH__.Corr_data.CS_CSxSS_AUTO;

inds_span_00_10 = (CS_inds_span >= 00) & (CS_inds_span < 09);
inds_span_10_15 = (CS_inds_span >= 09) & (CS_inds_span < 16);
inds_span_15_20 = (CS_inds_span >= 16) & (CS_inds_span < 50);

% inds_span_00_10 = (CS_inds_span >= 00) & (CS_inds_span < 13);
% inds_span_10_15 = (CS_inds_span >= 13) & (CS_inds_span < 17);
% inds_span_15_20 = (CS_inds_span >= 17) & (CS_inds_span < 50);

% inds_span_00_10 = (CS_inds_span >= 00) & (CS_inds_span < 15);
% inds_span_10_15 = (CS_inds_span >= 15) & (CS_inds_span < 18);
% inds_span_15_20 = (CS_inds_span >= 18) & (CS_inds_span < 50);

inds_SS_00_10 = sum(CS_CSxSS_AUTO(:,inds_span_00_10), 2) > 0 ;
inds_SS_10_15 = sum(CS_CSxSS_AUTO(:,inds_span_10_15), 2) > 0 ;
inds_SS_15_20 = sum(CS_CSxSS_AUTO(:,inds_span_15_20), 2) > 0 ;

inds_SS_15_20(inds_SS_10_15) = false;
inds_SS_15_20(inds_SS_00_10) = false;
inds_SS_10_15(inds_SS_00_10) = false;
fprintf(' --> Completed. \n')

%% Set up the movie.
% filter_mode = 'hi';
filter_mode = 'lo';
writerObj = VideoWriter(['CS_waveform_animation_' filter_mode '.avi']); % Name it.
writerObj.FrameRate = 1; % How many frames per second.
open(writerObj);

%% List of filter params
Wn_hipass_Hz_min_list = [1, 5, 10, 15, 20, 30, 40, 50, 70, 90, 110, 130, 150, 175, 200, 250, 300, 500, 750, 1000];
Wn_hipass_Hz_max_list = [10000, 6000, 3000, 2000, 1000, 800, 600, 500, 400, 300, 250, 200, 150, 100, 80, 60, 50, 40, 30, 20];
num_frames = length(Wn_hipass_Hz_min_list);

%% Plot empty fig
fig_num = 1;
fig_handle = figure(fig_num);
clf(fig_handle)
subplot(1,2,1)
hold on;
xlabel('Time (ms)')
ylabel('Voltage (uV)')
title('Raw trace')
subplot(1,2,2)
hold on
xlabel('Time (ms)')
title('Mean trace')

sgtitle('Filtering and CS Waveform', 'Interpreter', 'none');

ESN_Beautify_Plot;
hFig = fig_handle;
figure_size  = [10.0 6.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');
%% MAIN LOOP
for counter_frame = 1 : 1 : num_frames
%% Filter
fprintf(['Filter the signal ' ' ...'])
sample_freq = 30e3;

if contains(filter_mode, 'lo')
    Wn_hipass_Hz_min = Wn_hipass_Hz_min_list(counter_frame);
    Wn_hipass_Hz_max = 10000;
elseif contains(filter_mode, 'hi')
    Wn_hipass_Hz_min = 1;
    Wn_hipass_Hz_max = Wn_hipass_Hz_max_list(counter_frame);
else
    error('Incorrect filter_mode');
end


Wn_hipass_Hz = [Wn_hipass_Hz_min Wn_hipass_Hz_max];
Wn_hipass = Wn_hipass_Hz/(sample_freq/2);
[b_hipass_max,a_hipass_max] = butter(4, Wn_hipass(2), 'low');
[b_hipass_min,a_hipass_min] = butter(4, Wn_hipass(1), 'high');
CH__.CH_data.CH_hipass = filtfilt( b_hipass_max, a_hipass_max, CH__.CH_data.CH_data_cmn);
CH__.CH_data.CH_hipass = filtfilt( b_hipass_min, a_hipass_min, CH__.CH_data.CH_hipass);
fprintf(' --> Completed. \n')
%% Extract CS_waveform
fprintf(['Extract CS_waveform ' ' ...'])
inds_span = ((-(5*30)+1) : 1 : ((20*30)))';
CH__.CS_data.CS_inds = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
CH__.CS_data.CS_inds( CH__.CS_data.CS_inds < 1 ) = 1;
CH__.CS_data.CS_inds( CH__.CS_data.CS_inds > length( CH__.CH_data.ch_time ) ) = length( CH__.CH_data.ch_time );

CH__.CS_data.CS_waveform = CH__.CH_data.CH_hipass(CH__.CS_data.CS_inds);
fprintf(' --> Completed. \n')
%% Plot
fprintf(['Plot data ' ' ...'])
x_axis_ticks = inds_span(:) / 30;
% y_lim_ = [-1000; 800];
y_lim_ = [-500; 400];
x_lim_ = [min(x_axis_ticks) max(x_axis_ticks)];
fig_num = 1;
fig_handle = figure(fig_num);
clf(fig_handle)
subplot(1,2,1)
hold on;
line_width = 0.5;
line_alpha = 0.05;
line_color = [0.05 0.05 0.05];
y_axis_temp1 = CH__.CS_data.CS_waveform(inds_SS_15_20, :)';
x_axis_temp1 = repmat(x_axis_ticks, 1, size(y_axis_temp1, 2));
y_axis_temp2 = [y_axis_temp1; nan(1, size(y_axis_temp1, 2))];
x_axis_temp2 = [x_axis_temp1; nan(1, size(y_axis_temp1, 2))];
y_axis_temp3 = y_axis_temp2(:);
x_axis_temp3 = x_axis_temp2(:);
% plot(x_axis_temp3, y_axis_temp3, '-k');
patch(x_axis_temp3, y_axis_temp3,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',line_color);
line_color = [0.05 0.05 0.95];
y_axis_temp1 = CH__.CS_data.CS_waveform(inds_SS_10_15, :)';
x_axis_temp1 = repmat(x_axis_ticks, 1, size(y_axis_temp1, 2));
y_axis_temp2 = [y_axis_temp1; nan(1, size(y_axis_temp1, 2))];
x_axis_temp2 = [x_axis_temp1; nan(1, size(y_axis_temp1, 2))];
y_axis_temp3 = y_axis_temp2(:);
x_axis_temp3 = x_axis_temp2(:);
% plot(x_axis_temp3, y_axis_temp3, '-k');
patch(x_axis_temp3, y_axis_temp3,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',line_color);
line_color = [0.95 0.05 0.05];
y_axis_temp1 = CH__.CS_data.CS_waveform(inds_SS_00_10, :)';
x_axis_temp1 = repmat(x_axis_ticks, 1, size(y_axis_temp1, 2));
y_axis_temp2 = [y_axis_temp1; nan(1, size(y_axis_temp1, 2))];
x_axis_temp2 = [x_axis_temp1; nan(1, size(y_axis_temp1, 2))];
y_axis_temp3 = y_axis_temp2(:);
x_axis_temp3 = x_axis_temp2(:);
% plot(x_axis_temp3, y_axis_temp3, '-k');
patch(x_axis_temp3, y_axis_temp3,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',line_color);
ylim(y_lim_)
xlim(x_lim_)
xlabel('Time (ms)')
ylabel('Voltage (uV)')
title('Raw trace')

subplot(1,2,2)
plot(x_axis_ticks, mean(CH__.CS_data.CS_waveform), '-k', 'linewidth', 2)
ylim(y_lim_)
xlim(x_lim_)
xlabel('Time (ms)')
title('Mean trace')
% text(5, y_lim_(1)+500, {'Filter params', ['Lo-band: ' num2str(Wn_hipass_Hz(1))], ['Hi-band: ' num2str(Wn_hipass_Hz(2))]}, ...
%     'Interpreter', 'none', 'FontSize', 10)
% text(5, y_lim_(1)+200, {'E Sedaghat-Nejad', 'Shadmehr Lab, JHU'}, ...
%     'Interpreter', 'none', 'FontSize', 8)

text(5, y_lim_(1)+300, {'Filter params', ['Lo-band: ' num2str(Wn_hipass_Hz(1))], ['Hi-band: ' num2str(Wn_hipass_Hz(2))]}, ...
    'Interpreter', 'none', 'FontSize', 10)
text(5, y_lim_(1)+100, {'E Sedaghat-Nejad', 'Shadmehr Lab, JHU'}, ...
    'Interpreter', 'none', 'FontSize', 8)

sgtitle('Filtering and CS Waveform', 'Interpreter', 'none');

ESN_Beautify_Plot;
hFig = fig_handle;
figure_size  = [10.0 6.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');


fprintf(' --> Completed. \n')

pause(0.1);

%% Get the frame and store it
frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);
        
end

% Saves the movie.
close(writerObj); % Saves the movie.
fprintf('Saving the movie ... --> Completed. \n')

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
