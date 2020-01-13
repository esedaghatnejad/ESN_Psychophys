%% Load data
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

%% Clear vars
clc; close all;
clearvars -except CH_DATA BEHAVE

%% Extract variables
fprintf(['Extract variables to CH__ ' ' ...'])
CH__.CH_data.ch_time     = CH_DATA.CH_data.CH_time;
CH__.CH_data.CH_data_cmn = CH_DATA.CH_data.CH_data_cmn;
CH__.CS_data.CS_ind      = CH_DATA.CS_data.CS_ind;
fprintf(' --> Completed. \n')

%% Set up the movie.
writerObj = VideoWriter('CS_waveform_animation2.avi'); % Name it.
writerObj.FrameRate = 1; % How many frames per second.
open(writerObj);

%% List of filter params
Wn_hipass_Hz_min_list = [1, 5, 10, 15, 20, 30, 40, 50, 70, 90, 110, 130, 150, 175, 200, 250, 300, 500, 750, 1000];
Wn_hipass_Hz_max_list = [10000, 6000, 3000, 2000, 1000, 800, 600, 500, 400, 300, 250, 200, 150, 100, 80, 60, 50, 40, 30, 20];
num_frames = length(Wn_hipass_Hz_min_list);

%% MAIN LOOP
for counter_frame = 1 : 1 : num_frames
%% Filter
fprintf(['Filter the signal ' ' ...'])
sample_freq = 30e3;

% Wn_hipass_Hz_min = Wn_hipass_Hz_min_list(counter_frame);
% Wn_hipass_Hz_max = 10000;

Wn_hipass_Hz_max = Wn_hipass_Hz_max_list(counter_frame);
Wn_hipass_Hz_min = 1;

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
line_width = 0.5;
line_alpha = 0.05;
line_color = 0.05 * [1 1 1];
y_axis_temp1 = CH__.CS_data.CS_waveform';
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